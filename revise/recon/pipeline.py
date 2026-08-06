from __future__ import annotations

from pathlib import Path
from typing import Any
from typing import Dict
from revise.backend.contracts import EvaluationPolicy
from revise.backend.contracts import InputValidationPolicy
from revise.backend.contracts import LocalRefinementStrategy
from revise.svc import SVC
from revise.utils.provenance import completed_artifact


class UnifiedReconstructionPipeline:
    def __init__(
        self,
        strategy: LocalRefinementStrategy,
        validation_policy: InputValidationPolicy,
        evaluation_policy: EvaluationPolicy,
    ) -> None:
        self.strategy = strategy
        self.validation_policy = validation_policy
        self.evaluation_policy = evaluation_policy

    def run(self, ctx):
        # Template Method: this order is fixed across all modes/routes.
        # Strategy implementations can customize internals per stage but must
        # not change lifecycle topology.
        self._run_stage(ctx, "validate_inputs", self.validate_inputs)
        if ctx.dry_run:
            ctx.skip_pending_stages("dry_run")
            ctx.svc = SVC(
                expr=None,
                spatial=None,
                svc_kind=str(ctx.runtime.get("svc_kind", "sc")),
                provenance={"dry_run": True},
                artifacts={},
            )
            ctx.mark_run_succeeded()
            return ctx.svc

        self._run_stage(ctx, "global_anchoring", self.global_anchoring)
        self._run_stage(ctx, "local_refinement", self.local_refinement)
        self._run_stage(ctx, "finalize", self.finalize_svc)

        evaluation_skip = self._evaluation_skip_reason_or_fail(ctx)
        if evaluation_skip is None:
            self._run_stage(ctx, "evaluate", self.evaluate_if_needed)
        else:
            ctx.skip_stage("evaluate", evaluation_skip)
        ctx.mark_run_succeeded()
        return ctx.svc

    def _run_stage(self, ctx, stage: str, operation) -> None:
        ctx.start_stage(stage)
        ctx.logger.info("[pipeline] stage=%s", stage)
        try:
            operation(ctx)
        except KeyboardInterrupt as exc:
            try:
                ctx.terminate_stage(stage, exc)
            except BaseException as persistence_error:
                raise exc from persistence_error
            raise
        except Exception as exc:
            try:
                ctx.terminate_stage(stage, exc)
            except BaseException as persistence_error:
                raise exc from persistence_error
            raise
        else:
            ctx.succeed_stage(stage)

    def validate_inputs(self, ctx) -> None:
        self.validation_policy.validate(ctx)
        self.strategy.prepare_context(ctx)

    def global_anchoring(self, ctx) -> None:
        self.strategy.global_anchoring(ctx)

    def local_refinement(self, ctx) -> None:
        self.prepare_local_units(ctx)
        self.build_graph(ctx)
        self.build_ot_problem(ctx)
        self.solve_ot(ctx)
        self.update_expression(ctx)

    def prepare_local_units(self, ctx) -> None:
        self.strategy.prepare_local_units(ctx)

    def build_graph(self, ctx) -> None:
        self.strategy.build_graph(ctx)

    def build_ot_problem(self, ctx) -> None:
        self.strategy.build_ot_problem(ctx)

    def solve_ot(self, ctx) -> None:
        self.strategy.solve_ot(ctx)

    def update_expression(self, ctx) -> None:
        self.strategy.update_expression(ctx)

    def finalize_svc(self, ctx) -> None:
        ctx.svc = self.strategy.finalize_svc(ctx)
        self._persist_outputs(ctx)
        if ctx.finalize_callback is not None:
            ctx.finalize_callback(ctx)

    def evaluate_if_needed(self, ctx) -> None:
        from revise.analysis.metrics import compute_metric

        outputs = dict(ctx.svc.artifacts.get("outputs", {})) if ctx.svc else {}
        metrics_dir = Path(ctx.run_dir) / "metrics"
        metrics_dir.mkdir(parents=True, exist_ok=True)

        computed_metrics: Dict[str, Any] = {}
        for key, adata in outputs.items():
            common_index = adata.obs.index.intersection(ctx.real_st_adata.obs.index)
            if common_index.empty:
                ctx.logger.warning("[pipeline] no shared cells for %s; skip metric", key)
                continue

            pred = adata[common_index, :].copy()
            gt = ctx.real_st_adata[common_index, :].copy()
            metrics_df = compute_metric(
                gt,
                pred,
                ctx.logger,
                adata_process=False,
                gene_list=None,
                normalize=True,
            )
            out_file = metrics_dir / f"{key}_metrics_normalized.csv"
            metrics_df.to_csv(out_file)
            ctx.record_artifact(completed_artifact(f"metric:{key}", out_file))
            ctx.quality_metrics[key] = metrics_df
            computed_metrics[key] = metrics_df

        # Backward-compatible sink in compatibility mode. If a route exposes
        # multiple outputs, write the selected primary output explicitly instead
        # of depending on dictionary iteration order.
        if ctx.compatibility_mode and computed_metrics:
            primary_key = None
            if ctx.svc is not None:
                primary_key = ctx.svc.provenance.get("primary_output_key")
            if primary_key not in computed_metrics:
                primary_key = next(iter(computed_metrics.keys()))
            compatibility_path = Path(ctx.run_dir) / "metrics_normalized.csv"
            computed_metrics[primary_key].to_csv(compatibility_path)
            ctx.record_artifact(
                completed_artifact("metric:compatibility", compatibility_path)
            )

        if ctx.svc is not None:
            ctx.svc.quality_metrics = dict(ctx.quality_metrics)

    def _evaluation_skip_reason(self, ctx) -> str | None:
        if not self.evaluation_policy.should_evaluate(ctx):
            ctx.logger.info("[pipeline] evaluation skipped by policy")
            return "policy_disabled"
        outputs = dict(ctx.svc.artifacts.get("outputs", {})) if ctx.svc else {}
        if not outputs:
            ctx.logger.warning("[pipeline] no outputs available for evaluation")
            return "no_outputs"
        if ctx.real_st_adata is None:
            ctx.logger.warning(
                "[pipeline] real_st_adata is missing; benchmark evaluation skipped"
            )
            return "ground_truth_unavailable"
        if not any(
            not adata.obs.index.intersection(ctx.real_st_adata.obs.index).empty
            for adata in outputs.values()
        ):
            ctx.logger.warning("[pipeline] no aligned outputs available for evaluation")
            return "no_aligned_outputs"
        return None

    def _evaluation_skip_reason_or_fail(self, ctx) -> str | None:
        try:
            return self._evaluation_skip_reason(ctx)
        except KeyboardInterrupt as exc:
            self._terminate_pending_evaluation(ctx, exc)
            raise
        except Exception as exc:
            self._terminate_pending_evaluation(ctx, exc)
            raise

    @staticmethod
    def _terminate_pending_evaluation(
        ctx,
        error: BaseException,
    ) -> None:
        try:
            ctx.start_stage("evaluate")
            ctx.terminate_stage("evaluate", error)
        except BaseException as persistence_error:
            raise error from persistence_error

    def _persist_outputs(self, ctx) -> None:
        if ctx.svc is None:
            return
        if not bool(ctx.merged_config.get("io", {}).get("save_outputs", True)):
            return

        outputs = dict(ctx.svc.artifacts.get("outputs", {}))
        if not outputs:
            return

        benchmark_mode = str(ctx.runtime.get("mode")) == "benchmark"
        if not benchmark_mode:
            if ctx.compatibility_mode:
                raise RuntimeError(
                    "Application runs do not support compatibility output aliases"
                )
            self._persist_direct_application_outputs(ctx, outputs)
        elif not ctx.compatibility_mode:
            for key, adata in outputs.items():
                path = Path(ctx.run_dir) / f"{key}.h5ad"
                adata.write_h5ad(path)
                ctx.record_artifact(completed_artifact(f"output:{key}", path))

        if ctx.compatibility_mode:
            self._emit_compatibility_files(ctx, outputs)

    @staticmethod
    def _persist_direct_application_outputs(ctx, outputs: Dict[str, Any]) -> None:
        """Keep engine-level Application artifacts in the one run directory."""
        invalid_keys = [
            key
            for key in outputs
            if not isinstance(key, str)
            or not key
            or key in {".", ".."}
            or "/" in key
            or "\\" in key
        ]
        if invalid_keys:
            raise RuntimeError(
                f"Application output keys must be safe filename stems: {invalid_keys}"
            )

        run_dir = Path(ctx.run_dir).resolve()
        targets = {
            key: (run_dir / f"{key}.h5ad").resolve() for key in outputs
        }
        outside = [str(path) for path in targets.values() if path.parent != run_dir]
        if outside:
            raise RuntimeError(
                "Application engine outputs must be direct children of the run "
                f"directory: run_dir={run_dir}; invalid={outside}"
            )

        created: list[Path] = []
        artifact_identities = {
            (f"output:{key}", str(path)) for key, path in targets.items()
        }

        def cleanup() -> None:
            cleanup_errors = []
            for path in created:
                try:
                    path.unlink(missing_ok=True)
                except OSError as exc:
                    cleanup_errors.append(
                        {
                            "path": str(path),
                            "type": type(exc).__name__,
                            "message": str(exc),
                        }
                    )
            ctx.artifact_records[:] = [
                record
                for record in ctx.artifact_records
                if (record.get("role"), record.get("path"))
                not in artifact_identities
            ]
            if cleanup_errors:
                ctx.provenance.setdefault("output_cleanup_errors", []).extend(
                    cleanup_errors
                )
                ctx.logger.error(
                    "[pipeline] failed to clean %d direct output(s)",
                    len(cleanup_errors),
                )

        ctx.register_output_failure_cleanup(cleanup)
        try:
            for key, adata in outputs.items():
                path = targets[key]
                if path.exists():
                    raise FileExistsError(
                        f"Refusing to overwrite an existing run output: {path}"
                    )
                created.append(path)
                adata.write_h5ad(path)
                ctx.record_artifact(completed_artifact(f"output:{key}", path))
        except BaseException:
            ctx.cleanup_failed_outputs()
            raise

    def _emit_compatibility_files(self, ctx, outputs: Dict[str, Any]) -> None:
        compatibility_map = {
            "sp_svc": "sp_SVC.h5ad",
            "sc_svc_dec": "sc_SVC.h5ad",
            "sc_svc_dec_graphagg": "sc_SVC_graphagg.h5ad",
            "sc_svc_expr": "sc_SVC_expr.h5ad",
            "sc_svc_spatial": "sc_SVC_spatial.h5ad",
            "sc_svc_impute_in_panel": "sc_SVC_impute_in_panel.h5ad",
            "sc_svc_impute_all_panel": "sc_SVC_impute_all_panel.h5ad",
        }
        for key, filename in compatibility_map.items():
            if key in outputs:
                path = Path(ctx.run_dir) / filename
                outputs[key].write_h5ad(path)
                ctx.record_artifact(
                    completed_artifact(f"compatibility:{key}", path)
                )
