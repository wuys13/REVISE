from __future__ import annotations

import tempfile
from pathlib import Path

from revise.backend.contracts import EvaluationPolicy
from revise.backend.contracts import InputValidationPolicy
from revise.backend.ops.tacco_runtime import require_tacco
from revise.config.runner_conf import REQUIRED_IO_BY_MODE_TASK
from revise.config.runner_conf import resolve_input_specs
from revise.io import REVISEInputService
from revise.utils import completed_artifact, fingerprint_paths, write_json


class ModeValidationPolicy(InputValidationPolicy):
    def validate(self, ctx) -> None:
        runtime = ctx.runtime
        mode = runtime.get("mode")
        if mode not in {"application", "benchmark"}:
            raise ValueError(f"Unsupported mode: {mode}")
        task = runtime.get("task")
        if task in (None, ""):
            raise ValueError("runtime.task is required")

        io_cfg = ctx.io
        data_root = io_cfg.get("data_root")
        if not data_root:
            raise ValueError("io.data_root is required")
        if not Path(data_root).exists():
            raise FileNotFoundError(f"io.data_root does not exist: {data_root}")

        sample_name = io_cfg.get("sample_name")
        if not sample_name:
            raise ValueError("io.sample_name is required")
        output_root = io_cfg.get("output_root")
        if not output_root:
            raise ValueError("io.output_root is required")

        key = (mode, task)
        if key not in REQUIRED_IO_BY_MODE_TASK:
            supported = sorted(f"{m}:{t}" for (m, t) in REQUIRED_IO_BY_MODE_TASK.keys())
            raise ValueError(
                f"Unsupported mode/task combination: {mode}:{task}. "
                f"Supported combinations: {supported}"
            )

        required = REQUIRED_IO_BY_MODE_TASK[key]
        missing = sorted(k for k in required if io_cfg.get(k) in (None, ""))
        if missing:
            raise ValueError(
                f"Missing required io keys for {mode}:{task}: {missing}"
            )

        self._validate_solver_compatibility(ctx)
        if any(
            str(ctx.merged_config["ot"][phase]["solver"]) == "tacco"
            for phase in ("ga", "lr")
        ):
            require_tacco()

        specs = resolve_input_specs(runtime, io_cfg)
        ctx.input_specs = specs
        try:
            ctx.data_fingerprint = fingerprint_paths(specs)
        except (OSError, RuntimeError, TypeError, ValueError) as exc:
            ctx.data_fingerprint = None
            ctx.data_fingerprint_error = {
                "type": type(exc).__name__,
                "message": str(exc),
            }
            raise ValueError(
                "Invalid input content fingerprint: "
                f"{type(exc).__name__}: {exc}"
            ) from exc
        report = REVISEInputService.from_context(ctx).preflight(
            specs,
            runtime=runtime,
            columns=ctx.columns,
        )
        run_dir = Path(ctx.run_dir)
        run_dir.mkdir(parents=True, exist_ok=True)
        publication_dir = Path(output_root) / str(sample_name)
        publication_dir.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            dir=publication_dir,
            prefix=".preflight-write-",
        ):
            pass
        report["output"] = {
            "run_dir": str(run_dir),
            "publication_dir": str(publication_dir),
            "writable": True,
        }
        report_path = run_dir / "preflight.json"
        write_json(report_path, report)
        ctx.record_artifact(completed_artifact("preflight", report_path))

    @staticmethod
    def _validate_solver_compatibility(ctx) -> None:
        refinement = ctx.merged_config.get("local_refinement", {}) or {}
        guidance = str(refinement.get("guidance", "off"))
        if guidance == "off":
            return
        compatibility = refinement.get("compatibility", {}) or {}
        uses_reference = compatibility.get("mode") == "reference"
        runtime = ctx.runtime
        mode = str(runtime.get("mode"))
        task = str(runtime.get("task"))
        lr_solver = str(ctx.merged_config["ot"]["lr"]["solver"])

        if task == "sc_svc" and uses_reference:
            raise ValueError(
                "graph_edge local refinement does not support reference "
                "compatibility; use compatibility mode 'cost'"
            )
        if mode == "application" and uses_reference:
            raise ValueError(
                "application local refinement does not support reference "
                "compatibility; use compatibility mode 'cost'"
            )
        if uses_reference and lr_solver == "tacco":
            raise ValueError(
                "TACCO local refinement is incompatible with reference "
                "compatibility; use compatibility mode 'cost' or set "
                "ot.lr.solver=pot for a supported benchmark ablation"
            )
        if (
            mode == "benchmark"
            and task == "sc_svc_sr"
            and not bool(
                ctx.merged_config.get("sc", {}).get(
                    "sr_graph_agg_enabled",
                    False,
                )
            )
            and guidance == "require"
        ):
            callback = getattr(
                ctx,
                "assignment_guidance_callback",
                None,
            )
            if callback is not None:
                problem_key = "sr-benchmark:graph-branch"
                callback(
                    "start",
                    problem_key=problem_key,
                    route=str(ctx.route_key),
                    operator="virtual_cell_ot",
                    phase="lr",
                    mode="require",
                    applicability="applicable",
                    numerics={
                        "beta": compatibility.get("beta"),
                        "min_affinity": compatibility.get("min_affinity"),
                        "operator_strength": compatibility.get("strength"),
                    },
                    solver=lr_solver,
                )
                callback(
                    "terminal",
                    problem_key=problem_key,
                    outcome="failed",
                    availability="unavailable",
                )
            raise ValueError(
                "required local-refinement guidance is incompatible with the "
                "disabled SR graph branch"
            )


class ModeEvaluationPolicy(EvaluationPolicy):
    def should_evaluate(self, ctx) -> bool:
        runtime = ctx.runtime
        if runtime.get("mode") == "benchmark":
            return True
        return bool(ctx.merged_config.get("benchmark", {}).get("evaluate", False))
