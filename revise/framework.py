from __future__ import annotations

import copy
import signal
import threading
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, Iterable, Optional

from revise.backend import build_default_plugin_registry
from revise.backend import ModeEvaluationPolicy
from revise.backend import ModeValidationPolicy
from revise.backend import build_default_registry
from revise.config import infer_default_profile
from revise.config import load_raw_config
from revise.config import merge_unified_config
from revise.recon.context import PipelineContext
from revise.recon.pipeline import UnifiedReconstructionPipeline
from revise.svc import SVC
from revise.utils import (
    build_task_dir,
    build_run_dir,
    build_run_logger,
    canonical_config_projection,
    collect_package_versions,
    exclusive_run_directory,
    hash_jsonable,
    set_global_seed,
    write_json,
)


class _SigtermInterrupt(KeyboardInterrupt):
    def __init__(self, signum, previous_handler, frame):
        super().__init__("received SIGTERM")
        self.signum = signum
        self.previous_handler = previous_handler
        self.frame = frame


@contextmanager
def _temporary_sigterm_handler():
    """Translate SIGTERM into a handled interruption for one main-thread run."""
    if threading.current_thread() is not threading.main_thread():
        yield
        return

    previous_handler = signal.getsignal(signal.SIGTERM)
    if previous_handler is signal.SIG_IGN:
        yield
        return

    installed = True

    def restore_handler():
        nonlocal installed
        if installed:
            signal.signal(signal.SIGTERM, previous_handler)
            installed = False

    def handle_sigterm(signum, frame):
        restore_handler()
        raise _SigtermInterrupt(signum, previous_handler, frame)

    signal.signal(signal.SIGTERM, handle_sigterm)
    try:
        yield
    finally:
        restore_handler()


class REVISEPipeline:
    """Unified orchestration API for all REVISE tasks and modes."""

    def __init__(self, config_path: Optional[str] = None):
        if config_path is None:
            config_path = str(Path(__file__).with_name("revise.yaml"))
        self.config_path = str(self._resolve_config_path(config_path))
        self.raw_config = load_raw_config(self.config_path)
        self.registry = None
        self.plugin_registry = None

    @staticmethod
    def _resolve_config_path(config_path: str | Path) -> Path:
        path = Path(config_path)
        if path.exists():
            return path
        # Backward-compatible default used by README examples and root wrapper
        # scripts. In installed PyPI wheels, revise.yaml lives beside this file,
        # not under the caller's current working directory.
        if path.as_posix() == "revise/revise.yaml":
            packaged = Path(__file__).with_name("revise.yaml")
            if packaged.exists():
                return packaged
        return path

    def run(
        self,
        *,
        profile: Optional[str] = None,
        runtime_overrides: Optional[Dict[str, Any]] = None,
        io_overrides: Optional[Dict[str, Any]] = None,
        set_overrides: Optional[Iterable[str]] = None,
        dry_run: bool = False,
        finalize_callback=None,
    ):
        # 1) Resolve final runtime config from single YAML entry:
        # defaults -> profile -> CLI runtime/io overrides -> --set overrides.
        runtime_overrides = dict(runtime_overrides or {})
        io_overrides = dict(io_overrides or {})
        set_overrides = list(set_overrides or [])

        if profile is None:
            profile = infer_default_profile(self.raw_config, runtime_overrides)

        merged_config = merge_unified_config(
            raw_config=self.raw_config,
            profile=profile,
            runtime_overrides=runtime_overrides,
            io_overrides=io_overrides,
            set_overrides=set_overrides,
        )

        runtime = self._resolve_runtime_plugins(merged_config)
        config_hash = hash_jsonable(canonical_config_projection(merged_config))
        route_key = f"{runtime['platform']}:{runtime['confounding']}"
        output_root = merged_config["io"]["output_root"]
        sample_name = merged_config["io"]["sample_name"]
        run_dir = build_run_dir(
            output_root=output_root,
            sample_name=sample_name,
            route_key=route_key,
            io_cfg=merged_config["io"],
        )
        if route_key.startswith("sim2real:"):
            log_dir = build_task_dir(
                output_root=output_root,
                sample_name=sample_name,
                route_key=route_key,
                io_cfg=merged_config["io"],
            )
        else:
            log_dir = run_dir

        with exclusive_run_directory(run_dir):
            return self._run_in_directory(
                merged_config=merged_config,
                runtime=runtime,
                route_key=route_key,
                log_dir=log_dir,
                run_dir=run_dir,
                sample_name=sample_name,
                profile=profile,
                config_hash=config_hash,
                dry_run=dry_run,
                finalize_callback=finalize_callback,
            )

    def _run_in_directory(
        self,
        *,
        merged_config: Dict[str, Any],
        runtime: Dict[str, Any],
        route_key: str,
        log_dir: Path,
        run_dir: Path,
        sample_name: str,
        profile: Optional[str],
        config_hash: str,
        dry_run: bool,
        finalize_callback,
    ):
        logger_name = f"REVISEUnified::{sample_name}::{route_key}"
        if log_dir == run_dir:
            logger_name = f"{logger_name}::{Path(run_dir).name}"
        logger = build_run_logger(
            run_name=logger_name,
            run_dir=log_dir,
        )
        logger.info("[framework] start unified run route=%s strategy=%s", route_key, runtime["strategy"])

        set_global_seed(seed=runtime.get("seed"), deterministic=bool(runtime.get("deterministic", True)))

        ctx = PipelineContext(
            merged_config=merged_config,
            raw_config=self.raw_config,
            config_path=self.config_path,
            profile=profile,
            runtime=runtime,
            route_key=route_key,
            run_dir=run_dir,
            logger=logger,
            config_hash=config_hash,
            dry_run=bool(dry_run),
            finalize_callback=finalize_callback,
        )
        ctx.set_provenance_callback(self._write_final_metadata, notify=False)

        with _temporary_sigterm_handler():
            try:
                # Establish the run envelope before any validation or strategy
                # work so abrupt termination leaves parseable running truth.
                self._write_final_metadata(ctx)
                self._write_initial_metadata(ctx)

                if ctx.dry_run:
                    # Dry-run validates structural inputs without importing the
                    # heavy strategy registry. U8 extends this shared preflight.
                    self._run_dry_validation(ctx)
                    ctx.svc = SVC(
                        expr=None,
                        spatial=None,
                        svc_kind=str(runtime.get("svc_kind", "sc")),
                        provenance={"dry_run": True, "route": ctx.route},
                        artifacts={},
                    )
                    ctx.mark_run_succeeded()
                    logger.info("[framework] dry-run validated route=%s", route_key)
                    return ctx.svc

                if self.registry is None:
                    self.registry = build_default_registry()
                strategy = self.registry.get(runtime["strategy"])
                pipeline = UnifiedReconstructionPipeline(
                    strategy=strategy,
                    validation_policy=ModeValidationPolicy(),
                    evaluation_policy=ModeEvaluationPolicy(),
                )

                svc = pipeline.run(ctx)
                ctx.commit_pending_publication()
                logger.info("[framework] finished unified run route=%s", route_key)
                return svc
            except _SigtermInterrupt as exc:
                try:
                    ctx.rollback_pending_publication()
                    ctx.terminate_run(exc, interrupted=True)
                except BaseException as persistence_error:
                    raise KeyboardInterrupt("received SIGTERM") from persistence_error
                logger.warning("[framework] run interrupted by SIGTERM")
                previous_handler = exc.previous_handler
                if callable(previous_handler):
                    previous_handler(exc.signum, exc.frame)
                raise KeyboardInterrupt("received SIGTERM") from None
            except KeyboardInterrupt as exc:
                try:
                    ctx.rollback_pending_publication()
                    ctx.terminate_run(exc, interrupted=True)
                except BaseException as persistence_error:
                    raise exc from persistence_error
                logger.warning("[framework] run interrupted")
                raise
            except Exception as exc:
                try:
                    ctx.rollback_pending_publication()
                    ctx.terminate_run(exc)
                except BaseException as persistence_error:
                    raise exc from persistence_error
                logger.exception("[framework] run failed")
                raise

    def _run_dry_validation(self, ctx: PipelineContext) -> None:
        ctx.start_stage("validate_inputs")
        try:
            ModeValidationPolicy().validate(ctx)
        except KeyboardInterrupt as exc:
            try:
                ctx.terminate_stage("validate_inputs", exc, interrupted=True)
            except BaseException as persistence_error:
                raise exc from persistence_error
            raise
        except Exception as exc:
            try:
                ctx.terminate_stage("validate_inputs", exc)
            except BaseException as persistence_error:
                raise exc from persistence_error
            raise
        else:
            ctx.succeed_stage("validate_inputs")
        ctx.skip_pending_stages("dry_run")

    def _write_initial_metadata(self, ctx: PipelineContext) -> None:
        write_json(Path(ctx.run_dir) / "merged_config.json", self._export_merged_config(ctx))

    def _write_final_metadata(self, ctx: PipelineContext) -> None:
        provenance = {
            "schema_version": 2,
            "run": copy.deepcopy(
                getattr(
                    ctx,
                    "run_record",
                    {
                        "status": getattr(ctx, "run_status", "running"),
                        "dry_run": bool(getattr(ctx, "dry_run", False)),
                        "started_at": getattr(ctx, "run_started_at", None),
                        "ended_at": getattr(ctx, "run_ended_at", None),
                        "duration_seconds": getattr(
                            ctx, "run_duration_seconds", None
                        ),
                        "error": getattr(ctx, "run_error", None),
                    },
                )
            ),
            "config_path": ctx.config_path,
            "profile": ctx.profile,
            "route": ctx.route,
            "route_key": ctx.route_key,
            "run_dir": str(ctx.run_dir),
            "config_hash": getattr(ctx, "config_hash", None),
            "data_fingerprint": getattr(ctx, "data_fingerprint", None),
            "data_fingerprint_error": copy.deepcopy(
                getattr(ctx, "data_fingerprint_error", None)
            ),
            "packages": collect_package_versions(
                [
                    "revise-svc",
                    "scanpy",
                    "anndata",
                    "numpy",
                    "pandas",
                    "scipy",
                    "POT",
                    "tacco",
                    "leidenalg",
                ]
            ),
            "stages": copy.deepcopy(getattr(ctx, "stage_records", [])),
            "artifacts": copy.deepcopy(getattr(ctx, "artifact_records", [])),
            "quality_metric_keys": sorted(ctx.quality_metrics.keys()),
            "svc_summary": ctx.svc.summary() if ctx.svc else {},
            "ot_config": copy.deepcopy(ctx.merged_config["ot"]),
            "ot_events": copy.deepcopy(ctx.ot_events),
        }
        result = copy.deepcopy(getattr(ctx, "provenance", {}).get("result"))
        if result is not None:
            provenance["result"] = result

        write_json(Path(ctx.run_dir) / "provenance.json", provenance)

        if ctx.svc is not None:
            ctx.svc.provenance.update(provenance)

        if hasattr(ctx, "provenance"):
            ctx.provenance = copy.deepcopy(provenance)

    def _export_merged_config(self, ctx: PipelineContext) -> Dict[str, Any]:
        exported = copy.deepcopy(ctx.merged_config)
        return exported

    def _resolve_runtime_plugins(self, merged_config: Dict[str, Any]) -> Dict[str, Any]:
        """Resolve platform and confounding plugins before strategy instantiation."""
        if self.plugin_registry is None:
            self.plugin_registry = build_default_plugin_registry()

        runtime = merged_config["runtime"]
        payload: Dict[str, Any] = {
            "runtime": runtime,
            "merged_config": merged_config,
        }

        platform_adapter_id = runtime.get("platform_adapter") or runtime.get("platform") or "default"
        cf_strategy_id = runtime.get("cf_strategy") or runtime.get("confounding") or "default"
        payload = self.plugin_registry.get_platform_adapter(platform_adapter_id).adapt(payload)
        payload = self.plugin_registry.get_cf_strategy(cf_strategy_id).apply(payload)

        resolved_runtime = payload.get("runtime", runtime)
        merged_config["runtime"] = resolved_runtime
        return resolved_runtime
