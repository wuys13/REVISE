from __future__ import annotations

import copy
import signal
import threading
import warnings
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, Optional

from revise.backend import ModeEvaluationPolicy
from revise.backend import ModeValidationPolicy
from revise.backend import build_default_registry
from revise.config import ConfigError
from revise.config import load_raw_config
from revise.config import merge_unified_config
from revise.config import resolve_semantic_route
from revise.recon.context import PipelineContext
from revise.recon.pipeline import UnifiedReconstructionPipeline
from revise.svc import SVC
from revise.utils import (
    build_task_dir,
    build_run_dir,
    build_run_logger,
    canonical_config_projection,
    collect_software_versions,
    exclusive_run_directory,
    hash_jsonable,
    set_global_seed,
    sha256_file,
    write_json,
)
from revise.utils.logging import log_exception_to_run_file


_APPLICATION_CONFIG_PROVENANCE_KEYS = (
    "source_path",
    "source_sha256",
    "declared_root",
    "resolved_root",
    "cwd",
    "resolved_paths",
    "declared_action",
    "effective_action",
    "dry_run_override",
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


def _manifest_identity(path: Path) -> dict[str, int | str] | None:
    try:
        stat_result = path.stat()
        digest = sha256_file(path)
    except OSError:
        return None
    return {
        "mtime_ns": stat_result.st_mtime_ns,
        "size": stat_result.st_size,
        "sha256": digest,
    }


class REVISEPipeline:
    """Unified orchestration API for all REVISE tasks and modes."""

    def __init__(self, config_path: Optional[str] = None):
        if config_path is None:
            config_path = str(Path(__file__).with_name("revise.yaml"))
        self.config_path = str(self._resolve_config_path(config_path))
        self.raw_config = load_raw_config(self.config_path)
        self.registry = None

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
        svc_type: Optional[str] = None,
        cf: Optional[str] = None,
        runtime_overrides: Optional[Dict[str, Any]] = None,
        io_overrides: Optional[Dict[str, Any]] = None,
        dry_run: bool = False,
        finalize_callback=None,
    ):
        return self._execute_run(
            svc_type=svc_type,
            cf=cf,
            runtime_overrides=runtime_overrides,
            io_overrides=io_overrides,
            algorithm_overrides=None,
            dry_run=dry_run,
            finalize_callback=finalize_callback,
        )

    def _execute_run(
        self,
        *,
        svc_type: Optional[str] = None,
        cf: Optional[str] = None,
        runtime_overrides: Optional[Dict[str, Any]] = None,
        io_overrides: Optional[Dict[str, Any]] = None,
        algorithm_overrides: Optional[Dict[str, Any]] = None,
        dry_run: bool = False,
        finalize_callback=None,
        application_config_metadata: Optional[Dict[str, Any]] = None,
    ):
        # 1) Resolve final runtime config from single YAML entry:
        # defaults -> profile -> runtime/io overrides -> algorithm overrides.
        runtime_overrides = dict(runtime_overrides or {})
        io_overrides = dict(io_overrides or {})
        algorithm_overrides = dict(algorithm_overrides or {})

        route_identity_keys = {
            "platform",
            "application_route",
            "confounding",
            "mode",
            "task",
            "svc_kind",
            "strategy",
        }
        forbidden = sorted(route_identity_keys & set(runtime_overrides))
        if forbidden:
            raise ConfigError(
                "runtime_overrides cannot modify route identity: "
                + ", ".join(forbidden)
            )

        resolved_route = resolve_semantic_route(
            self.raw_config,
            svc_type=svc_type,
            cf=cf,
        )
        route_warning = resolved_route.pop("warning")
        profile = resolved_route.pop("profile")
        if route_warning:
            warnings.warn(route_warning, UserWarning, stacklevel=2)
        resolved_runtime = {**resolved_route, **runtime_overrides}

        merged_config = merge_unified_config(
            raw_config=self.raw_config,
            profile=profile,
            runtime_overrides=resolved_runtime,
            io_overrides=io_overrides,
            algorithm_overrides=algorithm_overrides,
        )

        runtime = merged_config["runtime"]
        config_hash = hash_jsonable(canonical_config_projection(merged_config))
        selector = (
            runtime["application_route"]
            if runtime["mode"] == "application"
            else runtime["confounding"]
        )
        route_key = f"{runtime['mode']}:{selector}"
        output_root = merged_config["io"]["output_root"]
        sample_name = merged_config["io"]["sample_name"]
        run_dir = build_run_dir(
            output_root=output_root,
            sample_name=sample_name,
            route_key=route_key,
            io_cfg=merged_config["io"],
            mode=runtime["mode"],
            cf=runtime.get("confounding"),
        )
        if runtime["mode"] == "benchmark":
            log_dir = build_task_dir(
                output_root=output_root,
                sample_name=sample_name,
                route_key=route_key,
                io_cfg=merged_config["io"],
                mode=runtime["mode"],
                cf=runtime.get("confounding"),
            )
        else:
            log_dir = run_dir

        with exclusive_run_directory(run_dir):
            manifest_path = run_dir / "provenance.json"
            manifest_before = _manifest_identity(manifest_path)
            try:
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
                    application_config_metadata=application_config_metadata,
                    route_warning=route_warning,
                )
            except BaseException as exc:
                manifest_after = _manifest_identity(manifest_path)
                if (
                    manifest_after is not None
                    and manifest_after != manifest_before
                ):
                    exc._revise_failure_context = {
                        "run_dir": str(run_dir),
                        "manifest_path": str(manifest_path),
                        "manifest_identity": manifest_after,
                    }
                raise

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
        application_config_metadata: Optional[Dict[str, Any]],
        route_warning: Optional[str],
    ):
        logger_name = f"REVISEUnified::{sample_name}::{route_key}"
        if log_dir == run_dir:
            logger_name = f"{logger_name}::{Path(run_dir).name}"
        logger = build_run_logger(
            run_name=logger_name,
            run_dir=log_dir,
        )
        if route_warning:
            logger.warning("[framework] %s", route_warning)
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
            application_config_metadata=copy.deepcopy(
                application_config_metadata or {}
            ),
            software_versions=collect_software_versions(merged_config),
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
                    ctx.terminate_run(exc)
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
                    ctx.terminate_run(exc)
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
                log_exception_to_run_file(
                    logger,
                    f"[framework] run failed: {type(exc).__name__}: {exc}",
                )
                raise

    def _run_dry_validation(self, ctx: PipelineContext) -> None:
        ctx.start_stage("validate_inputs")
        try:
            ModeValidationPolicy().validate(ctx)
        except KeyboardInterrupt as exc:
            try:
                ctx.terminate_stage("validate_inputs", exc)
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
            "runtime_seed": ctx.merged_config.get("runtime", {}).get("seed"),
            "config_hash": getattr(ctx, "config_hash", None),
            "input_identities": copy.deepcopy(
                getattr(ctx, "input_identities", [])
            ),
            "packages": copy.deepcopy(ctx.software_versions),
            "stages": copy.deepcopy(getattr(ctx, "stage_records", [])),
            "artifacts": copy.deepcopy(getattr(ctx, "artifact_records", [])),
            "quality_metric_keys": sorted(ctx.quality_metrics.keys()),
            "svc_summary": ctx.svc.summary() if ctx.svc else {},
            "ot_config": copy.deepcopy(ctx.merged_config["ot"]),
            "local_refinement": copy.deepcopy(
                ctx.local_refinement_record
            ),
            "sr_allocation": copy.deepcopy(
                getattr(ctx, "sr_allocation_records", [])
            ),
        }
        current_provenance = getattr(ctx, "provenance", {})
        application_config = copy.deepcopy(
            getattr(ctx, "application_config_metadata", {})
        )
        if application_config:
            provenance["application_config"] = {
                key: application_config[key]
                for key in _APPLICATION_CONFIG_PROVENANCE_KEYS
                if key in application_config
            }
        for result_key in ("result", "results"):
            result_value = copy.deepcopy(current_provenance.get(result_key))
            if result_value is not None:
                provenance[result_key] = result_value

        write_json(Path(ctx.run_dir) / "provenance.json", provenance)

        if ctx.svc is not None:
            ctx.svc.provenance.update(provenance)

        if hasattr(ctx, "provenance"):
            ctx.provenance = copy.deepcopy(provenance)

    def _export_merged_config(self, ctx: PipelineContext) -> Dict[str, Any]:
        exported = copy.deepcopy(ctx.merged_config)
        return exported
