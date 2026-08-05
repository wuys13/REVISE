"""Application route contract and reconstruction use case."""

from __future__ import annotations

from pathlib import Path
from typing import Any, NamedTuple, TYPE_CHECKING

from revise.framework import REVISEPipeline

if TYPE_CHECKING:
    from revise.application.request import ApplicationRequest


class ApplicationExecution(NamedTuple):
    """Normalized result returned by the application execution boundary."""

    status: str
    svc_type: str
    pipeline: dict[str, Any]
    result: Any = None
    output_path: Any = None
    preflight: Path | None = None


def _build_algorithm_overrides(request: ApplicationRequest) -> dict:
    overrides = {}
    if request.ot_method is not None:
        overrides["ot"] = {
            "ga": {"solver": request.ot_method},
            "lr": {"solver": request.ot_method},
        }
    if request.local_refinement_strength is not None:
        overrides["local_refinement"] = {
            "strength": request.local_refinement_strength,
        }
    columns = {"cell_type_col": request.broad_column}
    if request.subtype_column is not None:
        columns["sub_cell_type_col"] = request.subtype_column
    overrides["columns"] = columns
    if request.svc_type == "sc-SVC":
        overrides["sc"] = {"select_ct": request.select_cell_type}
    return overrides


def _build_io_overrides(request: ApplicationRequest) -> dict:
    overrides = {
        "output_root": str(request.output_root),
        "sample_name": request.sample_name,
        "patient_key": request.patient_key,
        "save_outputs": False,
        "input_format": request.st_format,
    }
    if request.input_mode == "direct":
        overrides.update(
            data_root="",
            st_file="",
            sc_ref_file="",
            st_path=str(request.st_path),
            sc_ref_path=str(request.reference_path),
        )
    else:
        overrides.update(
            data_root=str(request.data_root),
            st_path="",
            sc_ref_path="",
            st_file=request.st_file,
            sc_ref_file=request.reference_file,
        )
    if request.st_format in {"spatialdata", "auto"}:
        overrides["spatialdata_path"] = str(request.st_path)
        if request.spatialdata_table is not None:
            overrides["spatialdata_table"] = request.spatialdata_table
        if request.spatialdata_element is not None:
            overrides["spatialdata_spatial_element"] = request.spatialdata_element
    return overrides


def _run_pipeline(
    request: ApplicationRequest,
    *,
    finalize_callback=None,
):
    pipeline = REVISEPipeline(config_path=None)
    runtime_overrides = {}
    if request.seed is not None:
        runtime_overrides["seed"] = request.seed

    return pipeline._execute_run(
        svc_type=request.svc_type,
        cf=None,
        runtime_overrides=runtime_overrides,
        io_overrides=_build_io_overrides(request),
        algorithm_overrides=_build_algorithm_overrides(request),
        dry_run=request.effective_action == "preflight",
        finalize_callback=finalize_callback,
        application_config_metadata={
            "source_path": getattr(
                request,
                "source_path",
                str(getattr(request, "config_path", "")),
            ),
            "source_sha256": request.config_sha256,
            "declared_root": request.declared_root,
            "resolved_root": str(request.resolved_root),
            "cwd": str(request.cwd),
            "resolved_paths": request.resolved_paths,
            "declared_action": request.action,
            "effective_action": request.effective_action,
            "dry_run_override": request.dry_run_override,
        },
    )


def _pipeline_evidence(svc) -> dict[str, Any]:
    pipeline = svc.summary()
    route = svc.provenance.get("route", {})
    pipeline.update(
        profile=svc.provenance.get("profile"),
        task=route.get("task"),
        strategy=route.get("strategy"),
        route=route,
    )
    return pipeline


def execute_application(request: ApplicationRequest) -> ApplicationExecution:
    """Run or preflight one validated application request."""
    if request.effective_action == "preflight":
        svc = _run_pipeline(request)
        return ApplicationExecution(
            status="ready",
            svc_type=request.svc_type,
            preflight=Path(svc.provenance["run_dir"]) / "preflight.json",
            pipeline=_pipeline_evidence(svc),
        )

    published = {}

    def publish(ctx):
        from .publication import publish_pair, publish_single

        if request.svc_type == "sc-SVC":
            published["result"], published["path"] = publish_pair(request, ctx)
        else:
            published["result"], published["path"] = publish_single(request, ctx)

    svc = _run_pipeline(request, finalize_callback=publish)
    return ApplicationExecution(
        status="succeeded",
        svc_type=request.svc_type,
        result=published["result"],
        output_path=published["path"],
        pipeline=_pipeline_evidence(svc),
    )
