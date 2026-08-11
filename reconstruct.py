#!/usr/bin/env python3
"""Readable high-level REVISE Application entrypoint."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from anndata import AnnData

from revise.application.config import (
    ApplicationConfig,
    ApplicationConfigError,
    apply_application_overrides,
    compile_application_config,
    load_application_yaml,
)
from revise.application.preprocess import (
    filter_reference,
    preprocess_reference,
    preprocess_spatial,
)
from revise.application.publication import (
    application_metadata,
    output_paths,
    publish_outputs,
)
from revise.framework import REVISEPipeline


def _engine_overrides(
    config: ApplicationConfig,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    runtime = {"seed": config.seed} if config.seed is not None else {}
    io = {
        "st_path": str(config.st_path),
        "sc_ref_path": str(config.reference_path),
        "pm_on_cell_path": str(config.pm_on_cell_path) if config.pm_on_cell_path else "",
        "output_root": str(config.output_dir),
        "sample_name": config.svc_type,
        "patient_key": "",
        "save_outputs": False,
        "input_format": config.st_format,
        "data_root": "",
        "st_file": "",
        "sc_ref_file": "",
    }
    if config.st_format in {"spatialdata", "auto"}:
        io["spatialdata_path"] = str(config.st_path)
        if config.spatialdata_table is not None:
            io["spatialdata_table"] = config.spatialdata_table
        if config.spatialdata_element is not None:
            io["spatialdata_spatial_element"] = config.spatialdata_element

    algorithm: dict[str, Any] = {"columns": {"cell_type_col": config.broad_column}}
    if config.subtype_column is not None:
        algorithm["columns"]["sub_cell_type_col"] = config.subtype_column
    if config.select_cell_type is not None:
        algorithm["sc"] = {"select_ct": config.select_cell_type}
    if config.local_refinement_strength is not None:
        algorithm["local_refinement"] = {"strength": config.local_refinement_strength}
    if config.local_refinement_alpha is not None:
        algorithm["graph"] = {"alpha": config.local_refinement_alpha}
    if config.local_refinement_resolutions is not None:
        algorithm.setdefault("sc", {})["resolutions"] = list(
            config.local_refinement_resolutions
        )
    if config.local_refinement_graph_method is not None:
        algorithm["graph"] = {
            "method": config.local_refinement_graph_method,
            "alpha": config.local_refinement_graph_alpha,
            "n_neighbors": config.local_refinement_graph_n_neighbors,
            "exp_neighbors": config.local_refinement_graph_exp_neighbors,
            "spatial_neighbors": config.local_refinement_graph_spatial_neighbors,
        }
    if config.local_refinement_match_spot_sum is not None:
        algorithm.setdefault("sc", {})["match_spot_sum"] = (
            config.local_refinement_match_spot_sum
        )
    if config.ot_method is not None:
        algorithm["ot"] = {
            "ga": {"solver": config.ot_method},
            "lr": {"solver": config.ot_method},
        }
    return runtime, io, algorithm


def _pipeline_evidence(svc) -> dict[str, Any]:
    evidence = svc.summary()
    route = svc.provenance.get("route", {})
    evidence.update(
        profile=svc.provenance.get("profile"),
        task=route.get("task"),
        strategy=route.get("strategy"),
        route=route,
    )
    return evidence


def _result_shapes(svc_type: str, result) -> dict[str, list[int]]:
    if svc_type == "sc-SVC":
        return {
            "spatial": list(result[0].shape),
            "expression": list(result[1].shape),
        }
    return {"svc": list(result.shape)}


def _execute_application(
    config: str | Path,
    *,
    svc_type: str | None = None,
    root_dir: str | None = None,
    st_path: str | None = None,
    st_format: str | None = None,
    sc_ref_path: str | None = None,
    spatialdata_table: str | None = None,
    spatialdata_element: str | None = None,
    pm_on_cell_path: str | None = None,
    ot_method: str | None = None,
    cell_type_col: str | None = None,
    sub_cell_type_col: str | None = None,
    select_cell_type: str | None = None,
    local_refinement_strength: float | None = None,
    spatial_min_transcript_counts: int | None = None,
    spatial_min_counts: int | None = None,
    spatial_min_cell_counts: int | None = None,
    reference_min_transcript_counts: int | None = None,
    reference_min_genes: int | None = None,
    reference_min_cell_counts: int | None = None,
    output_dir: str | None = None,
    output_name: str | None = None,
    seed: int | None = None,
    dry_run: bool = False,
) -> tuple[AnnData | tuple[AnnData, AnnData] | None, dict[str, Any]]:
    source, document = load_application_yaml(config)
    overrides = apply_application_overrides(
        document,
        {
            "svc_type": svc_type,
            "root_dir": root_dir,
            "st_path": st_path,
            "st_format": st_format,
            "sc_ref_path": sc_ref_path,
            "spatialdata_table": spatialdata_table,
            "spatialdata_element": spatialdata_element,
            "pm_on_cell_path": pm_on_cell_path,
            "ot_method": ot_method,
            "cell_type_col": cell_type_col,
            "sub_cell_type_col": sub_cell_type_col,
            "select_cell_type": select_cell_type,
            "local_refinement_strength": local_refinement_strength,
            "spatial_min_transcript_counts": spatial_min_transcript_counts,
            "spatial_min_counts": spatial_min_counts,
            "spatial_min_cell_counts": spatial_min_cell_counts,
            "reference_min_transcript_counts": reference_min_transcript_counts,
            "reference_min_genes": reference_min_genes,
            "reference_min_cell_counts": reference_min_cell_counts,
            "output_dir": output_dir,
            "output_name": output_name,
            "seed": seed,
        },
    )
    effective = compile_application_config(document, source=source)
    paths = output_paths(effective)
    metadata = application_metadata(
        effective,
        cli_overrides=overrides,
        paths=paths,
        dry_run=dry_run,
    )
    pipeline_metadata = dict(
        metadata,
        reference_filter={
            "column": effective.reference_filter_column,
            "value": effective.reference_filter_value,
        },
    )
    runtime, io, algorithm = _engine_overrides(effective)
    result = None

    def application_preprocess(spatial, reference):
        reference = filter_reference(
            reference,
            effective.reference_filter_column,
            effective.reference_filter_value,
        )
        spatial = preprocess_spatial(
            spatial,
            effective.spatial_min_transcript_counts,
            effective.spatial_min_cell_counts,
            min_counts=effective.spatial_min_counts,
        )
        reference = preprocess_reference(
            reference,
            effective.reference_min_transcript_counts,
            effective.reference_min_cell_counts,
            min_genes=effective.reference_min_genes,
        )
        return spatial, reference

    def finalize(ctx):
        nonlocal result
        result = publish_outputs(effective, paths, ctx)

    svc = REVISEPipeline().run(
        svc_type=effective.svc_type,
        cf=None,
        runtime_overrides=runtime,
        io_overrides=io,
        algorithm_overrides=algorithm,
        dry_run=dry_run,
        application_preprocess_callback=application_preprocess,
        finalize_callback=finalize,
        application_config_metadata=pipeline_metadata,
    )
    report = {
        "status": "preflight_passed" if dry_run else "succeeded",
        "svc_type": effective.svc_type,
        "outputs": {key: str(path) for key, path in paths.items()},
        "application_config": metadata,
        "pipeline": _pipeline_evidence(svc),
    }
    if dry_run:
        report["preflight"] = str(Path(svc.provenance["run_dir"]) / "preflight.json")
        return None, report

    report["shapes"] = _result_shapes(effective.svc_type, result)
    return result, report


def run_application(
    config: str | Path,
    *,
    svc_type: str | None = None,
    root_dir: str | None = None,
    st_path: str | None = None,
    st_format: str | None = None,
    sc_ref_path: str | None = None,
    spatialdata_table: str | None = None,
    spatialdata_element: str | None = None,
    pm_on_cell_path: str | None = None,
    ot_method: str | None = None,
    cell_type_col: str | None = None,
    sub_cell_type_col: str | None = None,
    select_cell_type: str | None = None,
    local_refinement_strength: float | None = None,
    spatial_min_transcript_counts: int | None = None,
    spatial_min_counts: int | None = None,
    spatial_min_cell_counts: int | None = None,
    reference_min_transcript_counts: int | None = None,
    reference_min_genes: int | None = None,
    reference_min_cell_counts: int | None = None,
    output_dir: str | None = None,
    output_name: str | None = None,
    seed: int | None = None,
    dry_run: bool = False,
) -> AnnData | tuple[AnnData, AnnData] | None:
    """Run one Application route and return its published AnnData artifact(s)."""
    return _execute_application(
        config,
        svc_type=svc_type,
        root_dir=root_dir,
        st_path=st_path,
        st_format=st_format,
        sc_ref_path=sc_ref_path,
        spatialdata_table=spatialdata_table,
        spatialdata_element=spatialdata_element,
        pm_on_cell_path=pm_on_cell_path,
        ot_method=ot_method,
        cell_type_col=cell_type_col,
        sub_cell_type_col=sub_cell_type_col,
        select_cell_type=select_cell_type,
        local_refinement_strength=local_refinement_strength,
        spatial_min_transcript_counts=spatial_min_transcript_counts,
        spatial_min_counts=spatial_min_counts,
        spatial_min_cell_counts=spatial_min_cell_counts,
        reference_min_transcript_counts=reference_min_transcript_counts,
        reference_min_genes=reference_min_genes,
        reference_min_cell_counts=reference_min_cell_counts,
        output_dir=output_dir,
        output_name=output_name,
        seed=seed,
        dry_run=dry_run,
    )[0]


def main(argv: list[str] | None = None) -> None:
    from revise.application.cli import main as cli_main

    cli_main(argv)


__all__ = ["ApplicationConfigError", "run_application"]


if __name__ == "__main__":
    main()
