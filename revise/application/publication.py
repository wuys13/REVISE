"""Stable H5AD publication for application reconstruction runs."""

from __future__ import annotations

import copy
import json
import os
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData


def _cluster_labels(adata, cluster_col: str, source: str):
    if cluster_col not in adata.obs:
        raise KeyError(f"{source} is missing required obs column {cluster_col!r}")
    labels = adata.obs[cluster_col]
    if labels.isna().any():
        raise ValueError(f"{source} contains null {cluster_col!r} labels")
    return labels.to_numpy(dtype=object)


def _cluster_key(label) -> tuple[type, object]:
    return type(label), label


def _cluster_keys(labels) -> list[tuple[type, object]]:
    return [_cluster_key(label) for label in labels]


def _safe_output_component(value: object, *, field: str) -> str:
    component = str(value).strip().replace("/", "_").replace("\\", "_")
    if component in {"", ".", ".."} or "\x00" in component:
        raise ValueError(f"{field} cannot be used as an output directory: {value!r}")
    return component


def _validate_sc_pair_source(spatial, expression) -> None:
    spatial_labels = _cluster_labels(spatial, "SVC_cluster", "spatial SVC")
    expression_labels = _cluster_labels(expression, "SVC_cluster", "expression SVC")
    spatial_clusters = set(_cluster_keys(spatial_labels))
    expression_clusters = set(_cluster_keys(expression_labels))
    extra_expression_clusters = expression_clusters - spatial_clusters
    if extra_expression_clusters:
        raise ValueError(
            "expression SVC contains clusters absent from spatial SVC: "
            f"{sorted(map(repr, extra_expression_clusters))}"
        )
    if "spatial" not in spatial.obsm:
        raise KeyError("spatial SVC is missing required obsm['spatial'] coordinates")


def _validate_written_sc_result(
    path: Path,
    source,
    *,
    require_spatial: bool,
    require_cluster: bool = True,
) -> None:
    from anndata import read_h5ad

    written = read_h5ad(path, backed="r")
    try:
        if written.shape != source.shape:
            raise ValueError(
                "Published H5AD shape "
                f"{written.shape} does not match source {source.shape}"
            )
        if not written.obs_names.equals(source.obs_names):
            raise ValueError("Published H5AD observation names do not match the source")
        if not written.var_names.equals(source.var_names):
            raise ValueError("Published H5AD variable names do not match the source")
        if require_cluster and (
            "SVC_cluster" not in written.obs
            or written.obs["SVC_cluster"].isna().any()
        ):
            raise ValueError("Published H5AD must contain non-null SVC_cluster labels")
        if require_spatial and "spatial" not in written.obsm:
            raise ValueError("Published spatial H5AD is missing obsm['spatial']")
    finally:
        written.file.close()


def publish_single(request, ctx) -> tuple[AnnData, Path]:
    """Publish one route result and register its rollback transaction."""
    from revise.utils import completed_artifact

    svc = ctx.svc
    runtime = ctx.runtime
    if runtime.get("application_route") != request.svc_type:
        raise ValueError(
            "Publication request does not match the resolved application route"
        )
    seed = request.seed
    if seed is None:
        runtime = getattr(ctx, "runtime", None)
        if runtime is None:
            runtime = getattr(ctx, "merged_config", {}).get("runtime", {})
        seed = runtime.get("seed", 42)
    seed = int(seed)
    if svc.svc_kind != runtime.get("svc_kind"):
        raise ValueError(
            f"SVC type {request.svc_type!r} requires internal kind "
            f"{runtime.get('svc_kind')!r}; "
            f"strategy returned {svc.svc_kind!r}"
        )

    if request.svc_type == "sc-SVC":
        raise ValueError("standard sc-SVC must use the paired public-result publisher")

    output_key = svc.provenance.get("primary_output_key")
    outputs = dict(svc.artifacts.get("outputs", {}))
    if not output_key or output_key not in outputs:
        raise RuntimeError(
            f"{request.svc_type} pipeline did not return required output {output_key!r}; "
            f"available={sorted(outputs)}"
        )
    result = outputs[output_key].copy()
    result.uns["revise_reconstruction"] = {
        "svc_type": request.svc_type,
        "seed": seed,
    }

    output_dir = Path(request.output_root) / request.sample_name
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "SVC.h5ad"
    relative_run_dir = Path(os.path.relpath(ctx.run_dir, start=output_dir)).as_posix()
    provenance = result.uns.setdefault("revise_reconstruction", {})
    provenance.update(
        {
            "profile": ctx.profile,
            "run_dir": relative_run_dir,
            "run_manifest": (Path(relative_run_dir) / "provenance.json").as_posix(),
            "ot_method_override": request.ot_method,
            "ot_config": copy.deepcopy(ctx.merged_config["ot"]),
        }
    )

    temporary_path = None
    backup_path = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=output_dir,
            prefix=f".{output_path.stem}.",
            suffix=".h5ad",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
        result.write_h5ad(temporary_path)
        _validate_written_sc_result(
            temporary_path,
            result,
            require_spatial=False,
            require_cluster=False,
        )
        artifact = completed_artifact("public_result", temporary_path)
        artifact["path"] = str(output_path)
        result_record = {
            "filename": output_path.name,
            "type": request.svc_type,
        }
        had_previous_result = "result" in ctx.provenance
        previous_result = copy.deepcopy(ctx.provenance.get("result"))
        had_previous_svc_result = "result" in ctx.svc.provenance
        previous_svc_result = copy.deepcopy(ctx.svc.provenance.get("result"))

        if output_path.exists():
            with tempfile.NamedTemporaryFile(
                dir=output_dir,
                prefix=f".{output_path.stem}.previous.",
                suffix=output_path.suffix,
                delete=False,
            ) as handle:
                backup_path = Path(handle.name)
            backup_path.unlink()

        def commit():
            if backup_path is not None:
                backup_path.unlink(missing_ok=True)

        def rollback():
            if backup_path is None:
                output_path.unlink(missing_ok=True)
            elif backup_path.exists():
                os.replace(backup_path, output_path)

            if not had_previous_result:
                ctx.provenance.pop("result", None)
            else:
                ctx.provenance["result"] = previous_result
            if not had_previous_svc_result:
                ctx.svc.provenance.pop("result", None)
            else:
                ctx.svc.provenance["result"] = previous_svc_result

            for index in range(len(ctx.artifact_records) - 1, -1, -1):
                if ctx.artifact_records[index] == artifact:
                    del ctx.artifact_records[index]
                    break

        ctx.set_pending_publication(commit=commit, rollback=rollback)
        try:
            if backup_path is not None:
                os.replace(output_path, backup_path)
            os.replace(temporary_path, output_path)
            ctx.provenance["result"] = result_record
            ctx.record_artifact(artifact)
        except BaseException:
            ctx.rollback_pending_publication()
            raise
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
    return result, output_path


def publish_pair(request, ctx):
    """Publish the sc-SVC spatial/expression pair as one rollback unit."""
    from revise.utils import completed_artifact

    runtime = ctx.runtime
    if runtime.get("application_route") != "sc-SVC":
        raise ValueError("sc-SVC publication requires the sc-SVC application route")
    if runtime.get("svc_kind") != "sc" or ctx.svc.svc_kind != "sc":
        raise ValueError("sc-SVC publication requires an internal sc SVC result")

    outputs = dict(ctx.svc.artifacts.get("outputs", {}))
    required = {"sc_svc_spatial", "sc_svc_expr"}
    missing = sorted(required - outputs.keys())
    if missing:
        raise RuntimeError(
            f"sc-SVC pipeline did not return required outputs {missing}; "
            f"available={sorted(outputs)}"
        )

    spatial = outputs["sc_svc_spatial"]
    expression = outputs["sc_svc_expr"]
    _validate_sc_pair_source(spatial, expression)

    selected_cell_type = ctx.svc.provenance.get(
        "selected_cell_type",
        getattr(request, "select_cell_type", None),
    )
    if selected_cell_type in (None, ""):
        raise RuntimeError("sc-SVC result is missing its selected cell type")
    output_dir = (
        Path(request.output_root)
        / request.sample_name
        / "sc-SVC"
        / _safe_output_component(selected_cell_type, field="selected cell type")
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    published = {
        "spatial": (spatial, output_dir / "sc_SVC_spatial.h5ad", True),
        "expression": (expression, output_dir / "sc_SVC_expr.h5ad", False),
    }

    relative_run_dir = Path(os.path.relpath(ctx.run_dir, start=output_dir)).as_posix()
    pair_filenames = {
        role: path.name
        for role, (_adata, path, _require_spatial) in published.items()
        if role in {"spatial", "expression"}
    }
    tacco_parameters = copy.deepcopy(
        (ctx.merged_config.get("sc", {}) or {}).get("tacco_annotate")
    )
    for role, (adata, path, _require_spatial) in published.items():
        provenance = adata.uns.setdefault("revise_reconstruction", {})
        provenance.update(
            {
                "svc_type": request.svc_type,
                "output_contract": "sc_svc_notebook_pair_v1",
                "output_role": role,
                "paired_outputs": json.dumps(pair_filenames, sort_keys=True),
                "selected_cell_type": str(selected_cell_type),
                "profile": ctx.profile,
                "run_dir": relative_run_dir,
                "run_manifest": (Path(relative_run_dir) / "provenance.json").as_posix(),
                "ot_method_override": getattr(request, "ot_method", None),
                "ot_config": copy.deepcopy(ctx.merged_config["ot"]),
                "tacco_annotate": json.dumps(
                    tacco_parameters,
                    ensure_ascii=False,
                    sort_keys=True,
                ),
            }
        )

    temporary_paths = {}
    artifacts = {}
    backups = {}
    had_existing = {}
    previous_results = copy.deepcopy(ctx.provenance.get("results"))
    had_previous_results = "results" in ctx.provenance
    previous_svc_results = copy.deepcopy(ctx.svc.provenance.get("results"))
    had_previous_svc_results = "results" in ctx.svc.provenance

    try:
        for role, (adata, path, require_spatial) in published.items():
            with tempfile.NamedTemporaryFile(
                dir=output_dir,
                prefix=f".{path.stem}.",
                suffix=".h5ad",
                delete=False,
            ) as handle:
                temporary_path = Path(handle.name)
            temporary_paths[role] = temporary_path
            adata.write_h5ad(temporary_path)
            _validate_written_sc_result(
                temporary_path,
                adata,
                require_spatial=require_spatial,
            )
            artifact = completed_artifact("public_result", temporary_path)
            artifact.update(path=str(path), logical_role=role)
            artifacts[role] = artifact

        for role, (_adata, path, _require_spatial) in published.items():
            had_existing[role] = path.exists()
            if path.exists():
                with tempfile.NamedTemporaryFile(
                    dir=output_dir,
                    prefix=f".{path.stem}.previous.",
                    suffix=path.suffix,
                    delete=False,
                ) as handle:
                    backup = Path(handle.name)
                backup.unlink()
                backups[role] = backup
            else:
                backups[role] = None
        result_records = {
            role: {
                "filename": path.name,
                "path": str(path),
                "type": request.svc_type,
                "logical_role": role,
                "shape": list(adata.shape),
            }
            for role, (adata, path, _require_spatial) in published.items()
        }

        def commit():
            for backup in backups.values():
                if backup is not None:
                    backup.unlink(missing_ok=True)

        def rollback():
            for role, (_adata, path, _require_spatial) in published.items():
                backup = backups.get(role)
                if backup is not None and backup.exists():
                    path.unlink(missing_ok=True)
                    os.replace(backup, path)
                elif not had_existing.get(role, False):
                    path.unlink(missing_ok=True)
            if had_previous_results:
                ctx.provenance["results"] = previous_results
            else:
                ctx.provenance.pop("results", None)
            if had_previous_svc_results:
                ctx.svc.provenance["results"] = previous_svc_results
            else:
                ctx.svc.provenance.pop("results", None)

            artifact_records = getattr(ctx, "artifact_records", None)
            if artifact_records is not None:
                for artifact in artifacts.values():
                    for index in range(len(artifact_records) - 1, -1, -1):
                        if artifact_records[index] == artifact:
                            del artifact_records[index]
                            break

        ctx.set_pending_publication(commit=commit, rollback=rollback)
        try:
            for role, (_adata, path, _require_spatial) in published.items():
                backup = backups[role]
                if backup is not None:
                    os.replace(path, backup)
                os.replace(temporary_paths[role], path)
            ctx.provenance["results"] = copy.deepcopy(result_records)
            ctx.svc.provenance["results"] = copy.deepcopy(result_records)
            for role in published:
                ctx.record_artifact(artifacts[role])
        except BaseException:
            ctx.rollback_pending_publication()
            raise
    finally:
        for temporary_path in temporary_paths.values():
            temporary_path.unlink(missing_ok=True)

    return (
        {role: adata for role, (adata, _path, _require_spatial) in published.items()},
        {role: path for role, (_adata, path, _require_spatial) in published.items()},
    )
