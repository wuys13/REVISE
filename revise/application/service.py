"""Application route contract and reconstruction use case."""

from __future__ import annotations

import argparse
import copy
import os
import tempfile
from pathlib import Path
from typing import NamedTuple
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData


class ApplicationRoute(NamedTuple):
    route_id: str
    profile: str
    confounding: str
    output_key: str | None
    svc_kind: str


APPLICATION_ROUTES = {
    "hST-SVC": ApplicationRoute(
        route_id="sp_svc",
        profile="application_sp",
        confounding="bin2cell",
        output_key="sp_svc",
        svc_kind="sp",
    ),
    "iST-SVC": ApplicationRoute(
        route_id="sc_svc",
        profile="application_sc",
        confounding="segmentation",
        output_key=None,
        svc_kind="sc",
    ),
    "sST-SVC": ApplicationRoute(
        route_id="sc_svc_sr",
        profile="application_sc_sr",
        confounding="spot_size",
        output_key="sc_svc_dec",
        svc_kind="sc",
    ),
}

REVISEPipeline = None
RECONSTRUCTION_CONTRACT_KEYS = {
    "schema_version",
    "svc_type",
    "ist_mapping",
    "effective_seed",
    "expression_source",
    "donor_column",
    "donor_sha256",
}


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


def _build_algorithm_overrides(args: argparse.Namespace) -> dict:
    overrides = {}
    if args.ot_method is not None:
        overrides["ot"] = {
            "ga": {"solver": args.ot_method},
            "lr": {"solver": args.ot_method},
        }
    local_refinement_strength = getattr(args, "local_refinement_strength", None)
    if local_refinement_strength is not None:
        overrides["local_refinement"] = {
            "strength": local_refinement_strength,
        }
    columns = {}
    if args.cell_type_col is not None:
        columns["cell_type_col"] = args.cell_type_col
    if args.sub_cell_type_col is not None:
        columns["sub_cell_type_col"] = args.sub_cell_type_col
    if columns:
        overrides["columns"] = columns
    return overrides


def _validate_ist_mapping(args: argparse.Namespace) -> None:
    mapping = getattr(args, "ist_mapping", None)
    if args.svc_type == "iST-SVC":
        if mapping is None:
            args.ist_mapping = "mean"
        elif mapping not in {"mean", "random"}:
            raise ValueError(f"Unsupported iST mapping {mapping!r}")
    elif mapping is not None:
        raise ValueError("ist_mapping is only valid for iST-SVC")


def _run_pipeline(
    args: argparse.Namespace,
    *,
    dry_run: bool = False,
    finalize_callback=None,
):
    _validate_ist_mapping(args)
    route = APPLICATION_ROUTES[args.svc_type]
    pipeline_class = REVISEPipeline
    if pipeline_class is None:
        from revise.framework import REVISEPipeline as pipeline_class

    pipeline = pipeline_class(config_path=args.config)
    runtime_overrides = {
        "platform": route.route_id,
        "confounding": route.confounding,
    }
    if args.seed is not None:
        runtime_overrides["seed"] = args.seed

    svc = pipeline._run_with_algorithm_overrides(
        profile=route.profile,
        runtime_overrides=runtime_overrides,
        io_overrides={
            "data_root": args.data_root,
            "output_root": args.output_root,
            "sample_name": args.sample_name,
            "st_file": args.st_file,
            "sc_ref_file": args.sc_ref_file,
            "patient_key": args.patient_key,
            "save_outputs": False,
        },
        algorithm_overrides=_build_algorithm_overrides(args),
        dry_run=dry_run,
        finalize_callback=finalize_callback,
    )
    return route.profile, route.output_key, svc


def _effective_seed(args, ctx) -> int:
    seed = args.seed
    if seed is None:
        runtime = getattr(ctx, "runtime", None)
        if runtime is None:
            runtime = getattr(ctx, "merged_config", {}).get("runtime", {})
        seed = runtime.get("seed", 42)
    return int(seed)


def _is_finite(matrix) -> bool:
    import numpy as np
    from scipy import sparse

    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    return bool(np.isfinite(values).all())


def _copy_mapping(mapping) -> dict:
    return {key: value.copy() for key, value in mapping.items()}


def _merge_contract_metadata(existing, contract: dict) -> dict:
    if existing is None:
        merged = {}
    elif not isinstance(existing, dict):
        raise ValueError("revise_reconstruction namespace must be a mapping")
    else:
        merged = copy.deepcopy(existing)
    inapplicable = (merged.keys() & RECONSTRUCTION_CONTRACT_KEYS) - contract.keys()
    if inapplicable:
        raise ValueError(
            "revise_reconstruction contains inapplicable contract keys "
            f"{sorted(inapplicable)}"
        )
    for key, value in contract.items():
        if key in merged and merged[key] != value:
            raise ValueError(
                f"revise_reconstruction contract key {key!r} conflicts with "
                f"existing value {merged[key]!r}"
            )
        merged[key] = value
    return merged


def _ist_contract(mapping: str, seed: int, donor_ids: list[str] | None) -> dict:
    contract = {
        "schema_version": 2,
        "svc_type": "iST-SVC",
        "ist_mapping": mapping,
        "expression_source": "expression_carrier.X_as_is",
    }
    if donor_ids is not None:
        from revise.utils import hash_jsonable

        contract.update(
            effective_seed=seed,
            donor_column="revise_ist_donor_id",
            donor_sha256=hash_jsonable(donor_ids),
        )
    return contract


def _validate_ist_carriers(spatial, expression, *, random_mapping: bool):
    spatial_labels = _cluster_labels(spatial, "SVC_cluster", "spatial SVC")
    expression_labels = _cluster_labels(expression, "SVC_cluster", "expression SVC")
    spatial_keys = _cluster_keys(spatial_labels)
    expression_keys = _cluster_keys(expression_labels)
    if set(spatial_keys) != set(expression_keys):
        raise ValueError("spatial and expression SVC cluster sets must match exactly")
    if not expression.var_names.is_unique:
        raise ValueError("expression SVC var_names must be unique")
    if random_mapping:
        if not expression.obs_names.is_unique:
            raise ValueError("expression SVC obs_names must be unique for random mapping")
        if any(not str(name) for name in expression.obs_names):
            raise ValueError("expression SVC obs_names must be non-empty for random mapping")
    if not _is_finite(expression.X):
        raise ValueError("expression SVC X must contain only finite values")
    return spatial_keys, expression_keys


def _mean_expression(expression, spatial_keys, expression_keys):
    import numpy as np
    from scipy import sparse

    means = {}
    for cluster in set(spatial_keys):
        indices = [
            index for index, value in enumerate(expression_keys) if value == cluster
        ]
        block = expression.X[indices]
        mean = block.mean(axis=0)
        means[cluster] = (
            sparse.csr_matrix(mean)
            if sparse.issparse(expression.X)
            else np.asarray(mean).reshape(-1)
        )
    if sparse.issparse(expression.X):
        return sparse.vstack(
            [means[cluster] for cluster in spatial_keys],
            format="csr",
        )
    return np.vstack([means[cluster] for cluster in spatial_keys])


def _random_expression(expression, spatial_keys, expression_keys, seed: int):
    import numpy as np

    donors_by_cluster = {}
    for index, (name, cluster) in enumerate(
        zip(expression.obs_names, expression_keys, strict=True)
    ):
        donors_by_cluster.setdefault(cluster, []).append((str(name), index))
    for donors in donors_by_cluster.values():
        donors.sort(key=lambda item: item[0])

    rng = np.random.default_rng(seed)
    donor_ids = []
    donor_indices = []
    for cluster in spatial_keys:
        donors = donors_by_cluster[cluster]
        donor_name, donor_index = donors[int(rng.integers(len(donors)))]
        donor_ids.append(donor_name)
        donor_indices.append(donor_index)
    return expression.X[donor_indices].copy(), donor_ids


def _build_ist_result(args, ctx, seed: int):
    from anndata import AnnData

    outputs = dict(ctx.svc.artifacts.get("outputs", {}))
    required = {"sc_svc_spatial", "sc_svc_expr"}
    missing = sorted(required - outputs.keys())
    if missing:
        raise RuntimeError(
            f"iST-SVC pipeline did not return required outputs {missing}; "
            f"available={sorted(outputs)}"
        )

    spatial = outputs["sc_svc_spatial"]
    expression = outputs["sc_svc_expr"]
    random_mapping = args.ist_mapping == "random"
    spatial_keys, expression_keys = _validate_ist_carriers(
        spatial,
        expression,
        random_mapping=random_mapping,
    )
    donor_ids = None
    if random_mapping:
        output_x, donor_ids = _random_expression(
            expression,
            spatial_keys,
            expression_keys,
            seed,
        )
    else:
        output_x = _mean_expression(expression, spatial_keys, expression_keys)

    obs = spatial.obs.copy(deep=True)
    if donor_ids is not None:
        obs["revise_ist_donor_id"] = donor_ids
    existing = spatial.uns.get("revise_reconstruction")
    reconstruction = _merge_contract_metadata(
        existing,
        _ist_contract(args.ist_mapping, seed, donor_ids),
    )
    result = AnnData(
        X=output_x,
        obs=obs,
        var=expression.var.copy(deep=True),
        uns={"revise_reconstruction": reconstruction},
        obsm=_copy_mapping(spatial.obsm),
        varm=_copy_mapping(expression.varm),
        obsp=_copy_mapping(spatial.obsp),
        varp=_copy_mapping(expression.varp),
    )
    if not _is_finite(result.X):
        raise ValueError("iST-SVC output X must contain only finite values")
    return result


def _build_result(args, output_key, ctx, seed: int):
    if args.svc_type == "iST-SVC":
        return _build_ist_result(args, ctx, seed)

    outputs = dict(ctx.svc.artifacts.get("outputs", {}))
    if output_key not in outputs:
        raise RuntimeError(
            f"{args.svc_type} pipeline did not return required output {output_key!r}; "
            f"available={sorted(outputs)}"
        )
    result = outputs[output_key].copy()
    result.uns["revise_reconstruction"] = _merge_contract_metadata(
        result.uns.get("revise_reconstruction"),
        {"schema_version": 2, "svc_type": args.svc_type},
    )
    if not _is_finite(result.X):
        raise ValueError(f"{args.svc_type} output X must contain only finite values")
    return result


def _validate_written_result(path: Path, source) -> None:
    from anndata import read_h5ad

    written = read_h5ad(path, backed="r")
    try:
        if written.shape != source.shape:
            raise ValueError(
                f"Published H5AD shape {written.shape} does not match source "
                f"{source.shape}"
            )
        if not written.obs_names.equals(source.obs_names):
            raise ValueError("Published H5AD observation names do not match the source")
        if not written.var_names.equals(source.var_names):
            raise ValueError("Published H5AD variable names do not match the source")
        source_metadata = source.uns["revise_reconstruction"]
        written_metadata = written.uns.get("revise_reconstruction", {})
        expected = {
            key: value
            for key, value in source_metadata.items()
            if key in RECONSTRUCTION_CONTRACT_KEYS
        }
        actual = {
            key: value
            for key, value in written_metadata.items()
            if key in RECONSTRUCTION_CONTRACT_KEYS
        }
        if actual != expected:
            raise ValueError("Published H5AD reconstruction metadata is incomplete")
    finally:
        written.file.close()


def _build_public_result(args, output_key, ctx) -> tuple[AnnData, Path]:
    from revise.utils import completed_artifact

    route = APPLICATION_ROUTES[args.svc_type]
    if ctx.svc.svc_kind != route.svc_kind:
        raise ValueError(
            f"SVC type {args.svc_type!r} requires internal kind {route.svc_kind!r}; "
            f"strategy returned {ctx.svc.svc_kind!r}"
        )

    seed = _effective_seed(args, ctx)
    result = _build_result(args, output_key, ctx, seed)
    output_dir = Path(args.output_root) / args.sample_name / args.svc_type
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "SVC.h5ad"

    metadata = result.uns["revise_reconstruction"]
    result_record = {"filename": output_path.name, "type": args.svc_type}
    assembly_record = None
    if args.svc_type == "iST-SVC":
        random_mapping = args.ist_mapping == "random"
        assembly_record = {
            "ist_mapping": args.ist_mapping,
            "effective_seed": seed if random_mapping else None,
            "donor_column": metadata.get("donor_column"),
            "donor_sha256": metadata.get("donor_sha256"),
            "donor_count": result.n_obs if random_mapping else None,
        }

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
        result.copy().write_h5ad(temporary_path)
        _validate_written_result(temporary_path, result)
        artifact = completed_artifact("public_result", temporary_path)
        artifact["path"] = str(output_path)
        had_previous_result = "result" in ctx.provenance
        previous_result = copy.deepcopy(ctx.provenance.get("result"))
        had_previous_svc_result = "result" in ctx.svc.provenance
        previous_svc_result = copy.deepcopy(ctx.svc.provenance.get("result"))
        had_previous_assembly = "assembly" in ctx.provenance
        previous_assembly = copy.deepcopy(ctx.provenance.get("assembly"))
        had_previous_svc_assembly = "assembly" in ctx.svc.provenance
        previous_svc_assembly = copy.deepcopy(ctx.svc.provenance.get("assembly"))

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
                output_path.unlink(missing_ok=True)
                os.replace(backup_path, output_path)

            if had_previous_result:
                ctx.provenance["result"] = previous_result
            else:
                ctx.provenance.pop("result", None)
            if had_previous_svc_result:
                ctx.svc.provenance["result"] = previous_svc_result
            else:
                ctx.svc.provenance.pop("result", None)
            if had_previous_assembly:
                ctx.provenance["assembly"] = previous_assembly
            else:
                ctx.provenance.pop("assembly", None)
            if had_previous_svc_assembly:
                ctx.svc.provenance["assembly"] = previous_svc_assembly
            else:
                ctx.svc.provenance.pop("assembly", None)

            for index in range(len(ctx.artifact_records) - 1, -1, -1):
                if ctx.artifact_records[index] == artifact:
                    del ctx.artifact_records[index]
                    break

        ctx.set_pending_publication(commit=commit, rollback=rollback)
        try:
            if backup_path is not None:
                os.replace(output_path, backup_path)
            os.replace(temporary_path, output_path)
            ctx.provenance["result"] = copy.deepcopy(result_record)
            ctx.svc.provenance["result"] = copy.deepcopy(result_record)
            if assembly_record is not None:
                ctx.provenance["assembly"] = copy.deepcopy(assembly_record)
                ctx.svc.provenance["assembly"] = copy.deepcopy(assembly_record)
            ctx.record_artifact(artifact)
        except BaseException:
            ctx.rollback_pending_publication()
            raise
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
    return result, output_path


def reconstruct(args: argparse.Namespace):
    route = APPLICATION_ROUTES[args.svc_type]
    published = {}

    def publish(ctx):
        published["result"], published["path"] = _build_public_result(
            args,
            route.output_key,
            ctx,
        )

    profile, _, svc = _run_pipeline(args, finalize_callback=publish)
    summary = svc.summary()
    summary.update(profile=profile, route=svc.provenance.get("route_key"))
    return published["result"], published["path"], summary
