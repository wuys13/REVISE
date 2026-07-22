"""Application route contract and reconstruction use case."""

from __future__ import annotations

import argparse
import copy
import json
import os
import tempfile
from pathlib import Path
from typing import Iterable
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
    "sp-SVC": ApplicationRoute(
        route_id="sp_svc",
        profile="application_sp",
        confounding="bin2cell",
        output_key="sp_svc",
        svc_kind="sp",
    ),
    "sc-SVC": ApplicationRoute(
        route_id="sc_svc",
        profile="application_sc",
        confounding="segmentation",
        output_key=None,
        svc_kind="sc",
    ),
    "sc-SVC-sr": ApplicationRoute(
        route_id="sc_svc_sr",
        profile="application_sc_sr",
        confounding="spot_size",
        output_key="sc_svc_dec",
        svc_kind="sc",
    ),
}

REVISEPipeline = None
MANAGED_SET_KEYS = {
    "io.save_outputs",
    "runtime.seed",
}
SC_MANAGED_SET_KEYS = {
    "sc.select_ct",
}


def _copy_value(value):
    return value.copy() if hasattr(value, "copy") else copy.deepcopy(value)


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


def _matching_indices(keys: list[tuple[type, object]], cluster):
    import numpy as np

    return np.fromiter(
        (index for index, key in enumerate(keys) if key == cluster),
        dtype=np.int64,
    )


def _cluster_set_mismatch_message(spatial_labels, expr_labels) -> str | None:
    spatial_set = set(_cluster_keys(spatial_labels))
    expr_set = set(_cluster_keys(expr_labels))
    if spatial_set == expr_set:
        return None
    def sort_key(item):
        return item[0].__name__, repr(item[1])

    spatial_only = sorted(spatial_set - expr_set, key=sort_key)
    expr_only = sorted(expr_set - spatial_set, key=sort_key)
    return (
        "sc-SVC cluster sets do not match exactly: "
        f"spatial_only={spatial_only} (n={len(spatial_only)}), "
        f"expr_only={expr_only} (n={len(expr_only)})"
    )


def _mean_rows_by_cluster(matrix, labels, clusters):
    import numpy as np
    from scipy import sparse

    rows = []
    if sparse.issparse(matrix):
        matrix = matrix.tocsr()
        for cluster in clusters:
            mask = _matching_indices(labels, cluster)
            rows.append(matrix[mask].sum(axis=0) / float(mask.size))
        return sparse.vstack([sparse.csr_matrix(row) for row in rows], format="csr")

    dense = np.asarray(matrix)
    return np.vstack(
        [dense[_matching_indices(labels, cluster)].mean(axis=0) for cluster in clusters]
    )


def merge_sc_svc(
    spatial_adata,
    expr_adata,
    *,
    mode: str = "mean",
    seed: int = 42,
    cluster_col: str = "SVC_cluster",
):
    import numpy as np
    from anndata import AnnData
    from scipy import sparse

    if mode not in {"mean", "random"}:
        raise ValueError("sc-SVC mapping mode must be one of ['mean', 'random']")
    if not expr_adata.var_names.is_unique:
        duplicates = expr_adata.var_names[expr_adata.var_names.duplicated()].unique()
        raise ValueError(
            "expression SVC contains duplicate gene names: "
            f"{duplicates[:10].tolist()}"
        )

    spatial_labels = _cluster_labels(spatial_adata, cluster_col, "spatial SVC")
    expr_labels = _cluster_labels(expr_adata, cluster_col, "expression SVC")
    mismatch = _cluster_set_mismatch_message(spatial_labels, expr_labels)
    if mismatch:
        raise ValueError(mismatch)

    spatial_keys = _cluster_keys(spatial_labels)
    expr_keys = _cluster_keys(expr_labels)
    clusters = sorted(
        set(expr_keys), key=lambda item: (item[0].__name__, repr(item[1]))
    )
    if mode == "mean":
        cluster_means = _mean_rows_by_cluster(expr_adata.X, expr_keys, clusters)
        cluster_index = {cluster: index for index, cluster in enumerate(clusters)}
        X = cluster_means[[cluster_index[key] for key in spatial_keys]]
    else:
        rng = np.random.default_rng(seed)
        candidates = {
            cluster: _matching_indices(expr_keys, cluster) for cluster in clusters
        }
        selected = np.fromiter(
            (rng.choice(candidates[key]) for key in spatial_keys),
            dtype=np.int64,
            count=spatial_adata.n_obs,
        )
        X = expr_adata.X[selected].copy()

    values = X.data if sparse.issparse(X) else np.asarray(X)
    if not np.all(np.isfinite(values)):
        raise ValueError("mapped sc-SVC expression contains non-finite values")

    merged = AnnData(
        X=X,
        obs=spatial_adata.obs.copy(),
        var=expr_adata.var.copy(),
        uns=copy.deepcopy(spatial_adata.uns),
    )
    merged.obs_names = spatial_adata.obs_names.copy()
    merged.var_names = expr_adata.var_names.copy()
    for key, value in spatial_adata.obsm.items():
        merged.obsm[key] = _copy_value(value)
    for key, value in spatial_adata.obsp.items():
        merged.obsp[key] = _copy_value(value)
    for key, value in expr_adata.varm.items():
        merged.varm[key] = _copy_value(value)
    for key, value in expr_adata.varp.items():
        merged.varp[key] = _copy_value(value)

    merged.uns["revise_reconstruction"] = {
        "mapping_mode": mode,
        "mapping_seed": int(seed),
        "cluster_col": cluster_col,
        "spatial_n_obs": int(spatial_adata.n_obs),
        "expression_n_obs": int(expr_adata.n_obs),
        "n_vars": int(expr_adata.n_vars),
        "cluster_set_match": True,
    }
    return merged


def _override_keys(overrides: Iterable[str]) -> set[str]:
    return {
        item.split("=", 1)[0].strip()
        for item in overrides
        if "=" in item and item.split("=", 1)[0].strip()
    }


def _build_set_overrides(args: argparse.Namespace) -> list[str]:
    managed = set(MANAGED_SET_KEYS)
    values = ["io.save_outputs=false"]
    if args.ot_method is not None:
        managed.update({"ot.ga.solver", "ot.lr.solver"})
        values.extend(
            (
                f"ot.ga.solver={args.ot_method}",
                f"ot.lr.solver={args.ot_method}",
            )
        )
    if args.cell_type_col is not None:
        managed.add("columns.cell_type_col")
        values.append(f"columns.cell_type_col={args.cell_type_col}")
    if args.sub_cell_type_col is not None:
        managed.add("columns.sub_cell_type_col")
        values.append(f"columns.sub_cell_type_col={args.sub_cell_type_col}")
    if args.svc_type == "sc-SVC":
        managed.update(SC_MANAGED_SET_KEYS)
        values.append(f"sc.select_ct={args.select_ct}")
    conflicts = sorted(
        user_key
        for user_key in _override_keys(args.set_overrides)
        if any(
            user_key == managed_key
            or user_key.startswith(f"{managed_key}.")
            or managed_key.startswith(f"{user_key}.")
            for managed_key in managed
        )
    )
    if conflicts:
        raise ValueError(
            "Conflicting high-level CLI option and --set override for: "
            + ", ".join(conflicts)
        )
    return values + list(args.set_overrides)


def _run_pipeline(
    args: argparse.Namespace,
    *,
    dry_run: bool = False,
    finalize_callback=None,
):
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

    svc = pipeline.run(
        profile=route.profile,
        runtime_overrides=runtime_overrides,
        io_overrides={
            "data_root": args.data_root,
            "output_root": args.output_root,
            "sample_name": args.sample_name,
            "st_file": args.st_file,
            "sc_ref_file": args.sc_ref_file,
            "patient_key": args.patient_key,
        },
        set_overrides=_build_set_overrides(args),
        dry_run=dry_run,
        finalize_callback=finalize_callback,
    )
    return route.profile, route.output_key, svc


def _build_public_result(args, profile, output_key, ctx) -> tuple[AnnData, Path]:
    from revise.utils import completed_artifact

    route = APPLICATION_ROUTES[args.svc_type]
    svc = ctx.svc
    seed = args.seed
    if seed is None:
        runtime = getattr(ctx, "runtime", None)
        if runtime is None:
            runtime = getattr(ctx, "merged_config", {}).get("runtime", {})
        seed = runtime.get("seed", 42)
    seed = int(seed)
    if svc.svc_kind != route.svc_kind:
        raise ValueError(
            f"SVC type {args.svc_type!r} requires internal kind {route.svc_kind!r}; "
            f"strategy returned {svc.svc_kind!r}"
        )

    outputs = dict(svc.artifacts.get("outputs", {}))
    if args.svc_type == "sc-SVC":
        required = {"sc_svc_spatial", "sc_svc_expr"}
        missing = sorted(required - outputs.keys())
        if missing:
            raise RuntimeError(
                f"sc-SVC pipeline did not return required outputs {missing}; "
                f"available={sorted(outputs)}"
            )
        result = merge_sc_svc(
            outputs["sc_svc_spatial"],
            outputs["sc_svc_expr"],
            mode=args.sc_mapping,
            seed=seed,
        )
        result.uns["revise_reconstruction"]["svc_type"] = args.svc_type
    elif output_key not in outputs:
        raise RuntimeError(
            f"{args.svc_type} pipeline did not return required output {output_key!r}; "
            f"available={sorted(outputs)}"
        )
    else:
        result = outputs[output_key].copy()
        result.uns["revise_reconstruction"] = {
            "svc_type": args.svc_type,
            "seed": seed,
        }

    output_dir = Path(args.output_root) / args.sample_name
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "SVC.h5ad"
    relative_run_dir = Path(os.path.relpath(ctx.run_dir, start=output_dir)).as_posix()
    provenance = result.uns.setdefault("revise_reconstruction", {})
    provenance.update(
        {
            "profile": profile,
            "run_dir": relative_run_dir,
            "run_manifest": (Path(relative_run_dir) / "provenance.json").as_posix(),
            "ot_method_override": args.ot_method,
            "ot_config": copy.deepcopy(ctx.merged_config["ot"]),
            "ot_events": json.dumps(
                ctx.ot_events,
                ensure_ascii=False,
                sort_keys=True,
            ),
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
        artifact = completed_artifact("public_result", temporary_path)
        artifact["path"] = str(output_path)
        result_record = {
            "filename": output_path.name,
            "type": args.svc_type,
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


def reconstruct(args: argparse.Namespace):
    route = APPLICATION_ROUTES[args.svc_type]
    published = {}

    def publish(ctx):
        published["result"], published["path"] = _build_public_result(
            args,
            route.profile,
            route.output_key,
            ctx,
        )

    profile, _, svc = _run_pipeline(args, finalize_callback=publish)
    summary = svc.summary()
    summary.update(profile=profile, route=svc.provenance.get("route_key"))
    return published["result"], published["path"], summary
