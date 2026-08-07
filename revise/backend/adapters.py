from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, Iterable, List

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.spatial import cKDTree

from revise.backend.contracts import LocalRefinementStrategy
from revise.config.runner_conf import (
    ApplicationScConf,
    ApplicationScSrConf,
    ApplicationSpConf,
    BenchmarkImputeConf,
    BenchmarkSegConf,
    BenchmarkSrConf,
    resolved_input_path,
)
from revise.io import REVISEInputService
from revise.svc import SVC
from revise.utils import benchmark_case_leaf
from revise.utils.spot_sr_input import ensure_all_cells_in_spot


def _cfg_get(config: Dict[str, Any], *path: str, default: Any = None) -> Any:
    cur: Any = config
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def _ot_runner_kwargs(cfg: Dict[str, Any], *, impute: bool = False) -> Dict[str, Any]:
    """Translate the public GA/LR OT config into legacy runner field names."""
    ga = cfg["ot"]["ga"]
    lr = cfg["ot"]["lr"]
    kwargs = {
        "annotate_mode": str(ga["solver"]),
        "annotate_pot_reg": float(ga["pot"]["reg"]),
        "annotate_pot_reg_m": float(ga["pot"]["reg_m"]),
        "annotate_pot_reg_type": str(ga["pot"]["reg_type"]),
        "rec_ot_method": str(lr["solver"]),
    }
    if impute:
        numerics = cfg["ot"]["impute"]
        kwargs.update(
            rec_impute_pot_reg=float(numerics["reg"]),
            rec_impute_pot_reg_m=float(numerics["reg_m"]),
            rec_impute_pot_reg_type=str(numerics["reg_type"]),
        )
    else:
        numerics = lr["pot"]
        kwargs.update(
            rec_pot_reg=float(numerics["reg"]),
            rec_pot_reg_m=float(numerics["reg_m"]),
            rec_pot_reg_type=str(numerics["reg_type"]),
        )
    return kwargs


def _attach_local_refinement_strength(conf, cfg: Dict[str, Any]) -> None:
    refinement = cfg.get("local_refinement")
    if refinement is not None:
        conf.local_refinement_strength = float(refinement["strength"])


def _resolve_runtime_seed(ctx) -> int:
    runtime = getattr(ctx, "runtime", None) or ctx.merged_config.get("runtime", {})
    seed = runtime.get("seed", 42)
    if seed is None:
        return int(np.random.randint(0, np.iinfo(np.int32).max))
    return int(seed)


def _resolve_sample_seed(ctx) -> int:
    """Use parity-compatible seed when compatibility mode is requested."""
    if ctx.compatibility_mode:
        return 0
    return _resolve_runtime_seed(ctx)


def _subsample_obs(adata, n_obs: int, seed: int):
    """Deterministic obs-level subsampling helper.

    We avoid scanpy's in-place helper here so we can reuse the sampled index
    for paired objects (e.g. benchmark prediction/ground-truth adatas).
    """
    if n_obs <= 0 or n_obs >= adata.n_obs:
        return adata.copy(), adata.obs_names.to_numpy()
    rng = np.random.RandomState(seed)
    keep = rng.choice(adata.obs_names.to_numpy(), size=n_obs, replace=False)
    return adata[keep, :].copy(), keep


def _replace_slash_labels(adata, columns: Iterable[str]) -> None:
    for col in columns:
        if col not in adata.obs:
            continue
        # Force string dtype first so categorical columns do not keep stale
        # categories (e.g. "Mono/Macro") after replacement.
        series = adata.obs[col].astype(str)
        if series.str.contains("/", regex=False).any():
            normalized = series.str.replace("/", "_", regex=False)
            label_pairs = pd.DataFrame(
                {"original": series, "normalized": normalized}
            ).drop_duplicates()
            collisions = (
                label_pairs.groupby("normalized", sort=False)["original"].nunique()
            )
            if (collisions > 1).any():
                names = collisions[collisions > 1].index.tolist()
                raise ValueError(
                    f"Reference labels in {col!r} collide after slash normalization: "
                    f"{names[:5]}"
                )
            adata.obs[col] = normalized


def _ensure_transcript_counts(adata) -> None:
    if "transcript_counts" in adata.obs:
        return
    values = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["transcript_counts"] = values


def _input_service(ctx) -> REVISEInputService:
    return REVISEInputService.from_context(ctx)


def _input_path(ctx, role: str, fallback: str) -> str:
    return resolved_input_path(
        getattr(ctx, "input_specs", None),
        role,
        fallback,
    )


def _inject_sr_spatial_leakage_noise(
    adata,
    leak_ratio: float,
    k: int,
    weight_mode: str = "distance",
    preserve_total_counts: bool = True,
    seed: int = 42,
    logger=None,
):
    """Inject spot-level spatial transcript leakage by mixing neighbor compositions.

    The perturbation is applied on spot expression before global anchoring:
    each spot's gene composition is mixed with a weighted average of its
    spatial neighbors, optionally preserving the original per-spot total counts.
    """
    leak_ratio = float(leak_ratio)
    if leak_ratio <= 0:
        if logger is not None:
            logger.info("[sr-noise] leak_ratio<=0, skip noise injection")
        return adata
    if adata.n_obs < 2:
        if logger is not None:
            logger.info("[sr-noise] fewer than 2 spots, skip noise injection")
        return adata
    if "spatial" not in adata.obsm:
        raise KeyError("st_adata.obsm['spatial'] is required for SR spatial leakage noise")

    k = max(1, int(k))
    coords = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    if coords.ndim != 2 or coords.shape[1] < 2:
        raise ValueError("st_adata.obsm['spatial'] must have shape (n_spots, >=2)")
    coords = coords[:, :2].copy()

    # Tiny deterministic jitter avoids unstable tie ordering for duplicate coords.
    rng = np.random.RandomState(int(seed))
    coords += rng.normal(loc=0.0, scale=1e-9, size=coords.shape)

    X_in = adata.X
    X = X_in.toarray() if sparse.issparse(X_in) else np.asarray(X_in)
    X = np.asarray(X, dtype=np.float64)
    if X.ndim != 2:
        raise ValueError("st_adata.X must be a 2D matrix")

    tree = cKDTree(coords)
    query_k = min(adata.n_obs, k + 1)
    dists, nbrs = tree.query(coords, k=query_k)
    dists = np.asarray(dists)
    nbrs = np.asarray(nbrs)
    if dists.ndim == 1:
        dists = dists[:, None]
        nbrs = nbrs[:, None]

    row_sums = X.sum(axis=1, keepdims=True)
    safe_row_sums = np.maximum(row_sums, 1e-12)
    P = X / safe_row_sums
    P_nb = np.zeros_like(P)

    actual_k = 0
    for i in range(adata.n_obs):
        idx_row = nbrs[i]
        dist_row = dists[i]
        mask = idx_row != i
        idx = idx_row[mask][:k]
        dist = dist_row[mask][:k]

        if idx.size == 0:
            P_nb[i] = P[i]
            continue

        if weight_mode == "uniform":
            w = np.ones(idx.size, dtype=np.float64)
        elif weight_mode == "distance":
            # Adaptive local scale: use median non-zero neighbor distance.
            positive = dist[dist > 0]
            scale = float(np.median(positive)) if positive.size > 0 else 1.0
            scale = max(scale, 1e-8)
            w = np.exp(-dist / scale)
        else:
            raise ValueError(f"Unsupported sr_noise_weight='{weight_mode}', expected 'uniform' or 'distance'")

        w_sum = float(w.sum())
        if w_sum <= 0:
            w = np.full(idx.size, 1.0 / idx.size, dtype=np.float64)
        else:
            w = w / w_sum
        P_nb[i] = w @ P[idx]
        actual_k = max(actual_k, int(idx.size))

    if preserve_total_counts:
        P_noisy = (1.0 - leak_ratio) * P + leak_ratio * P_nb
        X_noisy = P_noisy * row_sums
    else:
        X_nb = P_nb * safe_row_sums
        X_noisy = (1.0 - leak_ratio) * X + leak_ratio * X_nb

    X_noisy = np.clip(X_noisy, 0.0, None)

    out = adata.copy()
    out.X = sparse.csr_matrix(X_noisy) if sparse.issparse(X_in) else X_noisy
    out.uns["sr_spatial_noise"] = {
        "enabled": True,
        "method": "spatial_neighbor_composition_mixing",
        "leak_ratio": leak_ratio,
        "k": int(k),
        "weight_mode": str(weight_mode),
        "preserve_total_counts": bool(preserve_total_counts),
        "seed": int(seed),
        "actual_max_neighbors_used": int(actual_k),
    }

    if logger is not None:
        logger.info(
            "[sr-noise] injected spatial leakage noise: lambda=%.4f, k=%d, weight=%s, preserve_total_counts=%s",
            leak_ratio,
            k,
            weight_mode,
            preserve_total_counts,
        )
    return out


def _extract_probs(adata, key: str) -> pd.DataFrame | None:
    if adata is None or key not in adata.obsm:
        return None
    values = adata.obsm[key]
    if isinstance(values, pd.DataFrame):
        return values.copy()
    if isinstance(values, np.ndarray):
        cols = [f"{key}_{i}" for i in range(values.shape[1])]
        return pd.DataFrame(values, index=adata.obs_names, columns=cols)
    return None


def _require_concrete_cell_type(value: Any) -> str:
    select_ct = value.strip() if isinstance(value, str) else ""
    if select_ct.lower() in {"", "all", "*", "__all__", "all_cell_types"}:
        raise ValueError(
            "route.select_cell_type must name one concrete broad cell type"
        )
    return select_ct


def _build_svc(
    ctx,
    outputs: Dict[str, Any],
    default_key: str | None,
    expr=None,
    spatial=None,
    extra_provenance: Dict[str, Any] | None = None,
) -> SVC:
    primary_key = None
    if default_key and default_key in outputs:
        primary_key = default_key
    elif default_key:
        raise RuntimeError(
            f"{ctx.runtime.get('application_route') or ctx.runtime.get('confounding')} "
            f"strategy did not return required output {default_key!r}; "
            f"available={sorted(outputs)}"
        )
    primary = outputs.get(primary_key) if primary_key is not None else None
    expr = primary if expr is None else expr
    spatial = primary if spatial is None else spatial

    cell_type_col = ctx.columns.get("cell_type_col", "Level1")
    confidence_col = ctx.columns.get("confidence_col", "Confidence")

    labels = None
    confidence = None
    probs = None
    if expr is not None:
        if cell_type_col in expr.obs:
            labels = expr.obs[cell_type_col].copy()
        if confidence_col in expr.obs:
            confidence = expr.obs[confidence_col].copy()
        elif "max_score" in expr.obs:
            confidence = expr.obs["max_score"].copy()
    if spatial is not None:
        probs = _extract_probs(spatial, cell_type_col)

    provenance = {
        "strategy": ctx.runtime.get("strategy"),
        "route": ctx.route,
        "output_keys": sorted(outputs.keys()),
        "primary_output_key": primary_key,
    }
    if extra_provenance:
        provenance.update(extra_provenance)

    return SVC(
        expr=expr,
        spatial=spatial,
        svc_kind=str(ctx.runtime.get("svc_kind", "sc")),
        cell_type_probs=probs,
        cell_type_label=labels,
        confidence=confidence,
        provenance=provenance,
        quality_metrics=dict(ctx.quality_metrics),
        artifacts={"outputs": outputs},
    )


def _install_safe_topology_patch() -> None:
    """Patch topology helper with dynamic PCA dims for tiny cohorts.

    The original helper hardcodes n_comps=50 (sc) / 30 (sp), which can crash on
    small sampled subsets. We patch the function reference used by
    `revise.backend.ops.meta.get_subcluster` to cap PCs to valid bounds.
    """
    import revise.backend.ops.meta as meta_mod
    import revise.backend.ops.topology as topo_mod

    if getattr(topo_mod, "_revise_unified_safe_patch", False):
        return

    original_fn = topo_mod.get_adjacency_graph

    def safe_get_adjacency_graph(
        adata,
        data_type,
        neighbors_method="pca",
        alpha=0.2,
        gene_neighbor_num=30,
        spatial_neighbor_num=30,
    ):
        if neighbors_method != "pca":
            return original_fn(
                adata,
                data_type=data_type,
                neighbors_method=neighbors_method,
                alpha=alpha,
                gene_neighbor_num=gene_neighbor_num,
                spatial_neighbor_num=spatial_neighbor_num,
            )

        if data_type == "sp":
            n_top_genes = 200
            base_n_comps = 30
        elif data_type == "sc":
            n_top_genes = 2000
            base_n_comps = 50
        elif data_type == "sc_app":
            n_top_genes = 100
            base_n_comps = 30
        else:
            return original_fn(
                adata,
                data_type=data_type,
                neighbors_method=neighbors_method,
                alpha=alpha,
                gene_neighbor_num=gene_neighbor_num,
                spatial_neighbor_num=spatial_neighbor_num,
            )

        ad = adata.copy()
        sc.pp.normalize_total(ad)
        sc.pp.log1p(ad)

        if ad.var.shape[0] > n_top_genes:
            sc.pp.highly_variable_genes(ad, n_top_genes=n_top_genes)
            ad = ad[:, ad.var["highly_variable"]]
        print(f"{data_type} Finial gene number: {ad.shape[1]}")

        if data_type == "sc":
            sc.pp.scale(ad, max_value=10)

        if ad.n_obs < 2 or ad.n_vars < 2:
            return sparse.eye(ad.n_obs, format="csr")

        max_pcs = min(base_n_comps, ad.n_obs - 1, ad.n_vars - 1)
        max_pcs = max(1, int(max_pcs))
        sc.pp.pca(ad, n_comps=max_pcs)
        sc.pp.neighbors(ad, n_pcs=max_pcs, use_rep="X_pca", n_neighbors=gene_neighbor_num)
        return ad.obsp["connectivities"]

    topo_mod.get_adjacency_graph = safe_get_adjacency_graph
    meta_mod.get_adjacency_graph = safe_get_adjacency_graph
    topo_mod._revise_unified_safe_patch = True


class RunnerBackedStrategy(LocalRefinementStrategy):
    strategy_id: str = "RunnerBackedStrategy"

    def global_anchoring(self, ctx) -> None:
        # Use backend kernel registry so all algorithm implementations are
        # centrally managed under revise.backend.kernels.
        from revise.backend.kernels import build_kernel

        if ctx.runtime.get("task") in {"sp_svc", "sc_svc_sr"}:
            ctx.runner_config.local_refinement_applied_callback = (
                lambda: ctx.record_local_refinement(True)
            )
        kernel = build_kernel("global_anchoring", config=ctx.runner_config, logger=ctx.logger)
        ctx.runner.st_adata = kernel.run(
            ctx.runner.st_adata,
            ctx.runner.sc_ref_adata,
            **ctx.runner_config.__dict__,
        )

    def solve_ot(self, ctx) -> None:
        ctx.runner.local_refinement()


class SpSvcApplicationStrategy(RunnerBackedStrategy):
    strategy_id = "SpSvcApplicationStrategy"

    def prepare_context(self, ctx) -> None:
        # Import lazily to keep route validation lightweight when a strategy
        # is not selected by current router resolution.
        from revise.backend.runners.sp_svc_application import SpSVC as SpAppRunner

        cfg = ctx.merged_config
        io_cfg = ctx.io
        columns = ctx.columns
        conf = ApplicationSpConf(
            sample_name=io_cfg["sample_name"],
            raw_data_path=io_cfg["data_root"],
            result_root_path=str(ctx.run_dir),
            cell_type_col=columns["cell_type_col"],
            confidence_col=columns["confidence_col"],
            unknown_key=columns["unknown_key"],
            st_file=io_cfg["st_file"],
            sc_ref_file=io_cfg["sc_ref_file"],
            prep_st_min_counts=int(_cfg_get(cfg, "preprocess", "st_min_counts", default=20)),
            prep_st_min_cells=int(_cfg_get(cfg, "preprocess", "st_min_cells", default=30)),
            prep_sc_min_counts=int(_cfg_get(cfg, "preprocess", "sc_min_counts", default=20)),
            prep_sc_min_cells=int(_cfg_get(cfg, "preprocess", "sc_min_cells", default=50)),
            plot_flag=bool(_cfg_get(cfg, "plot", "enabled", default=False)),
            plot_cluster_resolution=list(_cfg_get(cfg, "plot", "cluster_resolutions", default=[0.3, 0.5, 0.7])),
            plot_min_genes=int(_cfg_get(cfg, "plot", "min_genes", default=20)),
            plot_min_cells=int(_cfg_get(cfg, "plot", "min_cells", default=3)),
            plot_sample_size=int(_cfg_get(cfg, "plot", "sample_size", default=10000)),
            rec_graph_n_neighbors=int(_cfg_get(cfg, "graph", "n_neighbors", default=10)),
            rec_graph_exp_neighbor_num=int(_cfg_get(cfg, "graph", "exp_neighbors", default=10)),
            rec_graph_spatial_neighbor_num=int(_cfg_get(cfg, "graph", "spatial_neighbors", default=10)),
            rec_graph_method=str(_cfg_get(cfg, "graph", "method", default="joint")),
            rec_graph_alpha=float(_cfg_get(cfg, "graph", "alpha", default=0.5)),
            **_ot_runner_kwargs(cfg),
        )
        _attach_local_refinement_strength(conf, cfg)

        input_service = _input_service(ctx)
        adata_st = input_service.read_st_adata(
            _input_path(ctx, "st", conf.st_file_path)
        )
        sample_size = io_cfg.get("sample_size")
        if sample_size is not None:
            # Match original script behavior when compatibility_mode=true.
            sample_seed = _resolve_sample_seed(ctx)
            sc.pp.subsample(adata_st, n_obs=int(sample_size), random_state=sample_seed)

        # Preserve the established application sp-SVC preprocessing semantics.
        sc.pp.filter_cells(adata_st, min_counts=conf.prep_st_min_counts)
        sc.pp.filter_genes(adata_st, min_cells=conf.prep_st_min_cells)

        adata_sc = input_service.read_sc_ref_adata(
            _input_path(ctx, "sc_ref", conf.sc_ref_file_path)
        )
        # The original script uses min_genes for sc cells and min_cells for genes.
        sc.pp.filter_cells(adata_sc, min_genes=conf.prep_sc_min_counts)
        sc.pp.filter_genes(adata_sc, min_cells=conf.prep_sc_min_cells)
        _replace_slash_labels(
            adata_sc,
            [columns["cell_type_col"], columns["sub_cell_type_col"]],
        )
        ctx.runner_config = conf
        ctx.st_adata = adata_st
        ctx.sc_ref_adata = adata_sc
        ctx.runner = SpAppRunner(adata_st, adata_sc, conf, ctx.logger)

    def finalize_svc(self, ctx) -> SVC:
        outputs = dict(getattr(ctx.runner, "svc", {}))
        return _build_svc(ctx, outputs, default_key="sp_svc")


class ScSvcApplicationStrategy(RunnerBackedStrategy):
    strategy_id = "ScSvcApplicationStrategy"

    def prepare_context(self, ctx) -> None:
        from revise.backend.runners.sc_svc_application import ScSVC as ScAppRunner

        cfg = ctx.merged_config
        io_cfg = ctx.io
        columns = ctx.columns
        tacco_annotate_cfg = dict(
            (cfg.get("sc", {}) or {}).get("tacco_annotate", {}) or {}
        )

        conf = ApplicationScConf(
            sample_name=io_cfg["sample_name"],
            raw_data_path=io_cfg["data_root"],
            result_root_path=str(ctx.run_dir),
            cell_type_col=columns["cell_type_col"],
            confidence_col=columns["confidence_col"],
            unknown_key=columns["unknown_key"],
            st_file=io_cfg["st_file"],
            sc_ref_file=io_cfg["sc_ref_file"],
            prep_st_min_counts=int(_cfg_get(cfg, "preprocess", "st_min_transcripts", default=60)),
            prep_st_min_cells=int(_cfg_get(cfg, "preprocess", "st_min_cells", default=100)),
            prep_sc_min_cells=int(_cfg_get(cfg, "preprocess", "sc_min_cells", default=100)),
            rec_graph_n_neighbors=int(_cfg_get(cfg, "graph", "n_neighbors", default=10)),
            rec_graph_exp_neighbor_num=int(_cfg_get(cfg, "graph", "exp_neighbors", default=15)),
            rec_graph_spatial_neighbor_num=int(_cfg_get(cfg, "graph", "spatial_neighbors", default=6)),
            rec_graph_method=str(_cfg_get(cfg, "graph", "method", default="joint")),
            rec_graph_alpha=float(_cfg_get(cfg, "graph", "alpha", default=0.2)),
            rec_match_spot_sum=bool(_cfg_get(cfg, "sc", "match_spot_sum", default=False)),
            tacco_annotate_multi_center=tacco_annotate_cfg.get("multi_center"),
            tacco_annotate_lamb=tacco_annotate_cfg.get("lamb"),
            **_ot_runner_kwargs(cfg),
        )
        _attach_local_refinement_strength(conf, cfg)

        input_service = _input_service(ctx)
        adata_sp = input_service.read_st_adata(
            _input_path(ctx, "st", conf.st_file_path)
        )
        _ensure_transcript_counts(adata_sp)
        adata_sp = adata_sp[adata_sp.obs["transcript_counts"] >= conf.prep_st_min_counts, :].copy()
        sc.pp.filter_genes(adata_sp, min_cells=conf.prep_st_min_cells)

        adata_sc = input_service.read_sc_ref_adata(
            _input_path(ctx, "sc_ref", conf.sc_ref_file_path)
        )
        cell_type_col = columns.get("cell_type_col", "Level1")
        sub_cell_type_col = columns.get("sub_cell_type_col", "Level2")
        required_cols = list(dict.fromkeys([cell_type_col, sub_cell_type_col]))
        missing = [c for c in required_cols if c not in adata_sc.obs.columns]
        if missing:
            raise KeyError(f"Missing required columns in sc reference: {missing}")
        adata_sc.obs = adata_sc.obs.loc[:, required_cols].copy()
        sc.pp.filter_genes(adata_sc, min_cells=conf.prep_sc_min_cells)
        _replace_slash_labels(adata_sc, required_cols)

        overlap_genes = adata_sp.var_names.intersection(adata_sc.var_names)
        if overlap_genes.empty:
            raise ValueError("No overlapping genes between spatial and sc reference data")
        reference_overlap_mass = np.asarray(
            adata_sc[:, overlap_genes].X.sum(axis=1)
        ).ravel()
        zero_reference_cells = reference_overlap_mass == 0
        if np.any(zero_reference_cells):
            adata_sc = adata_sc[~zero_reference_cells, :].copy()
        adata_sp = adata_sp[:, overlap_genes].copy()

        ctx.runner_config = conf
        ctx.st_adata = adata_sp
        ctx.sc_ref_adata = adata_sc
        ctx.runner = ScAppRunner(adata_sp, adata_sc, conf, ctx.logger)

    def solve_ot(self, ctx) -> None:
        sc_cfg = ctx.merged_config.get("sc", {})
        sub_cell_type_col = ctx.columns.get("sub_cell_type_col", "Level2")

        select_ct = _require_concrete_cell_type(sc_cfg.get("select_ct"))

        resolutions = list(sc_cfg.get("resolutions", [0.6, 0.7, 0.8]))
        select_res = sc_cfg.get("select_resolution")
        sc_svc_spatial, sc_svc_expr = ctx.runner.local_refinement(
            select_ct,
            sub_cell_type_col,
            resolutions,
            select_res=select_res,
        )
        ctx.record_local_refinement(True)
        ctx.artifacts["outputs"] = {
            "sc_svc_spatial": sc_svc_spatial,
            "sc_svc_expr": sc_svc_expr,
        }
        ctx.artifacts["selected_cell_type"] = select_ct

    def finalize_svc(self, ctx) -> SVC:
        outputs = dict(ctx.artifacts.get("outputs", {}))
        return _build_svc(
            ctx,
            outputs,
            default_key="sc_svc_expr",
            expr=outputs.get("sc_svc_expr"),
            spatial=outputs.get("sc_svc_spatial"),
            extra_provenance={
                "selected_cell_type": ctx.artifacts.get("selected_cell_type"),
            },
        )


class ScSvcSrApplicationStrategy(RunnerBackedStrategy):
    strategy_id = "ScSvcSrApplicationStrategy"

    def prepare_context(self, ctx) -> None:
        from revise.backend.runners.sc_svc_sr_application import ScSVCSr as ScSrAppRunner

        cfg = ctx.merged_config
        io_cfg = ctx.io
        columns = ctx.columns

        conf = ApplicationScSrConf(
            sample_name=io_cfg["sample_name"],
            raw_data_path=io_cfg["data_root"],
            result_root_path=str(ctx.run_dir),
            cell_type_col=columns["cell_type_col"],
            confidence_col=columns["confidence_col"],
            unknown_key=columns["unknown_key"],
            st_file=io_cfg["st_file"],
            sc_ref_file=io_cfg["sc_ref_file"],
            prep_st_min_counts=int(_cfg_get(cfg, "preprocess", "st_min_transcripts", default=60)),
            prep_st_min_cells=int(_cfg_get(cfg, "preprocess", "st_min_cells", default=100)),
            prep_sc_min_cells=int(_cfg_get(cfg, "preprocess", "sc_min_cells", default=100)),
            rec_graph_n_neighbors=int(_cfg_get(cfg, "graph", "n_neighbors", default=20)),
            rec_graph_method=str(_cfg_get(cfg, "graph", "method", default="joint")),
            rec_graph_alpha=float(_cfg_get(cfg, "graph", "alpha", default=0.2)),
            rec_graph_exp_neighbor_num=int(_cfg_get(cfg, "graph", "exp_neighbors", default=10)),
            rec_graph_spatial_neighbor_num=int(_cfg_get(cfg, "graph", "spatial_neighbors", default=20)),
            rec_match_spot_sum=bool(_cfg_get(cfg, "sc", "match_spot_sum", default=False)),
            svc_completeness=_cfg_get(cfg, "sc", "svc_completeness"),
            sr_assignment_seed=_resolve_runtime_seed(ctx),
            **_ot_runner_kwargs(cfg),
        )
        _attach_local_refinement_strength(conf, cfg)
        conf.pm_on_cell = getattr(ctx, "pm_on_cell", None)

        input_service = _input_service(ctx)
        adata_st = input_service.read_st_adata(
            _input_path(ctx, "st", conf.st_file_path)
        )
        ensure_all_cells_in_spot(adata_st, logger=ctx.logger)
        adata_sc = input_service.read_sc_ref_adata(
            _input_path(ctx, "sc_ref", conf.sc_ref_file_path)
        )
        ctx.runner_config = conf
        ctx.st_adata = adata_st
        ctx.sc_ref_adata = adata_sc
        ctx.runner = ScSrAppRunner(adata_st, adata_sc, conf, ctx.logger)

    def solve_ot(self, ctx) -> None:
        ctx.runner_config.sr_allocation_callback = ctx.record_sr_allocation
        super().solve_ot(ctx)

    def finalize_svc(self, ctx) -> SVC:
        outputs = dict(getattr(ctx.runner, "svc", {}))
        default_key = "sc_svc_dec"
        if bool(getattr(ctx.runner_config, "rec_graph_agg_enabled", False)) and "sc_svc_dec_graphagg" in outputs:
            default_key = "sc_svc_dec_graphagg"
        return _build_svc(ctx, outputs, default_key=default_key)


class SpSvcBenchmarkSegStrategy(RunnerBackedStrategy):
    strategy_id = "SpSvcBenchmarkSegStrategy"

    def prepare_context(self, ctx) -> None:
        # Benchmark segmentation/bin2cell adapter:
        # route-level config is translated into runner dataclass contract.
        from revise.backend.runners.sp_svc_benchmark import SpSVC as SpBenchmarkRunner

        cfg = ctx.merged_config
        io_cfg = ctx.io
        columns = ctx.columns
        case_subdir = benchmark_case_leaf(ctx.runtime.get("confounding"), io_cfg)

        conf = BenchmarkSegConf(
            sample_name=io_cfg["sample_name"],
            raw_data_path=io_cfg["data_root"],
            result_root_path=str(io_cfg["output_root"]),
            cell_type_col=columns["cell_type_col"],
            confidence_col=columns["confidence_col"],
            unknown_key=columns["unknown_key"],
            st_file=io_cfg["st_file"],
            gt_svc_file=io_cfg["gt_svc_file"],
            sc_ref_file=io_cfg["sc_ref_file"],
            seg_method=io_cfg.get("seg_method", "seg_1"),
            case_subdir=case_subdir,
            rec_graph_n_neighbors=int(_cfg_get(cfg, "graph", "n_neighbors", default=50)),
            rec_graph_exp_neighbor_num=int(_cfg_get(cfg, "graph", "exp_neighbors", default=30)),
            rec_graph_spatial_neighbor_num=int(_cfg_get(cfg, "graph", "spatial_neighbors", default=30)),
            rec_graph_method=str(_cfg_get(cfg, "graph", "method", default="joint")),
            rec_graph_alpha=float(_cfg_get(cfg, "graph", "alpha", default=0.8)),
            rec_alpha=float(_cfg_get(cfg, "reconstruct", "alpha", default=1.0)),
            **_ot_runner_kwargs(cfg),
        )
        _attach_local_refinement_strength(conf, cfg)
        os.makedirs(conf.result_dir, exist_ok=True)

        input_service = _input_service(ctx)
        adata_st = input_service.read_st_adata(
            _input_path(ctx, "st", conf.st_file_path)
        )
        adata_real = input_service.read_real_adata(
            _input_path(ctx, "gt", conf.gt_svc_file_path)
        )
        adata_sc = input_service.read_sc_ref_adata(
            _input_path(ctx, "sc_ref", conf.sc_ref_file_path)
        )

        # Optional fast-compare mode: sample benchmark cells while keeping
        # prediction/ground-truth perfectly aligned on obs index.
        sample_size = io_cfg.get("sample_size")
        if sample_size is not None:
            sampled_st, keep = _subsample_obs(adata_st, int(sample_size), _resolve_sample_seed(ctx))
            keep_index = pd.Index(keep)
            common = keep_index.intersection(adata_real.obs_names)
            # Keep prediction and ground-truth obs perfectly aligned.
            adata_st = sampled_st[common, :].copy()
            adata_real = adata_real[common, :].copy()

        ctx.runner_config = conf
        ctx.st_adata = adata_st
        ctx.real_st_adata = adata_real
        ctx.sc_ref_adata = adata_sc
        ctx.runner = SpBenchmarkRunner(adata_st, adata_sc, conf, adata_real, ctx.logger)

    def finalize_svc(self, ctx) -> SVC:
        outputs = dict(getattr(ctx.runner, "svc", {}))
        return _build_svc(ctx, outputs, default_key="sp_svc")


class ScSvcSrBenchmarkStrategy(RunnerBackedStrategy):
    strategy_id = "ScSvcSrBenchmarkStrategy"

    def prepare_context(self, ctx) -> None:
        from revise.backend.runners.sc_svc_sr_benchmark import ScSVCSr as ScSrBenchmarkRunner

        cfg = ctx.merged_config
        io_cfg = ctx.io
        columns = ctx.columns
        case_subdir = benchmark_case_leaf(ctx.runtime.get("confounding"), io_cfg)

        conf = BenchmarkSrConf(
            sample_name=io_cfg["sample_name"],
            raw_data_path=io_cfg["data_root"],
            result_root_path=str(io_cfg["output_root"]),
            cell_type_col=columns["cell_type_col"],
            confidence_col=columns["confidence_col"],
            unknown_key=columns["unknown_key"],
            st_file=io_cfg["st_file"],
            gt_svc_file=io_cfg["gt_svc_file"],
            sc_ref_file=io_cfg["sc_ref_file"],
            spot_size=int(io_cfg.get("spot_size", 50)),
            case_subdir=case_subdir,
            svc_completeness=_cfg_get(cfg, "sc", "svc_completeness"),
            sr_assignment_seed=_resolve_runtime_seed(ctx),
            rec_graph_n_neighbors=int(_cfg_get(cfg, "graph", "n_neighbors", default=20)),
            rec_graph_exp_neighbor_num=int(_cfg_get(cfg, "graph", "exp_neighbors", default=10)),
            rec_graph_spatial_neighbor_num=int(_cfg_get(cfg, "graph", "spatial_neighbors", default=20)),
            rec_graph_method=str(_cfg_get(cfg, "graph", "method", default="joint")),
            rec_graph_alpha=float(_cfg_get(cfg, "graph", "alpha", default=0.2)),
            rec_alpha=float(_cfg_get(cfg, "reconstruct", "alpha", default=1.0)),
            rec_graph_agg_enabled=bool(_cfg_get(cfg, "sc", "sr_graph_agg_enabled", default=False)),
            rec_graph_agg_low_conf_only=bool(_cfg_get(cfg, "sc", "sr_graph_agg_low_conf_only", default=False)),
            rec_graph_agg_low_conf_quantile=float(
                _cfg_get(cfg, "sc", "sr_graph_agg_low_conf_quantile", default=0.2)
            ),
            rec_graph_agg_anchor_only=bool(_cfg_get(cfg, "sc", "sr_graph_agg_anchor_only", default=False)),
            rec_graph_agg_anchor_high_conf_quantile=float(
                _cfg_get(cfg, "sc", "sr_graph_agg_anchor_high_conf_quantile", default=0.8)
            ),
            rec_graph_agg_confidence_mode=str(
                _cfg_get(cfg, "sc", "sr_graph_agg_confidence_mode", default="auto")
            ),
            rec_graph_agg_conf_weighted_alpha=bool(
                _cfg_get(cfg, "sc", "sr_graph_agg_conf_weighted_alpha", default=False)
            ),
            rec_graph_agg_conf_alpha_min=float(
                _cfg_get(cfg, "sc", "sr_graph_agg_conf_alpha_min", default=0.0)
            ),
            rec_graph_agg_conf_alpha_max=float(
                _cfg_get(cfg, "sc", "sr_graph_agg_conf_alpha_max", default=-1.0)
            ),
            rec_graph_agg_conf_alpha_power=float(
                _cfg_get(cfg, "sc", "sr_graph_agg_conf_alpha_power", default=1.0)
            ),
            sr_noise_enabled=bool(_cfg_get(cfg, "sc", "sr_noise_enabled", default=False)),
            sr_noise_lambda=float(_cfg_get(cfg, "sc", "sr_noise_lambda", default=0.0)),
            sr_noise_k=int(_cfg_get(cfg, "sc", "sr_noise_k", default=4)),
            sr_noise_weight=str(_cfg_get(cfg, "sc", "sr_noise_weight", default="distance")),
            sr_noise_preserve_total_counts=bool(_cfg_get(cfg, "sc", "sr_noise_preserve_total_counts", default=True)),
            sr_noise_seed=int(_cfg_get(cfg, "sc", "sr_noise_seed", default=42)),
            **_ot_runner_kwargs(cfg),
        )
        _attach_local_refinement_strength(conf, cfg)
        conf.pm_on_cell = getattr(ctx, "pm_on_cell", None)
        os.makedirs(conf.result_dir, exist_ok=True)

        input_service = _input_service(ctx)
        adata_st = input_service.read_st_adata(
            _input_path(ctx, "st", conf.st_file_path)
        )
        adata_real = input_service.read_real_adata(
            _input_path(ctx, "gt", conf.gt_svc_file_path)
        )
        adata_sc = input_service.read_sc_ref_adata(
            _input_path(ctx, "sc_ref", conf.sc_ref_file_path)
        )
        ensure_all_cells_in_spot(adata_st, logger=ctx.logger, real_adata=adata_real)

        sample_size = io_cfg.get("sample_size")
        if sample_size is not None:
            # SR benchmark uses spot-level ST and cell-level GT with different
            # obs spaces; align GT via spot->cell mapping instead of obs names.
            sampled_st, _ = _subsample_obs(adata_st, int(sample_size), _resolve_sample_seed(ctx))
            adata_st = sampled_st
            mapping = adata_st.uns.get("all_cells_in_spot", {})
            keep_cells: List[str] = []
            for spot in adata_st.obs_names:
                keep_cells.extend(mapping.get(str(spot), []))
            if keep_cells:
                keep_index = pd.Index(np.asarray(keep_cells, dtype=str))
                if "cell_id" in adata_real.obs:
                    keep_mask = adata_real.obs["cell_id"].astype(str).isin(keep_index)
                    adata_real = adata_real[keep_mask.to_numpy(), :].copy()
                else:
                    real_index = keep_index.intersection(adata_real.obs_names)
                    adata_real = adata_real[real_index, :].copy()
                if isinstance(conf.pm_on_cell, pd.DataFrame):
                    conf.pm_on_cell = conf.pm_on_cell.loc[
                        keep_index.intersection(conf.pm_on_cell.index)
                    ].copy()

        ctx.runner_config = conf
        ctx.st_adata = adata_st
        ctx.real_st_adata = adata_real
        ctx.sc_ref_adata = adata_sc
        ctx.runner = ScSrBenchmarkRunner(adata_st, adata_sc, conf, adata_real, ctx.logger)

    def global_anchoring(self, ctx) -> None:
        conf = ctx.runner_config
        if bool(getattr(conf, "sr_noise_enabled", False)) and float(getattr(conf, "sr_noise_lambda", 0.0)) > 0:
            noisy_st = _inject_sr_spatial_leakage_noise(
                ctx.runner.st_adata,
                leak_ratio=float(conf.sr_noise_lambda),
                k=int(conf.sr_noise_k),
                weight_mode=str(conf.sr_noise_weight),
                preserve_total_counts=bool(conf.sr_noise_preserve_total_counts),
                seed=int(conf.sr_noise_seed),
                logger=ctx.logger,
            )
            ctx.runner.st_adata = noisy_st
            ctx.st_adata = noisy_st
            try:
                noisy_path = Path(ctx.run_dir) / "st_input_noisy.h5ad"
                noisy_path.parent.mkdir(parents=True, exist_ok=True)
                noisy_st.write_h5ad(noisy_path)
                ctx.artifacts["sr_noise"] = {
                    "enabled": True,
                    "artifact": str(noisy_path),
                    "lambda": float(conf.sr_noise_lambda),
                    "k": int(conf.sr_noise_k),
                    "weight": str(conf.sr_noise_weight),
                    "preserve_total_counts": bool(conf.sr_noise_preserve_total_counts),
                    "seed": int(conf.sr_noise_seed),
                }
                ctx.logger.info("[sr-noise] saved noisy ST input to %s", noisy_path)
            except Exception as exc:  # pragma: no cover - artifact persistence best effort
                ctx.logger.warning("[sr-noise] failed to save noisy ST artifact: %s", exc)
        else:
            ctx.logger.info("[sr-noise] disabled for SR benchmark run")

        super().global_anchoring(ctx)

    def solve_ot(self, ctx) -> None:
        ctx.runner_config.sr_allocation_callback = ctx.record_sr_allocation
        super().solve_ot(ctx)

    def finalize_svc(self, ctx) -> SVC:
        outputs = dict(getattr(ctx.runner, "svc", {}))
        default_key = "sc_svc_dec"
        if bool(getattr(ctx.runner_config, "rec_graph_agg_enabled", False)) and "sc_svc_dec_graphagg" in outputs:
            default_key = "sc_svc_dec_graphagg"
        return _build_svc(ctx, outputs, default_key=default_key)


class ScSvcImputeBenchmarkStrategy(RunnerBackedStrategy):
    strategy_id = "ScSvcImputeBenchmarkStrategy"

    def prepare_context(self, ctx) -> None:
        from revise.backend.runners.sc_svc_impute_benchmark import ScSVCImpute as ScImputeBenchmarkRunner

        _install_safe_topology_patch()

        cfg = ctx.merged_config
        io_cfg = ctx.io
        columns = ctx.columns
        case_subdir = benchmark_case_leaf(ctx.runtime.get("confounding"), io_cfg)

        conf = BenchmarkImputeConf(
            sample_name=io_cfg["sample_name"],
            raw_data_path=io_cfg["data_root"],
            result_root_path=str(io_cfg["output_root"]),
            cell_type_col=columns["cell_type_col"],
            confidence_col=columns["confidence_col"],
            unknown_key=columns["unknown_key"],
            st_file=io_cfg["st_file"],
            gt_svc_file=io_cfg["gt_svc_file"],
            sc_ref_file=io_cfg["sc_ref_file"],
            case_subdir=case_subdir,
            prep_min_cells=int(_cfg_get(cfg, "preprocess", "st_min_cells", default=30)),
            prep_min_counts=int(_cfg_get(cfg, "preprocess", "st_min_transcripts", default=60)),
            rec_graph_n_neighbors=int(_cfg_get(cfg, "graph", "n_neighbors", default=15)),
            rec_merge_subcluster_method=str(_cfg_get(cfg, "impute", "merge_subcluster_method", default="mean")),
            rec_subcluster_resolution=int(_cfg_get(cfg, "impute", "subcluster_resolution", default=3)),
            rec_in_panel_subcluster_resolution=(
                None
                if _cfg_get(cfg, "impute", "in_panel_subcluster_resolution", default=None) is None
                else int(_cfg_get(cfg, "impute", "in_panel_subcluster_resolution", default=None))
            ),
            rec_impute_prune_flag=bool(_cfg_get(cfg, "impute", "prune", default=True)),
            rec_impute_n_neighbors=int(_cfg_get(cfg, "impute", "n_neighbors", default=1)),
            rec_impute_method=str(_cfg_get(cfg, "impute", "method", default="mean")),
            **_ot_runner_kwargs(cfg, impute=True),
        )
        _attach_local_refinement_strength(conf, cfg)
        os.makedirs(conf.result_dir, exist_ok=True)

        input_service = _input_service(ctx)
        adata_st = input_service.read_st_adata(
            _input_path(ctx, "st", conf.st_file_path)
        )
        adata_real = input_service.read_real_adata(
            _input_path(ctx, "gt", conf.gt_svc_file_path)
        )
        adata_sc = input_service.read_sc_ref_adata(
            _input_path(ctx, "sc_ref", conf.sc_ref_file_path)
        )

        sample_size = io_cfg.get("sample_size")
        if sample_size is not None:
            sampled_st, keep = _subsample_obs(adata_st, int(sample_size), _resolve_sample_seed(ctx))
            keep_index = pd.Index(keep)
            common = keep_index.intersection(adata_real.obs_names)
            adata_st = sampled_st[common, :].copy()
            adata_real = adata_real[common, :].copy()
            # Fast benchmark mode: also bound sc reference size for imputation
            # routes so uncertainty/subcluster steps stay tractable.
            if adata_sc.n_obs > int(sample_size):
                adata_sc, _ = _subsample_obs(adata_sc, int(sample_size), _resolve_sample_seed(ctx))

        if columns["cell_type_col"] in adata_sc.obs:
            counts = adata_sc.obs[columns["cell_type_col"]].value_counts()
            valid_ct = counts[counts >= 2].index
            adata_sc = adata_sc[adata_sc.obs[columns["cell_type_col"]].isin(valid_ct), :].copy()

        # Gene-uncertainty builds per-cell-type PCA graphs; when sampled data
        # is small, cap n_pcs to avoid sklearn arpack dimension errors.
        if columns["cell_type_col"] in adata_sc.obs:
            ct_sizes = adata_sc.obs[columns["cell_type_col"]].value_counts()
            min_ct_size = int(ct_sizes.min()) if not ct_sizes.empty else int(adata_sc.n_obs)
        else:
            min_ct_size = int(adata_sc.n_obs)
        max_valid_pcs = min(int(conf.rec_graph_n_pcs), max(1, min_ct_size - 1), max(1, adata_sc.n_vars - 1))
        conf.rec_graph_n_pcs = max_valid_pcs
        ctx.logger.info("[adapter] impute rec_graph_n_pcs adjusted to %s", conf.rec_graph_n_pcs)

        ctx.runner_config = conf
        ctx.st_adata = adata_st
        ctx.real_st_adata = adata_real
        ctx.sc_ref_adata = adata_sc
        ctx.runner = ScImputeBenchmarkRunner(adata_st, adata_sc, conf, adata_real, ctx.logger)

    def finalize_svc(self, ctx) -> SVC:
        outputs = dict(getattr(ctx.runner, "svc", {}))
        preferred = "sc_svc_impute_in_panel" if "sc_svc_impute_in_panel" in outputs else None
        return _build_svc(ctx, outputs, default_key=preferred)
