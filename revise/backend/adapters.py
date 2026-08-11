from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, Iterable, List

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

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


def _resolve_runtime_seed(ctx) -> int:
    return int(ctx.runtime["seed"])


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


def _input_service(ctx) -> REVISEInputService:
    return REVISEInputService.from_context(ctx)


def _input_path(ctx, role: str) -> str:
    return resolved_input_path(
        ctx.input_specs,
        role,
    )


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

    cell_type_col = ctx.columns["cell_type_col"]
    confidence_col = ctx.columns["confidence_col"]

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
        svc_kind=str(ctx.runtime["svc_kind"]),
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
            plot_flag=bool(cfg["plot"]["enabled"]),
            plot_cluster_resolution=list(cfg["plot"]["cluster_resolutions"]),
            plot_min_genes=int(cfg["plot"]["min_genes"]),
            plot_min_cells=int(cfg["plot"]["min_cells"]),
            plot_sample_size=int(cfg["plot"]["sample_size"]),
            rec_graph_n_neighbors=int(cfg["graph"]["n_neighbors"]),
            rec_graph_exp_neighbor_num=int(cfg["graph"]["exp_neighbors"]),
            rec_graph_spatial_neighbor_num=int(cfg["graph"]["spatial_neighbors"]),
            rec_graph_method=str(cfg["graph"]["method"]),
            rec_graph_alpha=float(cfg["graph"]["alpha"]),
            rec_alpha=float(cfg["reconstruct"]["alpha"]),
            local_refinement_strength=float(
                cfg["local_refinement"]["strength"]
            ),
            **_ot_runner_kwargs(cfg),
        )

        input_service = _input_service(ctx)
        adata_st = input_service.read_st_adata(
            _input_path(ctx, "st")
        )
        sample_size = io_cfg.get("sample_size")
        if sample_size is not None:
            # Match original script behavior when compatibility_mode=true.
            sample_seed = _resolve_sample_seed(ctx)
            sc.pp.subsample(adata_st, n_obs=int(sample_size), random_state=sample_seed)

        adata_sc = input_service.read_sc_ref_adata(
            _input_path(ctx, "sc_ref")
        )
        application_preprocess = ctx.application_preprocess_callback
        if application_preprocess is None:
            raise RuntimeError("Application preprocessing callback is required")
        adata_st, adata_sc = application_preprocess(adata_st, adata_sc)
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
        tacco_annotate_cfg = cfg["sc"]["tacco_annotate"]

        conf = ApplicationScConf(
            sample_name=io_cfg["sample_name"],
            raw_data_path=io_cfg["data_root"],
            result_root_path=str(ctx.run_dir),
            cell_type_col=columns["cell_type_col"],
            confidence_col=columns["confidence_col"],
            unknown_key=columns["unknown_key"],
            st_file=io_cfg["st_file"],
            sc_ref_file=io_cfg["sc_ref_file"],
            rec_graph_n_neighbors=int(cfg["graph"]["n_neighbors"]),
            rec_graph_exp_neighbor_num=int(cfg["graph"]["exp_neighbors"]),
            rec_graph_spatial_neighbor_num=int(cfg["graph"]["spatial_neighbors"]),
            rec_graph_method=str(cfg["graph"]["method"]),
            rec_graph_alpha=float(cfg["graph"]["alpha"]),
            rec_random_state=int(cfg["graph"]["random_state"]),
            rec_alpha=float(cfg["reconstruct"]["alpha"]),
            rec_match_spot_sum=bool(cfg["sc"]["match_spot_sum"]),
            tacco_annotate_multi_center=tacco_annotate_cfg["multi_center"],
            tacco_annotate_lamb=tacco_annotate_cfg["lamb"],
            **_ot_runner_kwargs(cfg),
        )

        input_service = _input_service(ctx)
        adata_sp = input_service.read_st_adata(
            _input_path(ctx, "st")
        )
        adata_sc = input_service.read_sc_ref_adata(
            _input_path(ctx, "sc_ref")
        )
        application_preprocess = ctx.application_preprocess_callback
        if application_preprocess is None:
            raise RuntimeError("Application preprocessing callback is required")
        adata_sp, adata_sc = application_preprocess(adata_sp, adata_sc)
        cell_type_col = columns["cell_type_col"]
        sub_cell_type_col = columns["sub_cell_type_col"]
        required_cols = list(dict.fromkeys([cell_type_col, sub_cell_type_col]))
        missing = [c for c in required_cols if c not in adata_sc.obs.columns]
        if missing:
            raise KeyError(f"Missing required columns in sc reference: {missing}")
        adata_sc.obs = adata_sc.obs.loc[:, required_cols].copy()
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
        sc_cfg = ctx.merged_config["sc"]
        sub_cell_type_col = ctx.columns["sub_cell_type_col"]

        select_ct = _require_concrete_cell_type(sc_cfg.get("select_ct"))

        resolutions = list(sc_cfg["resolutions"])
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
            rec_graph_n_neighbors=int(cfg["graph"]["n_neighbors"]),
            rec_graph_method=str(cfg["graph"]["method"]),
            rec_graph_alpha=float(cfg["graph"]["alpha"]),
            rec_graph_exp_neighbor_num=int(cfg["graph"]["exp_neighbors"]),
            rec_graph_spatial_neighbor_num=int(cfg["graph"]["spatial_neighbors"]),
            rec_alpha=float(cfg["reconstruct"]["alpha"]),
            rec_match_spot_sum=bool(cfg["sc"]["match_spot_sum"]),
            rec_graph_agg_enabled=bool(cfg["sc"]["sr_graph_agg_enabled"]),
            svc_completeness=cfg["sc"]["svc_completeness"],
            sr_assignment_seed=_resolve_runtime_seed(ctx),
            local_refinement_strength=float(
                cfg["local_refinement"]["strength"]
            ),
            **_ot_runner_kwargs(cfg),
        )
        conf.pm_on_cell = getattr(ctx, "pm_on_cell", None)

        input_service = _input_service(ctx)
        adata_st = input_service.read_st_adata(
            _input_path(ctx, "st")
        )
        ensure_all_cells_in_spot(adata_st, logger=ctx.logger)
        adata_sc = input_service.read_sc_ref_adata(
            _input_path(ctx, "sc_ref")
        )
        application_preprocess = ctx.application_preprocess_callback
        if application_preprocess is None:
            raise RuntimeError("Application preprocessing callback is required")
        adata_st, adata_sc = application_preprocess(adata_st, adata_sc)
        _replace_slash_labels(
            adata_sc,
            [columns["cell_type_col"], columns["sub_cell_type_col"]],
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
        if ctx.runner_config.rec_graph_agg_enabled and "sc_svc_dec_graphagg" in outputs:
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
            seg_method=io_cfg["seg_method"],
            case_subdir=case_subdir,
            rec_graph_n_neighbors=int(cfg["graph"]["n_neighbors"]),
            rec_graph_exp_neighbor_num=int(cfg["graph"]["exp_neighbors"]),
            rec_graph_spatial_neighbor_num=int(cfg["graph"]["spatial_neighbors"]),
            rec_graph_method=str(cfg["graph"]["method"]),
            rec_graph_alpha=float(cfg["graph"]["alpha"]),
            rec_alpha=float(cfg["reconstruct"]["alpha"]),
            dropout_total_counts=int(cfg["benchmark"]["dropout_total_counts"]),
            swapping_total_counts=int(cfg["benchmark"]["swapping_total_counts"]),
            lower_ts=float(cfg["benchmark"]["lower_ts"]),
            upper_ts=float(cfg["benchmark"]["upper_ts"]),
            local_refinement_strength=float(
                cfg["local_refinement"]["strength"]
            ),
            **_ot_runner_kwargs(cfg),
        )
        os.makedirs(conf.result_dir, exist_ok=True)

        input_service = _input_service(ctx)
        adata_st = input_service.read_st_adata(
            _input_path(ctx, "st")
        )
        adata_real = input_service.read_real_adata(
            _input_path(ctx, "gt")
        )
        adata_sc = input_service.read_sc_ref_adata(
            _input_path(ctx, "sc_ref")
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
            spot_size=int(io_cfg["spot_size"]),
            case_subdir=case_subdir,
            svc_completeness=cfg["sc"]["svc_completeness"],
            sr_assignment_seed=_resolve_runtime_seed(ctx),
            rec_graph_n_neighbors=int(cfg["graph"]["n_neighbors"]),
            rec_graph_exp_neighbor_num=int(cfg["graph"]["exp_neighbors"]),
            rec_graph_spatial_neighbor_num=int(cfg["graph"]["spatial_neighbors"]),
            rec_graph_method=str(cfg["graph"]["method"]),
            rec_graph_alpha=float(cfg["graph"]["alpha"]),
            rec_alpha=float(cfg["reconstruct"]["alpha"]),
            rec_graph_agg_enabled=bool(cfg["sc"]["sr_graph_agg_enabled"]),
            rec_graph_agg_low_conf_only=bool(cfg["sc"]["sr_graph_agg_low_conf_only"]),
            rec_graph_agg_low_conf_quantile=float(cfg["sc"]["sr_graph_agg_low_conf_quantile"]),
            rec_graph_agg_anchor_only=bool(cfg["sc"]["sr_graph_agg_anchor_only"]),
            rec_graph_agg_anchor_high_conf_quantile=float(cfg["sc"]["sr_graph_agg_anchor_high_conf_quantile"]),
            rec_graph_agg_confidence_mode=str(cfg["sc"]["sr_graph_agg_confidence_mode"]),
            rec_graph_agg_conf_weighted_alpha=bool(cfg["sc"]["sr_graph_agg_conf_weighted_alpha"]),
            rec_graph_agg_conf_alpha_min=float(cfg["sc"]["sr_graph_agg_conf_alpha_min"]),
            rec_graph_agg_conf_alpha_max=float(cfg["sc"]["sr_graph_agg_conf_alpha_max"]),
            rec_graph_agg_conf_alpha_power=float(cfg["sc"]["sr_graph_agg_conf_alpha_power"]),
            local_refinement_strength=float(
                cfg["local_refinement"]["strength"]
            ),
            **_ot_runner_kwargs(cfg),
        )
        conf.pm_on_cell = getattr(ctx, "pm_on_cell", None)
        os.makedirs(conf.result_dir, exist_ok=True)

        input_service = _input_service(ctx)
        adata_st = input_service.read_st_adata(
            _input_path(ctx, "st")
        )
        adata_real = input_service.read_real_adata(
            _input_path(ctx, "gt")
        )
        adata_sc = input_service.read_sc_ref_adata(
            _input_path(ctx, "sc_ref")
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

    def solve_ot(self, ctx) -> None:
        ctx.runner_config.sr_allocation_callback = ctx.record_sr_allocation
        super().solve_ot(ctx)

    def finalize_svc(self, ctx) -> SVC:
        outputs = dict(getattr(ctx.runner, "svc", {}))
        default_key = "sc_svc_dec"
        if ctx.runner_config.rec_graph_agg_enabled and "sc_svc_dec_graphagg" in outputs:
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
            prep_min_cells=int(cfg["preprocess"]["st_min_cells"]),
            prep_min_counts=int(cfg["preprocess"]["st_min_transcripts"]),
            rec_graph_preprocess=bool(cfg["impute"]["graph_preprocess"]),
            rec_graph_n_pcs=int(cfg["impute"]["graph_n_pcs"]),
            rec_graph_n_neighbors=int(cfg["graph"]["n_neighbors"]),
            rec_merge_subcluster_method=str(cfg["impute"]["merge_subcluster_method"]),
            rec_subcluster_resolution=int(cfg["impute"]["subcluster_resolution"]),
            rec_in_panel_subcluster_resolution=(
                None
                if cfg["impute"]["in_panel_subcluster_resolution"] is None
                else int(cfg["impute"]["in_panel_subcluster_resolution"])
            ),
            rec_impute_prune_flag=bool(cfg["impute"]["prune"]),
            rec_impute_n_neighbors=int(cfg["impute"]["n_neighbors"]),
            rec_impute_method=str(cfg["impute"]["method"]),
            **_ot_runner_kwargs(cfg, impute=True),
        )
        os.makedirs(conf.result_dir, exist_ok=True)

        input_service = _input_service(ctx)
        adata_st = input_service.read_st_adata(
            _input_path(ctx, "st")
        )
        adata_real = input_service.read_real_adata(
            _input_path(ctx, "gt")
        )
        adata_sc = input_service.read_sc_ref_adata(
            _input_path(ctx, "sc_ref")
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
