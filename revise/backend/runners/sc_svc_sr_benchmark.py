import numpy as np
import pandas as pd
import scanpy as sc

from revise.backend.runners.benchmark_svc import BenchmarkSVC
from revise.backend.kernels import GraphAggregateKernel as GraphAggregate
from revise.backend.kernels import SpotSrKernel as SpotSr
from revise.backend.ops.distance import similarity_to_distance
from revise.backend.ops.local_ot import solve_local_ot, stabilize_local_ot_support
from revise.backend.ops.meta import construct_sc_ref
from revise.backend.ops.meta import get_sc_obs
from revise.backend.ops.meta import get_true_cell_type
from revise.backend.ops.meta import resolve_true_cell_type_key
from revise.backend.ops.sr_allocation import (
    condition_virtual_cell_ot_cost,
    mandatory_reference_allocation,
    project_spot_assignment_to_virtual_cells,
    record_mandatory_allocation,
    spot_global_assignment,
    subset_virtual_assignment,
)
from revise.backend.ops.topology import get_adjacency_graph


class ScSVCSr(BenchmarkSVC):
    """
    sc-SVC super-resolution for benchmark CFs: spot size/ batch effect.
    
    This class reconstructs single-cell resolution expression profiles
    from spatial transcriptomics data by redistributing spot-level
    expressions to virtual cells using cell type contributions.
    """
    def __init__(self, st_adata, sc_ref_adata, config, real_st_adata, logger):
        super().__init__(st_adata, sc_ref_adata, config, real_st_adata, logger)
        self._adata_validate()
        self._adata_validate_dec()
        self._adata_processing()
        self.svc_obs = self._get_svc_obs()
        self.spot_sr = SpotSr(self.config, self.logger)
        self.graph_aggregate = GraphAggregate(self.config, self.logger)
        self._graphagg_confidence_cache = None
        self._graphagg_confidence_source = None
        self._graphagg_alpha_weight_cache = None
        self._graphagg_posterior_source = None
        self.svc = {}

    def _adata_validate_dec(self):

        if "all_cells_in_spot" not in self.st_adata.uns:
            raise KeyError(
                "Invalid reconstruction input: "
                f"st_file_path={getattr(self.config, 'st_file_path', '<unavailable>')}; "
                "field=uns['all_cells_in_spot']; expected=present; actual=missing"
            )

    def _get_svc_obs(self):
        svc_obs = get_sc_obs(self.st_adata.obs.index, self.st_adata.uns['all_cells_in_spot'], self.st_adata.obsm["spatial"])
        configured_key = str(self.config.cell_type_col)
        if configured_key == "Level1" and configured_key not in self.real_st_adata.obs:
            ground_truth_label_source = resolve_true_cell_type_key(
                self.real_st_adata
            )
        else:
            ground_truth_label_source = resolve_true_cell_type_key(
                self.real_st_adata,
                configured_key,
            )
        self.ground_truth_label_source = ground_truth_label_source
        svc_obs = get_true_cell_type(
            svc_obs,
            self.real_st_adata,
            label_key=ground_truth_label_source,
        )
        return svc_obs

    def local_refinement(self, *args, **kwargs):
        """Reconstruct single-cell expression profiles from spot-level data.

        1. Assigns cell types to each virtual cell using SpotSr
        2. Constructs cell type reference profiles
        3. Calculates gene expression for each cell based on spot contributions
        4. Normalizes expressions to 10,000 counts per cell
        
        The reconstructed data is stored in self.svc["sc_svc_dec"].
        """
        graph_enabled = bool(
            getattr(self.config, "rec_graph_agg_enabled", False)
        )
        overlap_genes = list(self.st_adata.var_names.intersection(self.sc_ref_adata.var_names))
        key_type = self.config.cell_type_col
        assignment = spot_global_assignment(
            self.st_adata,
            broad_key=key_type,
            expected_categories=pd.Index(
                pd.unique(self.sc_ref_adata.obs[key_type])
            ),
        )
        try:
            st_adata_common = self.st_adata[:, overlap_genes].copy()
            sc.pp.normalize_total(st_adata_common, target_sum=1e4)
            cell_contributions = assignment.posterior.copy()
            sc.pp.normalize_total(self.st_adata, target_sum=1e4)
            sc.pp.normalize_total(self.sc_ref_adata, target_sum=1e4)
            self.spot_sr.run(self)
            type_list = sorted(
                list(self.sc_ref_adata.obs[key_type].unique().astype(str))
            )
            self.logger.info(
                f"There are {len(type_list)} cell types: {type_list}"
            )
            sc_ref_all = construct_sc_ref(
                self.sc_ref_adata,
                key_type=key_type,
                type_list=type_list,
            )
            sc_ref_all = sc_ref_all.loc[:, overlap_genes]
            type_list = list(sc_ref_all.index)
            norm_type_list = [str(t).replace("/", "_") for t in type_list]
            norm_type_to_idx = {
                name: idx for idx, name in enumerate(norm_type_list)
            }
            self.svc_obs["cell_type"] = (
                self.svc_obs["cell_type"]
                .astype(str)
                .str.replace("/", "_", regex=False)
            )
            spots = self.svc_obs["spot_name"].unique()
            spot_to_idx = {spot: idx for idx, spot in enumerate(spots)}
            self.logger.info("Using simple allocation method...")
            adata_spot = st_adata_common.copy()
            X = (
                adata_spot.X
                if type(adata_spot.X) is np.ndarray
                else adata_spot.X.toarray()
            )
            self.logger.info(
                "Using baseline closed-form reference allocation "
                "(posterior beta fixed at 1.0)"
            )
            Y = mandatory_reference_allocation(
                X,
                cell_contributions,
                sc_ref_all,
            )
            spot_indices = np.array(
                [spot_to_idx[spot] for spot in self.svc_obs["spot_name"]]
            )
            missing_types = sorted(
                set(self.svc_obs["cell_type"]) - set(norm_type_to_idx)
            )
            if missing_types:
                raise ValueError(
                    f"Missing cell types in reference index: {missing_types}"
                )
            type_indices = np.array(
                [norm_type_to_idx[t] for t in self.svc_obs["cell_type"]]
            )
            SVC_X = Y[spot_indices, :, type_indices]
        except Exception:
            record_mandatory_allocation(
                self.config,
                status="failed",
                broad_key=key_type,
                n_spots=int(self.st_adata.n_obs),
                n_virtual_cells=int(len(getattr(self, "svc_obs", ()))),
                allocation_method="posterior_reference_allocation",
                reason="allocation_failed",
            )
            raise
        record_mandatory_allocation(
            self.config,
            status="completed",
            broad_key=key_type,
            n_spots=int(self.st_adata.n_obs),
            n_virtual_cells=int(len(self.svc_obs)),
            allocation_method="posterior_reference_allocation",
        )
        self.logger.info("Extracted SVC expressions using simple allocation method")
        self._projected_assignment = (
            project_spot_assignment_to_virtual_cells(
                assignment,
                self.svc_obs,
            )
        )
        self._conditioning_strength = getattr(
            self.config,
            "local_refinement_strength",
            0.0,
        )
        self._refinement_applied = False
        self._graphagg_posterior_source = f"project(obsm[{key_type}])"

        SVC_X_raw = np.asarray(SVC_X, dtype=np.float64)
        SVC_X_graphagg = None
        if graph_enabled:
            self.logger.info("SR graph aggregation enabled: building additional graph-smoothed output")
            graphagg_target_mask = None
            graphagg_anchor_mask = None
            if bool(getattr(self.config, "rec_graph_agg_low_conf_only", False)):
                graphagg_target_mask = self._build_graphagg_low_conf_mask()
            if bool(getattr(self.config, "rec_graph_agg_anchor_only", False)):
                graphagg_anchor_mask = self._build_graphagg_high_conf_anchor_mask()
            SVC_X_graphagg = self._apply_graph_aggregation(
                SVC_X_raw,
                target_mask=graphagg_target_mask,
                anchor_mask=graphagg_anchor_mask,
            )
        else:
            self.logger.info("SR graph aggregation disabled: only raw output will be evaluated")

        SVC_X_raw = SVC_X_raw / (np.sum(SVC_X_raw, axis=1, keepdims=True) + 1e-10) * 1e4
        if SVC_X_graphagg is not None:
            SVC_X_graphagg = SVC_X_graphagg / (np.sum(SVC_X_graphagg, axis=1, keepdims=True) + 1e-10) * 1e4

        self.logger.info(f"Number of cells processed: {len(self.svc_obs)}")
        self.logger.info(f"Number of unique spots: {len(spots)}")
        self.logger.info(f"Shape of raw SVC_X: {SVC_X_raw.shape}")

        self.svc["sc_svc_dec"] = self._build_svc_adata(SVC_X_raw, st_adata_common.var_names)
        self.svc["sc_svc_dec"].uns["sr_allocation"] = {
            "broad_key": key_type,
            "posterior_key": key_type,
            "operator": "closed_form_reference_allocation",
            "beta": 1.0,
        }
        if SVC_X_graphagg is not None:
            self.svc["sc_svc_dec_graphagg"] = self._build_svc_adata(SVC_X_graphagg, st_adata_common.var_names)
            self.svc["sc_svc_dec_graphagg"].uns["graph_aggregation"] = dict(self._build_graphagg_meta())
            if "graphagg_alpha_weight" in self.svc_obs.columns:
                alpha_vals = pd.to_numeric(
                    self.svc_obs["graphagg_alpha_weight"], errors="coerce"
                ).to_numpy(dtype=np.float64)
                finite = np.isfinite(alpha_vals)
                if np.any(finite):
                    q_points = [0.0, 0.1, 0.5, 0.9, 1.0]
                    q_vals = np.quantile(alpha_vals[finite], q_points)
                    self.svc["sc_svc_dec_graphagg"].uns["graph_agg_conf_weighted_alpha"] = {
                        "enabled": True,
                        "n_finite": int(finite.sum()),
                        "n_total": int(alpha_vals.shape[0]),
                        "alpha_quantiles": {str(q): float(v) for q, v in zip(q_points, q_vals)},
                    }
        return self._refinement_applied

    def _build_svc_adata(self, X, var_names):
        svc_obs = self.svc_obs.copy()
        svc_obs["cell_id"] = svc_obs["cell_id"].astype(str)
        svc_obs.set_index("cell_id", inplace=True)
        svc_adata = sc.AnnData(X, obs=svc_obs)
        svc_adata.var_names = var_names
        return svc_adata

    def _build_graphagg_meta(self):
        return {
            "enabled": bool(getattr(self.config, "rec_graph_agg_enabled", False)),
            "low_conf_only": bool(getattr(self.config, "rec_graph_agg_low_conf_only", False)),
            "low_conf_quantile": float(getattr(self.config, "rec_graph_agg_low_conf_quantile", 0.2)),
            "anchor_only": bool(getattr(self.config, "rec_graph_agg_anchor_only", False)),
            "anchor_high_conf_quantile": float(
                getattr(self.config, "rec_graph_agg_anchor_high_conf_quantile", 0.8)
            ),
            "confidence_mode": str(getattr(self.config, "rec_graph_agg_confidence_mode", "auto")),
            "confidence_source": self._graphagg_confidence_source,
            "confidence_weighted_alpha": bool(
                getattr(self.config, "rec_graph_agg_conf_weighted_alpha", False)
            ),
            "confidence_alpha_min": float(getattr(self.config, "rec_graph_agg_conf_alpha_min", 0.0)),
            "confidence_alpha_max": float(getattr(self.config, "rec_graph_agg_conf_alpha_max", -1.0)),
            "confidence_alpha_power": float(getattr(self.config, "rec_graph_agg_conf_alpha_power", 1.0)),
            "posterior_source": self._graphagg_posterior_source,
        }

    def _infer_graphagg_assignment_confidence(self):
        """Per-cell assignment confidence for selective graph aggregation.

        Preference order:
        1) `pm_on_cell` assigned-type probability (cell-level)
        2) configured spot-level broad contribution for the assigned type
        3) spot-level global anchoring confidence (last fallback)
        """
        mode = str(getattr(self.config, "rec_graph_agg_confidence_mode", "auto") or "auto").lower()
        n_cells = len(self.svc_obs)
        assigned_types = self.svc_obs["cell_type"].astype(str).str.replace("/", "_", regex=False).to_numpy()
        cell_ids = self.svc_obs["cell_id"].astype(str).to_numpy()
        spot_names = self.svc_obs["spot_name"].astype(str).to_numpy()

        conf_pm = None
        conf_broad = None
        conf_spot = None

        pm_df = getattr(self.spot_sr, "pm_on_cell", None)
        if isinstance(pm_df, pd.DataFrame) and not pm_df.empty:
            conf_df = pm_df.copy()
            conf_df.index = conf_df.index.astype(str)
            conf_df.columns = [str(c).replace("/", "_") for c in conf_df.columns]
            if conf_df.columns.duplicated().any():
                conf_df = conf_df.T.groupby(level=0).sum().T
            conf_df = conf_df.reindex(cell_ids)
            values = conf_df.to_numpy(dtype=np.float64, copy=False)
            col_map = {str(c): i for i, c in enumerate(conf_df.columns)}
            type_idx = np.array([col_map.get(t, -1) for t in assigned_types], dtype=np.int32)
            conf = np.full(n_cells, np.nan, dtype=np.float64)
            valid = type_idx >= 0
            if np.any(valid):
                rows = np.arange(n_cells, dtype=np.int32)[valid]
                conf[valid] = values[rows, type_idx[valid]]
            conf_pm = conf

        broad_key = self.config.cell_type_col
        broad = self.st_adata.obsm.get(broad_key)
        if broad is not None:
            if isinstance(broad, pd.DataFrame):
                contrib = broad.copy()
            else:
                cols = [f"ct_{i}" for i in range(broad.shape[1])]
                contrib = pd.DataFrame(
                    broad,
                    index=self.st_adata.obs_names,
                    columns=cols,
                )
            contrib.index = contrib.index.astype(str)
            contrib.columns = [str(c).replace("/", "_") for c in contrib.columns]
            if contrib.columns.duplicated().any():
                contrib = contrib.T.groupby(level=0).sum().T

            contrib = contrib.reindex(spot_names)
            values = contrib.to_numpy(dtype=np.float64, copy=False)
            col_map = {str(c): i for i, c in enumerate(contrib.columns)}
            type_idx = np.array([col_map.get(t, -1) for t in assigned_types], dtype=np.int32)
            conf = np.full(n_cells, np.nan, dtype=np.float64)
            valid = type_idx >= 0
            if np.any(valid):
                rows = np.arange(n_cells, dtype=np.int32)[valid]
                conf[valid] = values[rows, type_idx[valid]]
            conf_broad = conf

        conf_col = getattr(self.config, "confidence_col", "Confidence")
        if conf_col in self.st_adata.obs:
            spot_conf = pd.to_numeric(self.st_adata.obs[conf_col], errors="coerce")
            spot_conf.index = self.st_adata.obs.index.astype(str)
            conf_spot = spot_conf.reindex(spot_names).to_numpy(dtype=np.float64, copy=False)

        def _combine_pm_and_spot(pm_vals, spot_vals, op: str):
            if pm_vals is None and spot_vals is None:
                return None
            if pm_vals is None:
                return np.asarray(spot_vals, dtype=np.float64)
            if spot_vals is None:
                return np.asarray(pm_vals, dtype=np.float64)
            pm_vals = np.asarray(pm_vals, dtype=np.float64)
            spot_vals = np.asarray(spot_vals, dtype=np.float64)
            out = np.full(pm_vals.shape[0], np.nan, dtype=np.float64)
            pm_finite = np.isfinite(pm_vals)
            spot_finite = np.isfinite(spot_vals)
            both = pm_finite & spot_finite
            if op == "product":
                out[both] = pm_vals[both] * spot_vals[both]
            else:
                out[both] = np.minimum(pm_vals[both], spot_vals[both])
            out[pm_finite & ~spot_finite] = pm_vals[pm_finite & ~spot_finite]
            out[~pm_finite & spot_finite] = spot_vals[~pm_finite & spot_finite]
            return out

        if mode in {"auto", "pm_on_cell", "pm"} and conf_pm is not None:
            return conf_pm, "pm_on_cell_assigned_prob"
        if mode in {"level1", "spot_contribution", "contribution"} and conf_broad is not None:
            return conf_broad, f"spot_{broad_key}_assigned_contribution"
        if mode in {"spot_confidence", "confidence"} and conf_spot is not None:
            return conf_spot, "spot_confidence"
        if mode in {"pm_x_spot", "product"}:
            combined = _combine_pm_and_spot(conf_pm, conf_spot, "product")
            if combined is not None:
                return combined, "pm_on_cell_x_spot_confidence"
        if mode in {"pm_min_spot", "min"}:
            combined = _combine_pm_and_spot(conf_pm, conf_spot, "min")
            if combined is not None:
                return combined, "min_pm_on_cell_spot_confidence"

        if conf_pm is not None:
            return conf_pm, "pm_on_cell_assigned_prob"
        if conf_broad is not None:
            return conf_broad, f"spot_{broad_key}_assigned_contribution"
        if conf_spot is not None:
            return conf_spot, "spot_confidence"
        return None, None

    def _get_graphagg_confidence(self):
        if self._graphagg_confidence_cache is None:
            conf, source = self._infer_graphagg_assignment_confidence()
            self._graphagg_confidence_cache = None if conf is None else np.asarray(conf, dtype=np.float64)
            self._graphagg_confidence_source = source
            if self._graphagg_confidence_cache is not None:
                self.svc_obs["graphagg_confidence"] = self._graphagg_confidence_cache
        return self._graphagg_confidence_cache, self._graphagg_confidence_source

    def _build_graphagg_conf_weighted_alpha_vector(self):
        """Low-confidence cells receive larger smoothing weights."""
        if self._graphagg_alpha_weight_cache is not None:
            return self._graphagg_alpha_weight_cache
        if not bool(getattr(self.config, "rec_graph_agg_conf_weighted_alpha", False)):
            return None

        conf, source = self._get_graphagg_confidence()
        if conf is None:
            self.logger.warning(
                "SR graph aggregation confidence-weighted alpha requested but no confidence source is available; "
                "falling back to fixed reconstruct alpha"
            )
            return None

        conf = np.asarray(conf, dtype=np.float64)
        finite = np.isfinite(conf)
        if not np.any(finite):
            self.logger.warning(
                "SR graph aggregation confidence-weighted alpha requested but confidence values are all non-finite; "
                "falling back to fixed reconstruct alpha"
            )
            return None

        alpha_min = float(getattr(self.config, "rec_graph_agg_conf_alpha_min", 0.0))
        alpha_max_cfg = float(getattr(self.config, "rec_graph_agg_conf_alpha_max", -1.0))
        alpha_max = float(self.config.rec_alpha) if alpha_max_cfg < 0 else alpha_max_cfg
        if alpha_max < alpha_min:
            alpha_min, alpha_max = alpha_max, alpha_min
        alpha_min = float(np.clip(alpha_min, 0.0, 1.0))
        alpha_max = float(np.clip(alpha_max, 0.0, 1.0))
        power = float(getattr(self.config, "rec_graph_agg_conf_alpha_power", 1.0))
        if power <= 0:
            power = 1.0

        out = np.full(conf.shape[0], float(self.config.rec_alpha), dtype=np.float64)
        cmin = float(np.min(conf[finite]))
        cmax = float(np.max(conf[finite]))
        if cmax - cmin <= 1e-12:
            out[finite] = float(np.clip(self.config.rec_alpha, alpha_min, alpha_max))
        else:
            conf_norm = (conf[finite] - cmin) / (cmax - cmin)
            strength = np.power(1.0 - conf_norm, power)
            out[finite] = alpha_min + (alpha_max - alpha_min) * strength
        out = np.clip(out, 0.0, 1.0)
        self._graphagg_alpha_weight_cache = out
        self.svc_obs["graphagg_alpha_weight"] = out

        q_points = [0.0, 0.1, 0.5, 0.9, 1.0]
        alpha_q = np.quantile(out[finite], q_points)
        self.logger.info(
            "SR graph aggregation confidence-weighted alpha enabled: source=%s, alpha_min=%.4f, alpha_max=%.4f, power=%.3f, conf_range=[%.6f, %.6f], alpha_quantiles=%s",
            source,
            alpha_min,
            alpha_max,
            power,
            cmin,
            cmax,
            {str(q): float(v) for q, v in zip(q_points, alpha_q)},
        )
        return out

    def _build_graphagg_low_conf_mask(self):
        conf, source = self._get_graphagg_confidence()
        if conf is None:
            self.logger.warning(
                "SR graph aggregation low-confidence-only requested but no confidence source is available; "
                "falling back to smoothing all cells"
            )
            return None

        finite = np.isfinite(conf)
        if not np.any(finite):
            self.logger.warning(
                "SR graph aggregation low-confidence-only requested but confidence values are all non-finite; "
                "falling back to smoothing all cells"
            )
            return None

        q = float(getattr(self.config, "rec_graph_agg_low_conf_quantile", 0.2))
        q = min(max(q, 0.0), 1.0)
        thr = float(np.quantile(conf[finite], q))
        mask = np.zeros(conf.shape[0], dtype=bool)
        mask[finite] = conf[finite] <= thr

        selected = int(mask.sum())
        total = int(mask.shape[0])
        q_points = [0.0, 0.1, 0.2, 0.5, 0.8, 0.9, 1.0]
        q_vals = np.quantile(conf[finite], q_points)
        self.logger.info(
            "SR graph aggregation low-confidence-only enabled: source=%s, quantile=%.3f, threshold=%.6f, selected=%d/%d (%.2f%%), finite_conf=%d, conf_quantiles=%s",
            source,
            q,
            thr,
            selected,
            total,
            100.0 * selected / max(total, 1),
            int(finite.sum()),
            {str(qp): float(qv) for qp, qv in zip(q_points, q_vals)},
        )
        if selected <= 0:
            self.logger.warning("No low-confidence cells selected for graph aggregation; smoothing all cells instead")
            return None
        return mask

    def _build_graphagg_high_conf_anchor_mask(self):
        conf, source = self._get_graphagg_confidence()
        if conf is None:
            self.logger.warning(
                "SR graph aggregation anchor-only requested but no confidence source is available; "
                "falling back to unrestricted donors"
            )
            return None
        finite = np.isfinite(conf)
        if not np.any(finite):
            self.logger.warning(
                "SR graph aggregation anchor-only requested but confidence values are all non-finite; "
                "falling back to unrestricted donors"
            )
            return None

        q = float(getattr(self.config, "rec_graph_agg_anchor_high_conf_quantile", 0.8))
        q = min(max(q, 0.0), 1.0)
        thr = float(np.quantile(conf[finite], q))
        mask = np.zeros(conf.shape[0], dtype=bool)
        mask[finite] = conf[finite] >= thr
        selected = int(mask.sum())
        total = int(mask.shape[0])
        self.svc_obs["graphagg_anchor"] = mask
        self.logger.info(
            "SR graph aggregation anchor-only enabled: source=%s, high_quantile=%.3f, threshold=%.6f, anchors=%d/%d (%.2f%%)",
            source,
            q,
            thr,
            selected,
            total,
            100.0 * selected / max(total, 1),
        )
        if selected <= 0:
            self.logger.warning("No high-confidence anchors selected; falling back to unrestricted donors")
            return None
        return mask

    def _apply_graph_aggregation(self, SVC_X, target_mask=None, anchor_mask=None):
        """Apply optional OT-based graph aggregation to SR virtual cells.

        The implementation mirrors the application-time SR graph smoothing but
        is kept optional in benchmark mode so raw vs. graph-aggregated metrics
        can be compared under the same noisy input.
        """
        n_cells = SVC_X.shape[0]
        if n_cells <= 1:
            self.logger.info("Skipping graph aggregation due to small cell count")
            return SVC_X.copy()
        if target_mask is not None:
            target_mask = np.asarray(target_mask, dtype=bool)
            if target_mask.shape[0] != n_cells:
                raise ValueError("target_mask length must match n_cells in _apply_graph_aggregation")
            self.logger.info(
                "Selective SR graph aggregation: updating %d/%d cells (%.2f%%)",
                int(target_mask.sum()),
                int(n_cells),
                100.0 * float(target_mask.sum()) / max(float(n_cells), 1.0),
            )
        if anchor_mask is not None:
            anchor_mask = np.asarray(anchor_mask, dtype=bool)
            if anchor_mask.shape[0] != n_cells:
                raise ValueError("anchor_mask length must match n_cells in _apply_graph_aggregation")
            self.logger.info(
                "Anchor-only donor mode: using %d/%d high-confidence cells as donors (%.2f%%)",
                int(anchor_mask.sum()),
                int(n_cells),
                100.0 * float(anchor_mask.sum()) / max(float(n_cells), 1.0),
            )
            if target_mask is not None:
                overlap = int(np.sum(anchor_mask & target_mask))
                self.logger.info("Anchor/target overlap count: %d", overlap)
        alpha_weight_global = self._build_graphagg_conf_weighted_alpha_vector()
        if alpha_weight_global is not None:
            self.logger.info(
                "Confidence-weighted alpha active for SR graph aggregation: full-graph smoothing strength varies by cell confidence"
            )

        cell_types = self.svc_obs["cell_type"].astype(str).to_numpy()
        unique_types = np.unique(cell_types)
        spatial_xy = self.svc_obs[["x", "y"]].to_numpy(dtype=np.float64)
        SVC_X_smoothed = SVC_X.copy()
        n_updated = 0

        for cell_type in unique_types:
            idx = np.where(cell_types == cell_type)[0]
            target_local = None if target_mask is None else target_mask[idx]
            anchor_local = None if anchor_mask is None else anchor_mask[idx]
            alpha_local = None if alpha_weight_global is None else alpha_weight_global[idx]
            if target_local is not None and not np.any(target_local):
                continue
            if anchor_local is not None and not np.any(anchor_local):
                if target_local is None or np.any(target_local):
                    self.logger.info(f"cell type: {cell_type}, no selected anchors, skip anchor-only graph aggregation")
                continue
            if idx.size < 50:
                self.logger.info(f"cell type: {cell_type}, has too few cells, skip graph aggregation")
                continue

            adata_cell = sc.AnnData(SVC_X[idx].copy())
            adata_cell.obsm["spatial"] = spatial_xy[idx]
            adjacent_matrix = get_adjacency_graph(
                adata_cell,
                data_type="sc_app",
                neighbors_method=self.config.rec_graph_method,
                alpha=self.config.rec_graph_alpha,
                gene_neighbor_num=self.config.rec_graph_exp_neighbor_num,
                spatial_neighbor_num=self.config.rec_graph_spatial_neighbor_num,
            )

            n_ct = idx.size
            K = min(int(self.config.rec_graph_n_neighbors), n_ct)
            if K <= 0:
                continue

            similarity_matrix = np.zeros((n_ct, K), dtype=np.float64)
            neighbor_idx_matrix = np.repeat(np.arange(n_ct, dtype=np.int32)[:, None], K, axis=1)
            valid_neighbor_mask = np.zeros((n_ct, K), dtype=bool)
            neighbor_margin_expr = np.zeros(K, dtype=np.float64)
            cell_gene_mean = np.asarray(adata_cell.X).mean(axis=1)

            for i in range(n_ct):
                row = adjacent_matrix.getrow(i)
                if row.nnz == 0:
                    continue
                data = row.data
                ridx = row.indices
                positive = data > 0
                data = data[positive]
                ridx = ridx[positive]
                use_anchor_donors = anchor_local is not None and (target_local is None or bool(target_local[i]))
                if use_anchor_donors:
                    donor_mask = anchor_local[ridx]
                    if np.any(donor_mask):
                        data = data[donor_mask]
                        ridx = ridx[donor_mask]
                    else:
                        # No anchor donor in the local graph neighborhood. Leave
                        # the row empty so the downstream kernel falls back to
                        # self-preserving behavior for this target.
                        continue
                take = min(K, data.size)
                if take <= 0:
                    continue
                if data.size > take:
                    top_idx = np.argpartition(-data, kth=take - 1)[:take]
                    top_idx = top_idx[np.argsort(-data[top_idx])]
                else:
                    top_idx = np.argsort(-data)
                sel_idx = ridx[top_idx]
                sel_data = data[top_idx]
                similarity_matrix[i, :take] = sel_data
                neighbor_idx_matrix[i, :take] = sel_idx.astype(np.int32)
                valid_neighbor_mask[i, :take] = True
                neighbor_margin_expr[:take] += cell_gene_mean[sel_idx]

            mu = np.ravel(np.asarray(adata_cell.X).sum(axis=1))
            nu = neighbor_margin_expr
            if not (np.any(mu) and np.any(nu)):
                self.logger.info(f"cell type: {cell_type}, skip graph aggregation due to empty marginals")
                continue

            source_idx, target_idx, active_support = stabilize_local_ot_support(
                nu,
                mu,
                valid_neighbor_mask.T,
            )
            if source_idx.size == 0 or target_idx.size == 0:
                self.logger.info(
                    f"cell type: {cell_type}, skip graph aggregation due to empty active support"
                )
                continue
            stable_support = np.zeros(valid_neighbor_mask.T.shape, dtype=bool)
            stable_support[np.ix_(source_idx, target_idx)] = active_support
            valid_neighbor_mask = stable_support.T
            distance_matrix = similarity_to_distance(
                similarity_matrix,
                valid_neighbor_mask,
            )
            group_assignment = subset_virtual_assignment(
                self._projected_assignment,
                self.svc_obs.iloc[idx]["cell_id"],
            )
            distance_matrix = condition_virtual_cell_ot_cost(
                distance_matrix,
                assignment=group_assignment,
                neighbor_indices=neighbor_idx_matrix,
                valid_support_mask=valid_neighbor_mask,
                strength=self._conditioning_strength,
            )
            distance_matrix[~valid_neighbor_mask] = np.inf
            T_transform = solve_local_ot(
                nu,
                mu,
                distance_matrix.T,
                method=self.config.rec_ot_method,
                pot_reg=self.config.rec_pot_reg,
                pot_reg_m=self.config.rec_pot_reg_m,
                pot_reg_type=self.config.rec_pot_reg_type,
                pot_verbose=False,
                pot_num_iter_max=5000,
                reference_measure=None,
                valid_support_mask=valid_neighbor_mask.T,
            )
            adata_cell = self.graph_aggregate.run(
                adata=adata_cell,
                neighbor_idx_matrix=neighbor_idx_matrix,
                coupling_matrix=T_transform,
                alpha_override=alpha_local,
                valid_neighbor_mask=valid_neighbor_mask,
            )
            callback = getattr(
                self.config,
                "local_refinement_applied_callback",
                None,
            )
            if callback is not None:
                callback()
            self._refinement_applied = True
            smoothed_block = np.asarray(adata_cell.X)
            if target_local is None:
                SVC_X_smoothed[idx] = smoothed_block
                n_updated += int(idx.size)
            else:
                SVC_X_smoothed[idx[target_local]] = smoothed_block[target_local]
                n_updated += int(np.sum(target_local))

        if target_mask is None:
            self.logger.info("SR graph aggregation updated all eligible cells")
        else:
            self.logger.info(
                "SR graph aggregation selectively updated %d cells and kept %d high-confidence cells as raw",
                int(n_updated),
                int(n_cells - n_updated),
            )
        return SVC_X_smoothed
