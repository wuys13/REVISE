import numpy as np
import pandas as pd
import scanpy as sc

from revise.backend.runners.application_svc import ApplicationSVC
from revise.backend.kernels import GraphAggregateKernel as GraphAggregate
from revise.backend.kernels import SpotSrKernel as SpotSr
from revise.backend.ops.distance import similarity_to_distance
from revise.backend.ops.meta import construct_sc_ref
from revise.backend.ops.meta import get_sc_obs
from revise.backend.ops.local_ot import solve_local_ot, stabilize_local_ot_support
from revise.backend.ops.sr_allocation import (
    condition_virtual_cell_ot_cost,
    mandatory_reference_allocation,
    project_spot_assignment_to_virtual_cells,
    record_mandatory_allocation,
    spot_global_assignment,
    subset_virtual_assignment,
)
from revise.backend.ops.topology import get_adjacency_graph
from revise.utils.spot_sr_input import CELL_LOCATIONS_KEY


class ScSVCSr(ApplicationSVC):
    """
    sc-SVC super-resolution for application usage.

    This class reconstructs single-cell resolution expression profiles
    from spatial transcriptomics data by redistributing spot-level
    expressions to virtual cells using cell type contributions.
    """
    def __init__(self, st_adata, sc_ref_adata, config, logger):
        super().__init__(st_adata, sc_ref_adata, config, None, logger)
        self._adata_validate_dec()
        self._adata_processing()
        self._adata_validate()
        self.svc_obs = self._get_svc_obs()
        self.spot_sr = SpotSr(self.config, self.logger)
        self.graph_aggregate = GraphAggregate(self.config, self.logger)
        self.svc = {}

    def _adata_validate_dec(self):
        if "all_cells_in_spot" not in self.st_adata.uns:
            raise KeyError(
                "Invalid reconstruction input: "
                f"st_file_path={getattr(self.config, 'st_file_path', '<unavailable>')}; "
                "field=uns['all_cells_in_spot']; expected=present; actual=missing"
            )

    def _get_svc_obs(self):
        svc_obs = get_sc_obs(
            self.st_adata.obs.index,
            self.st_adata.uns["all_cells_in_spot"],
            self.st_adata.obsm["spatial"],
            cell_locations=self.st_adata.uns.get(CELL_LOCATIONS_KEY),
        )
        svc_obs["cell_id"] = svc_obs["cell_id"].astype(str)
        # Keep SpotSr logging stable when ground truth is unavailable.
        svc_obs["true_cell_type"] = "Unknown"
        return svc_obs
    
    # TODO: out of spot imputation
    def local_refinement(self, *args, **kwargs):
        """Reconstruct single-cell expression profiles from spot-level data.

        1. Assigns cell types to each virtual cell using SpotSr
        2. Constructs cell type reference profiles
        3. Calculates gene expression for each cell based on spot contributions
        4. Normalizes expressions to 10,000 counts per cell

        The reconstructed data is stored in self.svc["sc_svc_dec"].
        """
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
        projected_assignment = project_spot_assignment_to_virtual_cells(
            assignment,
            self.svc_obs,
        )
        conditioning_strength = getattr(
            self.config,
            "local_refinement_strength",
            0.0,
        )
        conditioning_applied = False

        n_cells = SVC_X.shape[0]
        if n_cells > 1:
            self.logger.info("Applying OT-based neighbor enhancement among single cells")
            cell_types = self.svc_obs["cell_type"].astype(str).to_numpy()
            unique_types = np.unique(cell_types)
            spatial = self.svc_obs[["x", "y"]].to_numpy(dtype=np.float64)
            SVC_X_smoothed = SVC_X.copy()

            for cell_type in unique_types:
                idx = np.where(cell_types == cell_type)[0]
                if idx.size < 50:
                    self.logger.info(f"cell type: {cell_type}, has too few cells, skip OT smoothing")
                    SVC_X_smoothed[idx] = SVC_X[idx]
                    continue

                adata_cell = sc.AnnData(SVC_X[idx])
                adata_cell.obsm["spatial"] = spatial[idx]
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
                if np.any(mu) and np.any(nu):
                    source_idx, target_idx, active_support = stabilize_local_ot_support(
                        nu,
                        mu,
                        valid_neighbor_mask.T,
                    )
                    if source_idx.size == 0 or target_idx.size == 0:
                        self.logger.info(
                            f"cell type: {cell_type}, skip OT smoothing due to empty active support"
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
                        projected_assignment,
                        self.svc_obs.iloc[idx]["cell_id"],
                    )
                    distance_matrix = condition_virtual_cell_ot_cost(
                        distance_matrix,
                        assignment=group_assignment,
                        neighbor_indices=neighbor_idx_matrix,
                        strength=conditioning_strength,
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
                    if conditioning_strength != 0:
                        callback = getattr(
                            self.config,
                            "local_refinement_applied_callback",
                            None,
                        )
                        if callback is not None:
                            callback()
                        conditioning_applied = True
                    adata_cell = self.graph_aggregate.run(
                        adata=adata_cell,
                        neighbor_idx_matrix=neighbor_idx_matrix,
                        coupling_matrix=T_transform,
                        valid_neighbor_mask=valid_neighbor_mask,
                    )
                    SVC_X_smoothed[idx] = np.asarray(adata_cell.X)
                else:
                    self.logger.info(f"cell type: {cell_type}, skip OT smoothing due to empty marginals")

            SVC_X = SVC_X_smoothed
        else:
            self.logger.info("Skipping OT enhancement due to small cell count")

        if getattr(self.config, "rec_match_spot_sum", False):
            self.logger.info("Rescaling single-cell expressions to match spot totals")
            current_sum = np.zeros_like(X, dtype=np.float64)
            np.add.at(current_sum, spot_indices, SVC_X)
            ratio = X / (current_sum + 1e-10)
            SVC_X = SVC_X * ratio[spot_indices]
        else:
            SVC_X = SVC_X / (np.sum(SVC_X, axis=1, keepdims=True) + 1e-10) * 1e4

        self.logger.info(f"Number of cells processed: {len(self.svc_obs)}")
        self.logger.info(f"Number of unique spots: {len(spots)}")
        self.logger.info(f"Shape of SVC_X: {SVC_X.shape}")

        svc_obs = self.svc_obs.copy()
        svc_obs["cell_id"] = svc_obs["cell_id"].astype(str)
        svc_obs.set_index("cell_id", inplace=True)
        svc_adata = sc.AnnData(SVC_X, obs=svc_obs)
        svc_adata.var_names = st_adata_common.var_names
        self.svc["sc_svc_dec"] = svc_adata
        return conditioning_applied
