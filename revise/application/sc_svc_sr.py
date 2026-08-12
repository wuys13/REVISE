import numpy as np
import ot
import scanpy as sc

from revise.application.application_svc import ApplicationSVC
from revise.methods.graph_aggregate import GraphAggregate
from revise.methods.spot_sr import SpotSr
from revise.tools.meta import construct_sc_ref
from revise.tools.meta import get_sc_obs
from revise.tools.topology import get_adjacency_graph


class ScSVCSr(ApplicationSVC):
    """
    sc-SVC super-resolution for application usage.

    This class reconstructs single-cell resolution expression profiles
    from spatial transcriptomics data by redistributing spot-level
    expressions to virtual cells using cell type contributions.
    """
    def __init__(self, st_adata, sc_ref_adata, config, logger):
        super().__init__(st_adata, sc_ref_adata, config, None, logger)
        self._adata_validate()
        self._adata_validate_dec()
        self._adata_processing()
        self.svc_obs = self._get_svc_obs()
        self.spot_sr = SpotSr(self.config, self.logger)
        self.graph_aggregate = GraphAggregate(self.config, self.logger)
        self.svc = {}

    def _adata_validate_dec(self):
        assert "all_cells_in_spot" in self.st_adata.uns, "spot-sc mapping is not in st_adata.uns"

    def _get_svc_obs(self):
        svc_obs = get_sc_obs(self.st_adata.obs.index, self.st_adata.uns['all_cells_in_spot'], self.st_adata.obsm["spatial"])
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
        st_adata_common = self.st_adata[:, overlap_genes]
        sc.pp.normalize_total(st_adata_common, target_sum=1e4)
        cell_contributions = st_adata_common.obsm["Level1"].values if hasattr(
            st_adata_common.obsm["Level1"], 'values') else st_adata_common.obsm["Level1"]

        sc.pp.normalize_total(self.st_adata, target_sum=1e4)
        sc.pp.normalize_total(self.sc_ref_adata, target_sum=1e4)
        self.spot_sr.run(self)
        key_type = "clusters"
        if key_type not in self.sc_ref_adata.obs.columns:
            self.sc_ref_adata.obs[key_type] = self.sc_ref_adata.obs["Level1"].astype(str)
        type_list = sorted(list(self.sc_ref_adata.obs[key_type].unique().astype(str)))
        self.logger.info(f'There are {len(type_list)} cell types: {type_list}')
        sc_ref_all = construct_sc_ref(self.sc_ref_adata, key_type=key_type, type_list=type_list)
        sc_ref_all = sc_ref_all.loc[:, overlap_genes]
        type_list = sorted(list(sc_ref_all.index))
        spots = self.svc_obs['spot_name'].unique()

        spot_to_idx = {spot: idx for idx, spot in enumerate(spots)}
        self.logger.info("Using simple allocation method...")

        adata_spot = st_adata_common.copy()
        X = adata_spot.X if type(adata_spot.X) is np.ndarray else adata_spot.X.toarray()

        # cell_contributions shape: (13863, 14), sc_ref_all shape: (14, 407), Y shape: (13863, 407, 14)
        Y = cell_contributions[:, np.newaxis, :] * sc_ref_all.values.T
        Y = Y / (np.sum(Y, axis=2, keepdims=True) + 1e-10)
        Y = Y * X[:, :, np.newaxis]

        spot_indices = np.array([spot_to_idx[spot] for spot in self.svc_obs['spot_name']])
        type_indices = np.array([type_list.index(t) for t in self.svc_obs['cell_type']])

        # SVC_X shape: (14060, 407)
        SVC_X = Y[spot_indices, :, type_indices]
        self.logger.info("Extracted SVC expressions using simple allocation method")

        n_cells = SVC_X.shape[0]
        if n_cells > 1:
            self.logger.info("Applying OT-based neighbor enhancement among single cells")
            cell_types = self.svc_obs["cell_type"].astype(str).to_numpy()
            unique_types = np.unique(cell_types)
            spatial = self.svc_obs[["x", "y"]].to_numpy(dtype=np.float64)
            SVC_X_smoothed = SVC_X.copy()

            for cell_type in unique_types:
                idx = np.where(cell_types == cell_type)[0]
                if idx.size < 2:
                    self.logger.info(f"cell type: {cell_type}, has too few cells, skip OT smoothing")
                    SVC_X_smoothed[idx] = SVC_X[idx]
                    continue
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

                cost_matrix = np.zeros((n_ct, K), dtype=SVC_X.dtype)
                neighbor_idx_matrix = np.repeat(np.arange(n_ct, dtype=np.int32)[:, None], K, axis=1)
                neighbor_margin_expr = np.zeros(K, dtype=SVC_X.dtype)
                cell_gene_mean = np.asarray(adata_cell.X).mean(axis=1)

                for i in range(n_ct):
                    row = adjacent_matrix.getrow(i)
                    if row.nnz == 0:
                        continue
                    data = row.data
                    ridx = row.indices
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
                    cost_matrix[i, :take] = sel_data
                    neighbor_idx_matrix[i, :take] = sel_idx.astype(np.int32)
                    neighbor_margin_expr[:take] += cell_gene_mean[sel_idx]

                mu = np.ravel(np.asarray(adata_cell.X).sum(axis=1))
                nu = neighbor_margin_expr
                if np.any(mu) and np.any(nu):
                    cm_max = float(cost_matrix.max())
                    if cm_max <= 0:
                        cm_max = 1.0
                    T_transform = ot.unbalanced.sinkhorn_unbalanced(
                        nu,
                        mu,
                        cost_matrix.T / cm_max,
                        reg=self.config.rec_pot_reg,
                        reg_m=self.config.rec_pot_reg_m,
                        reg_type=self.config.rec_pot_reg_type,
                        verbose=False,
                        numItermax=5000
                    )
                    adata_cell = self.graph_aggregate.run(
                        adata=adata_cell,
                        neighbor_idx_matrix=neighbor_idx_matrix,
                        coupling_matrix=T_transform
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

        svc_adata = sc.AnnData(SVC_X)
        svc_adata.var_names = st_adata_common.var_names
        svc_adata.obs = self.svc_obs.copy()
        svc_adata.obs.set_index("cell_id", inplace=True)
        self.svc["sc_svc_dec"] = svc_adata
