import numpy as np
import ot
import scanpy as sc
import scipy
from tqdm import tqdm

from revise.benchmark.benchmark_svc import BenchmarkSVC
from revise.methods.seg_evaluate import SegEvaluate
from revise.tools.topology import get_adjacency_graph


class SpSVC(BenchmarkSVC):
    """
    sp-SVC class for benchmark CFs: segmentation/bin2cell.
    
    This class reconstructs single-cell resolution expression profiles
    from spatial transcriptomics data, with special handling for segmentation
    errors (diminishing, expanding, unchanged cells).
    """
    def __init__(self, st_adata, sc_ref_adata, config, real_st_adata, logger):
        super().__init__(st_adata, sc_ref_adata, config, real_st_adata, logger)
        self._adata_validate()
        self._adata_processing()
        self.seg_evaluate = SegEvaluate(self.config, self.logger)
        self.svc = {}

    def local_refinement(self):
        """Reconstruct expression profiles with segmentation-aware smoothing.

        1. Evaluate segmentation errors and flag cells that need correction.
        2. Split each cell type into ``replace`` and ``candidate`` groups.
        3. Use optimal transport between the two groups to obtain smoothed
           expressions for the ``replace`` cells.
        4. Merge corrected and unchanged cells to form ``self.svc["sp_svc"]``.
        """
        if "seg_error" in self.st_adata.obs.columns:
            self.st_adata = self.seg_evaluate.run(self.st_adata, self.logger)
        else:
            self.logger.warning("No 'seg_error' not in st_adata.obs, evaluation skip.")
        cell_type_adata_list = []
        for cell_type in tqdm(self.st_adata.obs[self.config.cell_type_col].unique().tolist(), desc="Reconstruting"):
            svc_adata_cell_type = self.st_adata[self.st_adata.obs[self.config.cell_type_col] == cell_type]
            svc_replace_adata = svc_adata_cell_type[~svc_adata_cell_type.obs["no_effect"]]
            svc_candidate_adata = svc_adata_cell_type[svc_adata_cell_type.obs["no_effect"]]
            if svc_replace_adata.shape[0] < 50:
                self.logger.info(f"cell type: {cell_type} has too few spots, skip OT smoothing")
                svc_replace_adata.layers["ot_smooth"] = svc_replace_adata.X.copy()
                cell_type_adata_list.append(svc_replace_adata)
            else:
                # Build adjacency on ordered data to align replace and candidate partitions
                svc_ordered = sc.concat([svc_replace_adata, svc_candidate_adata])
                adjacent_matrix_all = get_adjacency_graph(
                    svc_ordered,
                    data_type="sc",
                    neighbors_method=self.config.rec_graph_method,
                    alpha=self.config.rec_graph_alpha,
                    gene_neighbor_num=self.config.rec_graph_exp_neighbor_num,
                    spatial_neighbor_num=self.config.rec_graph_spatial_neighbor_num,
                )

                n_recon = svc_replace_adata.shape[0]
                n_cand = svc_candidate_adata.shape[0]
                cross_adj = adjacent_matrix_all[:n_recon, n_recon:n_recon + n_cand].tocsr()
                svc_replace_adata.obsm["cross_connectivities"] = cross_adj

                cost_matrix = np.zeros((n_recon, self.config.rec_graph_n_neighbors), dtype=svc_replace_adata.X.dtype)
                neighbor_idx_matrix = np.zeros((n_recon, self.config.rec_graph_n_neighbors), dtype=np.int32)

                nu_slots = np.zeros(self.config.rec_graph_n_neighbors, dtype=svc_replace_adata.X.dtype)
                cand_X_csr = svc_candidate_adata.X.tocsr()
                recon_X_csr = svc_replace_adata.X.tocsr()

                for i in tqdm(range(n_recon), desc="TopK expression"):
                    row = cross_adj.getrow(i).toarray().ravel()
                    if np.count_nonzero(row) == 0:
                        continue

                    take = min(self.config.rec_graph_n_neighbors, row.size)
                    idx = np.argpartition(-row, kth=take - 1)[:take]
                    idx = idx[np.argsort(-row[idx])]

                    cost_matrix[i, :take] = row[idx].copy()
                    neighbor_idx_matrix[i, :take] = idx.astype(np.int32)

                    slot_expr = cand_X_csr[idx].toarray().mean(axis=1).ravel()
                    nu_slots[:take] += slot_expr

                mu = np.ravel(recon_X_csr.mean(axis=1))
                nu = nu_slots

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
                    verbose=True,
                    numItermax=5000
                )
                alpha = float(self.config.rec_alpha)
                smoothed = scipy.sparse.lil_matrix(recon_X_csr.shape, dtype=recon_X_csr.dtype)

                for i in range(n_recon):
                    idx = neighbor_idx_matrix[i]
                    valid_mask = cost_matrix[i] > 0
                    if not np.any(valid_mask):
                        smoothed[i] = recon_X_csr.getrow(i)
                        continue

                    idx = idx[valid_mask]
                    w = T_transform[valid_mask, i]
                    w_sum = w.sum()
                    if w_sum > 0:
                        w = w / w_sum

                    neigh_expr = cand_X_csr[idx]
                    weighted = (neigh_expr.T @ w)
                    weighted = np.asarray(weighted).ravel()

                    base = recon_X_csr.getrow(i).toarray().ravel()
                    new_vec = (1.0 - alpha) * base + alpha * weighted

                    smoothed[i] = scipy.sparse.csr_matrix(new_vec)
                svc_replace_adata.layers["ot_smooth"] = smoothed.tocsr().copy()
                cell_type_adata_list.append(svc_replace_adata)

        svc_recon_adata = sc.concat(cell_type_adata_list)
        svc_recon_adata.X = svc_recon_adata.layers["ot_smooth"].copy()

        svc_no_effect = self.st_adata[self.st_adata.obs["no_effect"]]
        self.svc["sp_svc"] = sc.concat([svc_recon_adata, svc_no_effect])
