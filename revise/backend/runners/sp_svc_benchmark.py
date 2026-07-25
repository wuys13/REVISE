import numpy as np
import scanpy as sc
import scipy
from tqdm import tqdm

from revise.backend.runners.benchmark_svc import BenchmarkSVC
from revise.backend.kernels import SegEvaluateKernel as SegEvaluate
from revise.backend.ops.distance import similarity_to_distance
from revise.backend.ops.local_ot import solve_local_ot, stabilize_local_ot_support
from revise.backend.runners.sp_svc_assignment_guidance import (
    assignment_categories,
    guidance_mode,
    prepare_assignment_guidance,
    record_guidance_terminal,
    record_not_applicable,
)
from revise.backend.ops.topology import get_adjacency_graph


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
            guidance_problem_key = f"sp_svc_benchmark:{cell_type}"
            svc_adata_cell_type = self.st_adata[self.st_adata.obs[self.config.cell_type_col] == cell_type]
            svc_replace_adata = svc_adata_cell_type[~svc_adata_cell_type.obs["no_effect"]]
            svc_candidate_adata = svc_adata_cell_type[svc_adata_cell_type.obs["no_effect"]]
            if svc_replace_adata.shape[0] < 50:
                self.logger.info(f"cell type: {cell_type} has too few spots, skip OT smoothing")
                if getattr(
                    self.config,
                    "assignment_guidance_callback",
                    None,
                ) is not None:
                    record_not_applicable(
                        self.config,
                        route="sp_svc_benchmark",
                        operator="replacement_ot",
                        problem_key=guidance_problem_key,
                        reason="insufficient_observations",
                    )
                svc_replace_adata.layers["ot_smooth"] = svc_replace_adata.X.copy()
                cell_type_adata_list.append(svc_replace_adata)
            elif svc_candidate_adata.shape[0] == 0:
                self.logger.info(f"cell type: {cell_type} has no candidate cells, skip OT smoothing")
                if getattr(
                    self.config,
                    "assignment_guidance_callback",
                    None,
                ) is not None:
                    record_not_applicable(
                        self.config,
                        route="sp_svc_benchmark",
                        operator="replacement_ot",
                        problem_key=guidance_problem_key,
                        reason="missing_donors",
                    )
                svc_replace_adata.layers["ot_smooth"] = svc_replace_adata.X.copy()
                cell_type_adata_list.append(svc_replace_adata)
            else:
                assignment_category_labels = []
                if guidance_mode(self.config) != "off":
                    assignment_category_labels = assignment_categories(
                        self.sc_ref_adata,
                        self.st_adata,
                        key=self.config.cell_type_col,
                    )
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

                similarity_matrix = np.zeros(
                    (n_recon, self.config.rec_graph_n_neighbors),
                    dtype=np.float64,
                )
                neighbor_idx_matrix = np.zeros((n_recon, self.config.rec_graph_n_neighbors), dtype=np.int32)
                valid_neighbor_mask = np.zeros(
                    (n_recon, self.config.rec_graph_n_neighbors),
                    dtype=bool,
                )

                nu_slots = np.zeros(self.config.rec_graph_n_neighbors, dtype=np.float64)
                cand_X_csr = svc_candidate_adata.X.tocsr()
                recon_X_csr = svc_replace_adata.X.tocsr()

                for i in tqdm(range(n_recon), desc="TopK expression"):
                    row = cross_adj.getrow(i).toarray().ravel()
                    positive_idx = np.flatnonzero(row > 0)
                    if positive_idx.size == 0:
                        continue

                    take = min(self.config.rec_graph_n_neighbors, positive_idx.size)
                    if positive_idx.size > take:
                        local_idx = np.argpartition(
                            -row[positive_idx],
                            kth=take - 1,
                        )[:take]
                        idx = positive_idx[local_idx]
                    else:
                        idx = positive_idx
                    idx = idx[np.argsort(-row[idx])]

                    similarity_matrix[i, :take] = row[idx].copy()
                    neighbor_idx_matrix[i, :take] = idx.astype(np.int32)
                    valid_neighbor_mask[i, :take] = True

                    slot_expr = cand_X_csr[idx].toarray().mean(axis=1).ravel()
                    nu_slots[:take] += slot_expr

                mu = np.ravel(recon_X_csr.mean(axis=1))
                nu = nu_slots

                source_idx, target_idx, active_support = stabilize_local_ot_support(
                    nu,
                    mu,
                    valid_neighbor_mask.T,
                )
                if source_idx.size == 0 or target_idx.size == 0:
                    self.logger.info(
                        f"cell type: {cell_type}, skip OT smoothing due to empty active support"
                    )
                    if getattr(
                        self.config,
                        "assignment_guidance_callback",
                        None,
                    ) is not None:
                        record_not_applicable(
                            self.config,
                            route="sp_svc_benchmark",
                            operator="replacement_ot",
                            problem_key=guidance_problem_key,
                            reason="empty_active_support",
                        )
                    svc_replace_adata.layers["ot_smooth"] = svc_replace_adata.X.copy()
                    cell_type_adata_list.append(svc_replace_adata)
                    continue
                stable_support = np.zeros(valid_neighbor_mask.T.shape, dtype=bool)
                stable_support[np.ix_(source_idx, target_idx)] = active_support
                valid_neighbor_mask = stable_support.T
                distance_matrix = similarity_to_distance(
                    similarity_matrix,
                    valid_neighbor_mask,
                )
                (
                    distance_matrix,
                    reference_measure,
                    guidance_attempted,
                ) = prepare_assignment_guidance(
                    self.config,
                    route="sp_svc_benchmark",
                    operator="replacement_ot",
                    problem_key=guidance_problem_key,
                    left_adata=svc_replace_adata,
                    right_adata=svc_candidate_adata,
                    category_labels=assignment_category_labels,
                    support=neighbor_idx_matrix,
                    distance_matrix=distance_matrix,
                    source_mass=nu,
                    target_mass=mu,
                )
                distance_matrix[~valid_neighbor_mask] = np.inf
                try:
                    T_transform = solve_local_ot(
                        nu,
                        mu,
                        distance_matrix.T,
                        method=self.config.rec_ot_method,
                        pot_reg=self.config.rec_pot_reg,
                        pot_reg_m=self.config.rec_pot_reg_m,
                        pot_reg_type=self.config.rec_pot_reg_type,
                        pot_verbose=True,
                        pot_num_iter_max=5000,
                        reference_measure=reference_measure,
                        valid_support_mask=valid_neighbor_mask.T,
                        event_callback=getattr(self.config, "ot_event_callback", None),
                    )
                except Exception:
                    record_guidance_terminal(
                        self.config,
                        problem_key=guidance_problem_key,
                        attempted=guidance_attempted,
                        outcome="failed",
                        reason="solver_failed",
                    )
                    raise
                try:
                    alpha = float(self.config.rec_alpha)
                    smoothed = scipy.sparse.lil_matrix(
                        recon_X_csr.shape,
                        dtype=recon_X_csr.dtype,
                    )

                    for i in range(n_recon):
                        idx = neighbor_idx_matrix[i]
                        valid_mask = valid_neighbor_mask[i]
                        if not np.any(valid_mask):
                            smoothed[i] = recon_X_csr.getrow(i)
                            continue

                        idx = idx[valid_mask]
                        w = T_transform[valid_mask, i]
                        w_sum = w.sum()
                        if w_sum <= 0:
                            smoothed[i] = recon_X_csr.getrow(i)
                            continue
                        w = w / w_sum

                        neigh_expr = cand_X_csr[idx]
                        weighted = neigh_expr.T @ w
                        weighted = np.asarray(weighted).ravel()

                        base = recon_X_csr.getrow(i).toarray().ravel()
                        new_vec = (1.0 - alpha) * base + alpha * weighted

                        smoothed[i] = scipy.sparse.csr_matrix(new_vec)
                    svc_replace_adata.layers["ot_smooth"] = (
                        smoothed.tocsr().copy()
                    )
                    cell_type_adata_list.append(svc_replace_adata)
                except Exception:
                    record_guidance_terminal(
                        self.config,
                        problem_key=guidance_problem_key,
                        attempted=guidance_attempted,
                        outcome="failed",
                        reason="expression_update_failed",
                    )
                    raise
                record_guidance_terminal(
                    self.config,
                    problem_key=guidance_problem_key,
                    attempted=guidance_attempted,
                    outcome="applied",
                )

        svc_recon_adata = sc.concat(cell_type_adata_list)
        svc_recon_adata.X = svc_recon_adata.layers["ot_smooth"].copy()

        svc_no_effect = self.st_adata[self.st_adata.obs["no_effect"]]
        self.svc["sp_svc"] = sc.concat([svc_recon_adata, svc_no_effect])
