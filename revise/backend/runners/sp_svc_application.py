import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from tqdm import tqdm

from revise.backend.runners.application_svc import ApplicationSVC
from revise.backend.kernels import GraphAggregateKernel as GraphAggregate
from revise.backend.kernels.ot import OTKernel, stabilize_local_ot_support
from revise.backend.ops.distance import similarity_to_distance
from revise.backend.runners.sp_svc_assignment import (
    condition_sp_local_ot_cost,
    global_assignment_from_adata,
)
from revise.backend.ops.shaver import trim_sp_adata
from revise.backend.ops.topology import get_adjacency_graph


def _dense_stable_topk(adjacent_matrix, row_index, n_neighbors):
    row = adjacent_matrix[row_index].toarray().ravel()
    positive_idx = np.flatnonzero(row > 0)
    if positive_idx.size == 0:
        return None, None
    take = min(n_neighbors, positive_idx.size)
    if positive_idx.size > take:
        local_idx = np.argpartition(-row[positive_idx], kth=take - 1)[:take]
        idx = positive_idx[local_idx]
    else:
        idx = positive_idx
    idx = idx[np.argsort(-row[idx])]
    return idx.astype(np.int32, copy=False), row[idx].copy()


def _sparse_exact_topk(adjacent_csr, row_index, n_neighbors):
    start = adjacent_csr.indptr[row_index]
    end = adjacent_csr.indptr[row_index + 1]
    row_indices = adjacent_csr.indices[start:end]
    row_data = adjacent_csr.data[start:end]

    if row_data.size == 0:
        return None, None, False

    positive_mask = row_data > 0
    if not np.all(positive_mask):
        row_indices = row_indices[positive_mask]
        row_data = row_data[positive_mask]
        if row_data.size == 0:
            return None, None, False

    take = min(n_neighbors, row_data.size)
    if take > 0:
        local_idx = np.argpartition(-row_data, kth=take - 1)[:take]
        selected_data = row_data[local_idx]
        # Tied selected weights can change neighbor order versus the dense
        # np.argpartition path, so keep the original implementation for them.
        boundary_value = selected_data.min()
        boundary_is_unique = (
            np.count_nonzero(row_data == boundary_value)
            == np.count_nonzero(selected_data == boundary_value)
        )
        if boundary_is_unique and np.unique(selected_data).size == selected_data.size:
            order = np.argsort(-selected_data)
            local_idx = local_idx[order]
            return (
                row_indices[local_idx].astype(np.int32, copy=False),
                row_data[local_idx].copy(),
                False,
            )

    idx, values = _dense_stable_topk(adjacent_csr, row_index, n_neighbors)
    return idx, values, True


def _compute_topk_expression(adjacent_matrix, expression_matrix, n_neighbors, dtype, progress=True):
    adjacent_csr = adjacent_matrix.tocsr()
    n_obs = adjacent_csr.shape[0]

    similarity_matrix = np.zeros((n_obs, n_neighbors), dtype=np.float64)
    neighbor_margin_expr = np.zeros(n_neighbors, dtype=np.float64)
    neighbor_idx_matrix = np.zeros((n_obs, n_neighbors), dtype=np.int32)
    valid_neighbor_mask = np.zeros((n_obs, n_neighbors), dtype=bool)

    if sparse.issparse(expression_matrix):
        expression_csr = expression_matrix.tocsr()
        row_expr_mean = np.asarray(expression_csr.mean(axis=1)).ravel()
    else:
        row_expr_mean = np.mean(np.asarray(expression_matrix), axis=1)
    row_expr_mean = row_expr_mean.astype(np.float64, copy=False)

    dense_fallback_rows = 0
    iterator = range(n_obs)
    if progress:
        iterator = tqdm(iterator, desc="TopK expression")
    for i in iterator:
        idx, values, used_dense_fallback = _sparse_exact_topk(adjacent_csr, i, n_neighbors)
        if idx is None:
            continue
        if used_dense_fallback:
            dense_fallback_rows += 1

        take = idx.size
        neighbor_margin_expr[:take] += row_expr_mean[idx]
        similarity_matrix[i, :take] = values.astype(np.float64, copy=False)
        neighbor_idx_matrix[i, :take] = idx
        valid_neighbor_mask[i, :take] = True

    return (
        similarity_matrix,
        neighbor_margin_expr,
        neighbor_idx_matrix,
        valid_neighbor_mask,
        dense_fallback_rows,
    )


def _validate_expression_unchanged(current, original, *, cell_type):
    current_is_sparse = sparse.issparse(current)
    original_is_sparse = sparse.issparse(original)
    if current_is_sparse != original_is_sparse:
        current_representation = "sparse" if current_is_sparse else "dense"
        original_representation = "sparse" if original_is_sparse else "dense"
        raise RuntimeError(
            "Expression invariant failed before aggregation: "
            f"cell_type={cell_type}; field=X; expected=exactly unchanged; "
            f"actual=current_representation={current_representation}, "
            f"original_representation={original_representation}"
        )

    if current.shape != original.shape:
        raise RuntimeError(
            "Expression invariant failed before aggregation: "
            f"cell_type={cell_type}; field=X; expected=exactly unchanged; "
            f"actual=shape={current.shape}, original_shape={original.shape}"
        )

    if current_is_sparse:
        mismatch_count = (current != original).nnz
    else:
        current_array = np.asarray(current)
        original_array = np.asarray(original)
        mismatch_count = (
            0
            if np.array_equal(current_array, original_array)
            else int(np.count_nonzero(current_array != original_array))
        )

    if mismatch_count:
        raise RuntimeError(
            "Expression invariant failed before aggregation: "
            f"cell_type={cell_type}; field=X; expected=exactly unchanged; "
            f"actual=mismatches={mismatch_count}"
        )


class SpSVC(ApplicationSVC):
    """
    sp-SVC class for application usage.
    
    This class reconstructs single-cell resolution expression profiles
    from spatial transcriptomics data using optimal transport-based
    graph aggregation for each cell type.
    """
    def __init__(self, st_adata, sc_ref_adata, config, logger):
        super().__init__(st_adata, sc_ref_adata, config, None, logger)
        self._adata_validate()
        self.overlap_genes = list(self.st_adata.var_names.intersection(self.sc_ref_adata.var_names))
        self.st_adata = self.st_adata[:, self.overlap_genes]
        self.sc_ref_adata = self.sc_ref_adata[:, self.overlap_genes]
        self.svc = {}
        self.graph_aggregate = GraphAggregate(config, logger)

    def local_refinement(self):
        """
        Reconstruct single-cell resolution expression profiles.
        
        This method performs the following steps:
        1. Trims spatial data by removing low-expression genes
        2. For each cell type, constructs an adjacency graph
        3. Uses optimal transport to find neighbor relationships
        4. Aggregates neighbor expressions using graph-based smoothing
        5. Optionally generates UMAP plots for visualization
        
        The reconstructed data is stored in self.svc["sp_svc"].
        """
        if self.config.plot_flag:
            self.logger.info("Plotting Raw ...")
            self._umap_plot(self.st_adata, prefix="Raw")

        svc_recon_adata = self.st_adata.copy()
        self.logger.info(f"before trim: {svc_recon_adata.X.data.shape}")
        svc_recon_adata, celltype_genes = trim_sp_adata(
            svc_recon_adata,
            self.sc_ref_adata,
            self.config.cell_type_col,
        )
        self.logger.info(f"after trim: {svc_recon_adata.X.data.shape}")

        svc_recon_adata.obsm = self.st_adata.obsm.copy()
        # Global Anchoring stores its soft cell-type posterior Q in
        # obsm[cell_type_col]. Q carries the inferred mixture evidence used to
        # disambiguate mixed or mis-segmented spatial units and is propagated
        # unchanged into Local Refinement.
        assignment_categories = pd.Index(
            pd.unique(
                self.sc_ref_adata.obs[self.config.cell_type_col]
            )
        )
        assignment = global_assignment_from_adata(
            svc_recon_adata,
            key=self.config.cell_type_col,
            expected_categories=assignment_categories,
            unknown_key=self.config.unknown_key,
        )
        conditioning_strength = self.config.local_refinement_strength
        refinement_applied = False
        cell_type_adata_list = []
        for cell_type in tqdm(svc_recon_adata.obs[self.config.cell_type_col].unique().tolist(), desc="Reconstructing"):
            svc_recon_adata_cell_type = svc_recon_adata[svc_recon_adata.obs[self.config.cell_type_col] == cell_type]
            raw_st_adata_cell_type = svc_recon_adata_cell_type.copy()
            self.logger.info(f"begin OT smoothing for cell type: {cell_type}, adata shape: {svc_recon_adata_cell_type.shape}")
            # The graph path uses ARPACK PCA with up to 50 components, which
            # requires n_components < n_samples. Keep exactly-50 groups on the
            # established no-smoothing fallback instead of failing in PCA.
            if svc_recon_adata_cell_type.shape[0] <= 50:
                self.logger.info(f"cell type: {cell_type}, has too few spots, skip OT smoothing")
                cell_type_adata_list.append(svc_recon_adata_cell_type)
            else:
                adjacent_matrix = get_adjacency_graph(
                    svc_recon_adata_cell_type,
                    data_type="sc",
                    neighbors_method=self.config.rec_graph_method,
                    alpha=self.config.rec_graph_alpha,
                    gene_neighbor_num=self.config.rec_graph_exp_neighbor_num,
                    spatial_neighbor_num=self.config.rec_graph_spatial_neighbor_num,
                )
                svc_recon_adata_cell_type.obsp["joint_connectivities"] = adjacent_matrix

                (
                    similarity_matrix,
                    neighbor_margin_expr,
                    neighbor_idx_matrix,
                    valid_neighbor_mask,
                    dense_fallback_rows,
                ) = _compute_topk_expression(
                    adjacent_matrix=adjacent_matrix,
                    expression_matrix=svc_recon_adata_cell_type.X,
                    n_neighbors=self.config.rec_graph_n_neighbors,
                    dtype=svc_recon_adata_cell_type.X.dtype,
                    progress=True,
                )
                if dense_fallback_rows:
                    self.logger.info(
                        f"TopK expression used dense fallback for {dense_fallback_rows}/"
                        f"{adjacent_matrix.shape[0]} rows to preserve exact compatibility ordering."
                    )

                mu = np.ravel(svc_recon_adata_cell_type.X.sum(axis=1))
                nu = neighbor_margin_expr

                source_idx, target_idx, active_support = stabilize_local_ot_support(
                    nu,
                    mu,
                    valid_neighbor_mask.T,
                )
                if source_idx.size == 0 or target_idx.size == 0:
                    self.logger.info(
                        f"cell type: {cell_type}, skip OT smoothing due to empty active support"
                    )
                    cell_type_adata_list.append(svc_recon_adata_cell_type)
                    continue
                stable_support = np.zeros(valid_neighbor_mask.T.shape, dtype=bool)
                stable_support[np.ix_(source_idx, target_idx)] = active_support
                valid_neighbor_mask = stable_support.T
                distance_matrix = similarity_to_distance(
                    similarity_matrix,
                    valid_neighbor_mask,
                )
                # Local Refinement combines the propagated Q with spatially
                # informed neighbor support: posterior compatibility conditions
                # local OT costs, and the resulting coupling is applied gene-wise
                # to neighboring observed expression to reconstruct spatially
                # specific profiles.
                distance_matrix = condition_sp_local_ot_cost(
                    distance_matrix,
                    assignment=assignment,
                    left_observations=svc_recon_adata_cell_type.obs_names,
                    right_observations=svc_recon_adata_cell_type.obs_names,
                    neighbor_indices=neighbor_idx_matrix,
                    valid_support_mask=valid_neighbor_mask,
                    strength=conditioning_strength,
                )
                distance_matrix[~valid_neighbor_mask] = np.inf
                T_transform = OTKernel.couple(
                    nu,
                    mu,
                    distance_matrix.T,
                    method=self.config.rec_ot_method,
                    pot_reg=self.config.rec_pot_reg,
                    pot_reg_m=self.config.rec_pot_reg_m,
                    pot_reg_type=self.config.rec_pot_reg_type,
                    pot_verbose=True,
                    pot_num_iter_max=5000,
                    reference_measure=None,
                    valid_support_mask=valid_neighbor_mask.T,
                )
                # Ensure expressions are unchanged before aggregation
                _validate_expression_unchanged(
                    svc_recon_adata_cell_type.X,
                    raw_st_adata_cell_type.X,
                    cell_type=cell_type,
                )
                svc_recon_adata_cell_type = self.graph_aggregate.run(
                    adata=svc_recon_adata_cell_type,
                    neighbor_idx_matrix=neighbor_idx_matrix,
                    coupling_matrix=T_transform,
                    valid_neighbor_mask=valid_neighbor_mask,
                )
                callback = getattr(
                    self.config,
                    "local_refinement_applied_callback",
                    None,
                )
                if callback is not None:
                    callback()
                refinement_applied = True
                cell_type_adata_list.append(svc_recon_adata_cell_type)
        self.svc["sp_svc"] = sc.concat(cell_type_adata_list)
        self.svc["sp_svc"].X = sparse.csr_matrix(self.svc["sp_svc"].X)

        if self.config.plot_flag:
            self.logger.info("Plotting spSVC...")
            self._umap_plot(self.svc["sp_svc"], prefix="sp_SVC")
        return refinement_applied

    def _umap_plot(self, adata, prefix):
        """
        Generate UMAP visualization plots.
        
        Args:
            adata: AnnData object to plot
            prefix: Prefix string for output file names
            
        This method performs preprocessing (filtering, normalization, PCA),
        computes clustering at multiple resolutions, and generates UMAP
        and spatial scatter plots saved to the result directory.
        """
        from revise.analysis.metrics import compute_clustering_metrics

        adata = adata.copy()
        sc.pp.filter_cells(adata, min_genes=self.config.plot_min_genes)
        sc.pp.filter_genes(adata, min_cells=self.config.plot_min_cells)

        if self.config.plot_sample_size > 0:
            self.logger.info(f"Downsampling to {self.config.plot_sample_size} cells for plotting ...")
            np.random.seed(self.config.plot_sample_size)
            indices = np.random.choice(adata.shape[0], self.config.plot_sample_size, replace=False)
            adata = adata[indices, :].copy()

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        if adata.shape[1] > 2000:
            self.logger.info(f"Highly variable genes filtering for {adata.shape[1]} genes.")
            sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3", subset=True)
            adata = adata[:, adata.var.highly_variable].copy()

        sc.tl.pca(adata, svd_solver='arpack')
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

        for res in self.config.plot_cluster_resolution:
            sc.tl.leiden(adata, resolution=res, key_added=f"leiden_res_{res}")
            n_clusters = len(adata.obs[f"leiden_res_{res}"].cat.categories)
            self.logger.info(f"Number of clusters for leiden resolution {res}: {n_clusters}")
            ari, nmi = compute_clustering_metrics(adata, f"leiden_res_{res}", self.config.cell_type_col)
            self.logger.info(f"ari: {ari}, nmi: {nmi}")

        sc.tl.umap(adata)
        umap_resolution = [f"leiden_res_{res}" for res in self.config.plot_cluster_resolution]
        umap_resolution = [self.config.cell_type_col] + umap_resolution
        sc.pl.umap(adata, color=umap_resolution, show=False)

        plt.savefig(os.path.join(self.config.result_dir, f"{prefix}_umap.png"))
        plt.close()
        for res in self.config.plot_cluster_resolution:
            sc.pl.scatter(adata, x='x', y='y', color=f'leiden_res_{res}', show=False)
            plt.savefig(os.path.join(self.config.result_dir, f"{prefix}_resolution_{res}_scatter.png"))
            plt.close()
