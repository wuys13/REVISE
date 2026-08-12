import os

import matplotlib.pyplot as plt
import numpy as np
import ot
import scanpy as sc
from scipy import sparse
from tqdm import tqdm

from revise.application.application_svc import ApplicationSVC
from revise.methods.graph_aggregate import GraphAggregate
from revise.tools.metric import compute_clustering_metrics
from revise.tools.shaver import trim_sp_adata
from revise.tools.topology import get_adjacency_graph


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
        svc_recon_adata, celltype_genes = trim_sp_adata(svc_recon_adata, self.sc_ref_adata, "Level1")
        self.logger.info(f"after trim: {svc_recon_adata.X.data.shape}")

        svc_recon_adata.obsm = self.st_adata.obsm.copy()
        cell_type_adata_list = []
        for cell_type in tqdm(svc_recon_adata.obs[self.config.cell_type_col].unique().tolist(), desc="Reconstructing"):
            svc_recon_adata_cell_type = svc_recon_adata[svc_recon_adata.obs[self.config.cell_type_col] == cell_type]
            raw_st_adata_cell_type = svc_recon_adata_cell_type.copy()
            self.logger.info(f"begin OT smoothing for cell type: {cell_type}, adata shape: {svc_recon_adata_cell_type.shape}")
            if svc_recon_adata_cell_type.shape[0] < 50:
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

                cost_matrix = np.zeros((svc_recon_adata_cell_type.shape[0], self.config.rec_graph_n_neighbors),
                                       dtype=svc_recon_adata_cell_type.X.dtype)
                neighbor_margin_expr = np.zeros(self.config.rec_graph_n_neighbors, dtype=svc_recon_adata_cell_type.X.dtype)
                neighbor_idx_matrix = np.zeros((svc_recon_adata_cell_type.shape[0], self.config.rec_graph_n_neighbors),
                                               dtype=np.int32)
                for i in tqdm(range(adjacent_matrix.shape[0]), desc="TopK expression"):
                    row = adjacent_matrix[i].toarray().ravel()
                    if np.count_nonzero(row) == 0:
                        continue
                    idx = np.argpartition(-row, kth=min(self.config.rec_graph_n_neighbors, row.size) - 1)[
                          :self.config.rec_graph_n_neighbors]
                    idx = idx[np.argsort(-row[idx])]

                    neighbor_margin_expr += np.mean(svc_recon_adata_cell_type.X[idx].toarray(), axis=1).ravel()
                    cost_matrix[i] = row[idx].copy()
                    neighbor_idx_matrix[i] = idx

                mu = np.ravel(svc_recon_adata_cell_type.X.sum(axis=1))
                nu = neighbor_margin_expr
                # np.nan_to_num(cost_matrix, cost_matrix.max())

                cost_matrix = 1 - cost_matrix
                cost_matrix = 1 / cost_matrix

                T_transform = ot.unbalanced.sinkhorn_unbalanced(
                    nu,
                    mu,
                    cost_matrix.T / cost_matrix.max(),
                    reg=self.config.rec_pot_reg,
                    reg_m=self.config.rec_pot_reg_m,
                    reg_type=self.config.rec_pot_reg_type,
                    verbose=True,
                    numItermax=5000
                )

                # Ensure expressions are unchanged before aggregation
                if sparse.issparse(svc_recon_adata_cell_type.X) and sparse.issparse(raw_st_adata_cell_type.X):
                    assert (svc_recon_adata_cell_type.X != raw_st_adata_cell_type.X).nnz == 0
                else:
                    assert np.array_equal(
                        np.asarray(svc_recon_adata_cell_type.X),
                        np.asarray(raw_st_adata_cell_type.X)
                    )
                svc_recon_adata_cell_type = self.graph_aggregate.run(
                    adata=svc_recon_adata_cell_type,
                    neighbor_idx_matrix=neighbor_idx_matrix,
                    coupling_matrix=T_transform
                )
                cell_type_adata_list.append(svc_recon_adata_cell_type)
        self.svc["sp_svc"] = sc.concat(cell_type_adata_list)
        self.svc["sp_svc"].X = sparse.csr_matrix(self.svc["sp_svc"].X)

        if self.config.plot_flag:
            self.logger.info("Plotting spSVC...")
            self._umap_plot(self.svc["sp_svc"], prefix="sp_SVC")

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
