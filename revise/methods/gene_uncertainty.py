import os

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy.sparse import coo_matrix
from tqdm import tqdm

from revise.methods.base_method import BaseMethod
from revise.tools.coefficients import compute_autocorr_metrics


class GeneUncertainty(BaseMethod):
    """
    Gene uncertainty evaluation using spatial autocorrelation metrics.
    
    This class evaluates gene expression uncertainty by comparing
    spatial autocorrelation (Moran's I) under different HVG selection
    strategies (in-panel vs all-panel).
    """
    def __init__(self, config, logger):
        super().__init__(config, logger)

    def run(self, adata: AnnData, overlap_genes: list):
        """
        Compare gene uncertainty between in-panel and all-panel HVG selection.
        
        Args:
            adata: Single-cell AnnData object
            overlap_genes: List of genes to use as training genes for in-panel strategy
            
        Returns:
            pd.DataFrame: Comparison DataFrame containing:
                - moranI_in_panel: Moran's I using in-panel HVG selection
                - moranI_all_panel: Moran's I using all-panel HVG selection
                - delta_moranI: Difference between the two
                - winner: Which strategy performs better ("in_panel", "all_panel", or "tie")
                - test: Whether gene is in test set (not in overlap_genes)
                
        Results are saved to CSV files in the result directory.
        """
        self.logger.info("compare_uncertainty .......")

        # Compare: in_panel vs all_panel (at HVG selection level)
        results = self.compare_moran_in_vs_all_panel(
            adata=adata,
            train_genes=overlap_genes
        )

        # Save respective uncertainty results
        results["uncertainty_in_panel"].to_csv(os.path.join(self.config.result_dir, f"uncertainty_in_panel.csv"))
        results["uncertainty_all_panel"].to_csv(os.path.join(self.config.result_dir, "uncertainty_all_panel.csv"))

        # Save comparison table
        compare_df = results["compare"]
        # compare_df.to_csv(self.config.rec_gene_compare_file)

        # Statistics: which side wins
        win_counts = compare_df["winner"].value_counts()
        with open(os.path.join(self.config.result_dir, "summary.txt"), "w") as f:
            f.write("Winner counts (by moranI)\n")
            f.write(str(win_counts) + "\n")

        compare_df = compare_df[~compare_df['test']]
        win_counts = compare_df["winner"].value_counts()
        with open(os.path.join(self.config.result_dir, "summary.txt"), "a") as f:
            f.write("Winner counts (by moranI, test genes)\n")
            f.write(str(win_counts) + "\n")

        self.logger.info("Done. Results saved to:", self.config.result_dir)

        return compare_df

    def compare_moran_in_vs_all_panel(
            self,
            adata,
            train_genes
    ):
        """
        Compare MoranI under two HVG selection strategies (in_panel vs all_panel):
        - in_panel: Select HVG within train_genes candidate set
        - all_panel: Select HVG from all genes in adata, then intersect with train_genes

        Returns a DataFrame containing MoranI for both strategies and winner labels.
        """

        # Calculate uncertainty (MoranI, etc.)
        df_in = self.get_gene_uncertainty(
            adata=adata,
            train_genes=train_genes
        )
        df_all = self.get_gene_uncertainty(
            adata=adata,
            train_genes=adata.var_names.tolist()
        )

        # Only compare genes present in both sets
        common_genes = df_in.index.intersection(df_all.index)
        cmp = pd.DataFrame(index=common_genes)
        cmp["moranI_in_panel"] = df_in.loc[common_genes, "moranI"]
        cmp["moranI_all_panel"] = df_all.loc[common_genes, "moranI"]
        cmp["delta_moranI"] = cmp["moranI_in_panel"] - cmp["moranI_all_panel"]
        cmp["winner"] = np.where(cmp["delta_moranI"] > 0, "in_panel",
                                 np.where(cmp["delta_moranI"] < 0, "all_panel", "tie"))
        cmp["test"] = ~cmp.index.isin(train_genes)

        # Also include mean/variance for reference (using in_panel side's mean/variance)
        for col in ["mean", "variance"]:
            if col in df_in.columns:
                cmp[col] = df_in.loc[common_genes, col]

        return {
            "hvgs_in_panel": train_genes,
            "hvgs_all_panel": adata.var_names.tolist(),
            "uncertainty_in_panel": df_in,
            "uncertainty_all_panel": df_all,
            "compare": cmp,
        }

    def get_gene_uncertainty(self, adata, train_genes=None, sample_size=None):
        """
        Calculate gene expression uncertainty using spatial autocorrelation.
        
        This method computes Moran's I and Geary's C for each gene, which
        measure spatial autocorrelation. Higher values indicate more
        spatially coherent expression patterns (lower uncertainty).
        
        Args:
            adata: Single-cell AnnData object
            train_genes: List of genes to use for building kNN graph.
                If None, uses all genes in adata
            sample_size: If provided, randomly sample this many cells
                for faster computation
                
        Returns:
            pd.DataFrame: DataFrame with columns:
                - moranI: Moran's I statistic (higher = more spatial coherence)
                - gearyC: Geary's C statistic (lower = more spatial coherence)
                - mean: Mean expression level
                - variance: Expression variance
                - test: Whether gene is in test set (not in train_genes)
        """
        if sample_size is not None:
            print("sample size:", sample_size)
            adata = adata[np.random.choice(adata.n_obs, sample_size, replace=False), :]
        if self.config.cell_type_col is not None:
            cts = adata.obs[self.config.cell_type_col].unique()
        else:
            cts = None
        if train_genes is None:
            train_genes = adata.var_names

        # Keep normalized+log transformed matrix (without scale) for computing metrics
        adata_raw = adata.copy()

        if self.config.rec_graph_preprocess:
            sc.pp.normalize_total(adata_raw)
            sc.pp.log1p(adata_raw)

        X = adata_raw.X
        if hasattr(X, "toarray"):
            X = X.toarray()  # Convert to dense matrix
        adata.var['mean'] = X.mean(axis=0)
        adata.var['variance'] = X.var(axis=0)

        # Build graph
        if cts is not None:
            # Build graph by cell type and merge into global sparse matrix (block diagonal/block embedding)
            rows, cols, data = [], [], []
            for ct in tqdm(cts, desc="celltype"):
                ct_mask = adata.obs[self.config.cell_type_col] == ct
                ct_indices = np.where(ct_mask.values if hasattr(ct_mask, 'values') else ct_mask)[0]
                if ct_indices.size == 0:
                    continue
                ct_adata = adata[ct_mask].copy()
                ct_adata_tmp = self.build_knn_graph(ct_adata, train_genes)
                G_ct = ct_adata_tmp.obsp["connectivities"].tocoo()
                # Map to global indices
                rows.extend(ct_indices[G_ct.row])
                cols.extend(ct_indices[G_ct.col])
                data.extend(G_ct.data)
            G = coo_matrix((data, (rows, cols)), shape=(adata.n_obs, adata.n_obs)).tocsr()
        else:
            adata_tmp = self.build_knn_graph(adata, train_genes)
            G = adata_tmp.obsp["connectivities"]

        # Three metrics
        df = compute_autocorr_metrics(adata_raw, G)

        # Merge
        # df = pd.concat([df_auto, df_lap], axis=1)

        # Add mean and variance information
        df['mean'] = adata.var.loc[df.index, 'mean'].values
        df['variance'] = adata.var.loc[df.index, 'variance'].values

        # Mark test genes
        df["test"] = ~df.index.isin(train_genes)

        return df

    def build_knn_graph(self, adata, train_genes):
        """
        Build k-nearest neighbor graph based on selected genes.
        
        This method constructs a kNN graph using only train_genes for
        dimensionality reduction and neighbor computation. This is used
        to compute spatial autocorrelation metrics.
        
        Args:
            adata: AnnData object to process
            train_genes: List of genes to use for graph construction
            
        Returns:
            AnnData: Processed AnnData with kNN graph stored in obsp["connectivities"]
        """
        adata_tmp = adata.copy()

        # Normalize + log
        sc.pp.normalize_total(adata_tmp)
        sc.pp.log1p(adata_tmp)

        # Only keep train_genes for scale + PCA + kNN
        adata_tmp = adata_tmp[:, train_genes].copy()
        if adata_tmp.n_vars < 10000:
            n_hvgs = 200
            sc.pp.highly_variable_genes(adata_tmp, n_top_genes=n_hvgs)
            self.logger.info(f"Select {n_hvgs} HVGs from {adata_tmp.n_vars}")
        else:
            n_hvgs = 2000
            sc.pp.highly_variable_genes(adata_tmp, n_top_genes=n_hvgs)
            self.logger.info(f"Select {n_hvgs} HVGs from {adata_tmp.n_vars}")

        sc.pp.scale(adata_tmp, max_value=10)
        sc.tl.pca(adata_tmp, n_comps=self.config.rec_graph_n_pcs, svd_solver="arpack")
        sc.pp.neighbors(adata_tmp, n_neighbors=self.config.rec_graph_n_neighbors, n_pcs=self.config.rec_graph_n_pcs)

        return adata_tmp
