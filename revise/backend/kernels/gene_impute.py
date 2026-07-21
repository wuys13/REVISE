import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from tqdm import tqdm

from revise.backend.kernels.base import BaseKernel


class GeneImputeKernel(BaseKernel):
    """
    Gene imputation method using optimal transport neighbor weights.
    
    This class imputes missing gene expression values in spatial data
    by aggregating expressions from similar single cells based on
    optimal transport coupling weights.
    """
    def __init__(self, config, logger):
        super().__init__(config, logger)

    def run(self, adata_st, adata_sc, genes_to_predict, neighbor_weights):
        """
        Impute gene expression for spatial spots using single-cell reference.
        
        Args:
            adata_st: Spatial transcriptomics AnnData object
            adata_sc: Single-cell reference AnnData object
            genes_to_predict: List of gene names to impute, or None for all genes
            neighbor_weights: DataFrame or array of shape (n_spots, n_cells)
                containing optimal transport coupling weights
                
        Returns:
            AnnData: Imputed spatial data with predicted gene expressions
            
        Raises:
            ValueError: If neighbor_weights dimensions don't match input data,
                or if imputation method is not recognized
                
        The imputation method can be:
        - 'mean': Simple average of top-k neighbors
        - 'weighted': Weighted average using OT coupling weights
        - 'knn_weighted': KNN-weighted average with inverse weight transformation
        """
        if genes_to_predict is None:
            genes_to_predict = adata_sc.var_names
        else:
            genes_to_predict = [g for g in genes_to_predict if g in adata_sc.var_names]
        print(f"number of genes to predict: {len(genes_to_predict)}")

        sc_view = adata_sc[:, genes_to_predict]
        if sp.issparse(sc_view.X):
            sc_expr = sc_view.X.tocsr().astype(np.float32)
        else:
            sc_expr = sp.csr_matrix(np.asarray(sc_view.X, dtype=np.float32))

        n_spots, n_cells = neighbor_weights.shape
        if n_cells == 0 or n_spots == 0:
            raise ValueError("ot_mapping is empty after alignment, please check if spots / cells names match.")

        if n_spots != adata_st.n_obs:
            raise ValueError(f"ot_mapping rows ({n_spots}) != adata_st.n_obs ({adata_st.n_obs})")
        if n_cells != adata_sc.n_obs:
            raise ValueError(f"ot_mapping cols ({n_cells}) != adata_sc.n_obs ({adata_sc.n_obs})")

        print(f"OT mapping dimensions: {n_spots} spots × {n_cells} cells")
        k = int(min(max(self.config.rec_impute_n_neighbors, 1), n_cells))
        M = neighbor_weights.to_numpy(dtype=np.float32, copy=False)
        np.nan_to_num(M, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
        M[M < 0] = 0.0

        imputed = np.zeros((n_spots, len(genes_to_predict)), dtype=np.float32)
        global_mean = np.asarray(sc_expr.mean(axis=0)).ravel().astype(np.float32, copy=False)

        for i in tqdm(range(n_spots), desc="Imputing with OT"):
            row = M[i]
            topk = np.argpartition(row, -k)[-k:]
            topk = topk[row[topk] > 0]

            # OT may produce an all-zero row on degenerate local mappings.
            # Fall back to global reference mean to keep the pipeline stable.
            if topk.size == 0:
                imputed[i] = global_mean
                continue

            expr_subset = sc_expr[topk]

            if self.config.rec_impute_method == 'mean':
                imputed[i] = np.asarray(expr_subset.mean(axis=0)).ravel()
                weights = np.ones(topk.size)
            elif self.config.rec_impute_method == 'weighted':
                weights = row[topk]
                w_sum = float(weights.sum())
                if w_sum <= 0:
                    imputed[i] = np.asarray(expr_subset.mean(axis=0)).ravel()
                    continue
                weights = weights / w_sum
                imputed[i] = weights @ expr_subset
            elif self.config.rec_impute_method == 'knn_weighted':
                weights = row[topk]
                w_sum = float(weights.sum())
                if w_sum <= 0:
                    imputed[i] = np.asarray(expr_subset.mean(axis=0)).ravel()
                    continue
                w = 1 - weights / w_sum
                denom = max(1, len(w) - 1)
                weights = w / denom
                imputed[i] = weights @ expr_subset
            else:
                raise ValueError("mode must be 'weighted', 'mean', or 'knn_weighted'")

        adata_imputed = sc.AnnData(
            X=imputed,
            obs=adata_st.obs.copy(),
            var=pd.DataFrame(index=pd.Index(genes_to_predict))
        )
        adata_imputed.var_names = pd.Index(genes_to_predict)

        return adata_imputed
