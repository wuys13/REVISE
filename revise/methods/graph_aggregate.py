import numpy as np
from anndata import AnnData
from scipy import sparse

from revise.methods.base_method import BaseMethod


class GraphAggregate(BaseMethod):
    """
    Graph-based expression aggregation using optimal transport weights.
    
    This class aggregates gene expressions from neighboring spots/cells
    using optimal transport coupling weights to create a smoothed
    expression profile.
    """
    def __init__(self, config, logger):
        super().__init__(config, logger)

    def run(self, adata: AnnData, neighbor_idx_matrix, coupling_matrix):
        """
        Aggregate neighbor expressions using OT coupling weights.
        
        Args:
            adata: AnnData object to update (will be modified in place)
            neighbor_idx_matrix: Array of shape (n_spot, K) containing
                indices of K nearest neighbors for each spot
            coupling_matrix: Array of shape (K, n_spot) containing
                optimal transport coupling weights from neighbors to spots
                
        Returns:
            AnnData: Updated AnnData with smoothed expression in adata.X
            
        The method:
        1. Constructs a sparse weight matrix from coupling weights
        2. Normalizes weights so each spot's neighbor weights sum to 1
        3. Computes weighted average of neighbor expressions
        4. Fuses with original expression using config.rec_alpha
        """
        raw_X = adata.X.copy()
        n_spot = adata.n_obs
        K = neighbor_idx_matrix.shape[1]
        assert coupling_matrix.shape == (K, n_spot), "T_transform should be (K, n_spot)"

        # Construct sparse OT weight matrix W (n_spot x n_spot)
        rows = np.repeat(np.arange(n_spot, dtype=np.int32), K)
        cols = neighbor_idx_matrix.reshape(-1)
        data = coupling_matrix.T.reshape(-1).astype(np.float64)

        data = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)
        data[data < 0] = 0.0

        W = sparse.csr_matrix((data, (rows, cols)), shape=(n_spot, n_spot))

        # Row normalization (each spot's weights to neighbors sum to 1)
        row_sums = np.asarray(W.sum(axis=1)).ravel()

        # Handle all-zero rows: use uniform distribution as fallback
        zero_mask = row_sums == 0
        if np.any(zero_mask):
            zero_rows = np.where(zero_mask)[0]
            zr_rows_rep = np.repeat(zero_rows, K)
            zr_cols = neighbor_idx_matrix[zero_rows].reshape(-1)
            zr_data = np.full(zr_rows_rep.shape[0], 1.0 / K, dtype=float)

            W_backup = sparse.csr_matrix((zr_data, (zr_rows_rep, zr_cols)),
                                         shape=(n_spot, n_spot))
            W = W + W_backup
            row_sums = np.asarray(W.sum(axis=1)).ravel()

        d_inv = 1.0 / row_sums
        W = sparse.diags(d_inv) @ W

        assert np.allclose(np.asarray(W.sum(axis=1)).ravel(), 1.0), "Weight matrix not correctly normalized"

        # Compute neighbor weighted average, then fuse with original expression
        X = adata.X
        if sparse.issparse(X):
            X_nb = W @ X
            X_smooth = (1 - self.config.rec_alpha) * raw_X + self.config.rec_alpha * X_nb
        else:
            X_nb = W @ np.asarray(X)
            X_smooth = (1 - self.config.rec_alpha) * np.asarray(raw_X) + self.config.rec_alpha * X_nb

        adata.X = X_smooth

        return adata
