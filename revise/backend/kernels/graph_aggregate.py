import numpy as np
from anndata import AnnData
from scipy import sparse

from revise.backend.kernels.base import BaseKernel


class GraphAggregateKernel(BaseKernel):
    """
    Graph-based expression aggregation using optimal transport weights.
    
    This class aggregates gene expressions from neighboring spots/cells
    using optimal transport coupling weights to create a smoothed
    expression profile.
    """
    def __init__(self, config, logger):
        super().__init__(config, logger)

    def run(
        self,
        adata: AnnData,
        neighbor_idx_matrix,
        coupling_matrix,
        alpha_override=None,
        valid_neighbor_mask=None,
    ):
        """
        Aggregate neighbor expressions using OT coupling weights.
        
        Args:
            adata: AnnData object to update (will be modified in place)
            neighbor_idx_matrix: Array of shape (n_spot, K) containing
                indices of K nearest neighbors for each spot
            coupling_matrix: Array of shape (K, n_spot) containing
                optimal transport coupling weights from neighbors to spots
            valid_neighbor_mask: Optional boolean array of shape (n_spot, K)
                identifying semantic neighbor edges
                
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
        neighbor_idx = np.asarray(neighbor_idx_matrix)
        if neighbor_idx.ndim != 2 or neighbor_idx.shape[0] != n_spot:
            raise ValueError(
                "neighbor_idx_matrix must have shape (n_obs, K), got "
                f"{neighbor_idx.shape} for n_obs={n_spot}"
            )
        K = neighbor_idx.shape[1]

        coupling = np.asarray(coupling_matrix, dtype=np.float64)
        if coupling.shape != (K, n_spot):
            raise ValueError(
                "coupling_matrix must have shape (K, n_obs), got "
                f"{coupling.shape}; expected {(K, n_spot)}"
            )
        if not np.all(np.isfinite(coupling)):
            raise ValueError("coupling_matrix must contain only finite values")
        if np.any(coupling < 0):
            raise ValueError("coupling_matrix must be non-negative")

        if valid_neighbor_mask is None:
            valid = np.ones((n_spot, K), dtype=bool)
        else:
            valid = np.asarray(valid_neighbor_mask, dtype=bool)
            if valid.shape != (n_spot, K):
                raise ValueError(
                    "valid_neighbor_mask must have shape (n_obs, K), got "
                    f"{valid.shape}; expected {(n_spot, K)}"
                )

        coupling_by_obs = coupling.T.copy()
        invalid_mass = float(np.abs(coupling_by_obs[~valid]).sum())
        if invalid_mass > 1e-12:
            raise ValueError(
                "coupling_matrix has mass on invalid neighbor slots: "
                f"total absolute mass={invalid_mass:.3e}, tolerance=1e-12"
            )
        coupling_by_obs[~valid] = 0.0

        # Construct sparse OT weight matrix W (n_spot x n_spot)
        valid_flat = valid.reshape(-1)
        rows = np.repeat(np.arange(n_spot, dtype=np.int32), K)[valid_flat]
        cols = neighbor_idx.reshape(-1)[valid_flat]
        data = coupling_by_obs.reshape(-1)[valid_flat]

        W = sparse.csr_matrix((data, (rows, cols)), shape=(n_spot, n_spot))

        # Row normalization (each spot's weights to neighbors sum to 1)
        row_sums = np.asarray(W.sum(axis=1)).ravel()

        # Preserve expression exactly when a row has no valid coupling.
        zero_rows = np.flatnonzero(row_sums == 0)
        if zero_rows.size:
            W = W.tolil()
            W[zero_rows, zero_rows] = 1.0
            W = W.tocsr()
            row_sums = np.asarray(W.sum(axis=1)).ravel()

        d_inv = 1.0 / row_sums
        W = sparse.diags(d_inv) @ W

        assert np.allclose(np.asarray(W.sum(axis=1)).ravel(), 1.0), "Weight matrix not correctly normalized"

        # Compute neighbor weighted average, then fuse with original expression.
        # `alpha_override` can be a scalar or a per-row vector to support
        # confidence-weighted smoothing in benchmark experiments.
        X = adata.X
        alpha_row = None
        if alpha_override is not None:
            alpha_arr = np.asarray(alpha_override, dtype=np.float64)
            if alpha_arr.ndim == 0:
                alpha_arr = np.full((n_spot,), float(alpha_arr), dtype=np.float64)
            if alpha_arr.shape[0] != n_spot:
                raise ValueError("alpha_override length must match number of observations")
            alpha_arr = np.clip(alpha_arr, 0.0, 1.0)
            alpha_row = alpha_arr[:, None]

        if sparse.issparse(X):
            X_nb = W @ X
            if alpha_row is None:
                X_smooth = (1 - self.config.rec_alpha) * raw_X + self.config.rec_alpha * X_nb
            else:
                raw_dense = raw_X.toarray() if sparse.issparse(raw_X) else np.asarray(raw_X)
                nb_dense = X_nb.toarray() if sparse.issparse(X_nb) else np.asarray(X_nb)
                X_smooth = (1 - alpha_row) * raw_dense + alpha_row * nb_dense
        else:
            X_nb = W @ np.asarray(X)
            raw_arr = np.asarray(raw_X)
            if alpha_row is None:
                X_smooth = (1 - self.config.rec_alpha) * raw_arr + self.config.rec_alpha * X_nb
            else:
                X_smooth = (1 - alpha_row) * raw_arr + alpha_row * X_nb

        adata.X = X_smooth

        return adata
