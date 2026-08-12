import numba
import numpy as np
import pandas as pd
from numba import njit
from numba import prange
from scipy.sparse import csc_matrix
from scipy.sparse import csr_matrix
from scipy.sparse import issparse
from sklearn.utils.sparsefuncs import inplace_row_scale


def bhattacharyya_distance(A, B):
    """
    Calculate Bhattacharyya distance between two matrices.
    
    This distance metric measures the dissimilarity between probability
    distributions using the Bhattacharyya coefficient, transformed to
    a distance measure.
    
    Args:
        A: First matrix (n_samples_A, n_features)
        B: Second matrix (n_samples_B, n_features)
        
    Returns:
        np.ndarray: Distance matrix of shape (n_samples_A, n_samples_B)
            with values in [0, 1], where 0 means identical distributions
    """
    ABT = bhattacharyya_coefficient(A, B, parallel=True)

    return np.maximum(1 - ABT, 0)


def bhattacharyya_coefficient(A, B, parallel=True):
    """
    Calculate Bhattacharyya coefficient between two matrices.
    
    The Bhattacharyya coefficient measures the overlap between probability
    distributions. Rows are normalized as probabilities, then transformed
    to probability amplitudes (square root) before computing overlap.
    
    Args:
        A: First matrix (n_samples_A, n_features)
        B: Second matrix (n_samples_B, n_features)
        parallel: Whether to use parallel computation
        
    Returns:
        np.ndarray: Coefficient matrix of shape (n_samples_A, n_samples_B)
            with values in [0, 1], where 1 means identical distributions
    """
    Anorm = get_sum(A, axis=1)
    Bnorm = get_sum(B, axis=1)

    A = A.copy()
    B = B.copy()

    row_scale(A, 1 / Anorm)
    row_scale(B, 1 / Bnorm)

    # transform to probability amplitudes
    A = np.sqrt(A)
    B = np.sqrt(B)

    # get overlap
    ABT = gemmT(A, B, parallel=parallel)

    return ABT


def get_sum(
        X,
        axis,
        dtype=None,
):
    """
    Calculates the sum of a sparse matrix or array-like in a specified axis and
    returns a flattened result.

    Parameters
    ----------
    X
        A 2d :class:`~numpy.ndarray` or a `scipy` sparse matrix.
    axis
        The axis along which to calculate the sum
    dtype
        The dtype used in the accumulators and the result.

    Returns
    -------
    The flattened sums as 1d :class:`~numpy.ndarray`.

    """

    if issparse(X):
        result = X.sum(axis=axis, dtype=dtype).A.flatten()
    else:
        result = X.sum(axis=axis, dtype=dtype)

    import anndata._core.views
    if isinstance(result, anndata._core.views.ArrayView):
        result = result.toarray()

    return result


def row_scale(
        X,
        rescaling_factors,
        round=False,
):
    """
    Rescales rows of dense or sparse matrix inplace.

    Parameters
    ----------
    X
        A 2d :class:`~numpy.ndarray` array or a `scipy` sparse matrix.
    rescaling_factors
        A 1d :class:`~numpy.ndarray` containing the row-wise rescaling factors.
    round
        Whether to round the result to integer values

    Returns
    -------
    `None`. This is an inplace operation.

    """

    if hasattr(rescaling_factors, 'to_numpy'):
        rescaling_factors = rescaling_factors.to_numpy()
    rescaling_factors = rescaling_factors.flatten()

    if issparse(X):
        inplace_row_scale(X, rescaling_factors)
        if round:
            np.around(X.data, out=X.data)
    else:
        X *= rescaling_factors[:, None]
        if round:
            np.around(X, out=X)


def gemmT(
        A,
        B,
        parallel=True,
        sparse_result=False,
):
    """
    Perform a matrix-matrix multiplication A @ B.T for arbitrary sparseness of
    A and B in parallel. Uses `sparse_dot_mkl` if available.

    Parameters
    ----------
    A
        A 2d :class:`~numpy.ndarray` or a `scipy` sparse matrix.
    B
        A 2d :class:`~numpy.ndarray` or a `scipy` sparse matrix with the same
        second dimension as `A`.
    parallel
        Whether to run the multiplication in parallel.
    sparse_result
        Whether to return a sparse result when both inputs are sparse

    Returns
    -------
    Depending on `sparse_result` returns either a :class:`~numpy.ndarray` or a\
    scipy sparse matrix containing the result.

    """

    if not issparse(A) and not isinstance(A, np.ndarray):
        raise ValueError("`A` can only be a scipy sparse matrix or a numpy array!")
    if not issparse(B) and not isinstance(B, np.ndarray):
        raise ValueError("`B` can only be a scipy sparse matrix or a numpy array!")
    if issparse(A) and not isinstance(A, (csr_matrix, csc_matrix)):
        A = A.tocsr()
    if issparse(B) and not isinstance(B, (csr_matrix, csc_matrix)):
        B = B.tocsr()

    A, B = cast_down_common(A, B)

    numba_threads = numba.get_num_threads()
    if not parallel:
        numba.set_num_threads(1)

    result_shape = (A.shape[0], B.shape[0])
    inner_size = A.shape[1]
    if B.shape[1] != A.shape[1]:
        raise ValueError(
            "`A` and `B` dont have the same second dimension! The shapes are %s and %s!" % (A.shape, B.shape))

    # use numba/scipy implementation directly, without sparse_dot_mkl dependency
    if issparse(A) and issparse(B):
        if not sparse_result:
            A = A.tocsr()
            BT = B.tocsc()  # avoids intermediate during transposition: .tocsc() == .T.tocsr()
            return _csrcsr_gemm_dense(result_shape[0], result_shape[1], A.indptr, A.indices, A.data, BT.indptr, BT.indices, BT.data)
        else:
            return A @ (B.T)
    elif issparse(A) and not issparse(B):
        A = A.tocsc()
        return _cscdense_gemm_dense(result_shape[0], result_shape[1], A.indptr, A.indices, A.data, B.T)
    elif not issparse(A) and issparse(B):
        BT = B.tocsc()  # avoids intermediate during transposition: .tocsc() == .T.tocsr()
        return _densecsr_gemm_dense(result_shape[0], result_shape[1], A, BT.indptr, BT.indices, BT.data)
    else:
        return A @ (B.T)


def cast_down_common(A, B):
    if not pd.api.types.is_float_dtype(A):
        A = A.astype(np.float64)
    if not pd.api.types.is_float_dtype(B):
        B = B.astype(np.float64)

    if A.dtype == B.dtype:
        return A, B
    elif A.dtype == np.float64 and B.dtype == np.float32:
        return A.astype(np.float32), B
    elif A.dtype == np.float32 and B.dtype == np.float64:
        return A, B.astype(np.float32)
    else:
        return A, B


@njit(fastmath=True, parallel=True, cache=True)
def _csrcsr_gemm_dense(n_row, n_col, Ap, Aj, Ax, Bp, Bj, Bx):
    C = np.zeros((n_row, n_col))

    for i in prange(n_row):
        jj_start = Ap[i]
        jj_end = Ap[i + 1]
        for jj in range(jj_start, jj_end):
            j = Aj[jj]
            v = Ax[jj]

            kk_start = Bp[j]
            kk_end = Bp[j + 1]
            for kk in range(kk_start, kk_end):
                k = Bj[kk]

                C[i, k] += v * Bx[kk]

    return C


@njit(fastmath=True, parallel=True, cache=True)
def _cscdense_gemm_dense(n_row, n_col, Ap, Ai, Ax, B):
    nk = B.shape[0]

    C = np.zeros((n_col, n_row))

    for j in prange(n_col):
        for k in range(nk):
            v = B[k, j]

            ii_start = Ap[k]
            ii_end = Ap[k + 1]
            for ii in range(ii_start, ii_end):
                i = Ai[ii]

                C[j, i] += v * Ax[ii]

    return C.T


@njit(fastmath=True, parallel=True, cache=True)
def _densecsr_gemm_dense(n_row, n_col, A, Bp, Bj, Bx):
    nk = A.shape[1]

    C = np.zeros((n_row, n_col))

    for i in prange(n_row):
        for j in range(nk):
            v = A[i, j]

            kk_start = Bp[j]
            kk_end = Bp[j + 1]
            for kk in range(kk_start, kk_end):
                k = Bj[kk]

                C[i, k] += v * Bx[kk]

    return C
