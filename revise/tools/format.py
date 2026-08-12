import numpy as np
from scipy.sparse import issparse


def _sample_values(adata, max_samples=10000):
    X = adata.X
    if issparse(X):
        data = X.data
        if data.size > max_samples:
            data = data[:max_samples]
        return np.asarray(data)
    arr = np.asarray(X)
    if arr.size > max_samples:
        arr = arr.ravel()[:max_samples]
    else:
        arr = arr.ravel()
    return arr


def _looks_normalized(adata, target_sum=1e4, max_samples=200, tol=0.15, max_cv=0.1):
    n_obs = adata.n_obs
    if n_obs == 0:
        return False
    take = min(n_obs, max_samples)
    idx = np.linspace(0, n_obs - 1, num=take, dtype=int)
    X = adata.X
    if issparse(X):
        s = X[idx].sum(axis=1)
        if hasattr(s, "A1"):
            totals = s.A1
        else:
            totals = np.asarray(s).ravel()
    else:
        totals = np.asarray(X[idx].sum(axis=1)).ravel()
    mean = float(np.mean(totals))
    if mean <= 0:
        return False
    cv = float(np.std(totals) / mean)
    return (abs(mean - target_sum) / target_sum) <= tol and cv <= max_cv


def _looks_log1p(adata, max_val=30.0, atol=1e-6):
    if hasattr(adata, "uns") and isinstance(adata.uns, dict) and "log1p" in adata.uns:
        return True
    vals = _sample_values(adata)
    if vals.size == 0:
        return False
    integer_like = np.all(np.isclose(vals, np.round(vals), atol=atol))
    if integer_like:
        return False
    return float(np.max(vals)) <= max_val and float(np.min(vals)) >= 0.0


def warn_if_processed(adata, logger, name):
    normalized = _looks_normalized(adata)
    log1p = _looks_log1p(adata)
    if normalized or log1p:
        details = []
        if normalized:
            details.append("normalized")
        if log1p:
            details.append("log1p")
        detail_str = ", ".join(details)
        logger.warning(
            f"{name} appears to be {detail_str}; input is not raw counts. "
            "Results may be unreliable. Please provide raw counts if possible."
        )
