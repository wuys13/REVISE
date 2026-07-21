from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse


_VALID_MODES = {"off", "cost", "reference"}


def posterior_conditioning_mode(config: Any) -> str:
    """Return the configured posterior-conditioning mode."""
    enabled = bool(getattr(config, "posterior_conditioning_enabled", False))
    if not enabled:
        return "off"
    raw_mode = getattr(config, "posterior_conditioning_mode", "off")
    mode = "off" if raw_mode is False or raw_mode is None else str(raw_mode).strip().lower()
    if mode in {"", "false", "none", "null"}:
        mode = "off"
    if mode not in _VALID_MODES:
        raise ValueError(
            "posterior_conditioning_mode must be one of "
            f"{sorted(_VALID_MODES)}, got {mode!r}"
        )
    return mode


def posterior_conditioning_enabled(config: Any, kind: str) -> bool:
    """Check whether a conditioning component is active."""
    mode = posterior_conditioning_mode(config)
    if kind == "cost":
        return mode == "cost"
    if kind == "reference":
        return mode == "reference"
    raise ValueError(f"Unknown posterior-conditioning kind: {kind!r}")


def posterior_conditioning_strict(config: Any) -> bool:
    """Return whether requested posterior conditioning must fail instead of fallback."""
    return bool(getattr(config, "posterior_conditioning_strict", False))


def get_posterior_matrix(adata, key: str, eps: float = 1e-12) -> np.ndarray | None:
    """Extract a row-normalized posterior matrix from ``adata.obsm[key]``."""
    if adata is None or key not in adata.obsm:
        return None
    values = adata.obsm[key]
    if hasattr(values, "to_numpy"):
        q = values.to_numpy(dtype=np.float64)
    else:
        q = np.asarray(values, dtype=np.float64)
    if q.ndim != 2 or q.shape[0] != adata.n_obs:
        return None
    q = np.nan_to_num(q, nan=0.0, posinf=0.0, neginf=0.0)
    q[q < 0] = 0.0
    row_sums = q.sum(axis=1, keepdims=True)
    zero_rows = row_sums[:, 0] <= eps
    if np.any(zero_rows):
        q[zero_rows] = 1.0 / q.shape[1]
        row_sums = q.sum(axis=1, keepdims=True)
    return q / np.maximum(row_sums, eps)


def _normalized_category_axis(labels, name: str) -> pd.Index:
    labels = pd.Index(labels)
    if bool(pd.isna(labels).any()):
        raise ValueError(f"{name} contain null category labels")
    original_duplicates = labels[labels.duplicated()].tolist()
    if original_duplicates:
        raise ValueError(
            f"{name} contain duplicate category labels: {original_duplicates}"
        )
    normalized = pd.Index([str(label).replace("/", "_") for label in labels])
    normalized_duplicates = normalized[normalized.duplicated()].tolist()
    if normalized_duplicates:
        raise ValueError(
            f"{name} collide after '/' normalization: {normalized_duplicates}"
        )
    return normalized


def align_posterior_categories(
    posterior: pd.DataFrame,
    reference_labels,
    *,
    posterior_name: str = "posterior columns",
    reference_name: str = "reference categories",
) -> pd.DataFrame:
    """Validate and align posterior columns to a labeled category axis."""
    if not isinstance(posterior, pd.DataFrame):
        raise ValueError("posterior must be a pandas DataFrame with category columns")

    posterior_normalized = _normalized_category_axis(
        posterior.columns,
        posterior_name,
    )
    reference_normalized = _normalized_category_axis(
        reference_labels,
        reference_name,
    )
    posterior_set = set(posterior_normalized.tolist())
    reference_set = set(reference_normalized.tolist())
    missing = [label for label in reference_normalized if label not in posterior_set]
    extra = [label for label in posterior_normalized if label not in reference_set]
    if missing or extra:
        raise ValueError(
            "posterior and reference categories differ after '/' normalization: "
            f"missing={missing}, extra={extra}"
        )

    aligned = posterior.copy()
    aligned.columns = posterior_normalized
    return aligned.reindex(columns=reference_normalized)


def posterior_affinity(
    q_left: np.ndarray,
    q_right: np.ndarray | None = None,
    *,
    beta: float = 1.0,
    min_affinity: float = 1e-3,
) -> np.ndarray:
    """Compute posterior compatibility ``(q_left @ q_right.T)^beta``."""
    if q_right is None:
        q_right = q_left
    affinity = np.asarray(q_left, dtype=np.float64) @ np.asarray(q_right, dtype=np.float64).T
    affinity = np.nan_to_num(affinity, nan=0.0, posinf=0.0, neginf=0.0)
    affinity = np.clip(affinity, float(min_affinity), 1.0)
    beta = float(beta)
    if beta != 1.0:
        affinity = np.power(affinity, beta)
    return affinity


def neighbor_posterior_affinity(
    q_rows: np.ndarray,
    neighbor_idx_matrix: np.ndarray,
    *,
    q_neighbors: np.ndarray | None = None,
    beta: float = 1.0,
    min_affinity: float = 1e-3,
) -> np.ndarray:
    """Compute row-by-neighbor posterior compatibility for top-k OT costs."""
    if q_neighbors is None:
        q_neighbors = q_rows
    neighbor_idx_matrix = np.asarray(neighbor_idx_matrix, dtype=np.int64)
    q_rows = np.asarray(q_rows, dtype=np.float64)
    q_neighbors = np.asarray(q_neighbors, dtype=np.float64)
    right = q_neighbors[neighbor_idx_matrix]
    affinity = np.einsum("id,ikd->ik", q_rows, right)
    affinity = np.nan_to_num(affinity, nan=0.0, posinf=0.0, neginf=0.0)
    affinity = np.clip(affinity, float(min_affinity), 1.0)
    beta = float(beta)
    if beta != 1.0:
        affinity = np.power(affinity, beta)
    return affinity


def condition_cost_matrix(cost: np.ndarray, affinity: np.ndarray, strength: float) -> np.ndarray:
    """Apply the cost-side equivalent of posterior-conditioned entropy."""
    cost = np.asarray(cost, dtype=np.float64)
    affinity = np.asarray(affinity, dtype=np.float64)
    if cost.shape != affinity.shape:
        raise ValueError(f"cost and affinity shapes differ: {cost.shape} vs {affinity.shape}")
    penalty = -np.log(np.clip(affinity, 1e-12, 1.0))
    conditioned = cost + float(strength) * penalty
    finite = np.isfinite(conditioned)
    if np.any(finite):
        high = float(np.max(conditioned[finite]))
        conditioned = np.nan_to_num(conditioned, nan=0.0, posinf=high, neginf=0.0)
    else:
        conditioned = np.zeros_like(conditioned, dtype=np.float64)
    return conditioned


def reference_measure_from_marginals(
    source_mass: np.ndarray,
    target_mass: np.ndarray,
    affinity: np.ndarray,
    *,
    eps: float = 1e-12,
) -> np.ndarray:
    """Build ``c_ij proportional to source_i target_j affinity_ij`` for POT."""
    source_mass = np.asarray(source_mass, dtype=np.float64).ravel()
    target_mass = np.asarray(target_mass, dtype=np.float64).ravel()
    affinity = np.asarray(affinity, dtype=np.float64)
    product = np.outer(source_mass, target_mass)
    if product.shape != affinity.shape:
        raise ValueError(f"marginal product and affinity shapes differ: {product.shape} vs {affinity.shape}")
    reference = product * np.clip(affinity, eps, 1.0)
    product_sum = float(product.sum())
    reference_sum = float(reference.sum())
    if product_sum > eps and reference_sum > eps:
        reference *= product_sum / reference_sum
    reference = np.nan_to_num(reference, nan=eps, posinf=eps, neginf=eps)
    reference[reference <= 0] = eps
    return reference


def posterior_reference_allocation(
    spot_expression: np.ndarray,
    posterior: pd.DataFrame,
    reference_profiles: pd.DataFrame,
    *,
    beta: float = 1.0,
    eps: float = 1e-10,
) -> np.ndarray:
    """Closed-form Q-conditioned reference allocation for SR benchmark paths."""
    if not isinstance(posterior, pd.DataFrame):
        raise ValueError("posterior must be a pandas DataFrame with category columns")
    if not isinstance(reference_profiles, pd.DataFrame):
        raise ValueError(
            "reference_profiles must be a pandas DataFrame with category rows"
        )

    reference_labels = pd.Index(reference_profiles.index)
    posterior_aligned = align_posterior_categories(
        posterior,
        reference_labels,
        reference_name="reference profile index",
    )
    reference_aligned = reference_profiles.copy()
    reference_aligned.index = posterior_aligned.columns

    spot_expression = np.asarray(spot_expression, dtype=np.float64)
    posterior_values = posterior_aligned.to_numpy(dtype=np.float64, copy=False)
    reference_values = reference_aligned.to_numpy(dtype=np.float64, copy=False)
    if spot_expression.ndim != 2:
        raise ValueError(f"spot_expression must be 2D, got shape={spot_expression.shape}")
    if spot_expression.shape[0] != posterior_values.shape[0]:
        raise ValueError(
            "spot_expression and posterior must have the same number of rows: "
            f"{spot_expression.shape[0]} vs {posterior_values.shape[0]}"
        )
    if spot_expression.shape[1] != reference_values.shape[1]:
        raise ValueError(
            "spot_expression and reference_profiles must have the same number of genes: "
            f"{spot_expression.shape[1]} vs {reference_values.shape[1]}"
        )
    if posterior_values.shape[1] != reference_values.shape[0]:
        raise ValueError(
            "posterior columns must match reference profile rows: "
            f"{posterior_values.shape[1]} vs {reference_values.shape[0]}"
        )

    q = np.nan_to_num(posterior_values, nan=0.0, posinf=0.0, neginf=0.0)
    q[q < 0] = 0.0
    beta = float(beta)
    if beta != 1.0:
        q = np.power(q, beta)
    weights = q[:, np.newaxis, :] * reference_values.T[np.newaxis, :, :]
    weights = weights / (np.sum(weights, axis=2, keepdims=True) + float(eps))
    return weights * spot_expression[:, :, np.newaxis]


def condition_sparse_graph(
    graph,
    q_left: np.ndarray,
    q_right: np.ndarray | None = None,
    *,
    beta: float = 1.0,
    min_affinity: float = 1e-3,
):
    """Multiply sparse graph edges by posterior compatibility."""
    if q_right is None:
        q_right = q_left
    graph = graph.tocsr()
    coo = graph.tocoo()
    if coo.nnz == 0:
        return graph
    edge_affinity = np.einsum("ij,ij->i", q_left[coo.row], q_right[coo.col])
    edge_affinity = np.clip(edge_affinity, float(min_affinity), 1.0)
    beta = float(beta)
    if beta != 1.0:
        edge_affinity = np.power(edge_affinity, beta)
    data = coo.data * edge_affinity
    return sparse.csr_matrix((data, (coo.row, coo.col)), shape=graph.shape)
