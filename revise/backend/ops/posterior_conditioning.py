from __future__ import annotations

from numbers import Real

import numpy as np
import pandas as pd

from revise.backend.ops.assignment import (
    GlobalAssignment,
    GlobalAssignmentContractError,
    validate_global_assignment,
)


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


def condition_local_ot_cost(
    cost: np.ndarray,
    left_posterior: GlobalAssignment,
    neighbor_indices: np.ndarray,
    *,
    right_posterior: GlobalAssignment | None = None,
    strength: Real,
) -> np.ndarray:
    """Condition one ``(n_left, k)`` local-OT cost on fixed GA support."""
    if isinstance(strength, (bool, np.bool_)) or not isinstance(strength, Real):
        raise TypeError("strength must be a real number, excluding bool")
    strength_value = float(strength)
    if not np.isfinite(strength_value) or strength_value < 0:
        raise ValueError("strength must be finite and non-negative")
    if not isinstance(left_posterior, GlobalAssignment):
        raise GlobalAssignmentContractError(
            "left_posterior must be a GlobalAssignment"
        )
    if right_posterior is not None and not isinstance(
        right_posterior, GlobalAssignment
    ):
        raise GlobalAssignmentContractError(
            "right_posterior must be a GlobalAssignment"
        )
    if strength_value == 0.0:
        return cost

    if not isinstance(left_posterior.posterior, pd.DataFrame) or not isinstance(
        left_posterior.labels, pd.Series
    ):
        raise GlobalAssignmentContractError(
            "left_posterior must contain GA labels and posterior"
        )
    if right_posterior is not None and (
        not isinstance(right_posterior.posterior, pd.DataFrame)
        or not isinstance(right_posterior.labels, pd.Series)
    ):
        raise GlobalAssignmentContractError(
            "right_posterior must contain GA labels and posterior"
        )

    left = validate_global_assignment(
        left_posterior,
        expected_observations=left_posterior.posterior.index,
        expected_categories=left_posterior.posterior.columns,
    )
    if right_posterior is None or right_posterior is left_posterior:
        right = left
    else:
        right = validate_global_assignment(
            right_posterior,
            expected_observations=right_posterior.posterior.index,
            expected_categories=left.posterior.columns,
        )

    try:
        raw_cost = np.asarray(cost)
    except (TypeError, ValueError) as exc:
        raise ValueError("cost must use a real numeric dtype") from exc
    if not (
        np.issubdtype(raw_cost.dtype, np.integer)
        or np.issubdtype(raw_cost.dtype, np.floating)
    ):
        raise ValueError("cost must use a real numeric dtype")
    cost_values = raw_cost.astype(np.float64, copy=False)
    support = np.asarray(neighbor_indices)
    if cost_values.ndim != 2:
        raise ValueError("cost must have shape (n_left, k)")
    if support.ndim != 2:
        raise ValueError("neighbor_indices must have shape (n_left, k)")
    if cost_values.shape != support.shape:
        raise ValueError("cost and neighbor_indices shapes differ")
    if cost_values.shape[0] != left.posterior.shape[0]:
        raise ValueError("cost first axis does not match left posterior")
    if not np.all(np.isfinite(cost_values)):
        raise ValueError("cost values must be finite")
    if np.any(cost_values < 0):
        raise ValueError("cost values must be non-negative")
    if not np.issubdtype(support.dtype, np.integer):
        raise ValueError("neighbor_indices must use an integer dtype")
    if np.any(support < 0) or np.any(support >= right.posterior.shape[0]):
        raise ValueError("neighbor_indices contain an out-of-bounds index")

    left_values = left.posterior.to_numpy(dtype=np.float64, copy=False)
    right_values = right.posterior.to_numpy(dtype=np.float64, copy=False)
    affinity = np.einsum(
        "id,ikd->ik",
        left_values,
        right_values[support],
    )
    return cost_values + strength_value * -np.log(np.maximum(affinity, 1e-12))


def posterior_reference_allocation(
    spot_expression: np.ndarray,
    posterior: pd.DataFrame,
    reference_profiles: pd.DataFrame,
    *,
    beta: float = 1.0,
    eps: float = 1e-10,
) -> np.ndarray:
    """Closed-form Q-conditioned mandatory allocation for SR routes."""
    if not isinstance(posterior, pd.DataFrame):
        raise ValueError("posterior must be a pandas DataFrame with category columns")
    if not isinstance(reference_profiles, pd.DataFrame):
        raise ValueError(
            "reference_profiles must be a pandas DataFrame with category rows"
        )

    posterior_aligned = align_posterior_categories(
        posterior,
        pd.Index(reference_profiles.index),
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
