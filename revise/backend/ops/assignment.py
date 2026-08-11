from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import numpy as np
import pandas as pd


class GlobalAssignmentContractError(ValueError):
    """Raised when a global-anchoring assignment violates its output contract."""


@dataclass(eq=False)
class GlobalAssignment:
    """Mutable GA snapshot; validation returns independently owned values."""

    labels: pd.Series
    posterior: pd.DataFrame


def _strict_global_axis(labels: Iterable[Any], axis: str) -> pd.Index:
    index = pd.Index(labels)
    if index.empty:
        raise GlobalAssignmentContractError(f"{axis} axis must not be empty")
    if index.hasnans:
        raise GlobalAssignmentContractError(f"{axis} axis contains null values")
    if any(isinstance(value, str) and not value.strip() for value in index):
        raise GlobalAssignmentContractError(f"{axis} axis contains empty values")
    if not index.is_unique:
        raise GlobalAssignmentContractError(f"{axis} axis contains duplicate values")
    slash_normalized = [str(value).replace("/", "_") for value in index]
    if len(set(slash_normalized)) != len(slash_normalized):
        raise GlobalAssignmentContractError(
            f"{axis} axis values collide after '/' normalization"
        )
    return index


def _validate_global_assignment_axis(
    actual: Iterable[Any],
    expected: Iterable[Any],
    axis: str,
    *,
    require_order: bool = True,
) -> None:
    """Validate one GA axis by exact raw value and order."""
    actual_index = _strict_global_axis(actual, axis)
    expected_index = _strict_global_axis(expected, f"expected {axis}")
    if actual_index.equals(expected_index):
        return
    missing = expected_index.difference(actual_index, sort=False).tolist()
    extra = actual_index.difference(expected_index, sort=False).tolist()
    if missing or extra:
        raise GlobalAssignmentContractError(
            f"{axis} axis mismatch: missing={missing}, extra={extra}"
        )
    if not require_order:
        return
    raise GlobalAssignmentContractError(
        f"{axis} axis order does not match expected order"
    )


def validate_global_assignment(
    assignment: GlobalAssignment,
    *,
    expected_observations: Iterable[Any],
    expected_categories: Iterable[Any],
    require_category_order: bool = True,
    require_row_normalization: bool = True,
) -> GlobalAssignment:
    """Validate a GA posterior without reordering, normalizing, or repairing it."""
    if not isinstance(assignment, GlobalAssignment):
        raise GlobalAssignmentContractError(
            "assignment must be a GlobalAssignment"
        )
    if not isinstance(assignment.posterior, pd.DataFrame):
        raise GlobalAssignmentContractError(
            "posterior must be a pandas DataFrame"
        )
    if not isinstance(assignment.labels, pd.Series):
        raise GlobalAssignmentContractError("labels must be a pandas Series")

    _validate_global_assignment_axis(
        assignment.posterior.index,
        expected_observations,
        "observation",
    )
    _validate_global_assignment_axis(
        assignment.labels.index,
        expected_observations,
        "label observation",
    )
    _validate_global_assignment_axis(
        assignment.posterior.columns,
        expected_categories,
        "category",
        require_order=require_category_order,
    )

    try:
        values = assignment.posterior.to_numpy(dtype=np.float64, copy=True)
    except (TypeError, ValueError) as exc:
        raise GlobalAssignmentContractError(
            "posterior values must be numeric"
        ) from exc
    if not np.all(np.isfinite(values)):
        raise GlobalAssignmentContractError("posterior values must be finite")
    if np.any(values < 0):
        raise GlobalAssignmentContractError(
            "posterior values must be non-negative"
        )
    row_mass = values.sum(axis=1)
    if np.any(row_mass <= 0):
        raise GlobalAssignmentContractError(
            "posterior rows must have positive mass"
        )
    if require_row_normalization and not np.allclose(
        row_mass, 1.0, rtol=0.0, atol=1e-6
    ):
        raise GlobalAssignmentContractError(
            "posterior rows must be row-normalized within atol=1e-6"
        )

    if assignment.labels.isna().any():
        raise GlobalAssignmentContractError("labels must not contain null values")
    expected_labels = assignment.posterior.idxmax(axis=1)
    if not pd.Index(assignment.labels.to_numpy()).equals(
        pd.Index(expected_labels.to_numpy())
    ):
        raise GlobalAssignmentContractError(
            "labels must equal argmax(posterior) using pandas idxmax"
        )

    return GlobalAssignment(
        labels=assignment.labels.copy(deep=True),
        posterior=assignment.posterior.copy(deep=True),
    )
