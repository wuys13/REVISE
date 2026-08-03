"""Strict GA posterior conditioning for the sp-SVC routes."""

from __future__ import annotations

from numbers import Real

import numpy as np
import pandas as pd
from anndata import AnnData

from revise.backend.ops.assignment import (
    GlobalAssignment,
    GlobalAssignmentContractError,
    validate_global_assignment,
)
from revise.backend.ops.posterior_conditioning import condition_local_ot_cost


def global_assignment_from_adata(
    adata: AnnData,
    *,
    key: str,
    expected_categories: pd.Index,
) -> GlobalAssignment:
    """Load and strictly validate the GA labels and posterior from AnnData."""
    if key not in adata.obs:
        raise GlobalAssignmentContractError(f"missing obs[{key}] GA labels")
    if key not in adata.obsm:
        raise GlobalAssignmentContractError(f"missing obsm[{key}] GA posterior")

    return validate_global_assignment(
        GlobalAssignment(
            labels=adata.obs[key],
            posterior=adata.obsm[key],
        ),
        expected_observations=adata.obs_names,
        expected_categories=expected_categories,
    )


def _group_assignment(
    assignment: GlobalAssignment,
    observations: pd.Index,
) -> GlobalAssignment:
    return GlobalAssignment(
        labels=assignment.labels.loc[observations],
        posterior=assignment.posterior.loc[observations],
    )


def condition_sp_local_ot_cost(
    cost: np.ndarray,
    *,
    assignment: GlobalAssignment,
    left_observations: pd.Index,
    right_observations: pd.Index,
    neighbor_indices: np.ndarray,
    strength: Real,
) -> np.ndarray:
    """Condition one route-owned local OT cost from strict GA posteriors."""
    left = _group_assignment(
        assignment,
        left_observations,
    )
    right = (
        left
        if left_observations.equals(right_observations)
        else _group_assignment(assignment, right_observations)
    )
    return condition_local_ot_cost(
        cost,
        left,
        neighbor_indices,
        right_posterior=right,
        strength=strength,
    )
