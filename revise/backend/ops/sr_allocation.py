"""sc-SVC-sr allocation and projected-posterior conditioning seams."""

from __future__ import annotations

from numbers import Real

import numpy as np
import pandas as pd

from revise.backend.ops.assignment import (
    GlobalAssignment,
    GlobalAssignmentContractError,
    validate_global_assignment,
)
from revise.backend.ops.posterior_conditioning import (
    condition_local_ot_cost,
    posterior_reference_allocation,
)


def mandatory_reference_allocation(
    spot_expression: np.ndarray,
    composition: pd.DataFrame,
    reference_profiles: pd.DataFrame,
) -> np.ndarray:
    """Run the algorithm-defining SR allocation with its fixed beta."""
    return posterior_reference_allocation(
        spot_expression,
        composition,
        reference_profiles,
        beta=1.0,
    )


def spot_global_assignment(
    adata,
    *,
    broad_key: str,
    expected_categories: pd.Index,
) -> GlobalAssignment:
    """Load one strict GA assignment from the spot axis."""
    if broad_key not in adata.obs:
        raise GlobalAssignmentContractError(f"missing obs[{broad_key}] GA labels")
    if broad_key not in adata.obsm:
        raise GlobalAssignmentContractError(
            f"missing obsm[{broad_key}] GA posterior"
        )
    return validate_global_assignment(
        GlobalAssignment(
            labels=adata.obs[broad_key],
            posterior=adata.obsm[broad_key],
        ),
        expected_observations=adata.obs_names,
        expected_categories=expected_categories,
    )


def project_spot_assignment_to_virtual_cells(
    assignment: GlobalAssignment,
    svc_obs: pd.DataFrame,
) -> GlobalAssignment:
    """Project spot Q to virtual cells through an explicit stable mapping."""
    if not isinstance(svc_obs, pd.DataFrame) or not {
        "cell_id",
        "spot_name",
    } <= set(svc_obs.columns):
        raise GlobalAssignmentContractError(
            "projection mapping must contain cell_id and spot_name"
        )
    if svc_obs[["cell_id", "spot_name"]].isna().any(axis=None):
        raise GlobalAssignmentContractError(
            "projection mapping contains null cell or spot labels"
        )
    cell_ids = pd.Index(svc_obs["cell_id"].astype(str), name="cell_id")
    spot_names = pd.Index(svc_obs["spot_name"])
    if not cell_ids.is_unique:
        raise GlobalAssignmentContractError(
            "projection mapping cell IDs must be unique"
        )

    missing_spots = spot_names.difference(
        assignment.posterior.index,
        sort=False,
    ).tolist()
    if missing_spots:
        raise GlobalAssignmentContractError(
            "projection spot mapping references missing observations: "
            f"{missing_spots}"
        )

    posterior = assignment.posterior.loc[spot_names].copy()
    posterior.index = cell_ids
    labels = assignment.labels.loc[spot_names].copy()
    labels.index = cell_ids
    return validate_global_assignment(
        GlobalAssignment(labels=labels, posterior=posterior),
        expected_observations=cell_ids,
        expected_categories=assignment.posterior.columns,
    )


def subset_virtual_assignment(
    assignment: GlobalAssignment,
    observation_labels,
) -> GlobalAssignment:
    """Select one ordered virtual-cell block without repairing its axis."""
    requested = pd.Index(observation_labels)
    if requested.empty or not requested.is_unique or requested.hasnans:
        raise GlobalAssignmentContractError(
            "virtual-cell subset observations must be non-empty and unique"
        )
    missing = requested.difference(
        assignment.posterior.index,
        sort=False,
    ).tolist()
    if missing:
        raise GlobalAssignmentContractError(
            f"virtual-cell subset observations are missing: {missing}"
        )
    return GlobalAssignment(
        labels=assignment.labels.loc[requested],
        posterior=assignment.posterior.loc[requested],
    )


def condition_virtual_cell_ot_cost(
    cost: np.ndarray,
    *,
    assignment: GlobalAssignment,
    neighbor_indices: np.ndarray,
    strength: Real,
    valid_support_mask: np.ndarray | None = None,
) -> np.ndarray:
    """Condition one virtual-cell local OT cost on projected spot Q."""
    return condition_local_ot_cost(
        cost,
        assignment,
        neighbor_indices,
        strength=strength,
        valid_support_mask=valid_support_mask,
    )


def record_mandatory_allocation(
    config,
    *,
    status: str,
    broad_key: str,
    n_spots: int,
    n_virtual_cells: int,
    allocation_method: str,
    reason: str | None = None,
) -> None:
    """Record mandatory allocation independently from local OT smoothing."""
    callback = getattr(config, "sr_allocation_callback", None)
    if callback is None:
        return
    evidence = {
        "status": str(status),
        "broad_key": str(broad_key),
        "n_spots": int(n_spots),
        "n_virtual_cells": int(n_virtual_cells),
        "allocation_method": str(allocation_method),
    }
    if reason is not None:
        evidence["reason"] = str(reason)
    callback(evidence)
