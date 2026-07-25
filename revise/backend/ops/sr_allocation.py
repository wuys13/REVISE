"""sc-SVC SR mandatory allocation and optional guidance seams."""

from __future__ import annotations

from dataclasses import replace
from typing import Any, Callable

import numpy as np
import pandas as pd

from revise.backend.ops.assignment import (
    AssignmentState,
    AssignmentStateError,
    project_assignment,
    validate_assignment,
)
from revise.backend.ops.assignment_guidance import (
    assignment_guidance_mode,
    assignment_compatibility,
    ot_cost_guidance,
    resolve_assignment_guidance,
)
from revise.backend.ops.posterior_conditioning import (
    posterior_conditioning_mode,
    posterior_reference_allocation,
    reference_measure_from_marginals,
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


def spot_assignment_state(adata, *, broad_key: str) -> AssignmentState:
    """Build the configured broad soft assignment state on the spot axis."""
    if broad_key not in adata.obsm:
        raise AssignmentStateError("assignment_state_unavailable")
    posterior = adata.obsm[broad_key]
    if not isinstance(posterior, pd.DataFrame):
        raise AssignmentStateError("category_labels_missing")
    return validate_assignment(
        AssignmentState(
            values=posterior.to_numpy(dtype=np.float64, copy=True),
            observation_labels=posterior.index,
            category_labels=posterior.columns,
            source=f"obsm[{broad_key}]",
            level=str(broad_key),
            value_semantics="soft",
            lineage=[
                {
                    "operation": "source",
                    "axis": "spot",
                    "category_axis": str(broad_key),
                }
            ],
        )
    )


def projected_virtual_assignment(
    adata,
    svc_obs: pd.DataFrame,
    *,
    broad_key: str,
) -> AssignmentState:
    """Project the spot soft posterior through an explicit virtual-cell map."""
    if not isinstance(svc_obs, pd.DataFrame) or not {
        "cell_id",
        "spot_name",
    } <= set(svc_obs.columns):
        raise AssignmentStateError("projection_mapping_invalid")
    cell_ids = svc_obs["cell_id"].astype(str)
    spot_names = svc_obs["spot_name"].astype(str)
    if cell_ids.duplicated().any():
        raise AssignmentStateError("observation_labels_duplicate")
    state = spot_assignment_state(adata, broad_key=broad_key)
    projected = project_assignment(
        state,
        dict(zip(cell_ids, spot_names)),
        source=f"project(obsm[{broad_key}])",
        level="virtual_cell",
    )
    return validate_assignment(
        replace(
            projected,
            lineage=[
                *projected.lineage,
                {
                    "operation": "spot_to_virtual_projection",
                    "mapping": "svc_obs[cell_id->spot_name]",
                    "source_axis": "spot",
                    "target_axis": "virtual_cell",
                    "category_axis": str(broad_key),
                },
            ],
        )
    )


def subset_virtual_assignment(
    state: AssignmentState,
    observation_labels,
) -> AssignmentState:
    """Select an ordered virtual-cell block from a route-level state."""
    state = validate_assignment(state)
    requested = tuple(
        str(label).replace("/", "_") for label in observation_labels
    )
    if not requested or len(set(requested)) != len(requested):
        raise AssignmentStateError("observation_labels_invalid")
    positions = {
        label: index for index, label in enumerate(state.observation_labels)
    }
    if any(label not in positions for label in requested):
        raise AssignmentStateError("observation_labels_mismatch")
    order = [positions[label] for label in requested]
    return validate_assignment(
        replace(
            state,
            values=state.values[order],
            observation_labels=requested,
        )
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
    """Record allocation independently from optional-guidance outcomes."""
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


def guidance_mode(config) -> str:
    return assignment_guidance_mode(config)


def _start_guidance_event(
    config,
    *,
    problem_key: str,
    applicability: str,
):
    callback = getattr(config, "assignment_guidance_callback", None)
    if callback is not None:
        callback(
            "start",
            problem_key=problem_key,
            route=str(
                getattr(
                    config,
                    "assignment_guidance_route",
                    "sc_svc_sr",
                )
            ),
            operator="virtual_cell_ot",
            phase="lr",
            mode=guidance_mode(config),
            applicability=applicability,
            numerics={
                "beta": float(
                    getattr(config, "posterior_conditioning_beta", 1.0)
                ),
                "min_affinity": float(
                    getattr(
                        config,
                        "posterior_conditioning_min_affinity",
                        0.05,
                    )
                ),
                "operator_strength": float(
                    getattr(
                        config,
                        "posterior_conditioning_cost_strength",
                        0.2,
                    )
                ),
            },
            solver=str(config.rec_ot_method),
        )
    return callback


def record_virtual_cell_not_applicable(
    config,
    *,
    problem_key: str,
    reason: str,
) -> None:
    callback = _start_guidance_event(
        config,
        problem_key=problem_key,
        applicability="not_applicable",
    )
    if callback is not None:
        callback(
            "terminal",
            problem_key=problem_key,
            outcome="not_applicable",
            reason=reason,
        )


def record_virtual_cell_unavailable(
    config,
    *,
    problem_key: str,
    reason: str,
) -> None:
    """Record an applicable problem whose route capability is unavailable."""
    callback = _start_guidance_event(
        config,
        problem_key=problem_key,
        applicability="applicable",
    )
    mode = guidance_mode(config)
    if mode == "off":
        if callback is not None:
            callback(
                "terminal",
                problem_key=problem_key,
                outcome="off",
                reason="guidance_off",
            )
        return
    outcome = "fallback" if mode == "prefer" else "failed"
    if callback is not None:
        callback(
            "terminal",
            problem_key=problem_key,
            outcome=outcome,
            availability="unavailable",
            reason=reason,
        )
    if outcome == "failed":
        raise ValueError(
            f"required assignment guidance unavailable: {reason}"
        )


def prepare_virtual_cell_guidance(
    config,
    *,
    problem_key: str,
    state_loader: Callable[[], AssignmentState | None],
    neighbor_support: np.ndarray,
    distance_matrix: np.ndarray,
    source_mass: np.ndarray,
    target_mass: np.ndarray,
) -> tuple[np.ndarray, np.ndarray | None, bool]:
    """Resolve and inject guidance before the route-owned local solver."""
    callback = _start_guidance_event(
        config,
        problem_key=problem_key,
        applicability="applicable",
    )
    mode = guidance_mode(config)
    if mode == "off":
        if callback is not None:
            callback(
                "terminal",
                problem_key=problem_key,
                outcome="off",
                reason="guidance_off",
            )
        return distance_matrix, None, False

    resolution = resolve_assignment_guidance(mode, state_loader)
    if resolution.availability != "available":
        if callback is not None:
            callback(
                "terminal",
                problem_key=problem_key,
                outcome=resolution.outcome,
                availability=resolution.availability,
                reason=resolution.reason,
            )
        if resolution.outcome == "failed":
            raise ValueError(
                f"assignment guidance unavailable: {resolution.reason}"
            )
        return distance_matrix, None, False

    state = resolution.state
    assert state is not None
    affinity = assignment_compatibility(
        state,
        state,
        beta=getattr(config, "posterior_conditioning_beta", 1.0),
        min_affinity=getattr(
            config,
            "posterior_conditioning_min_affinity",
            0.05,
        ),
        support=neighbor_support,
    )
    compatibility_mode = posterior_conditioning_mode(config)
    reference_measure = None
    if compatibility_mode == "cost":
        distance_matrix = ot_cost_guidance(
            distance_matrix,
            affinity,
            getattr(
                config,
                "posterior_conditioning_cost_strength",
                0.2,
            ),
        )
    elif compatibility_mode == "reference":
        reference_measure = reference_measure_from_marginals(
            source_mass,
            target_mass,
            affinity.T,
        )
    else:  # pragma: no cover - resolved configuration rejects this
        raise ValueError(
            f"unsupported virtual-cell compatibility mode: {compatibility_mode}"
        )
    if callback is not None:
        callback(
            "attempt",
            problem_key=problem_key,
            availability="available",
            left_assignment=state,
            right_assignment=state,
        )
    return distance_matrix, reference_measure, True


def record_virtual_cell_guidance_terminal(
    config,
    *,
    problem_key: str,
    attempted: bool,
    outcome: str,
    reason: str | None = None,
) -> None:
    callback = getattr(config, "assignment_guidance_callback", None)
    if not attempted or callback is None:
        return
    fields: dict[str, Any] = {"outcome": outcome}
    if reason is not None:
        fields["reason"] = reason
    callback("terminal", problem_key=problem_key, **fields)
