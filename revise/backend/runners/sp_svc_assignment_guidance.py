"""sp-SVC route glue for the shared assignment-guidance contract."""

from __future__ import annotations

import numpy as np

from revise.backend.ops.assignment import (
    AssignmentState,
    AssignmentStateError,
    align_assignment_categories,
    align_assignment_observations,
    one_hot_assignment,
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
    reference_measure_from_marginals,
)


def guidance_mode(config):
    return assignment_guidance_mode(config)


def assignment_categories(*adatas, key):
    categories = []
    for adata in adatas:
        if adata is None or key not in adata.obs:
            continue
        for value in adata.obs[key].tolist():
            if value not in categories:
                categories.append(value)
    return categories


def assignment_state_from_adata(adata, *, key, category_labels):
    if key in adata.obsm:
        raw = adata.obsm[key]
        if hasattr(raw, "to_numpy") and hasattr(raw, "columns"):
            values = raw.to_numpy(dtype=np.float64)
            observations = raw.index
            categories = raw.columns
        else:
            values = np.asarray(raw, dtype=np.float64)
            observations = adata.obs_names
            if values.ndim != 2 or len(category_labels) != values.shape[1]:
                raise AssignmentStateError("category_labels_missing")
            categories = category_labels
        state = validate_assignment(
            AssignmentState(
                values=values,
                observation_labels=observations,
                category_labels=categories,
                source=f"obsm[{key}]",
                level=str(key),
                value_semantics="soft",
                lineage=[],
            )
        )
        return align_assignment_observations(state, adata.obs_names)
    if key not in adata.obs:
        return None
    return one_hot_assignment(
        adata.obs[key].tolist(),
        observation_labels=adata.obs_names,
        category_labels=category_labels,
        source=f"obs[{key}]",
        level=str(key),
    )


def _problem_start(
    config,
    *,
    route,
    operator,
    problem_key,
    applicability,
):
    callback = getattr(config, "assignment_guidance_callback", None)
    if callback is not None:
        callback(
            "start",
            problem_key=problem_key,
            route=str(
                getattr(config, "assignment_guidance_route", route)
            ),
            operator=operator,
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


def record_not_applicable(
    config,
    *,
    route,
    operator,
    problem_key,
    reason,
):
    callback = _problem_start(
        config,
        route=route,
        operator=operator,
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


def prepare_assignment_guidance(
    config,
    *,
    route,
    operator,
    problem_key,
    left_adata,
    right_adata,
    category_labels,
    support,
    distance_matrix,
    source_mass,
    target_mass,
):
    callback = _problem_start(
        config,
        route=route,
        operator=operator,
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

    loaded = {}

    def load_state():
        left = assignment_state_from_adata(
            left_adata,
            key=config.cell_type_col,
            category_labels=category_labels,
        )
        right = assignment_state_from_adata(
            right_adata,
            key=config.cell_type_col,
            category_labels=category_labels,
        )
        if left is None or right is None:
            raise KeyError(config.cell_type_col)
        loaded["left"] = left
        loaded["right"] = align_assignment_categories(
            right,
            left.category_labels,
        )
        return left

    resolution = resolve_assignment_guidance(mode, load_state)
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

    affinity = assignment_compatibility(
        loaded["left"],
        loaded["right"],
        beta=getattr(config, "posterior_conditioning_beta", 1.0),
        min_affinity=getattr(
            config,
            "posterior_conditioning_min_affinity",
            0.05,
        ),
        support=support,
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
    if callback is not None:
        callback(
            "attempt",
            problem_key=problem_key,
            availability="available",
            left_assignment=loaded["left"],
            right_assignment=loaded["right"],
        )
    return distance_matrix, reference_measure, True


def record_guidance_terminal(
    config,
    *,
    problem_key,
    attempted,
    outcome,
    reason=None,
):
    callback = getattr(config, "assignment_guidance_callback", None)
    if not attempted or callback is None:
        return
    fields = {"outcome": outcome}
    if outcome == "failed":
        fields["reason"] = reason
    callback("terminal", problem_key=problem_key, **fields)
