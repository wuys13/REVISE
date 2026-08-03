from __future__ import annotations

import copy
import hashlib
import json
import sys
import types
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from revise.backend.ops import posterior_conditioning
from revise.backend.ops.assignment import (
    AssignmentState,
    AssignmentStateError,
    GlobalAssignment,
    GlobalAssignmentContractError,
    aggregate_assignment,
    align_assignment_categories,
    align_assignment_observations,
    argmax_assignment,
    assignment_state_evidence,
    one_hot_assignment,
    project_assignment,
    validate_assignment,
)
from revise.backend.ops.assignment_guidance import (
    AssignmentGuidanceCollector,
    FallbackReason,
    NotApplicableReason,
    assignment_compatibility,
    assignment_guidance_mode,
    graph_guidance,
    ot_cost_guidance,
    resolve_assignment_guidance,
)


def _state(
    values,
    *,
    observations=("o1", "o2"),
    categories=("A", "B"),
    semantics="soft",
    lineage=None,
):
    return AssignmentState(
        values=np.asarray(values, dtype=float),
        observation_labels=observations,
        category_labels=categories,
        source="test",
        level="cell_type",
        value_semantics=semantics,
        lineage=list(lineage or []),
    )


def _global_assignment(
    values,
    *,
    observations=("o1", "o2"),
    categories=("A", "B"),
):
    posterior = pd.DataFrame(values, index=observations, columns=categories)
    return GlobalAssignment(labels=posterior.idxmax(axis=1), posterior=posterior)


@pytest.mark.parametrize("mode", ["off", "prefer", "require"])
def test_assignment_guidance_mode_requires_canonical_resolved_policy(mode):
    config = SimpleNamespace(
        assignment_guidance_policy=mode,
        posterior_conditioning_enabled=mode == "off",
        posterior_conditioning_mode="reference",
        posterior_conditioning_strict=mode != "require",
    )

    assert assignment_guidance_mode(config) == mode


def test_assignment_guidance_mode_rejects_legacy_only_runner_config():
    config = SimpleNamespace(
        posterior_conditioning_enabled=True,
        posterior_conditioning_mode="cost",
        posterior_conditioning_strict=False,
    )

    with pytest.raises(
        ValueError,
        match="missing resolved assignment_guidance_policy",
    ):
        assignment_guidance_mode(config)


def test_assignment_axes_align_by_label_and_argmax_is_one_hot():
    left = validate_assignment(_state([[0.8, 0.2], [0.1, 0.9]]))
    reversed_categories = validate_assignment(
        _state([[0.2, 0.8], [0.9, 0.1]], categories=("B", "A"))
    )
    reversed_observations = validate_assignment(
        _state(
            [[0.1, 0.9], [0.8, 0.2]],
            observations=("o2", "o1"),
        )
    )

    assert np.allclose(
        align_assignment_categories(reversed_categories, left.category_labels).values,
        left.values,
    )
    assert np.allclose(
        align_assignment_observations(
            reversed_observations, left.observation_labels
        ).values,
        left.values,
    )
    one_hot = argmax_assignment(left)
    explicit_one_hot = one_hot_assignment(
        ("A", "B"),
        observation_labels=("o1", "o2"),
        category_labels=("A", "B"),
        source="argmax_labels",
        level="cell_type",
    )
    assert one_hot.value_semantics == "one_hot"
    assert np.array_equal(one_hot.values, explicit_one_hot.values)
    assert one_hot.category_labels == left.category_labels
    projected_one_hot = project_assignment(
        explicit_one_hot,
        {"v1": "o1", "v2": "o2"},
        source="virtual_argmax",
        level="virtual",
    )
    assert projected_one_hot.value_semantics == "one_hot"
    assert projected_one_hot.lineage[-1]["operation"] == "project"


def test_assignment_evidence_summarizes_ordered_axes_without_storing_labels():
    state = validate_assignment(
        _state(
            [[0.8, 0.2], [0.1, 0.9]],
            observations=("细胞-1", "cell-2"),
            lineage=[{"operation": "aggregate"}],
        )
    )

    evidence = assignment_state_evidence(state)
    expected_obs = hashlib.sha256(
        json.dumps(
            ["细胞-1", "cell-2"],
            ensure_ascii=False,
            separators=(",", ":"),
        ).encode("utf-8")
    ).hexdigest()
    reversed_evidence = assignment_state_evidence(
        align_assignment_observations(state, ("cell-2", "细胞-1"))
    )

    assert evidence == {
        "source": "test",
        "level": "cell_type",
        "value_semantics": "soft",
        "lineage": [{"operation": "aggregate"}],
        "observation_axis": {
            "count": 2,
            "ordered_labels_sha256": expected_obs,
        },
        "category_axis": {
            "count": 2,
            "ordered_labels_sha256": hashlib.sha256(
                b'["A","B"]'
            ).hexdigest(),
        },
    }
    assert (
        reversed_evidence["observation_axis"]["ordered_labels_sha256"]
        != evidence["observation_axis"]["ordered_labels_sha256"]
    )
    assert "labels" not in evidence["observation_axis"]


def test_validate_assignment_returns_independently_owned_mutable_state():
    values = np.array([[0.8, 0.2], [0.1, 0.9]])
    lineage = [{"operation": "aggregate", "members": ["c1", "c2"]}]
    raw = AssignmentState(
        values=values,
        observation_labels=("o1", "o2"),
        category_labels=("A", "B"),
        source="raw",
        level="cell_type",
        value_semantics="soft",
        lineage=lineage,
    )

    validated = validate_assignment(raw)
    values[0, 0] = 0.0
    lineage[0]["members"].append("c3")
    validated.source = "owned"

    assert not np.shares_memory(validated.values, values)
    assert np.allclose(validated.values[0], [0.8, 0.2])
    assert validated.lineage == [
        {"operation": "aggregate", "members": ["c1", "c2"]}
    ]
    assert raw.source == "raw"


def test_validate_assignment_rejects_finite_values_with_overflowing_row_mass():
    state = _state([[1e308, 1e308], [0.5, 0.5]])

    with pytest.raises(AssignmentStateError) as exc_info:
        validate_assignment(state)

    assert exc_info.value.reason == "values_row_mass_nonfinite"


def test_event_records_bilateral_assignment_evidence_and_terminal_can_update_it():
    left = validate_assignment(
        _state(
            [[0.8, 0.2], [0.1, 0.9]],
            lineage=[{"operation": "aggregate"}],
        )
    )
    right = one_hot_assignment(
        ("A", "B"),
        observation_labels=("r1", "r2"),
        category_labels=("A", "B"),
        source="argmax",
        level="reference",
        lineage=[{"operation": "project"}],
    )
    collector = AssignmentGuidanceCollector()
    collector.callback(
        "start",
        problem_key="bilateral",
        route="sp_svc:bin2cell",
        operator="local_ot",
        phase="local_refinement",
        mode="prefer",
        applicability="applicable",
        numerics={},
        solver="pot",
        left_assignment=left,
        right_assignment=None,
    )
    collector.callback(
        "attempt",
        problem_key="bilateral",
        availability="available",
    )
    collector.callback(
        "terminal",
        problem_key="bilateral",
        outcome="applied",
        right_assignment=right,
    )

    event = collector.events[0]
    assert event["left_assignment"]["value_semantics"] == "soft"
    assert event["left_assignment"]["lineage"] == [
        {"operation": "aggregate"}
    ]
    assert event["right_assignment"]["value_semantics"] == "one_hot"
    assert event["right_assignment"]["lineage"] == [
        {"operation": "project"}
    ]


def test_aggregate_and_explicit_projection_preserve_axes_and_lineage():
    state = validate_assignment(
        _state(
            [[1.0, 0.0], [0.0, 1.0], [0.5, 0.5]],
            observations=("c1", "c2", "c3"),
        )
    )

    aggregated = aggregate_assignment(
        state,
        ("s1", "s1", "s2"),
        source="cell_to_spot",
        level="spot",
    )
    projected = project_assignment(
        aggregated,
        {"v2": "s2", "v1": "s1"},
        source="spot_to_virtual",
        level="virtual",
    )

    assert aggregated.observation_labels == ("s1", "s2")
    assert np.allclose(aggregated.values, [[0.5, 0.5], [0.5, 0.5]])
    assert projected.observation_labels == ("v2", "v1")
    assert projected.category_labels == ("A", "B")
    assert projected.value_semantics == "soft"
    assert [entry["operation"] for entry in projected.lineage] == [
        "aggregate",
        "project",
    ]
    repeated = project_assignment(
        aggregated,
        {"v1": "s1", "v2": "s1"},
        source="spot_to_virtual",
        level="virtual",
    )
    assert np.array_equal(repeated.values[0], repeated.values[1])


@pytest.mark.parametrize(
    ("state", "reason"),
    [
        (_state([[0.5, 0.5]], observations=("o1",), categories=()), "category_labels_empty"),
        (_state([[0.5, 0.5]], observations=("o1",), categories=("A", None)), "category_labels_null"),
        (_state([[0.5, 0.5]], observations=("o1",), categories=("A", "A")), "category_labels_duplicate"),
        (_state([[0.5, 0.5]], observations=("o1",), categories=("A/B", "A_B")), "category_labels_slash_collision"),
        (_state([[np.nan, 1.0]], observations=("o1",)), "values_nan"),
        (_state([[np.inf, 1.0]], observations=("o1",)), "values_infinite"),
        (_state([[-1.0, 2.0]], observations=("o1",)), "values_negative"),
        (_state([[0.0, 0.0]], observations=("o1",)), "values_zero_row"),
    ],
)
def test_invalid_assignment_has_stable_reason(state, reason):
    with pytest.raises(AssignmentStateError) as exc:
        validate_assignment(state)

    assert exc.value.reason == reason


@pytest.mark.parametrize(
    ("axis", "labels", "reason"),
    [
        ("observations", ("o1", "missing"), "observation_labels_mismatch"),
        ("observations", ("o1", "o1"), "observation_labels_duplicate"),
        ("categories", ("A", "missing"), "category_labels_mismatch"),
        ("categories", ("A", "A"), "category_labels_duplicate"),
    ],
)
def test_axis_alignment_rejects_missing_extra_or_duplicate(axis, labels, reason):
    state = validate_assignment(_state([[0.8, 0.2], [0.1, 0.9]]))

    with pytest.raises(AssignmentStateError) as exc:
        if axis == "observations":
            align_assignment_observations(state, labels)
        else:
            align_assignment_categories(state, labels)

    assert exc.value.reason == reason


def test_dense_and_supported_compatibility_share_penalty_semantics():
    left = validate_assignment(_state([[0.8, 0.2], [0.1, 0.9]]))
    right = validate_assignment(
        _state(
            [[0.2, 0.8], [0.9, 0.1], [0.5, 0.5]],
            observations=("r1", "r2", "r3"),
            categories=("B", "A"),
        )
    )
    dense = assignment_compatibility(left, right, beta=2.0, min_affinity=0.05)
    top_k = np.array([[2, 0], [1, 2]])
    supported = assignment_compatibility(
        left,
        right,
        beta=2.0,
        min_affinity=0.05,
        support=top_k,
    )
    edge_rows = np.array([0, 1, 1])
    edge_columns = np.array([2, 1, 0])
    edge_supported = assignment_compatibility(
        left,
        right,
        beta=2.0,
        min_affinity=0.05,
        support=(edge_rows, edge_columns),
    )

    assert np.allclose(supported, np.take_along_axis(dense, top_k, axis=1))
    assert np.allclose(edge_supported, dense[edge_rows, edge_columns])
    cost = np.ones_like(supported)
    weights = np.ones_like(supported)
    strength = 0.7
    assert np.allclose(
        ot_cost_guidance(cost, supported, strength) - cost,
        -strength * np.log(supported),
    )
    assert np.allclose(
        graph_guidance(weights, supported, strength),
        weights * np.power(supported, strength),
    )


def test_local_ot_cost_conditioning_uses_fixed_directed_top_k_formula():
    left = _global_assignment([[1.0, 0.0], [0.0, 1.0]])
    right = _global_assignment(
        [[1.0, 0.0], [0.0, 1.0], [0.5, 0.5]],
        observations=("r1", "r2", "r3"),
    )
    cost = np.array([[1.0, 2.0], [3.0, 4.0]])
    support = np.array([[1, 2], [0, 2]])

    conditioned = posterior_conditioning.condition_local_ot_cost(
        cost,
        left,
        support,
        right_posterior=right,
        strength=0.5,
    )
    expected_affinity = np.array([[1e-12, 0.5], [1e-12, 0.5]])

    assert np.allclose(
        conditioned,
        cost + 0.5 * -np.log(expected_affinity),
        rtol=0.0,
        atol=1e-14,
    )

    changed_left = _global_assignment([[0.75, 0.25], [0.25, 0.75]])
    changed = posterior_conditioning.condition_local_ot_cost(
        cost,
        changed_left,
        support,
        right_posterior=right,
        strength=0.5,
    )
    assert not np.array_equal(changed, conditioned)


def test_local_ot_cost_conditioning_does_not_clip_affinity_above_one():
    posterior = _global_assignment(
        [[1.0000005, 0.0]],
        observations=("o1",),
    )
    cost = np.array([[1.0]])

    conditioned = posterior_conditioning.condition_local_ot_cost(
        cost,
        posterior,
        np.array([[0]]),
        strength=1.0,
    )

    assert conditioned[0, 0] == pytest.approx(
        1.0 - np.log(1.0000005**2),
    )
    assert conditioned[0, 0] < cost[0, 0]


@pytest.mark.parametrize("right_mode", ["default", "same"])
def test_local_ot_cost_conditioning_validates_shared_posterior_once(
    monkeypatch,
    right_mode,
):
    posterior = _global_assignment([[1.0, 0.0], [0.0, 1.0]])
    original = posterior_conditioning.validate_global_assignment
    calls = []

    def counting_validator(*args, **kwargs):
        calls.append(args[0])
        return original(*args, **kwargs)

    monkeypatch.setattr(
        posterior_conditioning,
        "validate_global_assignment",
        counting_validator,
    )
    kwargs = {} if right_mode == "default" else {"right_posterior": posterior}

    posterior_conditioning.condition_local_ot_cost(
        np.ones((2, 1)),
        posterior,
        np.array([[0], [1]]),
        strength=0.2,
        **kwargs,
    )

    assert calls == [posterior]


def test_local_ot_zero_strength_returns_original_cost_unchanged():
    posterior = _global_assignment([[1.0, 0.0], [0.0, 1.0]])
    opaque_posterior = GlobalAssignment(labels=None, posterior=None)
    cost = np.array([[np.inf, np.nan], [2.0, 3.0]], dtype=np.float32)

    conditioned = posterior_conditioning.condition_local_ot_cost(
        cost,
        posterior,
        np.array([[0, 1], [1, 0]]),
        strength=0.0,
    )

    assert conditioned is cost
    assert conditioned.dtype == cost.dtype
    assert conditioned.shape == cost.shape
    assert np.array_equal(conditioned, cost, equal_nan=True)
    assert np.isposinf(conditioned[0, 0])
    assert (
        posterior_conditioning.condition_local_ot_cost(
            cost,
            opaque_posterior,
            np.array([[0, 1], [1, 0]]),
            right_posterior=opaque_posterior,
            strength=0.0,
        )
        is cost
    )
    with pytest.raises(GlobalAssignmentContractError, match="left_posterior"):
        posterior_conditioning.condition_local_ot_cost(
            cost,
            None,
            np.array([[0, 1], [1, 0]]),
            strength=0.0,
        )


@pytest.mark.parametrize(
    ("strength", "error"),
    [
        (True, TypeError),
        (False, TypeError),
        ("1", TypeError),
        (np.nan, ValueError),
        (np.inf, ValueError),
        (-0.1, ValueError),
    ],
)
def test_local_ot_cost_conditioning_rejects_invalid_strength(strength, error):
    posterior = _global_assignment([[1.0, 0.0], [0.0, 1.0]])

    with pytest.raises(error, match="strength"):
        posterior_conditioning.condition_local_ot_cost(
            np.ones((2, 1)),
            posterior,
            np.array([[0], [1]]),
            strength=strength,
        )


@pytest.mark.parametrize(
    ("cost", "support", "message"),
    [
        (np.ones((2, 2)), np.array([[0], [1]]), "shapes differ"),
        (np.ones((2, 1)), np.array([0, 1]), "neighbor_indices must have shape"),
        (np.ones((2, 1)), np.array([[0.0], [1.0]]), "integer dtype"),
        (np.ones((2, 1)), np.array([[0], [2]]), "out-of-bounds"),
        (np.array([[1.0], [np.inf]]), np.array([[0], [1]]), "finite"),
        (np.array([[1.0], [-1.0]]), np.array([[0], [1]]), "non-negative"),
    ],
)
def test_local_ot_cost_conditioning_rejects_invalid_cost_or_support(
    cost,
    support,
    message,
):
    posterior = _global_assignment([[1.0, 0.0], [0.0, 1.0]])

    with pytest.raises(ValueError, match=message):
        posterior_conditioning.condition_local_ot_cost(
            cost,
            posterior,
            support,
            strength=0.2,
        )


@pytest.mark.parametrize(
    "cost",
    [
        np.array([[True], [False]]),
        np.array([["1.0"], ["2.0"]]),
        np.array([[1.0], [2.0]], dtype=object),
        np.array([[1.0 + 1.0j], [2.0 + 0.0j]]),
    ],
)
def test_local_ot_cost_conditioning_rejects_non_real_numeric_cost_dtype(cost):
    posterior = _global_assignment([[1.0, 0.0], [0.0, 1.0]])

    with pytest.raises(ValueError, match="real numeric dtype"):
        posterior_conditioning.condition_local_ot_cost(
            cost,
            posterior,
            np.array([[0], [1]]),
            strength=0.2,
        )


def test_local_ot_cost_conditioning_rejects_misaligned_or_invalid_posterior():
    left = _global_assignment([[1.0, 0.0], [0.0, 1.0]])
    reversed_categories = _global_assignment(
        [[0.0, 1.0], [1.0, 0.0]],
        observations=("r1", "r2"),
        categories=("B", "A"),
    )
    invalid = _global_assignment([[0.6, 0.6], [0.0, 1.0]])

    with pytest.raises(GlobalAssignmentContractError, match="category"):
        posterior_conditioning.condition_local_ot_cost(
            np.ones((2, 1)),
            left,
            np.array([[0], [1]]),
            right_posterior=reversed_categories,
            strength=0.2,
        )
    with pytest.raises(GlobalAssignmentContractError, match="row-normalized"):
        posterior_conditioning.condition_local_ot_cost(
            np.ones((2, 1)),
            invalid,
            np.array([[0], [1]]),
            strength=0.2,
        )


@pytest.mark.parametrize(
    "support",
    [
        (np.array([[0, 1]]), np.array([0, 1])),
        (np.array([0, 1]), np.array([[0, 1]])),
        (np.array([0]), np.array([0, 1])),
        (np.array([0.0, 1.0]), np.array([0, 1])),
        (np.array([0, -1]), np.array([0, 1])),
        (np.array([0, 2]), np.array([0, 1])),
        (np.array([0, 1]), np.array([0, 3])),
        np.array([0, 1]),
        np.array([[0], [1], [2]]),
        np.array([[0.0], [1.5]]),
        np.array([[0], [-1]]),
        np.array([[0], [3]]),
    ],
)
def test_compatibility_rejects_invalid_support_indices_without_coercion(support):
    left = validate_assignment(_state([[0.8, 0.2], [0.1, 0.9]]))
    right = validate_assignment(
        _state(
            [[0.9, 0.1], [0.2, 0.8], [0.5, 0.5]],
            observations=("r1", "r2", "r3"),
        )
    )

    with pytest.raises(ValueError, match="support"):
        assignment_compatibility(left, right, support=support)


def test_policy_off_does_not_read_state_prefer_falls_back_and_require_raises():
    reads = []

    def invalid_loader():
        reads.append(True)
        raise AssignmentStateError("values_negative")

    off = resolve_assignment_guidance("off", invalid_loader)
    prefer = resolve_assignment_guidance("prefer", invalid_loader)
    with pytest.raises(AssignmentStateError, match="values_negative"):
        resolve_assignment_guidance("require", invalid_loader)

    assert reads == [True, True]
    assert (off.availability, off.outcome, off.reason) == (
        "not_checked",
        "off",
        None,
    )
    assert (prefer.availability, prefer.outcome, prefer.reason) == (
        "unavailable",
        "fallback",
        FallbackReason.ASSIGNMENT_INVALID,
    )
    assert prefer.reason_details == {"cause": "values_negative"}


def _event(collector, key, *, applicability="applicable", mode="prefer"):
    return collector.callback(
        "start",
        problem_key=key,
        route="sp_svc:bin2cell",
        operator="local_ot",
        phase="local_refinement",
        mode=mode,
        applicability=applicability,
        numerics={"strength": 0.7},
        solver="pot",
        left_assignment=None,
        right_assignment=None,
    )


def test_collector_records_pre_solver_availability_attempt_and_terminal_outcome():
    collector = AssignmentGuidanceCollector()
    _event(collector, "available-then-failed")
    collector.callback(
        "attempt",
        problem_key="available-then-failed",
        availability="available",
    )
    collector.callback(
        "terminal",
        problem_key="available-then-failed",
        outcome="failed",
    )
    _event(collector, "unavailable")
    collector.callback(
        "terminal",
        problem_key="unavailable",
        availability="unavailable",
        outcome="fallback",
        reason=FallbackReason.ASSIGNMENT_MISSING,
    )

    first, second = collector.events
    assert first["ordinal"] == 1
    assert first["requested"] is True
    assert first["availability"] == "available"
    assert first["attempted"] is True
    assert first["outcome"] == "failed"
    assert second["ordinal"] == 2
    assert second["availability"] == "unavailable"
    assert second["attempted"] is False
    assert second["outcome"] == "fallback"
    assert collector.summary() == "failed"


def test_collector_summary_distinguishes_zero_inapplicable_mixed_and_interrupted():
    empty = AssignmentGuidanceCollector()
    assert empty.summary() == "not_started"

    off = AssignmentGuidanceCollector()
    _event(off, "disabled", mode="off")
    off.callback(
        "terminal",
        problem_key="disabled",
        outcome="off",
    )
    assert off.summary() == "off"

    inapplicable = AssignmentGuidanceCollector()
    _event(inapplicable, "no-local-problem", applicability="not_applicable")
    inapplicable.callback(
        "terminal",
        problem_key="no-local-problem",
        outcome="not_applicable",
        reason=NotApplicableReason.INSUFFICIENT_UNITS,
    )
    assert inapplicable.summary() == "not_applicable"

    mixed = AssignmentGuidanceCollector()
    for key, outcome in (("one", "applied"), ("two", "fallback")):
        _event(mixed, key)
        if outcome == "applied":
            mixed.callback("attempt", problem_key=key, availability="available")
        mixed.callback(
            "terminal",
            problem_key=key,
            outcome=outcome,
            availability="unavailable" if outcome == "fallback" else "available",
            reason=(
                None
                if outcome == "applied"
                else FallbackReason.ASSIGNMENT_MISSING
            ),
        )
    assert mixed.summary() == "mixed"

    interrupted = AssignmentGuidanceCollector()
    _event(interrupted, "done")
    interrupted.callback("attempt", problem_key="done", availability="available")
    interrupted.callback("terminal", problem_key="done", outcome="applied")
    _event(interrupted, "stopped")
    interrupted.callback("attempt", problem_key="stopped", availability="available")
    interrupted.callback(
        "terminal",
        problem_key="stopped",
        outcome="interrupted",
    )
    assert [event["outcome"] for event in interrupted.events] == [
        "applied",
        "interrupted",
    ]
    assert interrupted.summary() == "interrupted"


@pytest.mark.parametrize("outcome", ["not_applicable", "fallback"])
def test_tolerated_terminal_outcomes_require_a_stable_reason(outcome):
    collector = AssignmentGuidanceCollector()
    if outcome == "not_applicable":
        _event(
            collector,
            "reason-required",
            applicability="not_applicable",
        )
        terminal_fields = {}
    else:
        _event(collector, "reason-required")
        terminal_fields = {"availability": "unavailable"}

    with pytest.raises(ValueError, match="Reason"):
        collector.callback(
            "terminal",
            problem_key="reason-required",
            outcome=outcome,
            **terminal_fields,
        )

    assert collector.events[0]["outcome"] == "not_started"


@pytest.mark.parametrize("outcome", ["off", "applied", "failed", "interrupted"])
def test_non_tolerated_terminal_outcomes_forbid_reason(outcome):
    collector = AssignmentGuidanceCollector()
    mode = "off" if outcome == "off" else "prefer"
    _event(collector, "reason-forbidden", mode=mode)
    fields = {}
    if outcome in {"applied", "failed", "interrupted"}:
        collector.callback(
            "attempt",
            problem_key="reason-forbidden",
            availability="available",
        )

    collector.callback(
        "terminal",
        problem_key="reason-forbidden",
        outcome=outcome,
        **fields,
    )
    assert collector.events[0]["reason"] is None
    assert collector.events[0]["reason_details"] == {}

    other = AssignmentGuidanceCollector()
    _event(other, "invalid-reason", mode=mode)
    if outcome in {"applied", "failed", "interrupted"}:
        other.callback(
            "attempt",
            problem_key="invalid-reason",
            availability="available",
        )
    with pytest.raises(ValueError, match="must not record a reason"):
        other.callback(
            "terminal",
            problem_key="invalid-reason",
            outcome=outcome,
            reason="duplicated_outcome_detail",
        )


@pytest.mark.parametrize(
    ("outcome", "reason"),
    [
        ("fallback", NotApplicableReason.EMPTY_SUPPORT),
        ("not_applicable", FallbackReason.ASSIGNMENT_MISSING),
        ("fallback", "new_unregistered_reason"),
        ("not_applicable", "new_unregistered_reason"),
    ],
)
def test_reason_vocabulary_is_closed_and_outcome_scoped(outcome, reason):
    collector = AssignmentGuidanceCollector()
    _event(
        collector,
        "closed-vocabulary",
        applicability=(
            "not_applicable"
            if outcome == "not_applicable"
            else "applicable"
        ),
    )
    fields = {"reason": reason}
    if outcome == "fallback":
        fields["availability"] = "unavailable"

    with pytest.raises(ValueError, match="Reason"):
        collector.callback(
            "terminal",
            problem_key="closed-vocabulary",
            outcome=outcome,
            **fields,
        )


@pytest.mark.parametrize(
    "case",
    [
        "off_attempt",
        "off_wrong_terminal",
        "off_wrong_availability",
        "not_applicable_attempt",
        "not_applicable_wrong_terminal",
        "not_applicable_wrong_availability",
        "attempt_unavailable",
        "applied_without_attempt",
        "fallback_without_unavailable",
        "fallback_after_attempt",
        "require_fallback",
        "prefer_unavailable_failed",
        "interrupted_without_attempt",
        "unknown_field",
    ],
)
def test_invalid_collector_transition_is_rejected_atomically(case):
    collector = AssignmentGuidanceCollector()
    if case.startswith("off_"):
        _event(collector, "problem", mode="off")
    elif case.startswith("not_applicable"):
        _event(collector, "problem", applicability="not_applicable")
    elif case in {"require_fallback"}:
        _event(collector, "problem", mode="require")
    else:
        _event(collector, "problem")

    if case in {"fallback_after_attempt"}:
        collector.callback("attempt", problem_key="problem", availability="available")

    before = copy.deepcopy(collector.events)
    action, fields = {
        "off_attempt": ("attempt", {"availability": "available"}),
        "off_wrong_terminal": (
            "terminal",
            {
                "outcome": "fallback",
                "availability": "unavailable",
                "reason": FallbackReason.ASSIGNMENT_MISSING,
            },
        ),
        "off_wrong_availability": (
            "terminal",
            {"outcome": "off", "availability": "unavailable"},
        ),
        "not_applicable_attempt": ("attempt", {"availability": "available"}),
        "not_applicable_wrong_terminal": (
            "terminal",
            {"outcome": "off"},
        ),
        "not_applicable_wrong_availability": (
            "terminal",
            {
                "outcome": "not_applicable",
                "availability": "available",
                "reason": NotApplicableReason.EMPTY_SUPPORT,
            },
        ),
        "attempt_unavailable": ("attempt", {"availability": "unavailable"}),
        "applied_without_attempt": ("terminal", {"outcome": "applied"}),
        "fallback_without_unavailable": (
            "terminal",
            {
                "outcome": "fallback",
                "reason": FallbackReason.ASSIGNMENT_MISSING,
            },
        ),
        "fallback_after_attempt": (
            "terminal",
            {
                "outcome": "fallback",
                "availability": "unavailable",
                "reason": FallbackReason.ASSIGNMENT_MISSING,
            },
        ),
        "require_fallback": (
            "terminal",
            {
                "outcome": "fallback",
                "availability": "unavailable",
                "reason": FallbackReason.ASSIGNMENT_MISSING,
            },
        ),
        "prefer_unavailable_failed": (
            "terminal",
            {"outcome": "failed", "availability": "unavailable"},
        ),
        "interrupted_without_attempt": (
            "terminal",
            {"outcome": "interrupted"},
        ),
        "unknown_field": (
            "terminal",
            {
                "outcome": "fallback",
                "availability": "unavailable",
                "reason": FallbackReason.ASSIGNMENT_MISSING,
                "unknown": True,
            },
        ),
    }[case]

    with pytest.raises((TypeError, ValueError)):
        collector.callback(action, problem_key="problem", **fields)

    assert collector.events == before


def test_require_unavailable_failed_and_attempted_solver_failure_are_valid():
    unavailable = AssignmentGuidanceCollector()
    _event(unavailable, "missing", mode="require")
    unavailable.callback(
        "terminal",
        problem_key="missing",
        outcome="failed",
        availability="unavailable",
    )
    solver = AssignmentGuidanceCollector()
    _event(solver, "solver")
    solver.callback("attempt", problem_key="solver", availability="available")
    solver.callback(
        "terminal",
        problem_key="solver",
        outcome="failed",
    )

    assert unavailable.events[0]["attempted"] is False
    assert unavailable.events[0]["outcome"] == "failed"
    assert solver.events[0]["availability"] == "available"
    assert solver.events[0]["outcome"] == "failed"


def test_run_termination_atomically_closes_pre_attempt_and_attempted_open_events():
    collector = AssignmentGuidanceCollector()
    _event(collector, "pre-attempt")
    _event(collector, "attempted", mode="require")
    collector.callback(
        "attempt",
        problem_key="attempted",
        availability="available",
    )

    collector.terminate_open(outcome="failed")

    assert [event["outcome"] for event in collector.events] == [
        "failed",
        "failed",
    ]
    assert [event["attempted"] for event in collector.events] == [False, True]
    assert [event["availability"] for event in collector.events] == [
        "not_checked",
        "available",
    ]
    assert {event["reason"] for event in collector.events} == {None}
    assert {tuple(event["reason_details"]) for event in collector.events} == {
        ()
    }


@pytest.mark.parametrize("outcome", ["applied", "fallback", "not_started"])
def test_invalid_run_termination_does_not_mutate_events(outcome):
    collector = AssignmentGuidanceCollector()
    _event(collector, "open")
    before = copy.deepcopy(collector.events)

    with pytest.raises(ValueError):
        collector.terminate_open(outcome=outcome)

    assert collector.events == before


def test_run_termination_validates_all_proposals_before_batch_commit():
    collector = AssignmentGuidanceCollector()
    _event(collector, "valid")
    _event(collector, "corrupt")
    collector.events[1]["mode"] = "invalid"
    before = copy.deepcopy(collector.events)

    with pytest.raises(ValueError, match="mode"):
        collector.terminate_open(
            outcome="interrupted",
        )

    assert collector.events == before


def test_runner_strategy_injects_context_callback_and_real_route_key(
    monkeypatch,
    request,
):
    backend_package = sys.modules["revise.backend"]
    missing = object()
    previous_module = sys.modules.pop("revise.backend.adapters", missing)
    previous_attribute = getattr(backend_package, "adapters", missing)
    if previous_attribute is not missing:
        delattr(backend_package, "adapters")

    def restore_adapters_module():
        sys.modules.pop("revise.backend.adapters", None)
        if previous_module is not missing:
            sys.modules["revise.backend.adapters"] = previous_module
        if previous_attribute is not missing:
            setattr(backend_package, "adapters", previous_attribute)
        elif hasattr(backend_package, "adapters"):
            delattr(backend_package, "adapters")

    request.addfinalizer(restore_adapters_module)
    scanpy = types.ModuleType("scanpy")
    scanpy.pp = SimpleNamespace()
    scanpy.pl = SimpleNamespace()
    scanpy.tl = SimpleNamespace()
    monkeypatch.setitem(sys.modules, "scanpy", scanpy)
    from revise.backend.adapters import RunnerBackedStrategy

    def callback(*args, **kwargs):
        return None

    observed = []
    ctx = SimpleNamespace(
        assignment_guidance_callback=callback,
        route_key="sp_svc:bin2cell",
        runner_config=SimpleNamespace(),
    )
    ctx.runner = SimpleNamespace(
        local_refinement=lambda: observed.append(
            (
                ctx.runner_config.assignment_guidance_callback,
                ctx.runner_config.assignment_guidance_route,
            )
        )
    )

    class Strategy(RunnerBackedStrategy):
        def prepare_context(self, ctx):
            raise NotImplementedError

        def finalize_svc(self, ctx):
            raise NotImplementedError

    Strategy().solve_ot(ctx)

    assert observed == [(callback, "sp_svc:bin2cell")]
