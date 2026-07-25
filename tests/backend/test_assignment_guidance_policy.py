from __future__ import annotations

import copy
import hashlib
import json
import sys
import types
from types import SimpleNamespace

import numpy as np
import pytest

from revise.backend.ops.assignment import (
    AssignmentState,
    AssignmentStateError,
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
    assignment_compatibility,
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


def test_policy_off_does_not_read_state_and_invalid_prefer_or_require_is_structured():
    reads = []

    def invalid_loader():
        reads.append(True)
        raise AssignmentStateError("values_negative")

    off = resolve_assignment_guidance("off", invalid_loader)
    prefer = resolve_assignment_guidance("prefer", invalid_loader)
    require = resolve_assignment_guidance("require", invalid_loader)

    assert reads == [True, True]
    assert (off.availability, off.outcome, off.reason) == (
        "not_checked",
        "off",
        "guidance_off",
    )
    assert (prefer.availability, prefer.outcome, prefer.reason) == (
        "unavailable",
        "fallback",
        "invalid_assignment_values",
    )
    assert (require.availability, require.outcome, require.reason) == (
        "unavailable",
        "failed",
        "invalid_assignment_values",
    )


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
        reason="solver_failed",
    )
    _event(collector, "unavailable")
    collector.callback(
        "terminal",
        problem_key="unavailable",
        availability="unavailable",
        outcome="fallback",
        reason="assignment_state_unavailable",
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
        reason="guidance_off",
    )
    assert off.summary() == "off"

    inapplicable = AssignmentGuidanceCollector()
    _event(inapplicable, "no-local-problem", applicability="not_applicable")
    inapplicable.callback(
        "terminal",
        problem_key="no-local-problem",
        outcome="not_applicable",
        reason="no_local_problem",
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
            reason=None if outcome == "applied" else "assignment_missing",
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
        reason="run_interrupted",
    )
    assert [event["outcome"] for event in interrupted.events] == [
        "applied",
        "interrupted",
    ]
    assert interrupted.summary() == "interrupted"


@pytest.mark.parametrize(
    "outcome",
    ["off", "not_applicable", "fallback", "failed", "interrupted"],
)
def test_non_applied_terminal_outcomes_require_a_stable_reason_code(outcome):
    collector = AssignmentGuidanceCollector()
    if outcome == "off":
        _event(collector, "reason-required", mode="off")
        terminal_fields = {}
    elif outcome == "not_applicable":
        _event(
            collector,
            "reason-required",
            applicability="not_applicable",
        )
        terminal_fields = {}
    elif outcome == "fallback":
        _event(collector, "reason-required")
        terminal_fields = {"availability": "unavailable"}
    elif outcome == "failed":
        _event(collector, "reason-required", mode="require")
        terminal_fields = {"availability": "unavailable"}
    else:
        _event(collector, "reason-required")
        collector.callback(
            "attempt",
            problem_key="reason-required",
            availability="available",
        )
        terminal_fields = {}

    with pytest.raises(ValueError, match="reason"):
        collector.callback(
            "terminal",
            problem_key="reason-required",
            outcome=outcome,
            **terminal_fields,
        )

    assert collector.events[0]["outcome"] == "not_started"


def test_applied_terminal_outcome_allows_no_reason():
    collector = AssignmentGuidanceCollector()
    _event(collector, "applied")
    collector.callback("attempt", problem_key="applied", availability="available")

    collector.callback(
        "terminal",
        problem_key="applied",
        outcome="applied",
    )

    assert collector.events[0]["reason"] is None


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
            {"outcome": "fallback", "availability": "unavailable", "reason": "assignment_missing"},
        ),
        "off_wrong_availability": (
            "terminal",
            {"outcome": "off", "availability": "unavailable", "reason": "guidance_off"},
        ),
        "not_applicable_attempt": ("attempt", {"availability": "available"}),
        "not_applicable_wrong_terminal": (
            "terminal",
            {"outcome": "off", "reason": "guidance_off"},
        ),
        "not_applicable_wrong_availability": (
            "terminal",
            {"outcome": "not_applicable", "availability": "available", "reason": "no_local_problem"},
        ),
        "attempt_unavailable": ("attempt", {"availability": "unavailable"}),
        "applied_without_attempt": ("terminal", {"outcome": "applied"}),
        "fallback_without_unavailable": (
            "terminal",
            {"outcome": "fallback", "reason": "assignment_missing"},
        ),
        "fallback_after_attempt": (
            "terminal",
            {"outcome": "fallback", "availability": "unavailable", "reason": "assignment_missing"},
        ),
        "require_fallback": (
            "terminal",
            {"outcome": "fallback", "availability": "unavailable", "reason": "assignment_missing"},
        ),
        "prefer_unavailable_failed": (
            "terminal",
            {"outcome": "failed", "availability": "unavailable", "reason": "assignment_missing"},
        ),
        "interrupted_without_attempt": (
            "terminal",
            {"outcome": "interrupted", "reason": "solver_interrupted"},
        ),
        "unknown_field": (
            "terminal",
            {"outcome": "fallback", "availability": "unavailable", "reason": "assignment_missing", "unknown": True},
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
        reason="assignment_missing",
    )
    solver = AssignmentGuidanceCollector()
    _event(solver, "solver")
    solver.callback("attempt", problem_key="solver", availability="available")
    solver.callback(
        "terminal",
        problem_key="solver",
        outcome="failed",
        reason="solver_failed",
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

    collector.terminate_open(outcome="failed", reason="upstream_failure")

    assert [event["outcome"] for event in collector.events] == [
        "failed",
        "failed",
    ]
    assert [event["attempted"] for event in collector.events] == [False, True]
    assert [event["availability"] for event in collector.events] == [
        "not_checked",
        "available",
    ]
    assert {event["reason"] for event in collector.events} == {
        "upstream_failure"
    }


@pytest.mark.parametrize(
    ("outcome", "reason"),
    [
        ("applied", "upstream_failure"),
        ("fallback", "upstream_failure"),
        ("not_started", "upstream_failure"),
        ("failed", None),
        ("interrupted", ""),
    ],
)
def test_invalid_run_termination_does_not_mutate_events(outcome, reason):
    collector = AssignmentGuidanceCollector()
    _event(collector, "open")
    before = copy.deepcopy(collector.events)

    with pytest.raises(ValueError):
        collector.terminate_open(outcome=outcome, reason=reason)

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
            reason="run_interrupted",
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
