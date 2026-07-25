from __future__ import annotations

import copy
import re
from dataclasses import dataclass
from typing import Any, Callable

import numpy as np

from revise.backend.ops.assignment import (
    AssignmentState,
    AssignmentStateError,
    align_assignment_categories,
    assignment_state_evidence,
    validate_assignment,
)
from revise.backend.ops.posterior_conditioning import (
    condition_cost_matrix,
    neighbor_posterior_affinity,
    posterior_affinity,
    posterior_conditioning_mode,
    posterior_conditioning_strict,
)

_REASON_CODE = re.compile(r"^[a-z][a-z0-9_]*$")


@dataclass(frozen=True)
class GuidanceResolution:
    state: AssignmentState | None
    availability: str
    outcome: str
    reason: str | None


def assignment_guidance_mode(config) -> str:
    """Read the canonical policy, with one legacy fallback for direct runners."""
    canonical = getattr(config, "assignment_guidance_policy", None)
    if canonical is not None:
        mode = str(canonical)
        if mode not in {"off", "prefer", "require"}:
            raise ValueError(
                "assignment guidance policy must be off, prefer, or require"
            )
        return mode
    if posterior_conditioning_mode(config) == "off":
        return "off"
    return "require" if posterior_conditioning_strict(config) else "prefer"


def resolve_assignment_guidance(
    mode: str,
    state_loader: Callable[[], AssignmentState | None],
) -> GuidanceResolution:
    if mode == "off":
        return GuidanceResolution(None, "not_checked", "off", "guidance_off")
    if mode not in {"prefer", "require"}:
        raise ValueError("assignment guidance mode must be off, prefer, or require")
    try:
        state = state_loader()
        if state is None:
            raise AssignmentStateError("assignment_state_unavailable")
        state = validate_assignment(state)
    except (AssignmentStateError, KeyError) as exc:
        detail_reason = str(getattr(exc, "reason", "invalid_assignment_values"))
        if detail_reason.startswith("observation_"):
            reason = "observation_axis_mismatch"
        elif detail_reason.startswith("category_"):
            reason = "category_axis_mismatch"
        elif isinstance(exc, KeyError) or detail_reason == "assignment_state_unavailable":
            reason = "assignment_missing"
        elif detail_reason.startswith("projection_"):
            reason = "projection_unavailable"
        else:
            reason = "invalid_assignment_values"
        return GuidanceResolution(
            None,
            "unavailable",
            "fallback" if mode == "prefer" else "failed",
            reason,
        )
    return GuidanceResolution(state, "available", "not_started", None)


class AssignmentGuidanceCollector:
    """One mutable invocation log; callers interact through ``callback`` only."""

    def __init__(
        self,
        *,
        request_evidence: dict[str, Any] | None = None,
        resolved_request: dict[str, Any] | None = None,
    ) -> None:
        self.events: list[dict[str, Any]] = []
        evidence = copy.deepcopy(request_evidence or {}).get(
            "assignment_guidance", {}
        )
        compatibility = copy.deepcopy(resolved_request or {}).get(
            "compatibility", {}
        )
        self.configured = {
            "guidance": evidence.get("configured_guidance"),
            "compatibility_mode": evidence.get(
                "configured_compatibility_mode"
            ),
            "source": evidence.get("resolution_source"),
            "deprecations": list(evidence.get("deprecations", [])),
        }
        self.resolved = {
            "guidance": (resolved_request or {}).get("guidance"),
            "compatibility_mode": compatibility.get("mode"),
            "beta": compatibility.get("beta"),
            "min_affinity": compatibility.get("min_affinity"),
            "operator_strength": compatibility.get("strength"),
        }

    def callback(
        self,
        action: str,
        *,
        problem_key: str,
        **fields: Any,
    ) -> dict[str, Any]:
        if action == "start":
            if any(event["problem_key"] == problem_key for event in self.events):
                raise ValueError(f"duplicate assignment guidance problem_key: {problem_key}")
            applicability = fields.pop("applicability")
            mode = str(fields.pop("mode"))
            if applicability not in {"applicable", "not_applicable"}:
                raise ValueError(f"invalid operator applicability: {applicability}")
            if mode not in {"off", "prefer", "require"}:
                raise ValueError(f"invalid assignment guidance mode: {mode}")
            event = {
                "ordinal": len(self.events) + 1,
                "problem_key": str(problem_key),
                "route": str(fields.pop("route")),
                "operator": str(fields.pop("operator")),
                "phase": str(fields.pop("phase")),
                "mode": mode,
                "applicability": applicability,
                "requested": mode != "off",
                "availability": "not_checked",
                "attempted": False,
                "outcome": "not_started",
                "numerics": copy.deepcopy(fields.pop("numerics", {})),
                "solver": fields.pop("solver", None),
                "reason": None,
                "left_assignment": self._assignment_evidence(
                    fields.pop("left_assignment", None)
                ),
                "right_assignment": self._assignment_evidence(
                    fields.pop("right_assignment", None)
                ),
            }
            if fields:
                raise TypeError(f"unknown assignment guidance fields: {sorted(fields)}")
            self.events.append(event)
            return copy.deepcopy(event)

        event_index = next(
            (
                index
                for index, item in enumerate(self.events)
                if item["problem_key"] == problem_key
            ),
            None,
        )
        if event_index is None:
            raise KeyError(f"unknown assignment guidance problem_key: {problem_key}")
        candidate = copy.deepcopy(self.events[event_index])
        if candidate["outcome"] != "not_started":
            raise ValueError(f"assignment guidance event {problem_key!r} is terminal")
        if action == "attempt":
            if candidate["attempted"]:
                raise ValueError(
                    f"assignment guidance event {problem_key!r} was already attempted"
                )
            if (
                candidate["applicability"] != "applicable"
                or candidate["mode"] not in {"prefer", "require"}
            ):
                raise ValueError(
                    "only applicable prefer/require guidance can be attempted"
                )
            availability = fields.pop("availability", "available")
            if availability != "available":
                raise ValueError("attempted assignment guidance must be available")
            candidate["availability"] = availability
            candidate["attempted"] = True
        elif action == "terminal":
            outcome = fields.pop("outcome")
            if outcome not in {
                "not_applicable",
                "off",
                "applied",
                "fallback",
                "failed",
                "interrupted",
            }:
                raise ValueError(f"invalid assignment guidance outcome: {outcome}")
            if "availability" in fields:
                availability = fields.pop("availability")
                if availability not in {
                    "not_checked",
                    "available",
                    "unavailable",
                }:
                    raise ValueError(
                        f"invalid assignment availability: {availability}"
                    )
                candidate["availability"] = availability
            reason = fields.get("reason", candidate["reason"])
            if outcome != "applied" and not self._valid_reason_code(reason):
                raise ValueError(
                    f"{outcome} assignment guidance requires a stable reason code"
                )
            if reason is not None and not self._valid_reason_code(reason):
                raise ValueError("assignment guidance reason must be a stable code")
            candidate["outcome"] = outcome
        else:
            raise ValueError(f"unknown assignment guidance action: {action}")
        for key in ("reason", "numerics", "solver"):
            if key in fields:
                candidate[key] = copy.deepcopy(fields.pop(key))
        for key in ("left_assignment", "right_assignment"):
            if key in fields:
                candidate[key] = self._assignment_evidence(fields.pop(key))
        if fields:
            raise TypeError(f"unknown assignment guidance fields: {sorted(fields)}")
        self._validate_transition(candidate, action=action)
        self.events[event_index] = candidate
        return copy.deepcopy(candidate)

    @staticmethod
    def _assignment_evidence(
        state: AssignmentState | None,
    ) -> dict[str, Any] | None:
        if state is None:
            return None
        return assignment_state_evidence(state)

    @staticmethod
    def _valid_reason_code(reason: Any) -> bool:
        return isinstance(reason, str) and bool(_REASON_CODE.fullmatch(reason))

    @staticmethod
    def _validate_transition(event: dict[str, Any], *, action: str) -> None:
        if action == "attempt":
            return
        outcome = event["outcome"]
        attempted = bool(event["attempted"])
        availability = event["availability"]
        mode = event["mode"]
        applicability = event["applicability"]

        if applicability == "not_applicable":
            if (
                attempted
                or availability != "not_checked"
                or outcome != "not_applicable"
            ):
                raise ValueError(
                    "not-applicable guidance must remain not_checked, cannot "
                    "be attempted, and must terminate as not_applicable"
                )
            return
        if mode == "off":
            if attempted or availability != "not_checked" or outcome != "off":
                raise ValueError(
                    "off guidance must remain not_checked, cannot be attempted, "
                    "and must terminate as off"
                )
            return
        if outcome == "applied":
            if not attempted or availability != "available":
                raise ValueError(
                    "applied guidance must be attempted and available"
                )
            return
        if outcome == "fallback":
            if mode != "prefer" or attempted or availability != "unavailable":
                raise ValueError(
                    "fallback requires unavailable, unattempted prefer guidance"
                )
            return
        if outcome == "failed":
            solver_failure = attempted and availability == "available"
            required_unavailable = (
                mode == "require"
                and not attempted
                and availability == "unavailable"
            )
            if not (solver_failure or required_unavailable):
                raise ValueError(
                    "failed guidance requires a solver failure or unavailable "
                    "require guidance"
                )
            return
        if outcome == "interrupted":
            if not attempted or availability != "available":
                raise ValueError(
                    "interrupted guidance must be attempted and available"
                )
            return
        raise ValueError(
            f"invalid terminal outcome for applicable guidance: {outcome}"
        )

    def summary(self) -> str:
        outcomes = [event["outcome"] for event in self.events]
        if not outcomes or set(outcomes) == {"not_started"}:
            outcome = "not_started"
        elif "failed" in outcomes:
            outcome = "failed"
        elif "interrupted" in outcomes:
            outcome = "interrupted"
        elif set(outcomes) == {"not_applicable"}:
            outcome = "not_applicable"
        elif "applied" in outcomes and "fallback" in outcomes:
            outcome = "mixed"
        elif len(set(outcomes)) == 1:
            outcome = outcomes[0]
        else:
            outcome = "mixed"
        return outcome

    def manifest(self) -> dict[str, Any]:
        return {
            "schema_version": 1,
            "configured": copy.deepcopy(self.configured),
            "resolved": copy.deepcopy(self.resolved),
            "events": copy.deepcopy(self.events),
            "summary": self.summary(),
        }

    def terminate_open(self, *, outcome: str, reason: str) -> None:
        if outcome not in {"failed", "interrupted"}:
            raise ValueError(
                "run termination can only close guidance as failed or interrupted"
            )
        if not self._valid_reason_code(reason):
            raise ValueError("assignment guidance reason must be a stable code")
        proposed = copy.deepcopy(self.events)
        for event in proposed:
            if event["outcome"] == "not_started":
                self._validate_event_vocabulary(event)
                event["outcome"] = outcome
                event["reason"] = reason
                self._validate_event_vocabulary(event)
        self.events = proposed

    @classmethod
    def _validate_event_vocabulary(cls, event: dict[str, Any]) -> None:
        if event.get("mode") not in {"off", "prefer", "require"}:
            raise ValueError("assignment guidance event has invalid mode")
        if event.get("applicability") not in {
            "applicable",
            "not_applicable",
        }:
            raise ValueError(
                "assignment guidance event has invalid applicability"
            )
        if event.get("availability") not in {
            "not_checked",
            "available",
            "unavailable",
        }:
            raise ValueError("assignment guidance event has invalid availability")
        if type(event.get("attempted")) is not bool:
            raise ValueError("assignment guidance event has invalid attempted flag")
        if event.get("outcome") not in {
            "not_started",
            "not_applicable",
            "off",
            "applied",
            "fallback",
            "failed",
            "interrupted",
        }:
            raise ValueError("assignment guidance event has invalid outcome")
        reason = event.get("reason")
        if reason is not None and not cls._valid_reason_code(reason):
            raise ValueError("assignment guidance event has invalid reason code")

    def snapshot(self) -> list[dict[str, Any]]:
        return copy.deepcopy(self.events)

    def restore(self, events: list[dict[str, Any]]) -> None:
        self.events = copy.deepcopy(events)


def assignment_compatibility(
    left: AssignmentState,
    right: AssignmentState,
    *,
    beta: float = 1.0,
    min_affinity: float = 1e-3,
    support: Any = None,
) -> np.ndarray:
    left = validate_assignment(left)
    right = align_assignment_categories(right, left.category_labels)
    if support is None:
        return posterior_affinity(
            left.values,
            right.values,
            beta=beta,
            min_affinity=min_affinity,
        )
    elif isinstance(support, tuple) and len(support) == 2:
        rows = _support_indices(
            support[0],
            name="edge row support",
            ndim=1,
            bound=left.values.shape[0],
        )
        columns = _support_indices(
            support[1],
            name="edge column support",
            ndim=1,
            bound=right.values.shape[0],
        )
        if rows.shape != columns.shape:
            raise ValueError("edge support rows and columns must have equal length")
        affinity = np.einsum(
            "ij,ij->i",
            left.values[rows],
            right.values[columns],
        )
    else:
        columns = _support_indices(
            support,
            name="top-k support",
            ndim=2,
            bound=right.values.shape[0],
        )
        if columns.shape[0] != left.values.shape[0]:
            raise ValueError("top-k support must have shape (n_left, k)")
        return neighbor_posterior_affinity(
            left.values,
            columns,
            q_neighbors=right.values,
            beta=beta,
            min_affinity=min_affinity,
        )
    return np.power(
        np.clip(affinity, float(min_affinity), 1.0),
        float(beta),
    )


def _support_indices(
    values: Any,
    *,
    name: str,
    ndim: int,
    bound: int,
) -> np.ndarray:
    indices = np.asarray(values)
    if indices.ndim != ndim:
        raise ValueError(f"{name} must be {ndim}D")
    if not np.issubdtype(indices.dtype, np.integer):
        raise ValueError(f"{name} must use an integer dtype")
    if np.any(indices < 0):
        raise ValueError(f"{name} cannot contain negative indices")
    if np.any(indices >= bound):
        raise ValueError(f"{name} index is out of bounds")
    return indices


def ot_cost_guidance(
    cost: np.ndarray,
    affinity: np.ndarray,
    strength: float,
) -> np.ndarray:
    cost = np.asarray(cost, dtype=np.float64)
    affinity = np.asarray(affinity, dtype=np.float64)
    if cost.shape != affinity.shape:
        raise ValueError("cost and affinity shapes differ")
    return condition_cost_matrix(cost, affinity, strength)


def graph_guidance(
    weights: np.ndarray,
    affinity: np.ndarray,
    strength: float,
) -> np.ndarray:
    weights = np.asarray(weights, dtype=np.float64)
    affinity = np.asarray(affinity, dtype=np.float64)
    if weights.shape != affinity.shape:
        raise ValueError("weights and affinity shapes differ")
    return weights * np.power(affinity, float(strength))
