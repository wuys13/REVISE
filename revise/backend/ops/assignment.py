from __future__ import annotations

import copy
import hashlib
import json
from dataclasses import dataclass, replace
from typing import Any, Iterable, Literal, Mapping, Sequence

import numpy as np
import pandas as pd


ValueSemantics = Literal["soft", "one_hot"]


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
    raise GlobalAssignmentContractError(
        f"{axis} axis order does not match expected order"
    )


def validate_global_assignment(
    assignment: GlobalAssignment,
    *,
    expected_observations: Iterable[Any],
    expected_categories: Iterable[Any],
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
    if not np.allclose(row_mass, 1.0, rtol=0.0, atol=1e-6):
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


class AssignmentStateError(ValueError):
    def __init__(self, reason: str, detail: str | None = None) -> None:
        self.reason = reason
        super().__init__(reason if detail is None else f"{reason}: {detail}")


@dataclass
class AssignmentState:
    """Legacy guidance carrier pending migration to ``GlobalAssignment``.

    ``validate_assignment`` returns independently owned values and lineage.
    Callers must not mutate or reuse a state after registering it as event
    evidence.
    """

    values: np.ndarray
    observation_labels: Sequence[Any]
    category_labels: Sequence[Any]
    source: str
    level: str
    value_semantics: ValueSemantics
    lineage: list[dict[str, Any]]


def assignment_state_evidence(state: AssignmentState) -> dict[str, Any]:
    """Serialize compact, order-sensitive Assignment State evidence.

    Axis digests hash canonical JSON arrays encoded as UTF-8 with
    ``ensure_ascii=False`` and ``separators=(",", ":")``.
    """
    state = validate_assignment(state)

    def axis_summary(labels: Sequence[Any]) -> dict[str, Any]:
        payload = json.dumps(
            list(labels),
            ensure_ascii=False,
            separators=(",", ":"),
        ).encode("utf-8")
        return {
            "count": len(labels),
            "ordered_labels_sha256": hashlib.sha256(payload).hexdigest(),
        }

    return {
        "source": state.source,
        "level": state.level,
        "value_semantics": state.value_semantics,
        "lineage": copy.deepcopy(state.lineage),
        "observation_axis": axis_summary(state.observation_labels),
        "category_axis": axis_summary(state.category_labels),
    }


def _normalize_axis(labels: Iterable[Any], axis: str) -> tuple[str, ...]:
    values = list(labels)
    if not values:
        raise AssignmentStateError(f"{axis}_labels_empty")
    if bool(pd.isna(pd.Index(values)).any()):
        raise AssignmentStateError(f"{axis}_labels_null")
    text = [str(value) for value in values]
    if any(not value.strip() for value in text):
        raise AssignmentStateError(f"{axis}_labels_empty")
    if len(set(text)) != len(text):
        raise AssignmentStateError(f"{axis}_labels_duplicate")
    normalized = [value.replace("/", "_") for value in text]
    if len(set(normalized)) != len(normalized):
        raise AssignmentStateError(f"{axis}_labels_slash_collision")
    return tuple(normalized)


def _normalize_group_labels(labels: Iterable[Any]) -> tuple[str, ...]:
    values = list(labels)
    if not values:
        raise AssignmentStateError("aggregation_mapping_empty")
    if bool(pd.isna(pd.Index(values)).any()):
        raise AssignmentStateError("aggregation_mapping_null")
    text = [str(value) for value in values]
    if any(not value.strip() for value in text):
        raise AssignmentStateError("aggregation_mapping_empty")
    normalized = [value.replace("/", "_") for value in text]
    originals_by_normalized: dict[str, set[str]] = {}
    for original, label in zip(text, normalized):
        originals_by_normalized.setdefault(label, set()).add(original)
    if any(len(originals) > 1 for originals in originals_by_normalized.values()):
        raise AssignmentStateError("aggregation_mapping_slash_collision")
    return tuple(normalized)


def validate_assignment(state: AssignmentState) -> AssignmentState:
    observations = _normalize_axis(state.observation_labels, "observation")
    categories = _normalize_axis(state.category_labels, "category")
    try:
        values = np.asarray(state.values, dtype=np.float64).copy()
    except (TypeError, ValueError) as exc:
        raise AssignmentStateError("values_not_numeric") from exc
    if values.ndim != 2 or values.shape != (len(observations), len(categories)):
        raise AssignmentStateError("values_shape")
    if np.isnan(values).any():
        raise AssignmentStateError("values_nan")
    if np.isinf(values).any():
        raise AssignmentStateError("values_infinite")
    if np.any(values < 0):
        raise AssignmentStateError("values_negative")
    with np.errstate(over="ignore", invalid="ignore"):
        row_mass = values.sum(axis=1)
    if not np.all(np.isfinite(row_mass)):
        raise AssignmentStateError("values_row_mass_nonfinite")
    if np.any(row_mass <= 0):
        raise AssignmentStateError("values_zero_row")
    if state.value_semantics not in {"soft", "one_hot"}:
        raise AssignmentStateError("value_semantics_invalid")
    if state.value_semantics == "one_hot":
        if not (
            np.all((values == 0.0) | (values == 1.0))
            and np.all(row_mass == 1.0)
        ):
            raise AssignmentStateError("values_not_one_hot")
        normalized_values = values.copy()
    else:
        normalized_values = values / row_mass[:, np.newaxis]
    return replace(
        state,
        values=normalized_values,
        observation_labels=observations,
        category_labels=categories,
        lineage=copy.deepcopy(state.lineage),
    )


def _align_axis(
    state: AssignmentState,
    labels: Iterable[Any],
    *,
    axis: Literal["observation", "category"],
) -> AssignmentState:
    state = validate_assignment(state)
    requested = _normalize_axis(labels, axis)
    current = (
        state.observation_labels if axis == "observation" else state.category_labels
    )
    if set(requested) != set(current):
        raise AssignmentStateError(f"{axis}_labels_mismatch")
    positions = {label: index for index, label in enumerate(current)}
    order = [positions[label] for label in requested]
    if axis == "observation":
        return replace(
            state,
            values=state.values[order],
            observation_labels=requested,
        )
    return replace(
        state,
        values=state.values[:, order],
        category_labels=requested,
    )


def align_assignment_observations(
    state: AssignmentState,
    observation_labels: Iterable[Any],
) -> AssignmentState:
    return _align_axis(state, observation_labels, axis="observation")


def align_assignment_categories(
    state: AssignmentState,
    category_labels: Iterable[Any],
) -> AssignmentState:
    return _align_axis(state, category_labels, axis="category")


def argmax_assignment(state: AssignmentState) -> AssignmentState:
    state = validate_assignment(state)
    values = np.zeros_like(state.values)
    values[np.arange(values.shape[0]), np.argmax(state.values, axis=1)] = 1.0
    return replace(state, values=values, value_semantics="one_hot")


def one_hot_assignment(
    labels: Iterable[Any],
    *,
    observation_labels: Iterable[Any],
    category_labels: Iterable[Any],
    source: str,
    level: str,
    lineage: list[dict[str, Any]] | None = None,
) -> AssignmentState:
    observations = _normalize_axis(observation_labels, "observation")
    categories = _normalize_axis(category_labels, "category")
    assigned = _normalize_group_labels(labels)
    if len(assigned) != len(observations):
        raise AssignmentStateError("values_shape")
    positions = {label: index for index, label in enumerate(categories)}
    if any(label not in positions for label in assigned):
        raise AssignmentStateError("category_axis_mismatch")
    values = np.zeros((len(observations), len(categories)), dtype=np.float64)
    values[
        np.arange(len(observations)),
        [positions[label] for label in assigned],
    ] = 1.0
    return validate_assignment(
        AssignmentState(
            values=values,
            observation_labels=observations,
            category_labels=categories,
            source=source,
            level=level,
            value_semantics="one_hot",
            lineage=list(lineage or []),
        )
    )


def aggregate_assignment(
    state: AssignmentState,
    group_labels: Iterable[Any],
    *,
    source: str,
    level: str,
) -> AssignmentState:
    state = validate_assignment(state)
    groups = _normalize_group_labels(group_labels)
    if len(groups) != len(state.observation_labels):
        raise AssignmentStateError("aggregation_mapping_shape")
    ordered = tuple(dict.fromkeys(groups))
    group_array = np.asarray(groups)
    values = np.vstack(
        [state.values[group_array == group].mean(axis=0) for group in ordered]
    )
    return validate_assignment(
        AssignmentState(
            values=values,
            observation_labels=ordered,
            category_labels=state.category_labels,
            source=source,
            level=level,
            value_semantics="soft",
            lineage=[
                *state.lineage,
                {
                    "operation": "aggregate",
                    "source": state.source,
                    "from_level": state.level,
                    "to_level": level,
                },
            ],
        )
    )


def project_assignment(
    state: AssignmentState,
    target_to_source: Mapping[Any, Any],
    *,
    source: str,
    level: str,
) -> AssignmentState:
    state = validate_assignment(state)
    targets = _normalize_axis(target_to_source.keys(), "observation")
    sources = _normalize_group_labels(target_to_source.values())
    positions = {label: index for index, label in enumerate(state.observation_labels)}
    if any(label not in positions for label in sources):
        raise AssignmentStateError("projection_source_missing")
    return validate_assignment(
        AssignmentState(
            values=state.values[[positions[label] for label in sources]],
            observation_labels=targets,
            category_labels=state.category_labels,
            source=source,
            level=level,
            value_semantics=state.value_semantics,
            lineage=[
                *state.lineage,
                {
                    "operation": "project",
                    "source": state.source,
                    "from_level": state.level,
                    "to_level": level,
                },
            ],
        )
    )
