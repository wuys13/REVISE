"""Portable contract for host-provided Global Anchoring runs."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Protocol

import numpy as np


ALLOWED_GA_METADATA = frozenset({
    "host_adapter",
    "backend",
    "route",
    "normalization",
    "reconstruction_config_sha256",
    "st_input_sha256",
    "reference_sha256",
    "effective_seed",
    "effective_solver",
    "effective_parameters_sha256",
    "st_axis_sha256",
    "cell_type_labels",
    "reference_n_obs",
    "reference_n_vars",
    "scoring_n_vars",
    "scoring_genes_sha256",
})


@dataclass(frozen=True)
class GlobalAnchoringResult:
    distribution: np.ndarray
    st_unit_ids: list[str]
    cell_type_labels: list[str]
    metadata: dict[str, Any]


@dataclass(frozen=True)
class GAResponse:
    """A host result plus the independently prepared ST-unit axis it must match."""

    result: GlobalAnchoringResult
    expected_st_unit_ids: list[str]


class GARunner(Protocol):
    def __call__(
        self,
        reference_path: Path,
        reconstruction_config: Path,
    ) -> GAResponse: ...


def _require_nonempty_string_axis(name: str, values: object) -> list[str]:
    if not isinstance(values, list) or not values:
        raise ValueError(f"Global Anchoring {name} must be a nonempty list")
    if any(not isinstance(value, str) or not value for value in values):
        raise ValueError(f"Global Anchoring {name} must contain only nonempty string values")
    if len(set(values)) != len(values):
        raise ValueError(f"Global Anchoring {name} must be unique")
    return values


def _controlled_metadata(metadata: object) -> dict[str, Any]:
    if not isinstance(metadata, dict):
        raise ValueError("Global Anchoring metadata must be a dictionary")
    controlled: dict[str, Any] = {}
    for key, value in metadata.items():
        # Existing callbacks may attach runtime objects.  They remain usable,
        # but only this explicit provenance vocabulary crosses into reports.
        if key not in ALLOWED_GA_METADATA:
            continue
        if not isinstance(key, str):
            raise ValueError("Global Anchoring metadata keys must be strings")
        if value is None or type(value) in (str, int, float, bool):
            if isinstance(value, float) and not np.isfinite(value):
                raise ValueError(f"Global Anchoring metadata {key} must be finite")
            controlled[key] = value
            continue
        if isinstance(value, list) and all(
            item is None or type(item) in (str, int, float, bool) for item in value
        ):
            if any(isinstance(item, float) and not np.isfinite(item) for item in value):
                raise ValueError(f"Global Anchoring metadata {key} must be finite")
            controlled[key] = list(value)
            continue
        raise ValueError(
            f"Global Anchoring metadata {key} must be a JSON scalar or flat scalar list"
        )
    controlled["normalization"] = "adapter_declared_probability"
    return controlled


def validate_global_anchoring_result(
    result: GlobalAnchoringResult,
    expected_st_unit_ids: list[str],
) -> GlobalAnchoringResult:
    """Validate a probability matrix without altering its scientific content."""
    if not isinstance(result, GlobalAnchoringResult):
        raise TypeError("Global Anchoring must return a GlobalAnchoringResult")

    result_ids = _require_nonempty_string_axis("ST-unit IDs", result.st_unit_ids)
    expected_ids = _require_nonempty_string_axis("expected ST-unit IDs", expected_st_unit_ids)
    labels = _require_nonempty_string_axis("cell-type labels", result.cell_type_labels)
    if len(labels) < 2:
        raise ValueError("Global Anchoring requires at least two cell-type columns")

    matrix = np.asarray(result.distribution, dtype=float)
    if matrix.ndim != 2:
        raise ValueError("Global Anchoring distribution must be two-dimensional")
    if matrix.shape != (len(result_ids), len(labels)):
        raise ValueError("Global Anchoring matrix shape does not match its axes")
    if result_ids != expected_ids:
        raise ValueError("Global Anchoring ST-unit IDs do not match the SP input order")
    if not np.isfinite(matrix).all() or (matrix < 0).any():
        raise ValueError("Global Anchoring distribution must be finite and non-negative")
    if not np.allclose(matrix.sum(axis=1), 1.0, rtol=1e-5, atol=1e-8):
        raise ValueError(
            "Global Anchoring distribution must be row-normalized; "
            "the outer layer does not silently repair it"
        )

    return replace(result, distribution=matrix, metadata=_controlled_metadata(result.metadata))


__all__ = [
    "ALLOWED_GA_METADATA",
    "GAResponse",
    "GARunner",
    "GlobalAnchoringResult",
    "validate_global_anchoring_result",
]
