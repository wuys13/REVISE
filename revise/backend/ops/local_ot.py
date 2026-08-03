from __future__ import annotations

import warnings
from typing import Any

import numpy as np

from revise.backend.ops.tacco_runtime import require_tacco


_TACCO_MARGINAL_CONTINUATION_MAX_CYCLES = 1000
_TACCO_FEASIBILITY_ATOL = 1e-10


class _TACCOMarginalContinuationWarning(UserWarning):
    """A verified TACCO coupling needed more marginal-scaling cycles."""


def _validated_inputs(source_mass, target_mass, cost_matrix):
    # Preserve the caller's dtype for exact POT compatibility. TACCO promotes
    # marginals explicitly when they are normalized in its dedicated branch.
    source = np.asarray(source_mass).ravel()
    target = np.asarray(target_mass).ravel()
    cost = np.asarray(cost_matrix)

    if cost.shape != (source.size, target.size):
        raise ValueError(
            "Local OT cost shape must equal (source size, target size): "
            f"{cost.shape} vs {(source.size, target.size)}"
        )
    for name, values in (("source mass", source), ("target mass", target)):
        if not np.all(np.isfinite(values)):
            raise ValueError(f"Local OT {name} must contain only finite values")
    if np.isnan(cost).any() or np.isneginf(cost).any():
        raise ValueError("Local OT cost must not contain NaN or negative infinity")
    if np.any(source < 0) or np.any(target < 0):
        raise ValueError("Local OT marginals must be non-negative")
    if float(source.sum()) <= 0 or float(target.sum()) <= 0:
        raise ValueError("Local OT marginals must each have positive total mass")
    return source, target, cost


def _validated_coupling(coupling: Any, expected_shape: tuple[int, int]) -> np.ndarray:
    result = np.asarray(coupling)
    if result.shape != expected_shape:
        raise ValueError(
            f"Local OT coupling has shape {result.shape}, expected {expected_shape}"
        )
    if not np.all(np.isfinite(result)):
        raise ValueError("Local OT coupling must contain only finite values")
    if np.any(result < 0):
        raise ValueError("Local OT coupling must be non-negative")
    return result


def _tacco_marginal_errors(
    coupling: np.ndarray,
    source: np.ndarray,
    target: np.ndarray,
) -> dict[str, float]:
    promoted = np.asarray(coupling, dtype=np.float64)
    source_error = np.abs(promoted.sum(axis=1) - source)
    target_error = np.abs(promoted.sum(axis=0) - target)
    errors = {
        "total": abs(float(promoted.sum()) - 1.0),
        "source_l1": float(source_error.sum()),
        "source_max": float(source_error.max()),
        "target_l1": float(target_error.sum()),
        "target_max": float(target_error.max()),
    }
    if not np.all(np.isfinite(list(errors.values()))):
        raise ValueError("TACCO marginal validation errors must be finite")
    return errors


def _tacco_marginal_contract_ok(errors: dict[str, float]) -> bool:
    return bool(
        errors["total"] <= 1e-6
        and errors["source_l1"] <= 1e-2
        and errors["source_max"] <= 5e-3
        and errors["target_l1"] <= 1e-6
        and errors["target_max"] <= 1e-6
    )


def _tacco_hard_support_is_feasible(
    source: np.ndarray,
    target: np.ndarray,
    support: np.ndarray,
) -> bool:
    """Numerically check normalized-marginal feasibility on a hard support."""
    if np.all(support):
        return True

    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix

    source_count, target_count = support.shape
    row, column = np.nonzero(support)
    edge_count = row.size
    variables = np.arange(edge_count)
    constraint_rows = np.concatenate((row, source_count + column))
    constraint_columns = np.concatenate((variables, variables))
    equality = coo_matrix(
        (
            np.ones(edge_count * 2, dtype=np.float64),
            (constraint_rows, constraint_columns),
        ),
        shape=(source_count + target_count, edge_count),
    ).tocsr()
    expected = np.concatenate((source, target))
    solved = linprog(
        np.zeros(edge_count, dtype=np.float64),
        A_eq=equality,
        b_eq=expected,
        bounds=(0.0, None),
        method="highs",
        options={"primal_feasibility_tolerance": _TACCO_FEASIBILITY_ATOL},
    )
    status = int(solved.status)
    if status == 2:
        return False
    if not solved.success:
        raise ValueError(
            "TACCO could not establish hard-support feasibility before "
            f"marginal continuation: status={status}, message={solved.message}"
        )

    witness = getattr(solved, "x", None)
    if (
        witness is None
        or np.asarray(witness).shape != (edge_count,)
        or not np.all(np.isfinite(witness))
    ):
        raise ValueError(
            "TACCO hard-support feasibility result failed validation: "
            "the LP witness is absent, malformed, or non-finite"
        )
    residual = np.abs(equality @ witness - expected)
    if not np.all(np.isfinite(residual)):
        raise ValueError(
            "TACCO hard-support feasibility result failed validation: "
            "the equality residual is non-finite"
        )
    residual_max = float(residual.max())
    minimum_mass = float(np.min(witness))
    if residual_max > _TACCO_FEASIBILITY_ATOL or minimum_mass < -1e-12:
        raise ValueError(
            "TACCO hard-support feasibility result failed validation before "
            "marginal continuation: "
            f"max_abs_residual={residual_max:.6g}, "
            f"minimum_edge_mass={minimum_mass:.6g}"
        )
    return True


def _continue_tacco_marginal_scaling(
    coupling: np.ndarray,
    source: np.ndarray,
    target: np.ndarray,
    support: np.ndarray,
) -> tuple[np.ndarray | None, int, str | None]:
    """Continue a feasible TACCO Sinkhorn orbit until its contract is met."""
    balanced = np.asarray(coupling, dtype=np.float64).copy()
    if np.any(balanced[support] <= 0):
        return None, 0, "supported_edge_loss"

    for cycle in range(1, _TACCO_MARGINAL_CONTINUATION_MAX_CYCLES + 1):
        row_sums = balanced.sum(axis=1)
        if not np.all(np.isfinite(row_sums)) or np.any(row_sums <= 0):
            return None, cycle - 1, "invalid_row_sum"
        balanced *= (source / row_sums)[:, None]
        if not np.all(np.isfinite(balanced)):
            return None, cycle, "nonfinite_row_scaling"
        if np.any(balanced[support] <= 0):
            return None, cycle, "supported_edge_loss"

        column_sums = balanced.sum(axis=0)
        if not np.all(np.isfinite(column_sums)) or np.any(column_sums <= 0):
            return None, cycle - 1, "invalid_column_sum"
        balanced *= (target / column_sums)[None, :]
        balanced[~support] = 0.0

        if not np.all(np.isfinite(balanced)):
            return None, cycle, "nonfinite_column_scaling"
        if np.any(balanced[support] <= 0):
            return None, cycle, "supported_edge_loss"
        if _tacco_marginal_contract_ok(
            _tacco_marginal_errors(balanced, source, target)
        ):
            return balanced, cycle, None

    return None, _TACCO_MARGINAL_CONTINUATION_MAX_CYCLES, "cycle_limit"


def stabilize_local_ot_support(
    source_mass: Any, target_mass: Any, valid_support_mask: Any
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Remove zero-mass and unsupported rows and columns from an OT support."""
    source = np.asarray(source_mass).ravel()
    target = np.asarray(target_mass).ravel()
    support = np.asarray(valid_support_mask, dtype=bool)
    expected_shape = (source.size, target.size)
    if support.shape != expected_shape:
        raise ValueError(
            f"valid_support_mask shape {support.shape} != {expected_shape}"
        )

    source_active = source > 0
    target_active = target > 0
    while True:
        current = support & source_active[:, None] & target_active[None, :]
        next_source = source_active & current.any(axis=1)
        next_target = target_active & current.any(axis=0)
        if np.array_equal(next_source, source_active) and np.array_equal(
            next_target, target_active
        ):
            break
        source_active = next_source
        target_active = next_target

    source_idx = np.flatnonzero(source_active)
    target_idx = np.flatnonzero(target_active)
    return source_idx, target_idx, support[np.ix_(source_idx, target_idx)]


def solve_local_ot(
    source_mass,
    target_mass,
    cost_matrix,
    *,
    method: str = "pot",
    pot_reg: float | None = None,
    pot_reg_m: float | None = None,
    pot_reg_type: str = "kl",
    pot_verbose: bool = False,
    pot_num_iter_max: int = 5000,
    reference_measure=None,
    valid_support_mask=None,
) -> np.ndarray:
    """Solve a local coupling with POT or TACCO using one matrix contract."""
    source, target, cost = _validated_inputs(source_mass, target_mass, cost_matrix)
    full_shape = cost.shape
    normalized_method = str(method).strip().lower()

    if normalized_method not in {"pot", "tacco"}:
        raise ValueError(
            f"Unsupported local OT method {method!r}; expected one of ['pot', 'tacco']"
        )
    if normalized_method == "pot" and (pot_reg is None or pot_reg_m is None):
        raise ValueError("POT local OT requires pot_reg and pot_reg_m")
    if normalized_method == "tacco" and reference_measure is not None:
        raise ValueError(
            "TACCO local OT does not support posterior reference-measure conditioning"
        )

    reference = None
    if reference_measure is not None:
        reference = np.asarray(reference_measure)
        if reference.shape != full_shape:
            raise ValueError(
                f"Local OT reference measure has shape {reference.shape}, "
                f"expected {full_shape}"
            )

    if valid_support_mask is None:
        support = np.ones(full_shape, dtype=bool)
        source_idx = np.arange(source.size)
        target_idx = np.arange(target.size)
    else:
        support = np.asarray(valid_support_mask, dtype=bool)
        source_idx, target_idx, _ = stabilize_local_ot_support(
            source, target, support
        )
        if source_idx.size == 0 or target_idx.size == 0:
            return np.zeros(full_shape, dtype=np.float64)

    active_source = source[source_idx]
    active_target = target[target_idx]
    active_cost = cost[np.ix_(source_idx, target_idx)]
    if normalized_method == "tacco" or not np.issubdtype(
        active_cost.dtype, np.floating
    ):
        active_cost = active_cost.astype(np.float64, copy=True)
    else:
        active_cost = active_cost.copy()
    active_support = support[np.ix_(source_idx, target_idx)]
    if not np.all(np.isfinite(active_cost[active_support])):
        raise ValueError("valid local OT costs must be finite")
    if np.any(active_cost[active_support] < 0):
        raise ValueError("valid local OT costs must be non-negative")
    active_cost[~active_support] = np.inf

    reference_active = None
    if reference is not None:
        reference_active = reference[np.ix_(source_idx, target_idx)]

    source_norm = None
    target_norm = None
    if normalized_method == "pot":
        import ot

        kwargs = {}
        if reference_active is not None:
            kwargs["c"] = reference_active
        active_coupling = ot.unbalanced.sinkhorn_unbalanced(
            active_source,
            active_target,
            active_cost,
            reg=pot_reg,
            reg_m=pot_reg_m,
            reg_type=pot_reg_type,
            verbose=pot_verbose,
            numItermax=pot_num_iter_max,
            **kwargs,
        )
    else:
        tacco = require_tacco()

        # Normalize in float64 before TACCO sees the marginals. Normalizing a
        # float32 array first can leave the promoted vectors with unequal totals.
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            source_float = np.asarray(active_source, dtype=np.float64)
            target_float = np.asarray(active_target, dtype=np.float64)
        if (
            not np.all(np.isfinite(source_float))
            or not np.all(np.isfinite(target_float))
            or np.any(source_float <= 0)
            or np.any(target_float <= 0)
        ):
            raise ValueError(
                "TACCO active marginals must remain finite and strictly positive "
                "after float64 conversion"
            )
        source_total = float(source_float.sum())
        target_total = float(target_float.sum())
        if (
            not np.isfinite(source_total)
            or not np.isfinite(target_total)
            or source_total <= 0
            or target_total <= 0
        ):
            raise ValueError(
                "TACCO active marginal totals must be finite and strictly positive"
            )
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            source_norm = source_float / source_total
            target_norm = target_float / target_total
        if (
            not np.all(np.isfinite(source_norm))
            or not np.all(np.isfinite(target_norm))
            or np.any(source_norm <= 0)
            or np.any(target_norm <= 0)
        ):
            raise ValueError(
                "TACCO normalized marginals must remain finite and strictly positive"
            )
        with warnings.catch_warnings(record=True) as caught_warnings:
            warnings.simplefilter("always")
            active_coupling = tacco.utils.solve_OT(
                source_norm,
                target_norm,
                active_cost,
            )
        numerical_warnings = [
            str(item.message)
            for item in caught_warnings
            if "Numerical errors at iteration" in str(item.message)
        ]
        if numerical_warnings:
            raise ValueError(
                f"TACCO local OT failed: {numerical_warnings[0]}"
            )
        runtime_warnings = [
            item
            for item in caught_warnings
            if issubclass(item.category, RuntimeWarning)
        ]
        if runtime_warnings:
            raise ValueError(
                "TACCO local OT failed with RuntimeWarning: "
                f"{runtime_warnings[0].message}"
            )
        for item in caught_warnings:
            warnings.warn(item.message, item.category, stacklevel=2)

    active_coupling = _validated_coupling(active_coupling, active_cost.shape)
    invalid_mass = float(np.abs(active_coupling[~active_support]).sum())
    if invalid_mass > 1e-12:
        raise ValueError(
            "Local OT coupling assigns mass to invalid support: "
            f"{invalid_mass:.6g} > 1e-12"
        )
    active_coupling = active_coupling.copy()
    active_coupling[~active_support] = 0.0

    if normalized_method == "pot":
        total_mass = float(active_coupling.sum())
        if total_mass <= 0:
            raise ValueError("POT local OT coupling must have positive total mass")
        row_mass = active_coupling.sum(axis=1)
        if np.any(row_mass[active_source > 0] <= 0):
            raise ValueError(
                "POT local OT coupling must have positive transported row mass "
                "for every positive source marginal"
            )
        column_mass = active_coupling.sum(axis=0)
        if np.any(column_mass[active_target > 0] <= 0):
            raise ValueError(
                "POT local OT coupling must have positive transported column mass "
                "for every positive target marginal"
            )

    if normalized_method == "tacco":
        errors = _tacco_marginal_errors(
            active_coupling,
            source_norm,
            target_norm,
        )
        initial_errors = errors.copy()
        source_failed = (
            errors["source_l1"] > 1e-2 or errors["source_max"] > 5e-3
        )
        target_passed = (
            errors["target_l1"] <= 1e-6 and errors["target_max"] <= 1e-6
        )
        continuation_failure = None
        continuation_candidate = bool(
            errors["total"] <= 1e-6 and source_failed and target_passed
        )
        if continuation_candidate and continuation_failure is None:
            if not _tacco_hard_support_is_feasible(
                source_norm,
                target_norm,
                active_support,
            ):
                raise ValueError(
                    "TACCO hard-support balanced OT is infeasible for the "
                    "normalized marginals: "
                    f"source_l1={errors['source_l1']:.6g}, "
                    f"source_max_abs={errors['source_max']:.6g}, "
                    f"support_edges={int(active_support.sum())}/"
                    f"{active_support.size}"
                )
            continued, cycles, continuation_failure = (
                _continue_tacco_marginal_scaling(
                    active_coupling,
                    source_norm,
                    target_norm,
                    active_support,
                )
            )
            if continued is not None:
                active_coupling = continued
                errors = _tacco_marginal_errors(
                    active_coupling,
                    source_norm,
                    target_norm,
                )
                warnings.warn(
                    "TACCO default coupling required "
                    f"{cycles} additional marginal-scaling cycles: "
                    f"source_l1={initial_errors['source_l1']:.6g}->"
                    f"{errors['source_l1']:.6g}, "
                    f"source_max={initial_errors['source_max']:.6g}->"
                    f"{errors['source_max']:.6g}, "
                    f"target_l1={errors['target_l1']:.6g}, "
                    f"support_edges={int(active_support.sum())}/"
                    f"{active_support.size}",
                    _TACCOMarginalContinuationWarning,
                    stacklevel=2,
                )

        support_count = f"{int(active_support.sum())}/{active_support.size}"
        if errors["total"] > 1e-6:
            raise ValueError(
                "TACCO coupling total-mass validation failed: "
                f"error={errors['total']:.6g} > 1e-6, "
                f"source_max_abs={errors['source_max']:.6g}, "
                f"target_max_abs={errors['target_max']:.6g}, "
                f"support_edges={support_count}"
            )
        if errors["source_l1"] > 1e-2 or errors["source_max"] > 5e-3:
            raise ValueError(
                "TACCO source marginal validation failed: "
                f"l1={errors['source_l1']:.6g}, "
                f"max_abs={errors['source_max']:.6g}, "
                f"support_edges={support_count}, "
                f"continuation={continuation_failure or 'not_attempted'}"
            )
        if errors["target_l1"] > 1e-6 or errors["target_max"] > 1e-6:
            raise ValueError(
                "TACCO target marginal validation failed: "
                f"l1={errors['target_l1']:.6g}, "
                f"max_abs={errors['target_max']:.6g}, "
                f"support_edges={support_count}"
            )

    coupling = np.zeros(full_shape, dtype=active_coupling.dtype)
    coupling[np.ix_(source_idx, target_idx)] = active_coupling
    return coupling
