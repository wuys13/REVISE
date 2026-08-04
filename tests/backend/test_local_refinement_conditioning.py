from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from revise.backend.ops.assignment import (
    GlobalAssignment,
    GlobalAssignmentContractError,
)
from revise.backend.ops.posterior_conditioning import condition_local_ot_cost


def _assignment(
    values,
    *,
    observations=("o1", "o2"),
    categories=("A", "B"),
):
    posterior = pd.DataFrame(values, index=observations, columns=categories)
    return GlobalAssignment(labels=posterior.idxmax(axis=1), posterior=posterior)


def test_fixed_directed_top_k_formula():
    left = _assignment([[1.0, 0.0], [0.0, 1.0]])
    right = _assignment(
        [[1.0, 0.0], [0.0, 1.0], [0.5, 0.5]],
        observations=("r1", "r2", "r3"),
    )
    cost = np.array([[1.0, 2.0], [3.0, 4.0]])
    support = np.array([[1, 2], [0, 2]])

    conditioned = condition_local_ot_cost(
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


def test_zero_strength_returns_original_cost_without_reading_payloads():
    opaque = GlobalAssignment(labels=None, posterior=None)
    cost = np.array([[np.inf, np.nan]], dtype=np.float32)

    assert condition_local_ot_cost(
        cost,
        opaque,
        np.array([[0, 1]]),
        right_posterior=opaque,
        strength=0.0,
        valid_support_mask=object(),
    ) is cost


def test_partial_support_conditions_only_valid_edges_and_restores_hard_support():
    assignment = _assignment([[1.0, 0.0], [0.0, 1.0]])
    cost = np.array([[0.5, 7.0], [8.0, 1.5]])
    support = np.array([[0, 999], [-1, 1]])
    valid_support = np.array([[True, False], [False, True]])
    original = cost.copy()

    conditioned = condition_local_ot_cost(
        cost,
        assignment,
        support,
        strength=0.2,
        valid_support_mask=valid_support,
    )

    np.testing.assert_array_equal(cost, original)
    np.testing.assert_allclose(conditioned[valid_support], [0.5, 1.5])
    assert np.isposinf(conditioned[~valid_support]).all()


@pytest.mark.parametrize("invalid_cost", [np.nan, np.inf])
def test_partial_support_still_rejects_nonfinite_valid_cost(invalid_cost):
    assignment = _assignment([[1.0, 0.0], [0.0, 1.0]])

    with pytest.raises(ValueError, match="valid cost values must be finite"):
        condition_local_ot_cost(
            np.array([[invalid_cost, np.inf], [1.0, np.inf]]),
            assignment,
            np.array([[0, -1], [1, 999]]),
            strength=0.2,
            valid_support_mask=np.array([[True, False], [True, False]]),
        )


def test_partial_support_still_rejects_valid_out_of_bounds_index():
    assignment = _assignment([[1.0, 0.0], [0.0, 1.0]])

    with pytest.raises(ValueError, match="out-of-bounds"):
        condition_local_ot_cost(
            np.array([[1.0, np.inf], [1.0, np.inf]]),
            assignment,
            np.array([[2, -1], [1, 999]]),
            strength=0.2,
            valid_support_mask=np.array([[True, False], [True, False]]),
        )


@pytest.mark.parametrize(
    ("valid_support_mask", "message"),
    [
        (np.ones((2, 1), dtype=bool), "shape"),
        (np.ones((2, 2), dtype=np.int8), "boolean dtype"),
    ],
)
def test_partial_support_mask_is_strict(valid_support_mask, message):
    assignment = _assignment([[1.0, 0.0], [0.0, 1.0]])

    with pytest.raises(ValueError, match=message):
        condition_local_ot_cost(
            np.ones((2, 2)),
            assignment,
            np.array([[0, 1], [0, 1]]),
            strength=0.2,
            valid_support_mask=valid_support_mask,
        )


@pytest.mark.parametrize(
    ("strength", "error"),
    [
        (True, TypeError),
        ("1", TypeError),
        (np.nan, ValueError),
        (np.inf, ValueError),
        (-0.1, ValueError),
    ],
)
def test_strength_is_finite_non_negative_real(strength, error):
    assignment = _assignment([[1.0, 0.0], [0.0, 1.0]])

    with pytest.raises(error, match="strength"):
        condition_local_ot_cost(
            np.ones((2, 1)),
            assignment,
            np.array([[0], [1]]),
            strength=strength,
        )


@pytest.mark.parametrize(
    ("cost", "support", "message"),
    [
        (np.ones((2, 2)), np.array([[0], [1]]), "shapes differ"),
        (np.ones((2, 1)), np.array([0, 1]), "must have shape"),
        (np.ones((2, 1)), np.array([[0.0], [1.0]]), "integer dtype"),
        (np.ones((2, 1)), np.array([[0], [2]]), "out-of-bounds"),
        (np.array([[1.0], [np.inf]]), np.array([[0], [1]]), "finite"),
        (np.array([[1.0], [-1.0]]), np.array([[0], [1]]), "non-negative"),
    ],
)
def test_cost_and_support_are_strict(cost, support, message):
    assignment = _assignment([[1.0, 0.0], [0.0, 1.0]])

    with pytest.raises(ValueError, match=message):
        condition_local_ot_cost(
            cost,
            assignment,
            support,
            strength=0.2,
        )


def test_misaligned_or_invalid_posterior_is_rejected():
    left = _assignment([[1.0, 0.0], [0.0, 1.0]])
    reversed_categories = _assignment(
        [[0.0, 1.0], [1.0, 0.0]],
        observations=("r1", "r2"),
        categories=("B", "A"),
    )
    invalid = _assignment([[0.6, 0.6], [0.0, 1.0]])

    with pytest.raises(GlobalAssignmentContractError, match="category"):
        condition_local_ot_cost(
            np.ones((2, 1)),
            left,
            np.array([[0], [1]]),
            right_posterior=reversed_categories,
            strength=0.2,
        )
    with pytest.raises(GlobalAssignmentContractError, match="row-normalized"):
        condition_local_ot_cost(
            np.ones((2, 1)),
            invalid,
            np.array([[0], [1]]),
            strength=0.2,
        )
