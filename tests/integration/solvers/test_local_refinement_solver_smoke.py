from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from revise.backend.ops.assignment import GlobalAssignment
from revise.backend.ops.local_ot import solve_local_ot
from revise.backend.ops.posterior_conditioning import condition_local_ot_cost


def _conditioned_problem() -> tuple[np.ndarray, np.ndarray]:
    posterior = pd.DataFrame(
        [[0.9, 0.1], [0.1, 0.9]],
        index=["left", "right"],
        columns=["A", "B"],
    )
    assignment = GlobalAssignment(
        labels=posterior.idxmax(axis=1),
        posterior=posterior,
    )
    valid_support = np.array([[True, False], [True, True]])
    conditioned = condition_local_ot_cost(
        np.array([[0.0, np.inf], [1.0, 0.0]]),
        assignment,
        np.array([[0, 999], [0, 1]]),
        strength=0.2,
        valid_support_mask=valid_support,
    )
    return conditioned.T, valid_support.T


def _solve_with_installed_solver(method: str) -> tuple[np.ndarray, np.ndarray]:
    cost, valid_support = _conditioned_problem()
    return (
        solve_local_ot(
            [0.5, 0.5],
            [0.5, 0.5],
            cost,
            method=method,
            pot_reg=0.1,
            pot_reg_m=1.0,
            pot_reg_type="kl",
            valid_support_mask=valid_support,
        ),
        valid_support,
    )


def test_real_pot_executes_conditioned_cost_candidate():
    pytest.importorskip("ot", reason="installed POT smoke requires POT")

    coupling, valid_support = _solve_with_installed_solver("pot")

    assert coupling.shape == (2, 2)
    assert np.isfinite(coupling).all()
    assert float(coupling.sum()) > 0.0
    assert np.count_nonzero(coupling[~valid_support]) == 0


def test_real_tacco_executes_conditioned_cost_candidate():
    pytest.importorskip("tacco", reason="optional TACCO smoke requires tacco")

    coupling, valid_support = _solve_with_installed_solver("tacco")

    assert coupling.shape == (2, 2)
    assert np.isfinite(coupling).all()
    assert float(coupling.sum()) > 0.0
    assert np.count_nonzero(coupling[~valid_support]) == 0
