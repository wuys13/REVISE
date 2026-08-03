from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from revise.backend.ops.assignment import GlobalAssignment
from revise.backend.ops.local_ot import solve_local_ot
from revise.backend.ops.posterior_conditioning import condition_local_ot_cost


def _conditioned_cost() -> np.ndarray:
    posterior = pd.DataFrame(
        [[0.9, 0.1], [0.1, 0.9]],
        index=["left", "right"],
        columns=["A", "B"],
    )
    assignment = GlobalAssignment(
        labels=posterior.idxmax(axis=1),
        posterior=posterior,
    )
    return condition_local_ot_cost(
        np.array([[0.0, 1.0], [1.0, 0.0]]),
        assignment,
        np.array([[0, 1], [0, 1]]),
        strength=0.2,
    )


def _solve_with_installed_solver(method: str) -> np.ndarray:
    return solve_local_ot(
        [0.5, 0.5],
        [0.5, 0.5],
        _conditioned_cost(),
        method=method,
        pot_reg=0.1,
        pot_reg_m=1.0,
        pot_reg_type="kl",
    )


def test_real_pot_executes_conditioned_cost_candidate():
    pytest.importorskip("ot", reason="installed POT smoke requires POT")

    coupling = _solve_with_installed_solver("pot")

    assert coupling.shape == (2, 2)
    assert np.isfinite(coupling).all()
    assert float(coupling.sum()) > 0.0


def test_real_tacco_executes_conditioned_cost_candidate():
    pytest.importorskip("tacco", reason="optional TACCO smoke requires tacco")

    coupling = _solve_with_installed_solver("tacco")

    assert coupling.shape == (2, 2)
    assert np.isfinite(coupling).all()
    assert float(coupling.sum()) > 0.0
