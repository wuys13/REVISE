from __future__ import annotations

import logging
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from revise.backend.kernels.spot_sr import SpotSrKernel


def _kernel(*, completeness=True):
    config = SimpleNamespace(
        pm_on_cell_file="/path/that/does/not/exist.csv",
        svc_completeness=completeness,
        sr_assignment_seed=42,
    )
    return SpotSrKernel(config, logging.getLogger("test-spot-sr-quota"))


def _svc_obs(counts):
    rows = []
    for spot, count in counts.items():
        rows.extend(
            {"spot_name": spot, "cell_id": f"{spot}-cell-{i}"}
            for i in range(count)
        )
    return pd.DataFrame(rows)


@pytest.mark.parametrize(
    ("values", "n_cells", "expected"),
    [
        ([0.5, 0.5], 1, [1, 0]),
        ([1 / 3, 1 / 2, 1 / 6], 1, [1, 0, 0]),
        ([1.0, 0.0], 3, [3, 0]),
        ([0.2, 0.2, 0.6], 2, [0, 0, 2]),
        ([1 / 6] * 6, 4, [0, 0, 1, 1, 1, 1]),
    ],
)
def test_quota_uses_numpy_round_and_exact_stable_repairs(values, n_cells, expected):
    columns = [f"type-{i}" for i in range(len(values))]
    contributions = pd.DataFrame([values], index=["spot-a"], columns=columns)

    result = _kernel().get_spot_cell_distribution(
        contributions,
        _svc_obs({"spot-a": n_cells}),
    )

    assert result.loc["spot-a"].tolist() == expected
    assert result.index.tolist() == ["spot-a"]
    assert result.columns.tolist() == columns
    assert all(dtype == np.dtype("int64") for dtype in result.dtypes)


def test_quota_preserves_order_and_integer_invariants_for_many_valid_rows():
    rng = np.random.default_rng(20260721)
    spots = ["spot-c", "spot-a", "spot-b"]
    columns = ["z", "a", "m", "q", "b"]
    values = rng.dirichlet(np.ones(len(columns)), size=len(spots))
    contributions = pd.DataFrame(values, index=spots, columns=columns)
    cell_counts = {"spot-a": 1, "spot-b": 17, "spot-c": 6}

    result = _kernel().get_spot_cell_distribution(
        contributions,
        _svc_obs(cell_counts),
    )

    assert result.index.tolist() == spots
    assert result.columns.tolist() == columns
    assert result.to_numpy().dtype == np.int64
    assert np.all(result.to_numpy() >= 0)
    assert result.sum(axis=1).to_dict() == cell_counts


@pytest.mark.parametrize(
    ("values", "expected"),
    [
        ([1.000009, 0.0], [200_000, 0]),
        ([0.999991, 0.0], [200_000, 0]),
    ],
)
def test_near_one_row_is_repaired_to_an_exact_quota(values, expected):
    n_cells = 200_000
    contributions = pd.DataFrame(
        [values],
        index=["spot-a"],
        columns=["A", "B"],
    )
    svc_obs = pd.DataFrame(
        {
            "spot_name": np.repeat("spot-a", n_cells),
            "cell_id": np.arange(n_cells),
        }
    )

    result = _kernel().get_spot_cell_distribution(contributions, svc_obs)

    assert result.loc["spot-a"].tolist() == expected


def test_quota_requires_svc_completeness_to_be_exactly_true():
    contributions = pd.DataFrame([[1.0]], index=["spot-a"], columns=["A"])

    with pytest.raises(ValueError, match="svc_completeness"):
        _kernel(completeness=1).get_spot_cell_distribution(
            contributions,
            _svc_obs({"spot-a": 1}),
        )


def test_quota_requires_a_dataframe_with_unique_axes_and_exact_spot_set():
    svc_obs = _svc_obs({"spot-a": 1})

    with pytest.raises(TypeError, match="DataFrame"):
        _kernel().get_spot_cell_distribution(np.array([[1.0]]), svc_obs)

    duplicate_spots = pd.DataFrame([[1.0], [1.0]], index=["spot-a", "spot-a"], columns=["A"])
    with pytest.raises(ValueError, match="unique spot"):
        _kernel().get_spot_cell_distribution(duplicate_spots, svc_obs)

    duplicate_classes = pd.DataFrame([[0.5, 0.5]], index=["spot-a"], columns=["A", "A"])
    with pytest.raises(ValueError, match="unique categor"):
        _kernel().get_spot_cell_distribution(duplicate_classes, svc_obs)

    wrong_spots = pd.DataFrame([[1.0]], index=["spot-b"], columns=["A"])
    with pytest.raises(ValueError, match="spot set"):
        _kernel().get_spot_cell_distribution(wrong_spots, svc_obs)


@pytest.mark.parametrize(
    "values",
    [
        [0.4, 0.4],
        [-0.1, 1.1],
        [np.nan, np.nan],
        [np.inf, 0.0],
        ["not-a-number", 1.0],
    ],
)
def test_quota_rejects_invalid_contribution_values(values):
    contributions = pd.DataFrame([values], index=["spot-a"], columns=["A", "B"])

    with pytest.raises((TypeError, ValueError)):
        _kernel().get_spot_cell_distribution(
            contributions,
            _svc_obs({"spot-a": 2}),
        )
