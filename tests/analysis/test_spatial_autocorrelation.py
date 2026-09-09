from __future__ import annotations

import importlib
import warnings

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise.analysis.basic.spatial_autocorrelation import (
    compute_groupwise_spatial_autocorrelation,
)


def _spatial_fixture() -> AnnData:
    x = np.arange(40, dtype=float).reshape(8, 5) + 1
    obs = pd.DataFrame(
        {"Level1": pd.Categorical(["A"] * 3 + ["B"] * 2 + ["C"] * 3)},
        index=[f"cell-{i}" for i in range(8)],
    )
    adata = AnnData(x, obs=obs, var=pd.DataFrame(index=[f"G{i}" for i in range(5)]))
    adata.obsm["spatial"] = np.array(
        [[0, 0], [1, 0], [0, 1], [10, 0], [11, 0], [20, 0], [21, 0], [20, 1]],
        dtype=float,
    )
    return adata


@pytest.mark.parametrize(
    ("mode", "uns_key", "statistic"),
    [("moran", "moranI", "I"), ("geary", "gearyC", "C")],
)
def test_groupwise_autocorrelation_copies_subsets_and_returns_stable_schema(
    monkeypatch, mode, uns_key, statistic
):
    from revise.analysis.basic import spatial_autocorrelation

    adata = _spatial_fixture()
    original_x = adata.X.copy()
    graph_calls = []
    autocorr_calls = []

    def fake_neighbors(work):
        graph_calls.append(tuple(work.obs_names))
        work.obsp["spatial_connectivities"] = np.eye(work.n_obs)

    def fake_autocorr(work, *, mode, genes, n_perms, n_jobs):
        autocorr_calls.append((tuple(work.obs_names), mode, tuple(genes), n_perms, n_jobs))
        work.uns[uns_key] = pd.DataFrame(
            {statistic: np.arange(work.n_vars, dtype=float) + work.n_obs},
            index=work.var_names,
        )

    monkeypatch.setattr(spatial_autocorrelation.sq.gr, "spatial_neighbors", fake_neighbors)
    monkeypatch.setattr(spatial_autocorrelation.sq.gr, "spatial_autocorr", fake_autocorr)

    result = compute_groupwise_spatial_autocorrelation(
        adata,
        groupby="Level1",
        mode=mode,
        min_obs=3,
    )

    np.testing.assert_array_equal(adata.X, original_x)
    assert "spatial_connectivities" not in adata.obsp
    assert list(result.index) == list(adata.var_names)
    assert list(result.columns) == ["All", "A", "C"]
    assert graph_calls == [
        tuple(adata.obs_names),
        tuple(adata.obs_names[:3]),
        tuple(adata.obs_names[5:]),
    ]
    assert [call[1:] for call in autocorr_calls] == [
        (mode, tuple(adata.var_names), 1, 1),
        (mode, tuple(adata.var_names), 1, 1),
        (mode, tuple(adata.var_names), 1, 1),
    ]
    assert result.loc["G0"].tolist() == [8.0, 3.0, 3.0]


def test_groupwise_autocorrelation_validates_before_squidpy(monkeypatch):
    from revise.analysis.basic import spatial_autocorrelation

    calls = []
    monkeypatch.setattr(
        spatial_autocorrelation.sq.gr,
        "spatial_neighbors",
        lambda *_args, **_kwargs: calls.append("called"),
    )
    adata = _spatial_fixture()

    with pytest.raises(ValueError, match="mode"):
        compute_groupwise_spatial_autocorrelation(
            adata, groupby="Level1", mode="invalid", min_obs=3
        )
    with pytest.raises(KeyError, match="missing"):
        compute_groupwise_spatial_autocorrelation(
            adata, groupby="missing", mode="moran", min_obs=3
        )
    del adata.obsm["spatial"]
    with pytest.raises(KeyError, match="spatial"):
        compute_groupwise_spatial_autocorrelation(
            adata, groupby="Level1", mode="moran", min_obs=3
        )

    assert calls == []


def test_import_does_not_change_global_warning_filters():
    import revise.analysis.basic.spatial_autocorrelation as module

    before = list(warnings.filters)
    importlib.reload(module)
    assert warnings.filters == before
