from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise.analysis.basic.unsupervised import compute_leiden_sweep


def _cluster_fixture() -> AnnData:
    rng = np.random.default_rng(20260812)
    x = rng.poisson(2, size=(120, 82)).astype(float) + 1
    x[:60, :12] += 8
    x[60:, 12:24] += 8
    var_names = ["MT-DROP", *[f"G{i}" for i in range(81)]]
    obs = pd.DataFrame(
        {"Level1": pd.Categorical(["A"] * 60 + ["B"] * 60)},
        index=[f"cell-{i}" for i in range(120)],
    )
    return AnnData(x, obs=obs, var=pd.DataFrame(index=var_names))


def test_compute_leiden_sweep_preserves_input_and_returns_metric_schema():
    adata = _cluster_fixture()
    original_x = adata.X.copy()

    processed, metrics = compute_leiden_sweep(
        adata,
        resolutions=[0.3, 0.5],
        reference_col="Level1",
        random_state=0,
    )

    np.testing.assert_array_equal(adata.X, original_x)
    assert "MT-DROP" in adata.var_names
    assert not any(column.startswith("leiden_res_") for column in adata.obs)

    assert "MT-DROP" not in processed.var_names
    assert "X_pca" in processed.obsm
    assert "connectivities" in processed.obsp
    assert "distances" in processed.obsp
    assert {"leiden_res_0.3", "leiden_res_0.5"} <= set(processed.obs)
    assert list(metrics.columns) == ["resolution", "ARI", "NMI", "cluster_num"]
    assert metrics["resolution"].tolist() == [0.3, 0.5]
    assert metrics["cluster_num"].tolist() == [
        processed.obs["leiden_res_0.3"].nunique(),
        processed.obs["leiden_res_0.5"].nunique(),
    ]


def test_compute_leiden_sweep_uses_the_existing_metric_authority(monkeypatch):
    from revise.analysis.basic import unsupervised

    calls = []

    def fake_metrics(adata, pred_label_key, true_label_key):
        calls.append((adata, pred_label_key, true_label_key))
        return 0.25, 0.75

    monkeypatch.setattr(unsupervised, "compute_clustering_metrics", fake_metrics)

    processed, metrics = compute_leiden_sweep(
        _cluster_fixture(),
        resolutions=[0.3],
        reference_col="Level1",
    )

    assert calls == [(processed, "leiden_res_0.3", "Level1")]
    assert metrics.loc[0, ["ARI", "NMI"]].tolist() == [0.25, 0.75]


@pytest.mark.parametrize("resolutions", [[], [0], [-0.1], [float("nan")]])
def test_compute_leiden_sweep_rejects_empty_or_invalid_resolutions(resolutions):
    with pytest.raises(ValueError, match="resolutions"):
        compute_leiden_sweep(
            _cluster_fixture(),
            resolutions=resolutions,
            reference_col="Level1",
        )


def test_compute_leiden_sweep_requires_the_reference_column():
    with pytest.raises(KeyError, match="missing"):
        compute_leiden_sweep(
            _cluster_fixture(),
            resolutions=[0.3],
            reference_col="missing",
        )
