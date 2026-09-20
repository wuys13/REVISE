from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

from biological_recovery_benchmark import (
    compute_conditional_moran_i,
    compute_global_moran_i,
    compute_identity_metrics,
    compute_tmp_mer,
    evaluate_adata,
    plot_metric_comparison,
    plot_moran_heatmap,
    plot_spatial_metric_comparison,
    save_evaluation_results,
)


def _formula_adata() -> AnnData:
    adata = AnnData(
        X=np.asarray(
            [
                [4.0, 2.0, 1.0, 1.0],
                [2.0, 0.0, 3.0, 1.0],
                [1.0, 1.0, 5.0, 3.0],
                [0.0, 2.0, 3.0, 1.0],
            ]
        ),
        obs=pd.DataFrame(
            {
                "cell_type": ["A", "A", "B", "B"],
                "cluster": ["0", "0", "1", "1"],
            },
            index=["c1", "c2", "c3", "c4"],
        ),
        var=pd.DataFrame(index=["A1", "A2", "B1", "B2"]),
    )
    adata.obsm["X_test"] = np.asarray([[0.0, 0.0], [0.1, 0.0], [5.0, 5.0], [5.1, 5.0]])
    return adata


def _spatial_adata() -> AnnData:
    rng = np.random.default_rng(42)
    n_a = 9
    n_b = 11
    coordinates_a = np.column_stack([np.arange(n_a), np.zeros(n_a)])
    coordinates_b = np.column_stack([np.arange(n_b), np.ones(n_b) * 5])
    coordinates = np.vstack([coordinates_a, coordinates_b]).astype(float)
    x = rng.poisson(3, size=(n_a + n_b, 4)).astype(float) + 1
    x[:n_a, :2] += 6
    x[n_a:, 2:] += 6
    adata = AnnData(
        X=x,
        obs=pd.DataFrame(
            {
                "cell_type": ["A"] * n_a + ["B"] * n_b,
                "cluster": ["0"] * n_a + ["1"] * n_b,
            },
            index=[f"cell-{index}" for index in range(n_a + n_b)],
        ),
        var=pd.DataFrame(index=["A1", "A2", "B1", "B2"]),
    )
    adata.obsm["spatial"] = coordinates
    adata.obsm["X_test"] = np.column_stack(
        [
            np.r_[np.linspace(0, 1, n_a), np.linspace(5, 6, n_b)],
            np.r_[np.linspace(1, 0, n_a), np.linspace(6, 5, n_b)],
        ]
    )
    return adata


def test_identity_and_tmp_mer_match_expected_values():
    adata = _formula_adata()
    identity = compute_identity_metrics(
        adata,
        cell_type_col="cell_type",
        pred_label_col="cluster",
        embedding_key="X_test",
    ).set_index("metric")
    assert identity.loc["ARI", "value"] == pytest.approx(1.0)
    assert identity.loc["AMI", "value"] == pytest.approx(1.0)
    assert identity.loc["NMI", "value"] == pytest.approx(1.0)
    assert identity.loc["ASW", "value"] > 0.9

    tmp_mer = compute_tmp_mer(
        adata,
        cell_type_col="cell_type",
        marker_map={"A": ["A1", "A2"], "B": ["B1", "B2"]},
        eps=0.5,
    ).set_index("target_cell_type")
    assert tmp_mer.loc["A", "TMP"] == pytest.approx(8.0 / 14.5)
    assert tmp_mer.loc["A", "MER"] == pytest.approx(2.0 / 2.0)
    assert tmp_mer.loc["B", "TMP"] == pytest.approx(12.0 / 16.5)
    assert tmp_mer.loc["B", "MER"] == pytest.approx(3.0 / 1.5)


def test_conditional_and_global_moran_match_manual_formula():
    adata = _formula_adata()
    graph = sparse.csr_matrix(
        np.asarray(
            [
                [0, 1, 1, 0],
                [1, 0, 0, 1],
                [1, 0, 0, 1],
                [0, 1, 1, 0],
            ],
            dtype=float,
        )
    )
    adata.obsp["spatial_connectivities"] = graph
    conditional = compute_conditional_moran_i(
        adata,
        cell_type_col="cell_type",
        genes=["B1"],
    )
    global_moran = compute_global_moran_i(adata, genes=["B1"])

    work = adata.copy()
    import scanpy as sc

    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    values = np.asarray(work[:, ["B1"]].X).reshape(-1)
    centered = values - values.mean()

    def manual(weights):
        return (
            len(values)
            / weights.sum()
            * (centered @ weights @ centered)
            / (centered @ centered)
        )

    labels = adata.obs["cell_type"].to_numpy()
    rows, cols = graph.nonzero()
    same = sparse.csr_matrix(
        (
            np.ones(np.sum(labels[rows] == labels[cols])),
            (rows[labels[rows] == labels[cols]], cols[labels[rows] == labels[cols]]),
        ),
        shape=graph.shape,
    )
    different = graph - same
    assert conditional.loc[0, "MISC"] == pytest.approx(manual(same))
    assert conditional.loc[0, "MIDC"] == pytest.approx(manual(different))
    assert global_moran.loc[0, "MoranI"] == pytest.approx(manual(graph))


def test_evaluate_save_and_plot_smoke(tmp_path):
    adata = _spatial_adata()
    original_x = adata.X.copy()
    results = evaluate_adata(
        adata,
        cell_type_col="cell_type",
        pred_label_col="cluster",
        embedding_key="X_test",
        marker_map={"A": ["A1", "A2"], "B": ["B1", "B2"]},
    )
    assert {
        "identity_metrics",
        "tmp_mer",
        "tmp_mer_summary",
        "conditional_moran",
        "conditional_moran_summary",
        "global_moran",
        "cell_type_moran",
        "run_metadata",
    } == set(results)
    assert set(results["cell_type_moran"]["group"]) == {"B"}
    np.testing.assert_array_equal(adata.X, original_x)
    assert "spatial_connectivities" not in adata.obsp

    save_evaluation_results(
        results,
        tmp_path,
        metadata={"dataset": "fixture", "method": "REVISE"},
    )
    expected = {
        "identity_metrics.csv",
        "tmp_mer_by_cell_type.csv",
        "tmp_mer_summary.csv",
        "misc_midc_by_gene.csv",
        "misc_midc_summary.csv",
        "moran_all_by_gene.csv",
        "moran_by_cell_type_and_gene.csv",
        "run_metadata.json",
    }
    assert {path.name for path in tmp_path.iterdir()} == expected
    metadata = json.loads((tmp_path / "run_metadata.json").read_text())
    assert metadata["min_cell_type_size"] == 10
    assert metadata["dataset"] == "fixture"
    assert metadata["skipped_cell_types"] == ["A"]

    identity_plot = results["identity_metrics"].assign(method="REVISE")
    conditional_plot = results["conditional_moran"].assign(method="REVISE")
    moran_plot = pd.concat(
        [results["global_moran"], results["cell_type_moran"]],
        ignore_index=True,
    ).assign(method="REVISE")
    figures = [
        plot_metric_comparison(identity_plot)[0],
        plot_spatial_metric_comparison(conditional_plot, metric="MISC")[0],
        plot_moran_heatmap(moran_plot)[0],
    ]
    assert all(figure.axes for figure in figures)
    for figure in figures:
        plt.close(figure)
