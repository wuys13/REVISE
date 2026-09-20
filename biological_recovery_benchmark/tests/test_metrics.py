from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import squidpy as sq
from anndata import AnnData
from scipy import sparse

from biological_recovery_benchmark import (
    compute_cell_type_moran_i,
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


def _irregular_connectivity(n_obs: int) -> sparse.csr_matrix:
    edges = [
        (0, 1),
        (0, 2),
        (0, 9),
        (1, 2),
        (1, 10),
        (2, 3),
        (3, 4),
        (4, 5),
        (4, 10),
        (5, 6),
        (6, 7),
        (7, 8),
        (9, 10),
        (9, 11),
        (9, 12),
        (10, 11),
        (10, 13),
        (11, 12),
        (12, 13),
        (13, 14),
        (13, 15),
        (14, 15),
        (15, 16),
        (16, 17),
        (16, 19),
        (17, 18),
        (18, 19),
    ]
    rows = [start for edge in edges for start in edge]
    cols = [end for edge in edges for end in edge[::-1]]
    return sparse.csr_matrix(
        (np.ones(len(rows)), (rows, cols)),
        shape=(n_obs, n_obs),
    )


def _direct_squidpy_moran(
    adata: AnnData,
    connectivity: sparse.csr_matrix,
    genes: list[str],
) -> np.ndarray:
    work = adata.copy()
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    work = work[:, genes].copy()
    work.obsp["spatial_connectivities"] = connectivity
    result = sq.gr.spatial_autocorr(
        work,
        connectivity_key="spatial_connectivities",
        genes=genes,
        mode="moran",
        transformation=True,
        n_perms=None,
        corr_method=None,
        use_raw=False,
        copy=True,
    )
    return result.loc[genes, "I"].to_numpy()


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


def test_moran_wrappers_match_squidpy_on_irregular_graph_and_regraph_by_type():
    adata = _spatial_adata()
    graph = _irregular_connectivity(adata.n_obs)
    adata.obsp["spatial_connectivities"] = graph
    genes = ["A1", "B1"]
    conditional = compute_conditional_moran_i(
        adata,
        cell_type_col="cell_type",
        genes=genes,
    )
    global_moran = compute_global_moran_i(adata, genes=genes)

    labels = adata.obs["cell_type"].astype(str).to_numpy()
    graph_coo = graph.tocoo()
    same_mask = labels[graph_coo.row] == labels[graph_coo.col]
    same = sparse.coo_matrix(
        (
            graph_coo.data[same_mask],
            (graph_coo.row[same_mask], graph_coo.col[same_mask]),
        ),
        shape=graph.shape,
    ).tocsr()
    different = sparse.coo_matrix(
        (
            graph_coo.data[~same_mask],
            (graph_coo.row[~same_mask], graph_coo.col[~same_mask]),
        ),
        shape=graph.shape,
    ).tocsr()
    assert np.unique(np.diff(graph.indptr)).size > 1
    np.testing.assert_allclose(
        conditional["MISC"],
        _direct_squidpy_moran(adata, same, genes),
    )
    np.testing.assert_allclose(
        conditional["MIDC"],
        _direct_squidpy_moran(adata, different, genes),
    )
    np.testing.assert_allclose(
        global_moran["MoranI"],
        _direct_squidpy_moran(adata, graph, genes),
    )

    by_type = compute_cell_type_moran_i(
        adata,
        cell_type_col="cell_type",
        genes=genes,
        min_cell_type_size=9,
    )
    for cell_type in ["A", "B"]:
        selected = adata.obs["cell_type"].astype(str) == cell_type
        subset = adata[selected].copy()
        inherited_n_edges = int(subset.obsp["spatial_connectivities"].nnz)
        del subset.obsp["spatial_connectivities"]
        sq.gr.spatial_neighbors(
            subset,
            spatial_key="spatial",
            key_added="spatial",
        )
        rebuilt = subset.obsp["spatial_connectivities"].tocsr()
        actual = by_type[by_type["group"] == cell_type].set_index("Gene")
        expected = _direct_squidpy_moran(subset, rebuilt, genes)
        np.testing.assert_allclose(actual.loc[genes, "MoranI"], expected)
        assert actual["n_edges"].iloc[0] == rebuilt.nnz
        assert rebuilt.nnz != inherited_n_edges


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
    saved_tmp_mer = pd.read_csv(tmp_path / "tmp_mer_by_cell_type.csv")
    saved_conditional = pd.read_csv(tmp_path / "misc_midc_by_gene.csv")
    assert list(saved_tmp_mer["target_cell_type"]) == ["A", "B"]
    assert list(saved_conditional["Gene"]) == list(adata.var_names)

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
