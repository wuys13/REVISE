from __future__ import annotations

from types import SimpleNamespace
import warnings

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy.sparse import SparseEfficiencyWarning, csr_matrix

from revise.analysis.basic.partition_change import (
    align_observation_pairs,
    compare_partitions,
    filter_paired_sp_svc_inputs,
    leiden_labels,
    prepare_leiden_graph,
    select_raw_hvg_feature_names,
    select_resolution_for_target_cluster_count,
    select_shared_feature_names,
    select_level1_resolution,
    summarize_change_by_level1,
)


def test_leiden_labels_uses_only_the_precomputed_graph(monkeypatch):
    graph = AnnData(
        X=np.ones((4, 10)),
        obs=pd.DataFrame(index=["u0", "u1", "u2", "u3"]),
    )
    graph.obsp["connectivities"] = csr_matrix(np.eye(4))
    captured = {}

    def fake_leiden(adata, **kwargs):
        captured["n_vars"] = adata.n_vars
        captured["adjacency"] = kwargs["adjacency"]
        adata.obs[kwargs["key_added"]] = ["0", "0", "1", "1"]

    monkeypatch.setitem(__import__("sys").modules, "scanpy", SimpleNamespace(tl=SimpleNamespace(leiden=fake_leiden)))

    labels = leiden_labels(graph, resolution=0.5)

    assert labels.tolist() == ["0", "0", "1", "1"]
    assert captured["n_vars"] == 0
    assert captured["adjacency"] is graph.obsp["connectivities"]


def test_prepare_leiden_graph_suppresses_expected_zero_feature_warning(monkeypatch):
    adata = AnnData(
        X=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
        obs=pd.DataFrame(index=["u0", "u1", "u2"]),
        var=pd.DataFrame(index=["g0", "g1"]),
    )

    def fake_normalize(*_args, **_kwargs):
        warnings.warn("Some cells have zero counts", UserWarning)

    fake_scanpy = SimpleNamespace(
        pp=SimpleNamespace(
            normalize_total=fake_normalize,
            log1p=lambda *_args, **_kwargs: None,
            neighbors=lambda work, **_kwargs: work.obsp.__setitem__("connectivities", csr_matrix(np.eye(work.n_obs))),
        ),
        tl=SimpleNamespace(pca=lambda *_args, **_kwargs: None),
    )
    monkeypatch.setitem(__import__("sys").modules, "scanpy", fake_scanpy)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        prepare_leiden_graph(adata, feature_names=["g0", "g1"])

    assert not [item for item in caught if "zero counts" in str(item.message)]


def test_select_level1_resolution_uses_highest_ari_then_lower_resolution():
    sweep = pd.DataFrame(
        {
            "resolution": [0.3, 0.5, 0.8],
            "ARI": [0.42, 0.61, 0.61],
            "NMI": [0.5, 0.7, 0.7],
        }
    )

    selected = select_level1_resolution(sweep)

    assert selected["resolution"] == 0.5
    assert selected["ARI"] == 0.61


def test_select_level1_resolution_rejects_missing_or_nonfinite_ari():
    with pytest.raises(ValueError, match="ARI"):
        select_level1_resolution(pd.DataFrame({"resolution": [0.3], "ARI": [np.nan]}))


def test_target_cluster_resolution_prefers_exact_k_then_lower_resolution():
    sweep = pd.DataFrame(
        {
            "resolution": [0.30, 0.42, 0.50, 0.62],
            "n_clusters": [5, 6, 6, 7],
        }
    )

    selected = select_resolution_for_target_cluster_count(sweep, target_cluster_count=6)

    assert selected["resolution"] == pytest.approx(0.42)
    assert selected["n_clusters"] == 6
    assert selected["cluster_count_gap"] == 0
    assert selected["status"] == "ok"


def test_target_cluster_resolution_marks_unmatched_when_gap_exceeds_one():
    sweep = pd.DataFrame(
        {"resolution": [0.2, 0.4, 0.6], "n_clusters": [2, 3, 4]}
    )

    selected = select_resolution_for_target_cluster_count(sweep, target_cluster_count=6)

    assert selected["resolution"] == pytest.approx(0.6)
    assert selected["cluster_count_gap"] == 2
    assert selected["status"] == "unmatched_cluster_complexity"


def test_change_by_level1_uses_global_assignments_and_wilson_intervals():
    assignments = pd.DataFrame(
        {"unit_changed": [True, False, True, False]}, index=["u0", "u1", "u2", "u3"]
    )
    level1 = pd.Series(["Tumor", "Tumor", "T", "T"], index=assignments.index)

    summary = summarize_change_by_level1(assignments, level1, min_report_n=3).set_index("level1")

    assert summary.loc["Tumor", "total_units"] == 2
    assert summary.loc["Tumor", "changed_units"] == 1
    assert summary.loc["Tumor", "change_fraction"] == pytest.approx(0.5)
    assert summary.loc["Tumor", "low_sample_size"]
    assert summary.loc["Overall", "wilson_ci_lower"] < 0.5 < summary.loc["Overall", "wilson_ci_upper"]


def test_compare_partitions_is_invariant_to_label_permutation_and_reports_complements():
    raw = pd.Series(["a", "a", "b", "b"], index=["u0", "u1", "u2", "u3"])
    recon = pd.Series(["x", "x", "y", "y"], index=raw.index)

    result = compare_partitions(raw, recon, comparison_edge="raw_to_recon")

    summary = result.summary.iloc[0]
    assert summary["comparison_edge"] == "raw_to_recon"
    assert summary["matched_accuracy"] == 1.0
    assert summary["st_unit_change_fraction"] == 0.0
    assert summary["matched_macro_f1"] == 1.0
    assert summary["balanced_cluster_change"] == 0.0
    assert summary["ARI"] == 1.0
    assert summary["AMI"] == 1.0
    assert summary["VI"] == 0.0
    assert summary["matched_accuracy"] + summary["st_unit_change_fraction"] == pytest.approx(1.0, abs=1e-12)
    assert summary["matched_macro_f1"] + summary["balanced_cluster_change"] == pytest.approx(1.0, abs=1e-12)
    assert result.assignments["unit_changed"].sum() == 0


def test_compare_partitions_marks_unmatched_clusters_changed_and_keeps_rare_cluster():
    raw = pd.Series(["a", "a", "b", "rare"], index=["u0", "u1", "u2", "u3"])
    recon = pd.Series(["x", "x", "y", "y"], index=raw.index)

    result = compare_partitions(raw, recon)

    assert set(result.mapping["raw_cluster"].dropna()) == {"a", "b", "rare"}
    assert result.mapping["matched"].sum() == 2
    assert bool(result.assignments.loc["u3", "unit_changed"])
    assert 0 < result.summary.loc[0, "st_unit_change_fraction"] < 1


def test_compare_partitions_requires_identical_unit_ids():
    raw = pd.Series(["a", "b"], index=["u0", "u1"])
    recon = pd.Series(["x", "y"], index=["u1", "u2"])

    with pytest.raises(ValueError, match="same unit IDs"):
        compare_partitions(raw, recon)


def test_compare_partitions_orders_numeric_cluster_labels_naturally():
    raw = pd.Series(["10", "2", "1", "0"], index=["u0", "u1", "u2", "u3"])
    recon = pd.Series(["10", "2", "1", "0"], index=raw.index)

    result = compare_partitions(raw, recon)

    assert result.contingency.index.tolist() == ["0", "1", "2", "10"]
    assert result.contingency.columns.tolist() == ["0", "1", "2", "10"]


def test_align_observation_pairs_reorders_reconstruction_without_mutating_inputs():
    raw = AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame(index=["u0", "u1"]),
        var=pd.DataFrame(index=["g0", "g1"]),
    )
    recon = AnnData(
        X=np.array([[2.0, 2.0], [3.0, 3.0]]),
        obs=pd.DataFrame(index=["u1", "u0"]),
        var=pd.DataFrame(index=["g1", "g0"]),
    )

    aligned_raw, aligned_recon, audit = align_observation_pairs(raw, recon)

    assert aligned_raw.obs_names.tolist() == ["u0", "u1"]
    assert aligned_recon.obs_names.tolist() == ["u0", "u1"]
    assert aligned_raw.var_names.tolist() == aligned_recon.var_names.tolist() == ["g0", "g1"]
    assert audit == {"n_units": 2, "n_shared_genes": 2, "reordered_reconstruction": True}
    assert recon.obs_names.tolist() == ["u1", "u0"]
    assert recon.var_names.tolist() == ["g1", "g0"]


def test_align_observation_pairs_rejects_missing_ids_and_duplicate_ids():
    raw = AnnData(X=np.ones((2, 1)), obs=pd.DataFrame(index=["u0", "u1"]))
    missing = AnnData(X=np.ones((2, 1)), obs=pd.DataFrame(index=["u0", "u2"]))
    duplicate = AnnData(X=np.ones((2, 1)), obs=pd.DataFrame(index=["u0", "u0"]))

    with pytest.raises(ValueError, match="same observation IDs"):
        align_observation_pairs(raw, missing)
    with pytest.raises(ValueError, match="unique"):
        align_observation_pairs(raw, duplicate)


def test_select_shared_feature_names_is_deterministic_and_preserves_raw_order():
    raw = AnnData(
        X=np.array([[1.0, 0.0, 1.0], [5.0, 0.0, 1.0], [9.0, 0.0, 1.0]]),
        var=pd.DataFrame(index=["g0", "g1", "g2"]),
    )
    recon = AnnData(
        X=np.array([[0.0, 3.0, 1.0], [0.0, 8.0, 1.0], [0.0, 13.0, 1.0]]),
        var=pd.DataFrame(index=["g0", "g1", "g2"]),
    )

    selected = select_shared_feature_names(raw, recon, n_top_genes=2)

    assert selected == ["g0", "g1"]
    assert raw.var_names.tolist() == ["g0", "g1", "g2"]


def test_sp_svc_qc_is_defined_on_raw_and_applied_to_both_carriers():
    raw = AnnData(
        X=csr_matrix(
            [
                [1, 1, 0, 1, 1],
                [1, 0, 0, 1, 0],
                [1, 1, 0, 1, 1],
                [0, 0, 1, 1, 1],
                [1, 1, 0, 1, 1],
            ]
        ),
        obs=pd.DataFrame(index=["u0", "u1", "u2", "u3", "u4"]),
        var=pd.DataFrame(index=["g0", "g1", "g2", "MT-g3", "MT-g4"]),
    )
    recon = AnnData(
        X=np.full((5, 5), 9.0),
        obs=pd.DataFrame(index=["u2", "u0", "u4", "u3", "u1"]),
        var=pd.DataFrame(index=["g0", "g1", "g2", "MT-g3", "MT-g4"]),
    )

    filtered_raw, filtered_recon, audit = filter_paired_sp_svc_inputs(
        raw,
        recon,
        min_genes=3,
        min_cells=3,
    )

    assert filtered_raw.obs_names.tolist() == ["u0", "u2", "u4"]
    assert filtered_recon.obs_names.tolist() == ["u0", "u2", "u4"]
    assert filtered_raw.var_names.tolist() == filtered_recon.var_names.tolist() == ["g0", "g1"]
    assert audit["input_units"] == 5
    assert audit["excluded_raw_qc_units"] == 2
    assert audit["excluded_post_gene_filter_units"] == 1
    assert audit["n_units"] == 3
    assert audit["excluded_raw_qc_genes"] == 3


def test_raw_hvg_selection_uses_canonical_raw_preprocessing_and_preserves_gene_order(monkeypatch):
    raw = AnnData(
        X=np.array([[1.0, 2.0, 3.0], [2.0, 3.0, 4.0], [3.0, 4.0, 5.0]]),
        var=pd.DataFrame(index=["g0", "g1", "g2"]),
    )
    captured = {"calls": []}

    def fake_normalize(adata, **kwargs):
        captured["calls"].append("normalize")
        captured.update(kwargs)
        adata.X = np.asarray(adata.X) * 2

    def fake_log1p(adata):
        captured["calls"].append("log1p")
        adata.X = np.log1p(np.asarray(adata.X))

    def fake_hvg(adata, **kwargs):
        captured["calls"].append("hvg")
        captured["matrix"] = np.asarray(adata.X).copy()
        captured.update(kwargs)
        warnings.warn("sparse structure", SparseEfficiencyWarning)
        adata.var["highly_variable"] = [False, True, True]

    monkeypatch.setitem(
        __import__("sys").modules,
        "scanpy",
        SimpleNamespace(
            pp=SimpleNamespace(
                normalize_total=fake_normalize,
                log1p=fake_log1p,
                highly_variable_genes=fake_hvg,
            )
        ),
    )

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        selected = select_raw_hvg_feature_names(raw, n_top_genes=2)

    assert selected == ["g1", "g2"]
    assert captured["calls"] == ["normalize", "log1p", "hvg"]
    assert not np.array_equal(captured["matrix"], np.asarray(raw.X))
    assert captured["target_sum"] == 1e4
    assert captured["flavor"] == "seurat_v3"
    assert captured["n_top_genes"] == 2
    assert captured["check_values"] is False
    assert not [item for item in caught if issubclass(item.category, SparseEfficiencyWarning)]
    assert np.array_equal(raw.X, np.array([[1.0, 2.0, 3.0], [2.0, 3.0, 4.0], [3.0, 4.0, 5.0]]))
