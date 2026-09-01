from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise.analysis.basic.partition_change import (
    align_observation_pairs,
    compare_partitions,
    select_resolution_for_target_cluster_count,
    select_shared_feature_names,
    select_level1_resolution,
    summarize_change_by_level1,
)


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
