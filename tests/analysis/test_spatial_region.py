from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from revise.analysis.basic.spatial_region import (
    assign_anatomy_candidates,
    assign_square_windows,
    assign_anatomy_regions,
    compute_window_diversity,
    compute_rarefied_window_diversity,
    convert_coordinates_to_microns,
    effective_number,
    select_region_threshold,
    select_window_scale,
    select_min_parent_support,
    summarize_anatomy_context,
    summarize_cluster_change_by_anatomy,
    summarize_diversity_by_anatomy,
    summarize_region,
    summarize_region_by_anatomy,
    summarize_region_extent_by_anatomy,
)


def test_effective_number_handles_pure_balanced_and_empty_labels():
    assert effective_number(["a", "a"]) == pytest.approx(1.0)
    assert effective_number(["a", "b"]) == pytest.approx(2.0)
    assert effective_number([]) != effective_number([])


def test_assign_square_windows_is_non_overlapping_and_assigns_boundaries_once():
    coordinates = pd.DataFrame(
        {"x": [0.0, 39.999, 40.0], "y": [0.0, 39.999, 0.0]},
        index=["u0", "u1", "u2"],
    )

    assigned = assign_square_windows(coordinates, window_side_length=40.0)

    assert assigned.loc["u0", "window_id"] == "0_0"
    assert assigned.loc["u1", "window_id"] == "0_0"
    assert assigned.loc["u2", "window_id"] == "1_0"
    assert assigned.index.is_unique
    assert assigned["window_id"].notna().all()


def test_coordinate_conversion_uses_explicit_platform_scale_without_mutation():
    coordinates = pd.DataFrame({"x": [0.0, 37.6470588], "y": [0.0, 188.2352941]}, index=["a", "b"])

    converted = convert_coordinates_to_microns(coordinates, microns_per_coordinate=0.2125)

    assert converted.loc["b", "x"] == pytest.approx(8.0)
    assert converted.loc["b", "y"] == pytest.approx(40.0)
    assert coordinates.loc["b", "x"] == pytest.approx(37.6470588)


def test_window_scale_knee_uses_occupancy_only_and_lower_tie_break():
    coordinates = pd.DataFrame(
        {"x": [0, 1, 2, 3, 20, 21, 22, 23], "y": [0] * 8},
        index=[f"u{i}" for i in range(8)],
    )

    selected, sensitivity = select_window_scale(
        coordinates,
        candidate_window_sides=[4.0, 8.0, 16.0],
        min_parent_units=4,
        origin=(0.0, 0.0),
    )

    assert selected["status"] == "ok"
    assert selected["window_side_length"] in {4.0, 8.0, 16.0}
    assert sensitivity["retained_parent_unit_fraction"].between(0, 1).all()
    assert set(sensitivity["window_side_length"]) == {4.0, 8.0, 16.0}


def test_paired_rarefaction_is_deterministic_and_keeps_uniform_level1_at_one():
    windows = pd.DataFrame(
        {"window_id": ["a"] * 6 + ["b"] * 4, "x": range(10), "y": [0] * 10},
        index=[f"u{i}" for i in range(10)],
    )
    raw = pd.Series(["Fibroblast"] * 10, index=windows.index)
    recon = pd.Series(["a", "a", "b", "b", "c", "c", "a", "b", "a", "b"], index=windows.index)

    first = compute_rarefied_window_diversity(
        windows, raw, recon, min_parent_units=4, n_draws=25, random_state=42
    )
    second = compute_rarefied_window_diversity(
        windows, raw, recon, min_parent_units=4, n_draws=25, random_state=42
    )

    pd.testing.assert_frame_equal(first, second)
    assert (first.loc[first["valid_window"], "neff_raw"] == 1.0).all()
    assert (first.loc[first["valid_window"], "neff_recon"] >= 1.0).all()
    assert (first.loc[first["valid_window"], "k_obs_raw"] == 1.0).all()
    assert (first.loc[first["valid_window"], "evenness_raw"] == 1.0).all()


def test_paired_rarefaction_reports_richness_effective_diversity_and_evenness():
    windows = pd.DataFrame(
        {"window_id": ["a"] * 4, "x": range(4), "y": [0] * 4},
        index=[f"u{i}" for i in range(4)],
    )
    raw = pd.Series(["a", "a", "a", "b"], index=windows.index)
    recon = pd.Series(["a", "b", "c", "d"], index=windows.index)

    metrics = compute_rarefied_window_diversity(
        windows, raw, recon, min_parent_units=4, n_draws=3, random_state=42
    ).iloc[0]

    assert metrics["k_obs_raw"] == pytest.approx(2.0)
    assert metrics["neff_raw"] == pytest.approx(np.exp(-0.75 * np.log(0.75) - 0.25 * np.log(0.25)))
    assert metrics["evenness_raw"] == pytest.approx(metrics["neff_raw"] / 2.0)
    assert metrics["k_obs_recon"] == pytest.approx(4.0)
    assert metrics["neff_recon"] == pytest.approx(4.0)
    assert metrics["evenness_recon"] == pytest.approx(1.0)
    assert metrics["delta_k_obs"] == pytest.approx(2.0)


def test_paired_rarefaction_uses_the_same_draws_for_raw_leiden_level2_and_recon():
    windows = pd.DataFrame(
        {"window_id": ["a"] * 5, "x": range(5), "y": [0] * 5},
        index=[f"u{i}" for i in range(5)],
    )
    raw_leiden = pd.Series(["a", "a", "b", "b", "c"], index=windows.index)
    raw_level2 = pd.Series(["x", "x", "y", "y", "z"], index=windows.index)
    recon = raw_leiden.copy()

    metrics = compute_rarefied_window_diversity(
        windows,
        raw_leiden,
        recon,
        raw_level2_labels=raw_level2,
        min_parent_units=4,
        n_draws=25,
        random_state=42,
    ).iloc[0]

    assert metrics["k_obs_level2"] == pytest.approx(metrics["k_obs_raw"])
    assert metrics["neff_level2"] == pytest.approx(metrics["neff_raw"])
    assert metrics["evenness_level2"] == pytest.approx(metrics["evenness_raw"])
    assert metrics["delta_k_obs_vs_raw_leiden"] == pytest.approx(0.0)
    assert metrics["delta_neff_vs_raw_leiden"] == pytest.approx(0.0)
    assert metrics["delta_evenness_vs_raw_leiden"] == pytest.approx(0.0)
    assert metrics["delta_k_obs_vs_raw_level2"] == pytest.approx(0.0)
    assert metrics["delta_neff_vs_raw_level2"] == pytest.approx(0.0)
    assert metrics["delta_evenness_vs_raw_level2"] == pytest.approx(0.0)


def test_region_threshold_refuses_small_or_unstable_window_sets():
    selected, bootstrap = select_region_threshold(np.linspace(1.0, 3.0, 50), n_bootstrap=20)

    assert selected["status"] == "no_stable_threshold"
    assert selected["threshold"] is None
    assert bootstrap.empty


def test_parent_support_knee_uses_lower_tie_and_never_returns_less_than_two():
    assignments = pd.DataFrame(
        {"window_id": ["a"] * 10 + ["b"] * 8 + ["c"] * 5 + ["d"] * 2},
        index=[f"u{i}" for i in range(25)],
    )

    selected, sensitivity = select_min_parent_support(assignments)

    assert selected["status"] == "ok"
    assert selected["min_parent_units"] >= 2
    assert sensitivity["retained_unit_fraction"].iloc[0] == pytest.approx(1.0)
    assert sensitivity["valid_window_fraction"].iloc[0] == pytest.approx(1.0)


def test_parent_support_marks_sparse_cohort_insufficient():
    assignments = pd.DataFrame({"window_id": ["a", "b", "c"]}, index=["u0", "u1", "u2"])

    selected, sensitivity = select_min_parent_support(assignments)

    assert selected == {"status": "insufficient_support", "min_parent_units": None, "p95_occupancy": 1}
    assert sensitivity.empty


def test_compute_window_diversity_preserves_invalid_windows_as_nan():
    windows = pd.DataFrame(
        {"window_id": ["0_0", "0_0", "1_0"], "x": [0.0, 1.0, 50.0], "y": [0.0, 1.0, 0.0]},
        index=["u0", "u1", "u2"],
    )
    raw = pd.Series(["a", "a", "b"], index=windows.index)
    recon = pd.Series(["a", "b", "b"], index=windows.index)

    metrics = compute_window_diversity(windows, raw, recon, min_units_per_window=2)

    valid = metrics.set_index("window_id").loc["0_0"]
    invalid = metrics.set_index("window_id").loc["1_0"]
    assert valid["neff_raw"] == pytest.approx(1.0)
    assert valid["neff_recon"] == pytest.approx(2.0)
    assert valid["delta_neff"] == pytest.approx(1.0)
    assert bool(valid["valid_window"])
    assert not bool(invalid["valid_window"])
    assert np.isnan(invalid["neff_raw"])
    assert np.isnan(invalid["neff_recon"])


def test_anatomy_regions_map_configured_epithelial_source_to_normal():
    windows = pd.DataFrame(
        {"window_id": ["0_0", "0_0", "1_0", "2_0"], "x": [0, 1, 50, 90], "y": [0, 1, 0, 0]},
        index=["u0", "u1", "u2", "u3"],
    )
    labels = pd.Series(["Tumor", "Intestinal Epithelial", "Tumor", "Other"], index=windows.index)
    anatomy = assign_anatomy_regions(windows, labels)
    assert anatomy.set_index("window_id").loc["0_0", "region_label"] == "Interface"
    assert anatomy.set_index("window_id").loc["1_0", "region_label"] == "Tumor"
    assert anatomy.set_index("window_id").loc["0_0", "normal_units"] == 1
    assert anatomy.set_index("window_id").loc["2_0", "region_label"] == "Other"

    normal_window = windows.loc[["u1"]].copy()
    normal_anatomy = assign_anatomy_regions(normal_window, labels.loc[["u1"]])
    assert normal_anatomy.loc[0, "region_label"] == "Normal"

    metrics = pd.DataFrame(
        {
            "window_id": ["0_0", "1_0", "2_0"],
            "valid_window": [True, True, False],
            "neff_recon": [2.1, 1.9, np.nan],
            "n_units": [2, 1, 1],
        }
    )
    summary, flagged = summarize_region(metrics, threshold=2.0, window_side_length=40.0)
    assert flagged.set_index("window_id").loc["0_0", "in_region"]
    assert summary.loc[0, "region_area"] == pytest.approx(1600.0)
    assert summary.loc[0, "region_area_fraction"] == pytest.approx(0.5)
    assert summary.loc[0, "region_unit_fraction"] == pytest.approx(2 / 3)


def test_anatomy_candidates_keep_binary_sources_union_and_context_stratification():
    windows = pd.DataFrame(
        {"window_id": ["0_0", "0_0", "1_0", "2_0"], "x": [0, 1, 50, 90], "y": [0, 1, 0, 0]},
        index=["u0", "u1", "u2", "u3"],
    )
    labels = pd.Series(["Tumor", "Intestinal Epithelial", "Tumor", "Other"], index=windows.index)

    anatomy = assign_anatomy_candidates(windows, labels).set_index("window_id")

    assert anatomy.loc["0_0", ["tumor_candidate", "normal_candidate", "union_candidate"]].tolist() == [1, 1, 2]
    assert anatomy.loc["0_0", "level1_region"] == "Interface"
    assert anatomy.loc["1_0", "level1_region"] == "Tumor"
    assert anatomy.loc["2_0", "level1_region"] == "Other"

    metrics = pd.DataFrame(
        {
            "window_id": ["0_0", "1_0", "2_0"],
            "valid_window": [True, True, False],
            "in_region": [True, False, False],
            "n_units": [2, 1, 1],
            "neff_recon": [2.2, 1.2, np.nan],
            "delta_neff": [0.5, -0.1, np.nan],
        }
    )
    summary = summarize_region_by_anatomy(metrics, anatomy.reset_index(), window_side_length=40.0)
    interface = summary.set_index("level1_region").loc["Interface"]
    assert interface["region_area"] == pytest.approx(1600.0)
    assert interface["region_unit_fraction"] == pytest.approx(1.0)


def test_anatomy_summary_accepts_metrics_already_annotated_with_context():
    anatomy = pd.DataFrame({"window_id": ["a"], "level1_region": ["Tumor"]})
    metrics = pd.DataFrame(
        {
            "window_id": ["a"],
            "valid_window": [True],
            "in_region": [True],
            "n_units": [3],
            "neff_recon": [2.2],
            "delta_neff": [0.4],
            "level1_region": ["Tumor"],
        }
    )

    summary = summarize_region_by_anatomy(metrics, anatomy, window_side_length=40.0)

    assert summary.loc[0, "level1_region"] == "Tumor"


def test_anatomy_context_summary_uses_all_tissue_windows_as_area_denominator():
    assignments = pd.DataFrame(
        {"window_id": ["a", "a", "b", "c", "d"]},
        index=["u0", "u1", "u2", "u3", "u4"],
    )
    anatomy = pd.DataFrame(
        {
            "window_id": ["a", "b", "c", "d"],
            "level1_region": ["Tumor", "Normal", "Interface", "Other"],
        }
    )

    summary = summarize_anatomy_context(
        assignments,
        anatomy,
        window_side_length=40.0,
    ).set_index("level1_region")

    assert summary.loc["Tumor", "full_level1_units"] == 2
    assert summary.loc["Tumor", "tissue_windows"] == 1
    assert summary.loc["Tumor", "area_um2"] == pytest.approx(1600.0)
    assert summary.loc["Tumor", "area_mm2"] == pytest.approx(0.0016)
    assert summary.loc["Tumor", "area_fraction"] == pytest.approx(0.25)
    assert summary.loc["Other", "area_fraction"] == pytest.approx(0.25)


def test_cluster_change_summary_uses_existing_global_mapping_and_valid_window_distribution():
    units = pd.DataFrame(
        {
            "window_id": ["a", "a", "b", "c", "d"],
            "level1_region": ["Tumor", "Tumor", "Normal", "Interface", "Other"],
            "unit_changed": [True, False, True, False, True],
        },
        index=["u0", "u1", "u2", "u3", "u4"],
    )
    windows = pd.DataFrame(
        {
            "window_id": ["a", "b", "c", "d"],
            "level1_region": ["Tumor", "Normal", "Interface", "Other"],
            "valid_window": [True, True, False, True],
            "unit_change_fraction": [0.5, 1.0, 0.0, 1.0],
        }
    )

    summary = summarize_cluster_change_by_anatomy(units, windows).set_index("level1_region")

    assert summary.loc["Overall", "paired_units"] == 5
    assert summary.loc["Overall", "changed_units"] == 3
    assert summary.loc["Overall", "change_fraction"] == pytest.approx(0.6)
    assert summary.loc["Tumor", "change_fraction"] == pytest.approx(0.5)
    assert summary.loc["Interface", "n_valid_windows"] == 0
    assert np.isnan(summary.loc["Interface", "median_window_change"])


def test_diversity_summary_reports_median_iqr_and_excludes_invalid_windows():
    metrics = pd.DataFrame(
        {
            "level1_region": ["Tumor", "Tumor", "Normal", "Other"],
            "valid_window": [True, True, True, False],
            "neff_raw": [1.0, 3.0, 2.0, 99.0],
            "neff_recon": [2.0, 4.0, 3.0, 99.0],
            "delta_neff": [1.0, 1.0, 1.0, 0.0],
        }
    )

    summary = summarize_diversity_by_anatomy(metrics).set_index("level1_region")

    assert summary.loc["Overall", "n_valid_windows"] == 3
    assert summary.loc["Overall", "median_neff_raw"] == pytest.approx(2.0)
    assert summary.loc["Overall", "q1_neff_raw"] == pytest.approx(1.5)
    assert summary.loc["Overall", "q3_neff_raw"] == pytest.approx(2.5)
    assert summary.loc["Tumor", "median_neff_recon"] == pytest.approx(3.0)
    assert summary.loc["Normal", "median_delta_neff"] == pytest.approx(1.0)


def test_region_extent_summary_reports_exact_denominators_and_area_units():
    metrics = pd.DataFrame(
        {
            "level1_region": ["Tumor", "Tumor", "Normal", "Other"],
            "valid_window": [True, True, True, False],
            "in_region": [True, False, True, False],
            "n_units": [4, 6, 5, 100],
        }
    )

    summary = summarize_region_extent_by_anatomy(
        metrics,
        window_side_length=40.0,
    ).set_index("level1_region")

    assert summary.loc["Overall", "valid_windows"] == 3
    assert summary.loc["Overall", "region_windows"] == 2
    assert summary.loc["Overall", "region_area_um2"] == pytest.approx(3200.0)
    assert summary.loc["Overall", "region_area_mm2"] == pytest.approx(0.0032)
    assert summary.loc["Overall", "area_fraction"] == pytest.approx(2 / 3)
    assert summary.loc["Overall", "valid_units"] == 15
    assert summary.loc["Overall", "region_units"] == 9
    assert summary.loc["Overall", "unit_fraction"] == pytest.approx(0.6)
