from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from revise.analysis.basic.spatial_region import (
    assign_anatomy_candidates,
    assign_square_windows,
    assign_anatomy_regions,
    compute_window_diversity,
    convert_coordinates_to_microns,
    effective_number,
    select_min_parent_support,
    summarize_region,
    summarize_region_by_anatomy,
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
