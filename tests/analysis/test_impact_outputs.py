"""Contract tests for the reconstruction-impact long-table writer."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from revise.analysis.basic.partition_change import PartitionComparison
from revise.analysis.reconstruction_impact import (
    AnatomyRegionAnalysis,
    PartitionAnalysis,
    RawLevel2Mapping,
    SpatialImpactAnalysis,
    compute_spatial_impact,
)
from revise.analysis.impact_outputs import ImpactScope, write_impact_outputs


def _partition() -> PartitionAnalysis:
    def comparison(edge: str, changed: bool) -> PartitionComparison:
        return PartitionComparison(
            summary=pd.DataFrame(
                [
                    {
                        "comparison_edge": edge,
                        "n_units": 2,
                        "n_raw_clusters": 2,
                        "n_recon_clusters": 2,
                        "matched_accuracy": 0.5 if changed else 1.0,
                        "st_unit_change_fraction": 0.5 if changed else 0.0,
                        "matched_macro_f1": 0.5 if changed else 1.0,
                        "balanced_cluster_change": 0.5 if changed else 0.0,
                        "ARI": 0.0 if changed else 1.0,
                        "AMI": 0.0 if changed else 1.0,
                        "NMI": 0.0 if changed else 1.0,
                        "VI": 1.0 if changed else 0.0,
                        "status": "ok",
                    }
                ]
            ),
            mapping=pd.DataFrame(
                {
                    "raw_cluster": ["a", "b"],
                    "recon_cluster": ["x", "y"],
                    "matched": [True, True],
                    "overlap_n": [1, 1],
                    "raw_size": [1, 1],
                    "recon_size": [1, 1],
                    "precision": [1.0, 1.0],
                    "recall": [1.0, 1.0],
                    "f1": [1.0, 1.0],
                }
            ),
            contingency=pd.DataFrame(
                [[1, 0], [0, 1]], index=pd.Index(["a", "b"], name="raw_cluster"),
                columns=pd.Index(["x", "y"], name="recon_cluster"),
            ),
            assignments=pd.DataFrame(
                {
                    "raw_cluster": ["a", "b"],
                    "matched_raw_cluster": ["a", "b"],
                    "recon_cluster": ["x", "x" if changed else "y"],
                    "unit_changed": [False, changed],
                    "comparison_edge": edge,
                },
                index=pd.Index(["u1", "u2"], name="unit_id"),
            ),
        )

    return PartitionAnalysis(
        resolution=0.42,
        resolution_source="matched_cluster_count",
        sweep=pd.DataFrame(
            {"resolution": [0.3, 0.42], "n_clusters": [1, 2], "target_cluster_count": [2, 2]}
        ),
        audit={"n_units": 2, "raw_qc_status": "ok"},
        comparisons={
            "raw_to_recon_expression": comparison("raw_to_recon_expression", True),
        },
        representation_audit={"expression_difference_nnz": 1},
        feature_names=["G1", "G2"],
        complexity_comparisons={
            "raw_to_final_svc": comparison("raw_to_final_svc", False),
        },
        complexity_sweep=pd.DataFrame(
            {"resolution": [0.5], "ARI": [0.25], "n_clusters": [2]}
        ),
        complexity_resolution=0.5,
        change_by_level1=pd.DataFrame(
            {
                "level1": ["Overall"],
                "total_units": [2],
                "changed_units": [1],
                "change_fraction": [0.5],
            }
        ),
        matched_cluster_status="ok",
    )


def _spatial() -> SpatialImpactAnalysis:
    return SpatialImpactAnalysis(
        scale_audit={"main_window_side_um": 12.345678901234567, "min_parent_units": 2},
        support_selection={"status": "ok", "window_side_length": 12.345678901234567},
        support_sensitivity=pd.DataFrame(
            {"window_side_length": [8.0], "n_valid_windows": [1]}
        ),
        unit_assignments=pd.DataFrame(
            {
                "x": [1.0, 2.0],
                "y": [3.0, 4.0],
                "window_id": ["w1", "w1"],
                "raw_cluster": ["a", "b"],
                "reconstructed_cluster": ["x", "y"],
                "unit_changed": [False, True],
            },
            index=pd.Index(["u1", "u2"], name="unit_id"),
        ),
        anatomy_unit_assignments=pd.DataFrame(
            {"window_id": ["a1", "a1"], "level1_region": ["Tumor", "Normal"]},
            index=pd.Index(["u1", "u2"], name="unit_id"),
        ),
        window_metrics=pd.DataFrame(
            {
                "window_id": ["w1"],
                "window_x": [1.2345678901234567],
                "window_y": [2.0],
                "n_units": [2],
                "valid_window": [True],
                "level1_region": ["Tumor"],
                "neff_raw": [1.25],
                "neff_recon": [2.5],
                "delta_neff_vs_raw_leiden": [1.25],
                "neff_level2": [1.75],
                "delta_neff_vs_raw_level2": [0.75],
                "neff_sd_raw": [0.1],
                "neff_sd_recon": [0.2],
                "delta_neff_vs_raw_leiden_sd": [0.3],
                "neff_sd_level2": [0.4],
                "delta_neff_vs_raw_level2_sd": [0.5],
            }
        ),
        anatomy_windows=pd.DataFrame(
            {"window_id": ["a1"], "level1_region": ["Tumor"], "area_um2": [4.0]}
        ),
        anatomy_context_summary=pd.DataFrame(
            {"level1_region": ["Tumor"], "tissue_windows": [1]}
        ),
        cluster_change_by_anatomy=pd.DataFrame(
            {"level1_region": ["Tumor"], "changed_units": [1], "total_units": [2]}
        ),
        diversity_by_anatomy=pd.DataFrame(
            {
                "level1_region": ["Tumor"],
                "n_valid_windows": [1],
                "median_neff_raw": [1.25],
                "median_neff_recon": [2.5],
                "median_delta_neff": [1.25],
                "median_delta_neff_vs_raw_level2": [0.75],
            }
        ),
        region_extent_by_anatomy=pd.DataFrame(
            {
                "level1_region": ["Tumor"],
                "n_region_windows": [1],
                "area_um2": [1.1234567890123457],
                "area_fraction": [0.25],
            }
        ),
        scale_sensitivity=pd.DataFrame(
            {
                "window_side_length": [8.0],
                "median_delta_neff_vs_raw_leiden": [1.25],
                "median_delta_neff_vs_raw_level2": [0.75],
            }
        ),
        threshold_sensitivity=pd.DataFrame(
            [
                {"region_type": "state", "status": "no_stable_threshold", "threshold": np.nan},
                {"region_type": "gain", "status": "ok", "threshold": 0.5},
            ]
        ),
        gain_region_extent_by_anatomy=pd.DataFrame(
            {
                "level1_region": ["Tumor"],
                "n_region_windows": [pd.NA],
                "area_um2": [pd.NA],
                "area_fraction": [pd.NA],
            }
        ),
        state_threshold={"status": "no_stable_threshold", "threshold": None, "n_windows": 1},
        gain_threshold={"status": "ok", "threshold": 0.5, "n_windows": 1},
        state_threshold_bootstrap=pd.DataFrame(columns=["bootstrap_threshold"]),
        gain_threshold_bootstrap=pd.DataFrame({"bootstrap_threshold": [0.4, 0.6]}),
    )


def _anatomy() -> AnatomyRegionAnalysis:
    return AnatomyRegionAnalysis(
        scale_audit={"main_window_side_um": 8.0, "origin_x_um": 0.0},
        support_selection={"status": "ok", "window_side_length": 8.0},
        support_sensitivity=pd.DataFrame({"window_side_length": [8.0], "n_valid_windows": [1]}),
        anatomy_unit_assignments=pd.DataFrame(
            {"window_id": ["a1"], "level1_region": ["Tumor"]},
            index=pd.Index(["u1"], name="unit_id"),
        ),
        anatomy_windows=pd.DataFrame({"window_id": ["a1"], "level1_region": ["Tumor"]}),
        anatomy_context_summary=pd.DataFrame({"level1_region": ["Tumor"], "tissue_windows": [1]}),
    )


def _raw_level2() -> RawLevel2Mapping:
    return RawLevel2Mapping(
        labels=pd.Series(["L2a", "L2b"], index=pd.Index(["u1", "u2"], name="unit_id")),
        posterior=pd.DataFrame(
            [[0.12345678901234567, 0.8765432109876543], [0.9, 0.1]],
            index=pd.Index(["u1", "u2"], name="unit_id"),
            columns=pd.Index(["L2a", "L2b"], name="raw_level2"),
        ),
        assignments=pd.DataFrame(
            {
                "unit_id": ["u1", "u2"],
                "raw_level2": ["L2a", "L2b"],
                "confidence": [0.8765432109876543, 0.9],
                "mapping_method": ["knn", "knn"],
            }
        ),
        audit={"method": "knn", "n_missing": 0},
    )


def _scope() -> ImpactScope:
    return ImpactScope(
        comparison_id="cmp-main",
        sample_id="sample-1",
        task_cell_type="T-cell",
        scope="parent-T-cell",
        comparison={
            "comparison_edge": "raw_to_recon_expression",
            "raw_view": "raw_leiden",
            "reconstruction_view": "reconstructed",
            "label_source": "observed",
            "observation_basis": "paired_units",
            "gene_rule": "shared_hvg",
            "normalization": "log1p",
        },
        partition=_partition(),
        spatial=_spatial(),
        anatomy=_anatomy(),
        raw_level2=_raw_level2(),
    )


def test_writer_emits_registered_long_tables_and_deterministic_gzip(tmp_path: Path):
    first = write_impact_outputs(
        tmp_path / "first", scopes=[_scope()], audit={"schema_version": 1},
        calculation={"input_view": "paired", "parameters": {}, "comparison_basis": "impact"},
    )
    second = write_impact_outputs(
        tmp_path / "second", scopes=[_scope()], audit={"schema_version": 1},
        calculation={"input_view": "paired", "parameters": {}, "comparison_basis": "impact"},
    )

    first_paths = {item["path"] for item in first["artifacts"].values()}
    assert first_paths
    assert all((tmp_path / "first" / path).is_file() for path in first_paths)
    assert set(first["artifacts"]) == set(second["artifacts"])
    for role, metadata in first["artifacts"].items():
        assert metadata["description"]
        assert metadata["path"] == second["artifacts"][role]["path"]
        assert (tmp_path / "first" / metadata["path"]).read_bytes() == (
            tmp_path / "second" / second["artifacts"][role]["path"]
        ).read_bytes()

    assignments = pd.read_csv(tmp_path / "first" / "partition/assignments.csv.gz")
    assert {"comparison_id", "scope", "comparison_edge", "unit_id"} <= set(assignments)
    assert not any(column.startswith("Unnamed") for column in assignments)
    assert assignments[["comparison_id", "comparison_edge", "unit_id"]].duplicated().sum() == 0

    summary = pd.read_csv(tmp_path / "first" / "partition/summary.csv")
    assert set(summary["comparison_edge"]) == {"raw_to_recon_expression", "raw_to_final_svc"}

    comparisons = pd.read_csv(tmp_path / "first" / "comparisons.csv")
    assert {
        "comparison_id", "sample_id", "task_cell_type", "scope", "comparison_edge",
        "raw_view", "reconstruction_view", "label_source", "observation_basis",
        "gene_rule", "normalization",
    } <= set(comparisons)


def test_writer_keeps_spatial_baselines_deltas_and_state_gain_missing(tmp_path: Path):
    write_impact_outputs(tmp_path, scopes=[_scope()], audit={"schema_version": 1})
    metrics = pd.read_csv(tmp_path / "spatial/window_metrics.csv")
    neff = metrics.loc[metrics["metric"].eq("neff")].set_index("baseline")
    assert set(neff.index) == {"raw_leiden", "raw_level2"}
    assert neff.loc["raw_leiden", "delta"] == pytest.approx(1.25)
    assert neff.loc["raw_level2", "delta"] == pytest.approx(0.75)
    assert neff.loc["raw_level2", "reconstruction"] == pytest.approx(2.5)

    sensitivity = pd.read_csv(tmp_path / "spatial/sensitivity.csv")
    assert set(sensitivity.loc[sensitivity["sensitivity_type"].eq("scale"), "baseline"]) == {
        "raw_leiden", "raw_level2"
    }

    extent = pd.read_csv(tmp_path / "spatial/region_extent_by_anatomy.csv")
    assert set(extent["region_type"]) == {"state", "gain"}
    gain = extent.loc[extent["region_type"].eq("gain")]
    assert gain["area_um2"].isna().all()

    bootstrap = pd.read_csv(tmp_path / "spatial/threshold_bootstrap.csv")
    assert set(bootstrap["region_type"]) == {"state", "gain"}
    audit = json.loads((tmp_path / "spatial/audit.json").read_text())
    assert audit["state_threshold"]["threshold"] is None
    assert "NaN" not in (tmp_path / "spatial/audit.json").read_text()


def test_writer_rejects_duplicate_scope_ids(tmp_path: Path):
    scope = _scope()
    duplicate = ImpactScope(
        comparison_id=scope.comparison_id,
        sample_id="sample-2",
        task_cell_type=None,
        scope="other",
        comparison=scope.comparison,
    )
    with pytest.raises(ValueError, match="unique comparison_id"):
        write_impact_outputs(tmp_path, scopes=[scope, duplicate], audit={})


def _real_spatial(*, with_level2: bool = True) -> SpatialImpactAnalysis:
    ids = pd.Index([f"u{i}" for i in range(8)])
    coordinates = pd.DataFrame({"x": np.arange(8), "y": np.zeros(8)}, index=ids)
    return compute_spatial_impact(
        full_coordinates=coordinates,
        full_level1_labels=pd.Series(["Tumor"] * 8, index=ids),
        paired_coordinates=coordinates,
        raw_labels=pd.Series(["a", "a", "b", "b", "a", "a", "b", "b"], index=ids),
        raw_level2_labels=(
            pd.Series(["l2a", "l2a", "l2b", "l2b", "l2a", "l2a", "l2b", "l2b"], index=ids)
            if with_level2
            else None
        ),
        reconstructed_labels=pd.Series(["r0", "r0", "r1", "r1", "r0", "r0", "r1", "r1"], index=ids),
        unit_changed=pd.Series([False, True, False, True, False, False, True, False], index=ids),
        microns_per_coordinate=1.0,
        candidate_window_sides_um=[16.0],
        min_parent_units=4,
        rarefaction_draws=5,
        threshold_bootstraps=5,
        random_state=17,
    )


def test_writer_assigns_edge_specific_comparison_ids_and_headline_fields(tmp_path: Path):
    write_impact_outputs(tmp_path, scopes=[_scope()], audit={})
    summary = pd.read_csv(tmp_path / "partition/summary.csv")
    assert {
        "matched_cluster_status", "headline_eligible", "resolution", "scope_comparison_id"
    } <= set(summary)
    assert summary.set_index("comparison_edge").loc["raw_to_recon_expression", "headline_eligible"]
    assert not summary.set_index("comparison_edge").loc["raw_to_final_svc", "headline_eligible"]
    expected = "cmp-main::partition::raw_to_recon_expression"
    assert summary.set_index("comparison_edge").loc["raw_to_recon_expression", "comparison_id"] == expected

    comparisons = pd.read_csv(tmp_path / "comparisons.csv")
    assert "multiple" not in set(comparisons["comparison_edge"])
    assert {expected, "cmp-main"} <= set(comparisons["comparison_id"])
    assert comparisons.loc[comparisons["comparison_id"].eq("cmp-main"), "comparison_edge"].iloc[0] == "spatial_context"


def test_writer_preserves_window_flags_draw_count_and_missing_level2_status(tmp_path: Path):
    scope = ImpactScope(
        comparison_id="cmp-real",
        sample_id="sample-1",
        task_cell_type="T-cell",
        scope="parent-T-cell",
        comparison={},
        spatial=_real_spatial(with_level2=False),
    )
    write_impact_outputs(tmp_path, scopes=[scope], audit={})
    metrics = pd.read_csv(tmp_path / "spatial/window_metrics.csv")
    assert {"in_state_region", "in_gain_region", "unit_change_fraction", "n_changed_units", "n_draws"} <= set(metrics)
    assert (metrics["n_draws"] == 5).all()
    level2 = metrics.loc[metrics["baseline"].eq("raw_level2")]
    assert not level2.empty
    assert (level2["status"] == "unmeasured").all()
    assert not (level2["status"] == "computed").any()


def test_writer_keeps_empty_published_tables_with_explicit_headers(tmp_path: Path):
    empty_mapping = RawLevel2Mapping(
        labels=pd.Series(dtype=object),
        posterior=pd.DataFrame(index=pd.Index([], name="unit_id")),
        assignments=pd.DataFrame(columns=["unit_id", "raw_level2"]),
        audit={"status": "unmeasured"},
    )
    scope = ImpactScope(
        comparison_id="cmp-empty",
        sample_id="sample-1",
        task_cell_type=None,
        scope="global",
        raw_level2=empty_mapping,
    )
    write_impact_outputs(tmp_path, scopes=[scope], audit={})
    posterior = pd.read_csv(tmp_path / "raw_level2/posterior.csv.gz")
    assignments = pd.read_csv(tmp_path / "raw_level2/assignments.csv.gz")
    assert {"comparison_id", "scope_comparison_id", "unit_id", "raw_level2", "posterior"} <= set(posterior)
    assert {"comparison_id", "scope_comparison_id", "unit_id", "raw_level2"} <= set(assignments)


def test_real_spatial_diversity_has_all_metric_baseline_quantiles(tmp_path: Path):
    scope = ImpactScope(
        comparison_id="cmp-real",
        sample_id="sample-1",
        task_cell_type=None,
        scope="global",
        spatial=_real_spatial(with_level2=True),
    )
    write_impact_outputs(tmp_path, scopes=[scope], audit={})
    diversity = pd.read_csv(tmp_path / "spatial/diversity_by_anatomy.csv")
    assert {"baseline", "metric", "raw", "reconstruction", "delta"} <= set(diversity)
    assert {"raw_leiden", "raw_level2"} <= set(diversity["baseline"])
    assert {"k_obs", "entropy", "neff", "evenness"} <= set(diversity["metric"])
    source = _real_spatial(with_level2=True).diversity_by_anatomy
    assert {
        "median_neff_level2", "q1_neff_level2", "q3_neff_level2",
        "median_delta_neff_vs_raw_level2", "q1_delta_neff_vs_raw_level2",
        "q3_delta_neff_vs_raw_level2",
    } <= set(source)
