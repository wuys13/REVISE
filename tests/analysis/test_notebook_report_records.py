"""Focused contracts for the saved-output reconstruction-impact report records."""

from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

from reproduce.case.reconstruction_impact.report_records import (
    QUESTIONS,
    SCOPES,
    _enrich_saved_conclusions,
    build_records,
    format_notebook_records,
    render_report,
    save_records,
)


def _write_csv(root: Path, relative: str, rows: list[dict]) -> None:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False, compression="gzip" if path.suffix == ".gz" else None)


def _write_json(root: Path, relative: str, value) -> None:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2), encoding="utf-8")


def _output(tmp_path: Path, *, hd: bool = True, no_threshold: bool = False) -> Path:
    sample = "P1CRC_VisiumHD" if hd else "P2CRC_Xenium_ThreeParents"
    route = "sp_svc" if hd else "sc_svc"
    output = tmp_path / sample
    analysis = output / "analysis"
    _write_json(
        analysis,
        "audit.json",
        {
            "sample_id": sample,
            "route": route,
            "seed": 42,
            "cohort": "Raw-defined paired scope",
            "expression": "complete gene space per side",
            "region": "high reconstructed diversity",
            "inputs": [{"role": "raw", "sha256": "a" * 64}],
        },
    )
    _write_json(
        analysis,
        "parameters.json",
        {
            "sample": sample,
            "route": route,
            "seed": 42,
            "partition": {"parent_values": ["Fibroblast", "Mono_Macro", "T"]},
            "spatial": {"min_parent_units": 4, "rarefaction_draws": 200, "threshold_bootstraps": 20},
        },
    )
    scopes = ["All", *SCOPES] if hd else list(SCOPES)
    summary_rows = []
    moran_rows = []
    distribution_rows = []
    pathway_rows = []
    availability_rows = []
    local_rows = []
    extent_rows = []
    for scope in scopes:
        matched = scope == "Fibroblast"
        status = "ok" if matched else "unmatched_cluster_complexity"
        summary_rows.append(
            {
                "scope": scope,
                "paired units": 10,
                "Raw K": 2,
                "Recon K": 3 if not matched else 2,
                "ST-unit change": 0.25,
                "balanced change": 0.5,
                "ARI": 0.5,
                "status": status,
                "headline_eligible": matched,
                "union_gene_n": 100,
                "raw_nvalid": 90,
                "raw_median_moran": 0.1,
                "reconstruction_nvalid": 80,
                "reconstruction_median_moran": 0.2,
                "matched_gene_n": 70,
                "shared_raw_nvalid": 70,
                "shared_raw_median": 0.11,
                "shared_reconstruction_nvalid": 70,
                "shared_reconstruction_median": 0.21,
                "paired_delta_n": 70,
                "paired_delta_median_moran": 0.1,
                "raw_unmeasured_n": 0,
                "raw_not_computable_n": 10,
                "raw_insufficient_support_n": 0,
                "reconstruction_unmeasured_n": 5,
                "reconstruction_not_computable_n": 1,
                "reconstruction_insufficient_support_n": 0,
                "pathway": "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                "raw_status": "computed",
                "reconstruction_status": "computed",
                "comparison_status": "computed",
                "raw_n_valid": 10,
                "raw_median_aucell": 0.1,
                "reconstruction_n_valid": 10,
                "reconstruction_median_aucell": 0.2,
                "paired_n": 10,
                "paired_delta_median_aucell": 0.1,
                "scale_um": 40,
                "total_windows": 20,
                "n_valid_windows": 10,
                "valid_units": 15,
                "min_parent_units": 4,
                "rarefaction_draws": 200,
                "local_k_obs_raw_leiden_baseline_median": 2,
                "local_k_obs_raw_leiden_delta_median": -0.5,
                "local_k_obs_raw_leiden_delta_n": 10,
                "local_k_obs_raw_level2_baseline_median": 1.5,
                "local_k_obs_raw_level2_delta_median": 0.1,
                "local_k_obs_raw_level2_delta_n": 10,
                "local_k_obs_reconstruction_median": 1.5,
                "local_neff_raw_leiden_baseline_median": 1.8,
                "local_neff_raw_leiden_delta_median": -0.2,
                "local_neff_raw_leiden_delta_n": 10,
                "local_neff_raw_level2_baseline_median": 1.5,
                "local_neff_raw_level2_delta_median": 0.1,
                "local_neff_raw_level2_delta_n": 10,
                "local_neff_reconstruction_median": 1.6,
                "local_evenness_raw_leiden_baseline_median": 0.9,
                "local_evenness_raw_leiden_delta_median": 0.0,
                "local_evenness_raw_leiden_delta_n": 10,
                "local_evenness_raw_level2_baseline_median": 0.9,
                "local_evenness_raw_level2_delta_median": 0.01,
                "local_evenness_raw_level2_delta_n": 10,
                "local_evenness_reconstruction_median": 0.91,
                "state_threshold_status": "no_stable_threshold" if no_threshold else "ok",
                "state_threshold": None if no_threshold else 1.5,
                "state_n_windows": 10,
                "state_n_valid_bootstrap": 0 if no_threshold else 20,
                "state_region_overall_region_available": False if no_threshold else True,
                "state_region_overall_valid_windows": 10,
                "state_region_overall_region_windows": None if no_threshold else 3,
                "state_region_overall_area_fraction": None if no_threshold else 0.3,
                "state_region_overall_valid_units": 15,
                "state_region_overall_region_units": None if no_threshold else 4,
                "state_region_overall_unit_fraction": None if no_threshold else 0.27,
                "state_region_overall_threshold_status": "no_stable_threshold" if no_threshold else "ok",
            }
        )
        moran_rows.extend(
            [
                {"scope": scope, "gene_set": "all_valid", "side": "Raw", "n_valid": 90, "median": 0.1, "q1": 0.0, "q75": 0.2},
                {"scope": scope, "gene_set": "all_valid", "side": "Reconstruction", "n_valid": 80, "median": 0.2, "q1": 0.1, "q75": 0.3},
                {"scope": scope, "gene_set": "shared_valid", "side": "Raw", "n_valid": 70, "median": 0.11, "q1": 0.01, "q75": 0.21},
                {"scope": scope, "gene_set": "shared_valid", "side": "Reconstruction", "n_valid": 70, "median": 0.21, "q1": 0.11, "q75": 0.31},
            ]
        )
        pathway_rows.append({"scope": scope, "pathway": "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "raw_status": "computed", "reconstruction_status": "computed", "comparison_status": "computed", "raw_n_valid": 10, "raw_median": 0.1, "reconstruction_n_valid": 10, "reconstruction_median": 0.2, "paired_n": 10, "paired_delta_median": 0.1})
        availability_rows.extend(
            [
                {"scope": scope, "side": "raw", "pathway": "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "resource_gene_count": 200, "available_gene_count": 20, "coverage": 0.1, "status": "computed"},
                {"scope": scope, "side": "reconstruction", "pathway": "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "resource_gene_count": 200, "available_gene_count": 180, "coverage": 0.9, "status": "computed"},
            ]
        )
        for metric in ("k_obs", "neff", "evenness"):
            for baseline in ("raw_leiden", "raw_level2"):
                local_rows.append({"scope": scope, "scale_um": 40, "total_windows": 20, "n_valid_windows": 10, "valid_units": 15, "metric": metric, "baseline": baseline, "baseline_n_observations": 10, "baseline_median": 1.0, "reconstruction_n_observations": 10, "reconstruction_median": 0.8, "delta_n_observations": 10, "delta_median": -0.2})
        extent_rows.append({"scope": scope, "level1_region": "Overall", "valid_windows": 10, "region_windows": None if no_threshold else 3, "region_area_um2": None if no_threshold else 300, "area_fraction": None if no_threshold else 0.3, "valid_units": 15, "region_units": None if no_threshold else 4, "unit_fraction": None if no_threshold else 0.27, "region_available": False if no_threshold else True, "threshold_status": "no_stable_threshold" if no_threshold else "ok", "scale_um": 40})

    _write_csv(analysis, "cross_parent_summary.csv", summary_rows)
    if hd:
        complexity_rows = [{"scope": scope, "Raw K": 2, "Recon K at Raw resolution": 3 if scope != "Fibroblast" else 2, "ARI": 0.4} for scope in scopes]
    else:
        complexity_rows = [{"scope": scope, "Raw K at reference resolution": 2, "Fixed final-cluster K": 3 if scope != "Fibroblast" else 2, "ARI": 0.4} for scope in scopes]
    _write_csv(analysis, "reconstruction_impact/partition/complexity_summary.csv", complexity_rows)
    _write_csv(
        analysis,
        "reconstruction_impact/partition/matched_k_summary.csv",
        [{"scope": scope, "paired units": 10, "Raw K": 2, "Recon K": 3 if scope != "Fibroblast" else 2, "ST-unit change": 0.25, "balanced change": 0.5, "ARI": 0.5, "status": "ok" if scope == "Fibroblast" else "unmatched_cluster_complexity", "headline_eligible": scope == "Fibroblast"} for scope in scopes],
    )
    _write_csv(analysis, "moran/moran_summary.csv", summary_rows)
    _write_csv(analysis, "moran/moran_distribution_summary.csv", moran_rows)
    _write_csv(analysis, "pathway_activity/summary.csv", pathway_rows)
    _write_csv(analysis, "pathway_activity/availability.csv", availability_rows)
    _write_csv(analysis, "pathway_activity/cutoff_audit.csv", [])
    _write_csv(analysis, "cross_parent_local_diversity.csv", local_rows)
    _write_csv(analysis, "cross_parent_state_region_extent.csv", extent_rows)
    if hd:
        _write_csv(analysis, "reconstruction_impact/partition/level1_summary.csv", [{"level1": "Fibroblast", "total_units": 5, "changed_units": 1, "change_fraction": 0.2, "scope": "All"}, {"level1": "Overall", "total_units": 10, "changed_units": 3, "change_fraction": 0.3, "scope": "All"}])
    _write_csv(analysis, "anatomy/anatomy_context_summary.csv", [{"level1_region": "Tumor", "full_level1_units": 20, "tissue_windows": 10, "area_um2": 1000, "area_fraction": 1.0}])
    _write_csv(analysis, "anatomy/support_sensitivity.csv", [{"window_side_length": 40, "n_tissue_windows": 20, "n_valid_windows": 10, "valid_window_fraction": 0.5, "retained_parent_units": 15, "retained_parent_unit_fraction": 0.75}])
    for scope in scopes:
        _write_json(analysis, f"reconstruction_impact/partition/{scope}_audit.json", {"audit": {"n_units": 10, "n_shared_genes": 100, "input_shared_genes": 100, "n_features": 20, "feature_selection": "test"}, "matched_cluster_status": "ok" if scope == "Fibroblast" else "unmatched_cluster_complexity"})
        _write_csv(analysis, f"{scope}/inputs/observations.csv.gz", [{"unit_id": f"{scope}-0", "x": 1.0, "y": 2.0, "included": True, "reason": "included"}])
        if scope == "All":
            continue
        _write_csv(analysis, f"local_state/{scope}/change_by_anatomy.csv", [{"level1_region": "Overall", "paired_units": 10, "changed_units": 3, "change_fraction": 0.3, "n_valid_windows": 10, "scale_um": 40, "matched_cluster_status": "ok", "headline_eligible": True}])
        _write_csv(analysis, f"local_state/{scope}/support_sensitivity.csv", [{"window_side_length": 40, "n_tissue_windows": 20, "n_valid_windows": 10, "valid_window_fraction": 0.5, "retained_parent_units": 15, "retained_parent_unit_fraction": 0.75}])
        _write_csv(analysis, f"local_state/{scope}/scale_sensitivity.csv", [{"parent": scope, "window_side_length": 40, "n_valid_windows": 10, "median_delta_neff_vs_raw_level2": 0.1}])
        _write_json(analysis, f"local_state/{scope}/window_decision.json", {"selection": {"status": "ok", "window_side_length": 40, "min_parent_units": 4}, "scale": {"rarefaction_draws": 200}})
        _write_json(analysis, f"local_state/{scope}/state_threshold.json", [{"status": "no_stable_threshold", "threshold": None, "ci_lower": 1.4, "ci_upper": 1.6, "n_windows": 10, "n_valid_bootstrap": 20}] if no_threshold else [{"status": "ok", "threshold": 1.5, "ci_lower": 1.4, "ci_upper": 1.6, "n_windows": 10, "n_valid_bootstrap": 20}])
        _write_json(analysis, f"local_state/{scope}/gain_threshold_audit.json", [{"status": "ok", "threshold": 0.5, "ci_lower": 0.4, "ci_upper": 0.6, "n_windows": 5, "n_valid_bootstrap": 20}])
        _write_csv(analysis, f"local_state/{scope}/window_metrics.csv", [{"valid_window": True, "neff_recon": 1.0}, {"valid_window": True, "neff_recon": 2.0}])
    figures = output / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    for name in ("matched_k_contingency.png", "anatomy_regions.png", "changed_units.png", "moran_shared_gene_scatter.png", "aucell_fibroblast_distribution_and_delta.png", "fibroblast_kobs_matrix.png"):
        (figures / name).write_bytes(b"PNG")
    return output


def _records(package: dict, question: str, scope: str):
    return next(r for r in package["records"] if r["identity"]["question"] == question and r["identity"]["scope"] == scope)


def test_factual_conclusions_use_saved_results_and_keep_region_na():
    records = [
        {
            "identity": {"question": "window_support", "scope": "Fibroblast"},
            "results": {
                "main_window_side_um": 40,
                "min_parent_units": 4,
                "rarefaction_draws": 200,
                "support_sensitivity": [
                    {"window_side_length": 16, "n_tissue_windows": 30, "n_valid_windows": 5, "valid_window_fraction": 1 / 6, "retained_parent_units": 6},
                    {"window_side_length": 40, "n_tissue_windows": 20, "n_valid_windows": 10, "valid_window_fraction": 0.5, "retained_parent_units": 15},
                ],
            },
            "conclusion": {"zh": "template", "en": "template"},
        },
        {
            "identity": {"question": "scale_sensitivity", "scope": "Fibroblast"},
            "results": {
                "main_scale_um": 40,
                "rows": [
                    {"window_side_length": 16, "median_delta_neff_vs_raw_leiden": -0.3, "median_delta_neff_vs_raw_level2": 0.0},
                    {"window_side_length": 40, "median_delta_neff_vs_raw_leiden": -0.2, "median_delta_neff_vs_raw_level2": 0.1},
                ],
            },
            "conclusion": {"zh": "template", "en": "template"},
        },
        {
            "identity": {"question": "region_extent", "scope": "Mono_Macro"},
            "results": {
                "threshold_status": "no_stable_threshold",
                "overall": {
                    "valid_windows": 10,
                    "region_windows": None,
                    "region_area_mm2": None,
                    "region_area_fraction": None,
                    "valid_units": 15,
                    "region_units": None,
                    "region_unit_fraction": None,
                    "region_available": False,
                },
            },
            "conclusion": {"zh": "template", "en": "template"},
        },
        {
            "identity": {"question": "changed_units", "scope": "All"},
            "results": {
                "stratification": "level1",
                "overall": {"level1": "Overall", "total_units": 10, "changed_units": 3, "change_fraction": 0.3},
            },
            "conclusion": {"zh": "template", "en": "template"},
        },
    ]

    returned = _enrich_saved_conclusions(records)

    assert returned is records
    assert "10/20" in records[0]["conclusion"]["en"]
    assert "50%" in records[0]["conclusion"]["en"]
    assert "16–40" in records[1]["conclusion"]["en"]
    assert "-0.3" in records[1]["conclusion"]["en"]
    assert "Region extent remains NA" in records[2]["conclusion"]["en"]
    assert "0 mm²" not in records[2]["conclusion"]["en"]
    assert "0%" not in records[2]["conclusion"]["en"]
    assert "3/10" in records[3]["conclusion"]["en"]
    assert "30%" in records[3]["conclusion"]["en"]


def test_foundation_uses_included_scope_cohort_not_observation_carrier_rows(tmp_path: Path):
    output = _output(tmp_path, hd=False)
    _write_csv(
        output / "analysis",
        "Fibroblast/inputs/observations.csv.gz",
        [
            {"unit_id": "Fibroblast-0", "x": 1.0, "y": 2.0, "included": True, "reason": "included"},
            {"unit_id": "Fibroblast-1", "x": 2.0, "y": 3.0, "included": True, "reason": "included"},
            {"unit_id": "Fibroblast-2", "x": 3.0, "y": 4.0, "included": False, "reason": "not_sampled"},
        ],
    )

    package = build_records(output, "P2CRC_Xenium_ThreeParents", "xenium")
    foundation = _records(package, "foundation", "Fibroblast")

    assert foundation["results"]["input_units"] == 2
    assert foundation["results"]["excluded_units"]["observation_carrier"] == 1
    assert foundation["results"]["input_units"] != 3
    assert "input_units=2" in foundation["conclusion"]["en"]


def test_build_records_has_fixed_order_and_explicit_scope_records(tmp_path: Path):
    package = build_records(_output(tmp_path), "P1CRC_VisiumHD", "visiumhd")

    global_questions = ["foundation", "complexity", "matched_k", "moran_all_valid", "moran_shared_valid", "changed_units"]
    assert [r["identity"]["question"] for r in package["records"]] == [
        *global_questions, *(question for scope in SCOPES for question in QUESTIONS)
    ]
    assert len(package["records"]) == 54
    assert package["review"]["status"] == "pending"
    foundation = _records(package, "foundation", "All")
    assert foundation["results"]["gene_space"] == "full"
    assert foundation["results"]["union_gene_n"] == 100
    assert foundation["identity"]["spec_version"] == package["schema_version"]
    assert foundation["identity"]["rule_id"].endswith("/foundation")
    assert foundation["evidence_condition"]["status"] == "valid"


def test_missing_and_invalid_values_map_to_na_without_zero_coercion(tmp_path: Path):
    package = build_records(_output(tmp_path, no_threshold=True), "P1CRC_VisiumHD", "sp_svc")

    region = _records(package, "region_extent", "Fibroblast")
    threshold = _records(package, "threshold_reliability", "Fibroblast")
    assert region["results"]["overall"]["region_area_um2"] is None
    assert region["results"]["overall"]["region_area_fraction"] is None
    assert threshold["results"]["state"]["threshold"] is None
    assert threshold["results"]["state"]["relative_ci_width"] > 0
    assert threshold["performance_judgment"] != "reliable"


def test_unmatched_matched_k_is_not_headline_and_has_downstream_limitation(tmp_path: Path):
    package = build_records(_output(tmp_path), "P1CRC_VisiumHD", "visiumhd")

    record = _records(package, "matched_k", "Mono_Macro")
    assert record["performance"]["headline_eligible"] is False
    assert "downstream" in " ".join(record["limitations"]["en"]).lower()
    assert record["results"]["raw_k"] == 2
    assert record["results"]["reconstruction_k"] == 3


def test_moran_all_valid_is_marginal_and_shared_valid_is_paired(tmp_path: Path):
    package = build_records(_output(tmp_path), "P1CRC_VisiumHD", "visiumhd")

    all_valid = _records(package, "moran_all_valid", "Fibroblast")
    shared_valid = _records(package, "moran_shared_valid", "Fibroblast")
    assert all_valid["results"]["comparison_kind"] == "marginal_side_specific"
    assert all_valid["results"]["paired_delta"] == {"n": None, "median": None}
    assert all_valid["performance_judgment"] == "side-specific distributions"
    assert shared_valid["results"]["comparison_kind"] == "paired_shared_valid"
    assert shared_valid["results"]["paired_delta"]["n"] == 70
    assert shared_valid["results"]["paired_delta"]["median"] == 0.1


def test_xenium_complexity_and_expression_carrier_scope_are_explicit(tmp_path: Path):
    output = _output(tmp_path, hd=False)
    _write_csv(
        output / "analysis",
        "inputs/carriers.csv",
        [
            {"role": "Raw full context", "n_units": 100},
            {"role": "Fibroblast expression carrier", "n_units": 10},
            {"role": "Mono_Macro expression carrier", "n_units": 20},
            {"role": "T expression carrier", "n_units": 30},
        ],
    )
    package = build_records(output, "P2CRC_Xenium_ThreeParents", "xenium")
    complexity = _records(package, "complexity", "T")
    score = _records(package, "emt_score", "T")
    assert complexity["results"]["comparison_kind"] == "raw_reference_vs_fixed_final_clusters"
    assert complexity["results"]["fixed_final_cluster_k"] == 3
    assert complexity["results"]["reconstruction_k_at_raw_resolution"] is None
    assert [row["role"] for row in score["results"]["expression_carriers"]] == ["Raw full context", "T expression carrier"]
    mono = _records(package, "emt_score", "Mono_Macro")
    assert [row["role"] for row in mono["results"]["expression_carriers"]] == ["Raw full context", "Mono_Macro expression carrier"]
    assert score["results"]["carrier_metadata"]["raw_native"] is True
    assert score["results"]["carrier_metadata"]["reconstruction_cluster_mean_projection"] is True


def test_complexity_and_local_baselines_are_separate_records(tmp_path: Path):
    output = _output(tmp_path)
    _write_csv(
        output / "analysis",
        "local_state/Fibroblast/diversity_by_anatomy.csv",
        [{
            "scope": "Fibroblast",
            "level1_region": "Overall",
            "n_valid_windows": 7,
            "scale_um": 40,
            "median_k_obs_raw": 1,
            "median_k_obs_recon": 2,
            "median_delta_k_obs_vs_raw_leiden": -0.1,
            "median_k_obs_level2": 3,
            "median_delta_k_obs_vs_raw_level2": 0.2,
            "median_neff_raw": 4,
            "median_neff_recon": 5,
            "median_delta_neff_vs_raw_leiden": -0.3,
            "median_neff_level2": 6,
            "median_delta_neff_vs_raw_level2": 0.4,
            "median_evenness_raw": 0.7,
            "median_evenness_recon": 0.8,
            "median_delta_evenness_vs_raw_leiden": 0.1,
            "median_evenness_level2": 0.9,
            "median_delta_evenness_vs_raw_level2": 0.2,
        }],
    )
    package = build_records(output, "P1CRC_VisiumHD", "visiumhd")

    complexity = _records(package, "complexity", "Mono_Macro")
    matched = _records(package, "matched_k", "Mono_Macro")
    local_leiden = _records(package, "local_vs_raw_leiden", "Fibroblast")
    local_level2 = _records(package, "local_vs_raw_level2", "Fibroblast")
    assert complexity["results"]["raw_k"] == 2
    assert complexity["results"]["reconstruction_k_at_raw_resolution"] == 3
    assert matched["results"]["raw_k"] == 2
    assert local_leiden["results"]["baseline"] == "raw_leiden"
    assert local_level2["results"]["baseline"] == "raw_level2"
    assert local_leiden["identity"]["task_cell_type"] == "Fibroblast"
    assert local_leiden["identity"]["baseline"] == "raw_leiden"
    assert complexity["identity"]["baseline"] is None
    assert local_leiden["results"]["metrics"]["Neff"]["delta_median"] == -0.2
    assert local_leiden["performance_judgment"] == "supports concentrated state direction"
    assert local_level2["performance_judgment"] == "does not support requested direction"
    assert local_leiden["results"]["anatomy_summary"] == [{
        "level1_region": "Overall",
        "n_valid_windows": 7,
        "scale_um": 40,
        "median_k_obs_raw": 1,
        "median_k_obs_recon": 2,
        "median_delta_k_obs_vs_raw_leiden": -0.1,
        "median_neff_raw": 4,
        "median_neff_recon": 5,
        "median_delta_neff_vs_raw_leiden": -0.3,
        "median_evenness_raw": 0.7,
        "median_evenness_recon": 0.8,
        "median_delta_evenness_vs_raw_leiden": 0.1,
    }]
    assert local_level2["results"]["anatomy_summary"] == [{
        "level1_region": "Overall",
        "n_valid_windows": 7,
        "scale_um": 40,
        "median_k_obs_level2": 3,
        "median_k_obs_recon": 2,
        "median_delta_k_obs_vs_raw_level2": 0.2,
        "median_neff_level2": 6,
        "median_neff_recon": 5,
        "median_delta_neff_vs_raw_level2": 0.4,
        "median_evenness_level2": 0.9,
        "median_evenness_recon": 0.8,
        "median_delta_evenness_vs_raw_level2": 0.2,
    }]
    assert "analysis/local_state/Fibroblast/diversity_by_anatomy.csv" in [
        item["path"] for item in local_leiden["evidence"]["tables"]
    ]


def test_threshold_keeps_unstable_ci_diagnostics_but_extent_is_na(tmp_path: Path):
    package = build_records(_output(tmp_path, no_threshold=True), "P1CRC_VisiumHD", "sp_svc")

    threshold = _records(package, "threshold_reliability", "Fibroblast")
    extent = _records(package, "region_extent", "Fibroblast")
    assert threshold["results"]["state"]["status"] == "no_stable_threshold"
    assert threshold["results"]["state"]["relative_ci_width"] > 0
    assert extent["results"]["overall"]["region_area_um2"] is None
    assert extent["results"]["overall"]["region_area_fraction"] is None


def test_hd_all_changed_units_uses_level1_strata_and_formatter_is_english(tmp_path: Path):
    package = build_records(_output(tmp_path), "P1CRC_VisiumHD", "visiumhd")

    changed = _records(package, "changed_units", "All")
    assert changed["results"]["stratification"] == "level1"
    assert changed["results"]["overall"]["level1"] == "Overall"
    assert changed["evidence"]["tables"][0]["path"].endswith("partition/level1_summary.csv")
    assert "Overview" in format_notebook_records(package)
    assert "Partition complexity" in format_notebook_records(package, layer=2)


def test_save_records_requires_matching_explicit_review_digest(tmp_path: Path):
    package = build_records(_output(tmp_path), "P1CRC_VisiumHD", "visiumhd")
    destination = tmp_path / "records"
    save_records(package, destination)
    saved = json.loads((destination / "report_records.json").read_text())
    review = json.loads((destination / "review.json").read_text())
    assert saved["record_digest"] == package["record_digest"]
    assert (destination / "report_records.csv").is_file()
    assert review["status"] == "pending"

    approved = {"status": "reviewed", "record_digest": package["record_digest"], "reviewer": "root"}
    save_records(package, destination, review_manifest=approved)
    assert json.loads((destination / "review.json").read_text())["status"] == "reviewed"
    save_records(package, destination, review_manifest={"status": "reviewed", "record_digest": "bad"})
    assert json.loads((destination / "review.json").read_text())["status"] == "pending"


def test_render_report_orders_sections_and_links_each_figure_once(tmp_path: Path):
    package = build_records(_output(tmp_path), "P1CRC_VisiumHD", "visiumhd")
    report_path = render_report(package, tmp_path / "report.html")
    html = report_path.read_text(encoding="utf-8")

    assert html.index('id="overview"') < html.index('id="layer-foundation"') < html.index('id="layer-overall"') < html.index('id="layer-localization"') < html.index('id="layer-region"')
    assert html.index('id="overview"') < html.index("整体对照矩阵")
    assert len(re.findall(r'<img[^>]+src="[^"]*matched_k_contingency\.png"', html)) == 1
    assert "Fibroblast" in html and "Mono_Macro" in html and ">T<" in html
    assert "All" in html
    assert "阈值可靠性" in html or "Threshold reliability" in html


def test_foundation_exposes_saved_carrier_denominator_gene_and_pairing_facts(tmp_path: Path):
    output = _output(tmp_path, hd=False)
    analysis = output / "analysis"
    _write_csv(
        analysis,
        "Fibroblast/inputs/observations.csv.gz",
        [
            {"unit_id": "in", "x": 1.0, "y": 2.0, "raw_level1": "Fibroblast", "reconstruction_covered": True, "included": True, "reason": "included"},
            {"unit_id": "notcovered", "x": 2.0, "y": 3.0, "raw_level1": "Fibroblast", "reconstruction_covered": False, "included": False, "reason": "not_reconstructed"},
            {"unit_id": "notsampled", "x": 3.0, "y": 4.0, "raw_level1": "Fibroblast", "reconstruction_covered": True, "included": False, "reason": "not_sampled"},
            {"unit_id": "outside", "x": 4.0, "y": 5.0, "raw_level1": "Other", "reconstruction_covered": True, "included": True, "reason": "outside_raw_level1"},
        ],
    )
    _write_csv(
        analysis,
        "inputs/carriers.csv",
        [
            {"role": "Raw spatial", "path": "raw.h5ad", "n_units": 2},
            {"role": "Reconstruction expression", "path": "reconstruction.h5ad", "n_units": 2},
        ],
    )
    _write_csv(analysis, "inputs/source_files.csv", [{"role": "raw", "path": "raw.h5ad", "sha256": "a" * 64}])
    _write_csv(
        analysis,
        "inputs/gene_availability.csv",
        [
            {"scope": "Fibroblast", "gene_id": "G1", "side": "raw", "status": "provided", "provided": True},
            {"scope": "Fibroblast", "gene_id": "G2", "side": "raw", "status": "unmeasured", "provided": False},
            {"scope": "Fibroblast", "gene_id": "G1", "side": "reconstruction", "status": "provided", "provided": True},
        ],
    )
    _write_csv(
        analysis,
        "inputs/pairing_audit.csv",
        [{
            "scope": "Fibroblast",
            "coordinate_source": "saved_observations",
            "coordinate_unit": "um",
            "checked_basis": "full_spatial_carrier_before_scope_sampling",
            "rtol": 1e-6,
            "atol": 1e-8,
            "n_checked": 2,
            "id_unique": True,
            "id_subset": True,
            "coordinates_match": True,
            "coordinate_finite": True,
        }],
    )

    package = build_records(output, "P2CRC_Xenium_ThreeParents", "xenium")
    foundation = _records(package, "foundation", "Fibroblast")
    results = foundation["results"]

    assert results["carrier"]["status"] == "valid"
    assert results["carrier"]["rows"][0]["role"] == "Raw spatial"
    assert results["source_files"]["rows"][0]["path"] == "raw.h5ad"
    assert results["cohort"] == "Raw-defined paired scope"
    assert results["cohort_audit"]["label"] == "Raw-defined paired scope"
    assert results["denominator"]["input_units"] == results["input_units"]
    assert results["denominator"]["moran_units"] is None
    assert results["denominator"]["emt_raw_units"] is None
    assert results["selection"]["raw_parent_n"] == 3
    assert results["selection"]["reconstruction_covered_parent_n"] == 2
    assert results["selection"]["eligible_n"] == 2
    assert results["selection"]["sampled_n"] == 1
    assert results["selection"]["stage_exclusions"] == {"outside_raw_level1": 1, "not_reconstructed": 1, "not_sampled": 1}
    assert results["gene"]["gene_space"] == "full"
    assert results["gene"]["by_side"]["raw"]["status_counts"] == {"provided": 1, "unmeasured": 1}
    assert results["pairing_audit"]["coordinate_unit"] == "um"
    assert results["pairing_audit"]["checked_basis"] == "full_spatial_carrier_before_scope_sampling"
    assert results["pairing_audit"]["n_checked"] == 2
    assert results["pairing_audit"]["coordinates_match"] is True


def test_foundation_marks_optional_saved_facts_missing_without_guessing(tmp_path: Path):
    package = build_records(_output(tmp_path, hd=False), "P2CRC_Xenium_ThreeParents", "xenium")
    results = _records(package, "foundation", "Fibroblast")["results"]

    assert results["carrier"]["status"] == "missing"
    assert results["source_files"]["status"] == "missing"
    assert results["gene"]["status"] == "missing"
    assert results["pairing_audit"]["status"] == "missing"
    assert results["pairing_audit"]["coordinate_unit"] is None
    assert results["pairing_audit"]["coordinates_match"] is None


def test_section_summaries_start_pending_and_are_digest_bound(tmp_path: Path):
    package = build_records(_output(tmp_path), "P1CRC_VisiumHD", "visiumhd")
    summaries = package["section_summaries"]

    assert [item["node_id"] for item in summaries] == ["1.6", "2.6", "3.9", "4.5", "summary"]
    assert all(item["conclusion"] == {"zh": "", "en": ""} for item in summaries)
    assert all(item["review"] == {"status": "pending", "record_digest": package["record_digest"]} for item in summaries)
    assert all(set(item["record_ids"]).issubset({record["identity"]["record_id"] for record in package["records"]}) for item in summaries)


def test_stale_section_summary_review_returns_to_pending_on_save(tmp_path: Path):
    package = build_records(_output(tmp_path), "P1CRC_VisiumHD", "visiumhd")
    summary = package["section_summaries"][0]
    summary["conclusion"] = {"zh": "已审阅", "en": "Reviewed"}
    summary["review"] = {"status": "reviewed", "record_digest": "stale"}

    destination = tmp_path / "records"
    save_records(package, destination)
    saved = json.loads((destination / "report_records.json").read_text())

    assert saved["section_summaries"][0]["review"]["status"] == "pending"
    assert saved["section_summaries"][0]["review"]["record_digest"] == saved["record_digest"]


def test_package_review_propagates_to_nonempty_section_summaries(tmp_path: Path):
    package = build_records(_output(tmp_path), "P1CRC_VisiumHD", "visiumhd")
    package["section_summaries"][0]["conclusion"] = {"zh": "已审阅", "en": "Reviewed"}

    destination = tmp_path / "records"
    save_records(package, destination)
    approved_digest = package["record_digest"]
    save_records(package, destination, review_manifest={"status": "reviewed", "record_digest": approved_digest, "reviewer": "root"})
    saved = json.loads((destination / "report_records.json").read_text())

    assert saved["section_summaries"][0]["review"] == {"status": "reviewed", "record_digest": approved_digest, "reviewer": "root"}


def test_foundation_preprocessing_reads_saved_graph_and_emt_cutoff_fields(tmp_path: Path):
    output = _output(tmp_path, hd=False)
    analysis = output / "analysis"
    _write_csv(
        analysis,
        "moran/moran_graph_audit.csv",
        [{
            "scope": "Fibroblast",
            "graph_id": "graph-fibroblast",
            "graph_status": "computed",
            "coordinate_source": "saved_raw_spatial",
            "coordinate_unit": "um",
            "n_units": 2,
            "n_edges_row_normalized": 4,
            "min_units": 51,
            "n_neighbors_actual": 6,
        }],
    )
    _write_csv(
        analysis,
        "pathway_activity/cutoff_audit.csv",
        [
            {"scope": "Fibroblast", "side": "raw", "n_observations": 2, "n_genes": 10, "effective_rank_length": 3, "provider_auc_threshold": 0.3, "cutoff_formula": "saved", "status": "computed"},
            {"scope": "Fibroblast", "side": "reconstruction", "n_observations": 2, "n_genes": 12, "effective_rank_length": 4, "provider_auc_threshold": 0.4, "cutoff_formula": "saved", "status": "computed"},
        ],
    )
    _write_json(analysis, "pathway_activity/resource.json", {"pathway": "EMT", "scorer": "AUCell", "cutoff_policy": "saved", "provider": {"name": "saved", "version": "1"}})

    package = build_records(output, "P2CRC_Xenium_ThreeParents", "xenium")
    preprocessing = _records(package, "foundation", "Fibroblast")["results"]["preprocessing"]

    assert preprocessing["moran"]["status"] == "valid"
    assert preprocessing["moran"]["fields"]["min_units"] == 51
    assert preprocessing["moran"]["fields"]["n_edges_row_normalized"] == 4
    assert preprocessing["emt"]["status"] == "valid"
    assert preprocessing["emt"]["fields"]["cutoff_by_side"]["raw"]["effective_rank_length"] == 3
    assert preprocessing["emt"]["fields"]["cutoff_by_side"]["reconstruction"]["effective_rank_length"] == 4


def test_real_notebook_attachment_routes_each_node_through_record_wrapper(tmp_path: Path):
    from reproduce.case.reconstruction_impact.run_route_notebook import attach_conclusions
    output = _output(tmp_path)
    package = build_records(output, 'P1CRC_VisiumHD', 'sp_svc')
    notebook = {'cells': [
        {'id': 'carrier', 'cell_type': 'markdown', 'metadata': {'impact_node': '1.1'}, 'source': 'placeholder'},
        {'id': 'moran', 'cell_type': 'markdown', 'metadata': {'impact_node': '2.4'}, 'source': 'placeholder'},
        {'id': 'code', 'cell_type': 'code', 'source': 'unchanged = True', 'outputs': []},
    ]}
    target = tmp_path / 'executed.ipynb'
    target.write_text(json.dumps(notebook))
    attach_conclusions(target, package)
    updated = json.loads(target.read_text())
    assert updated['cells'][0]['source'] != updated['cells'][1]['source']
    assert all('Saved conclusions: Overview' not in c['source'] for c in updated['cells'][:2])
    assert updated['cells'][2] == notebook['cells'][2]
    assert 'Observation units' in format_notebook_records(package, node_id='1.1')
