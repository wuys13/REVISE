"""Contracts for the bilingual saved-record interpretation candidates."""

from __future__ import annotations

import json

from reproduce.case.reconstruction_impact.report_interpretations import (
    NARRATIVE_NODE_IDS,
    build_node_narratives,
)
from reproduce.case.reconstruction_impact.report_records import save_records


def _record(scope: str, question: str, results: dict) -> dict:
    return {
        "identity": {
            "record_id": f"sample::route::{scope}::{question}",
            "scope": scope,
            "question": question,
        },
        "results": results,
    }


def _records() -> list[dict]:
    return [
        _record("All", "complexity", {"raw_k": 99, "reconstruction_k_at_raw_resolution": 100, "ari": 0.01}),
        _record("Fibroblast", "complexity", {"raw_k": 2, "reconstruction_k_at_raw_resolution": 3, "ari": 0.4, "status": "available"}),
        _record("Fibroblast", "matched_k", {"raw_k": 2, "reconstruction_k": 2, "matched": True, "status": "ok", "headline_eligible": True, "unit_change_fraction": 0.25, "balanced_change": 0.3, "ari": 0.5}),
        _record("Fibroblast", "moran_all_valid", {"raw": {"n_valid": 10, "median": 0.1}, "reconstruction": {"n_valid": 12, "median": 0.2}, "median_difference": 0.1}),
        _record("Fibroblast", "moran_shared_valid", {"raw": {"n_valid": 8, "median": 0.11}, "reconstruction": {"n_valid": 8, "median": 0.19}, "paired_delta": {"n": 8, "median": 0.07}}),
        _record("Fibroblast", "emt_coverage", {"raw": {"coverage": 0.1, "available_gene_count": 10, "resource_gene_count": 100}, "reconstruction": {"coverage": 0.9, "available_gene_count": 90, "resource_gene_count": 100}}),
        _record("Fibroblast", "emt_score", {"raw": {"median": 0.2}, "reconstruction": {"median": 0.4}, "paired": {"median_delta": 0.15}, "cutoff_audit": {"raw": {"effective_rank_length": 2, "provider_auc_threshold": 0.2}, "reconstruction": {"effective_rank_length": 9, "provider_auc_threshold": 0.9}}, "carrier_condition": "cluster_mean_projection"}),
        _record("Fibroblast", "anatomy", {"n_regions": 1, "regions": [{"level1_region": "Tumor"}], "support_scales": [{"window_side_length": 16}]}),
        _record("Fibroblast", "changed_units", {"overall": {"changed_units": 3, "paired_units": 10, "change_fraction": 0.3, "scale_um": 16, "n_valid_windows": 4, "matched_cluster_status": "ok", "headline_eligible": True}, "stratification": "anatomy"}),
        _record("Fibroblast", "emt_spatial_fields", {"available": True, "n_units": 10, "features": ["EMT"], "expression_layers": ["score"], "coordinate_units": ["um"], "carrier_condition": "cluster_mean_projection", "carrier_metadata": {"raw_native": True, "reconstruction_cluster_mean_projection": True}}),
        _record("Fibroblast", "window_support", {"main_window_side_um": 16, "min_parent_units": 4, "rarefaction_draws": 20, "support_sensitivity": [{"window_side_length": 16, "n_tissue_windows": 5, "n_valid_windows": 4, "valid_window_fraction": 0.8, "retained_parent_units": 10}]}),
        _record("Fibroblast", "local_vs_raw_leiden", {"scale_um": 16, "n_valid_windows": 4, "metrics": {"Kobs": {"baseline_median": 2, "reconstruction_median": 1.5, "delta_median": -0.5, "delta_n_observations": 4}, "Neff": {"baseline_median": 2, "reconstruction_median": 1.5, "delta_median": -0.5, "delta_n_observations": 4}, "evenness": {"baseline_median": 0.8, "reconstruction_median": 0.8, "delta_median": 0, "delta_n_observations": 4}}, "anatomy_summary": [{"level1_region": "Tumor", "median_delta_k_obs_vs_raw_leiden": -0.5, "median_delta_neff_vs_raw_leiden": -0.5, "median_delta_evenness_vs_raw_leiden": 0, "n_valid_windows": 4}]}),
        _record("Fibroblast", "local_vs_raw_level2", {"scale_um": 16, "n_valid_windows": 4, "metrics": {"Kobs": {"baseline_median": 1, "reconstruction_median": 1.5, "delta_median": 0.5, "delta_n_observations": 4}, "Neff": {"baseline_median": 1, "reconstruction_median": 1.5, "delta_median": 0.5, "delta_n_observations": 4}, "evenness": {"baseline_median": 1, "reconstruction_median": 0.8, "delta_median": -0.2, "delta_n_observations": 4}}, "anatomy_summary": [{"level1_region": "Tumor", "median_delta_k_obs_vs_raw_level2": 0.5, "median_delta_neff_vs_raw_level2": 0.5, "median_delta_evenness_vs_raw_level2": -0.2, "n_valid_windows": 4}]}),
        _record("Fibroblast", "threshold_reliability", {"state": {"status": "ok", "analysis_status": "stable", "threshold": 1.5, "ci_lower": 1.4, "ci_upper": 1.6, "relative_ci_width": 0.2, "n_valid_bootstrap": 20, "n_total_bootstrap": 20, "n_windows": 4}}),
        _record("Fibroblast", "region_extent", {"threshold_status": "ok", "overall": {"region_available": True, "valid_windows": 4, "region_windows": 2, "region_area_mm2": 0.2, "region_area_fraction": 0.5, "valid_units": 10, "region_units": 5, "region_unit_fraction": 0.5}, "by_region": [{"level1_region": "Tumor", "region_windows": 2, "valid_windows": 4, "region_area_mm2": 0.2, "region_area_fraction": 0.5, "region_units": 5, "valid_units": 10, "region_unit_fraction": 0.5, "threshold_status": "ok"}]}),
        _record("Fibroblast", "scale_sensitivity", {"main_scale_um": 16, "rows": [{"window_side_length": 8, "n_valid_windows": 3, "median_delta_neff_vs_raw_leiden": -0.4, "median_delta_neff_vs_raw_level2": 0.3}, {"window_side_length": 16, "n_valid_windows": 4, "median_delta_neff_vs_raw_leiden": -0.5, "median_delta_neff_vs_raw_level2": 0.5}]}),
    ]


def test_narratives_have_fixed_nodes_and_keep_hd_global_out_of_parent_reading():
    narratives = build_node_narratives(_records())

    assert [item["node_id"] for item in narratives] == list(NARRATIVE_NODE_IDS)
    assert all(set(item) == {"node_id", "title", "lead", "interpretations", "next", "record_ids", "review"} for item in narratives)
    assert all(item["review"] == {"status": "pending"} for item in narratives)
    assert all("All" not in " ".join(item["interpretations"]["en"]) for item in narratives)


def test_narratives_use_saved_paired_deltas_and_critical_boundaries():
    by_node = {item["node_id"]: item for item in build_node_narratives(_records())}

    moran = " ".join(by_node["2.4"]["interpretations"]["en"])
    assert "paired median delta=0.07" in moran
    assert "paired median delta=0.08" not in moran

    emt = " ".join(by_node["2.5"]["interpretations"]["en"])
    assert "paired median delta=0.15" in emt
    assert "effective rank length is Raw=2" in emt

    spatial = " ".join(by_node["3.3"]["interpretations"]["en"])
    assert "no region-level quantitative interpretation or enrichment claim" in spatial

    changed = " ".join(by_node["3.2"]["interpretations"]["en"])
    assert "count=3/10" in changed and "fraction=0.3 (30%)" in changed and "valid windows=4" in changed

    local = " ".join(by_node["3.6"]["interpretations"]["en"])
    assert "Raw Leiden baseline 2 to Recon 1.5 (delta=-0.5" in local
    assert "Raw Level2 baseline 1 to Recon 1.5 (delta=0.5" in local

    state = " ".join(by_node["4.2"]["interpretations"]["en"])
    assert "not an EMT hotspot" in state

    scale = " ".join(by_node["4.4"]["interpretations"]["en"])
    assert "8→16 μm" in scale and "-0.4 to -0.5" in scale and "0.3 to 0.5" in scale


def test_save_records_binds_narrative_review_to_package_digest(tmp_path):
    records = _records()
    package = {
        "schema_version": "reconstruction-impact-records/v1",
        "sample_id": "sample",
        "route_kind": "sp_svc",
        "run_identity": {},
        "records": records,
        "node_narratives": build_node_narratives(records),
    }
    destination = tmp_path / "records"
    save_records(package, destination)
    saved = json.loads((destination / "report_records.json").read_text(encoding="utf-8"))
    assert saved["review"]["status"] == "pending"
    assert all(item["review"] == {"status": "pending", "record_digest": saved["record_digest"]} for item in saved["node_narratives"])

    save_records(saved, destination, review_manifest={"status": "reviewed", "record_digest": saved["record_digest"], "reviewer": "test"})
    reviewed = json.loads((destination / "report_records.json").read_text(encoding="utf-8"))
    assert reviewed["review"]["status"] == "reviewed"
    assert all(item["review"]["status"] == "reviewed" and item["review"]["record_digest"] == reviewed["record_digest"] for item in reviewed["node_narratives"])
