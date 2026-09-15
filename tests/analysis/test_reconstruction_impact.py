from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise.analysis.reconstruction_impact import (
    compute_anatomy_regions,
    compute_spatial_impact,
    file_sha256,
    load_reconstruction_impact_config,
    map_raw_level2_labels,
    run_partition_analysis,
    select_raw_level1_parent_cohort,
    write_anatomy_artifacts,
    write_raw_level2_artifacts,
)


def test_config_loader_rejects_expression_side_as_spatial_comparison_input(tmp_path: Path):
    config_path = tmp_path / "invalid.yaml"
    config_path.write_text(
        """
schema_version: 1
sample: {id: example}
partition_change:
  comparisons:
    - name: Fibroblast
      raw_h5ad: raw.h5ad
      reconstructed_spatial_h5ad: results/Fibroblast/expr.h5ad
""".strip(),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="expression-side"):
        load_reconstruction_impact_config(config_path)


def test_file_sha256_is_content_based(tmp_path: Path):
    path = tmp_path / "input.txt"
    path.write_text("paired input\n", encoding="utf-8")

    assert file_sha256(path) == "0e597ce4292cf36f6d5fab492f7566e32a9052d5c1192fd4f52fbfde9fce0ae9"


def test_config_loader_resolves_defaults_and_preserves_explicit_output(tmp_path: Path):
    config_path = tmp_path / "valid.yaml"
    config_path.write_text(
        """
schema_version: 1
sample: {id: example}
partition_change: {comparisons: []}
output: {dir: output/reconstruction_impact/example}
""".strip(),
        encoding="utf-8",
    )

    config = load_reconstruction_impact_config(config_path)

    assert config["partition_change"]["level1_resolution_candidates"] == [0.3, 0.5, 0.8]
    assert config["partition_change"]["within_level1_resolution"] == 0.5
    assert config["spatial_region"]["candidate_window_sides_um"] == [16.0, 24.0, 32.0, 40.0, 56.0, 80.0]
    assert config["spatial_region"]["min_parent_units"] == 4
    assert config["spatial_region"]["anatomy_region"]["tumor_label"] == "Tumor"
    assert config["output"]["dir"] == "output/reconstruction_impact/example"


def test_raw_level2_mapping_filters_reference_and_uses_only_raw_parent(monkeypatch):
    from revise.analysis import reconstruction_impact

    raw = AnnData(
        X=np.array([[2.0, 0.0], [0.0, 2.0]]),
        obs=pd.DataFrame({"Level1": ["Mono/Macro", "Mono/Macro"]}, index=["u0", "u1"]),
        var=pd.DataFrame(index=["g0", "g1"]),
    )
    reference = AnnData(
        X=np.ones((4, 2)),
        obs=pd.DataFrame(
            {
                "Patient": ["P2CRC", "P2CRC", "P1CRC", "P2CRC"],
                "Level1": ["Mono/Macro", "Mono/Macro", "Mono/Macro", "T"],
                "Level2": ["Mono", "Macro_C1QC", "Macro_SPP1", "CD4T"],
            },
            index=["r0", "r1", "r2", "r3"],
        ),
        var=pd.DataFrame(index=["g0", "g1"]),
    )
    captured = {}

    def fake_annotate(target, selected_reference, **kwargs):
        captured["target_X"] = target.X.copy()
        captured["reference_obs"] = selected_reference.obs.copy()
        captured["kwargs"] = kwargs
        result = target.copy()
        result.obs["Level2"] = ["Mono", "Macro_C1QC"]
        result.obs["Level2_confidence"] = [0.8, 0.9]
        result.obsm["Level2"] = pd.DataFrame(
            [[0.8, 0.2], [0.1, 0.9]],
            index=result.obs_names,
            columns=["Mono", "Macro_C1QC"],
        )
        return result

    monkeypatch.setattr(reconstruction_impact.OTKernel, "annotate", fake_annotate)

    mapping = map_raw_level2_labels(
        raw,
        reference,
        parent_value="Mono_Macro",
        method="tacco",
        reference_filter_column="Patient",
        reference_filter_value="P2CRC",
        tacco_multi_center=1,
        tacco_lamb=0.001,
    )

    np.testing.assert_array_equal(captured["target_X"], raw.X)
    assert captured["reference_obs"]["Patient"].eq("P2CRC").all()
    assert captured["reference_obs"]["Level1"].str.replace("/", "_", regex=False).eq("Mono_Macro").all()
    assert captured["kwargs"]["method"] == "tacco"
    assert mapping.labels.to_dict() == {"u0": "Mono", "u1": "Macro_C1QC"}
    assert mapping.audit["source"] == "raw_expression_reference_level2"
    assert mapping.audit["n_reference_level2"] == 2


def test_raw_level1_parent_cohort_explicitly_audits_carrier_mismatches():
    raw_level1 = pd.Series(
        ["Fibroblast", "T", "Fibroblast"], index=["u0", "u1", "u2"]
    )

    selected, audit = select_raw_level1_parent_cohort(
        raw_level1, pd.Index(["u0", "u1", "u2"]), parent_value="Fibroblast"
    )

    assert selected.tolist() == ["u0", "u2"]
    assert audit == {
        "parent": "Fibroblast",
        "carrier_units": 3,
        "raw_level1_parent_units": 2,
        "excluded_raw_level1_mismatch": 1,
    }


def test_partition_analysis_keeps_level1_ari_complexity_diagnostic_separate_from_matched_k(monkeypatch):
    from revise.analysis import reconstruction_impact

    raw = AnnData(
        X=np.ones((4, 2)),
        obs=pd.DataFrame({"Level1": ["A", "A", "B", "B"]}, index=["u0", "u1", "u2", "u3"]),
        var=pd.DataFrame(index=["g0", "g1"]),
    )
    recon = AnnData(
        X=np.ones((4, 2)),
        obs=pd.DataFrame({"Level1": ["A", "A", "B", "B"]}, index=["u0", "u1", "u2", "u3"]),
        var=pd.DataFrame(index=["g0", "g1"]),
    )
    monkeypatch.setattr(reconstruction_impact, "prepare_leiden_graph", lambda adata, **_: adata)

    def fake_labels(adata, *, resolution, random_state):
        if resolution == 0.3:
            return pd.Series(["0", "1", "0", "1"], index=adata.obs_names)
        return pd.Series(["0", "0", "1", "1"], index=adata.obs_names)

    monkeypatch.setattr(reconstruction_impact, "leiden_labels", fake_labels)

    result = run_partition_analysis(
        raw,
        recon,
        level1_col="Level1",
        resolution_candidates=[0.3, 0.5],
        raw_qc_min_genes=1,
        raw_qc_min_cells=1,
    )

    assert result.complexity_resolution == 0.5
    assert result.sweep["target_cluster_count"].iloc[0] == 2
    assert result.resolution_source == "matched_cluster_count"
    assert result.sweep["n_clusters"].iloc[0] == 2
    assert result.matched_k_sweep is result.sweep
    assert {"resolution", "ARI"} <= set(result.raw_complexity_sweep.columns)
    assert list(result.complexity_comparisons) == ["raw_to_recon_same_resolution"]
    assert list(result.comparisons) == ["raw_to_recon_expression"]


def test_partition_analysis_accepts_preselected_shared_features(monkeypatch):
    from revise.analysis import reconstruction_impact

    raw = AnnData(
        X=np.ones((4, 2)),
        obs=pd.DataFrame({"Level1": ["A", "A", "B", "B"]}, index=["u0", "u1", "u2", "u3"]),
        var=pd.DataFrame(index=["g0", "g1"]),
    )
    recon = raw.copy()
    monkeypatch.setattr(reconstruction_impact, "prepare_leiden_graph", lambda adata, **_: adata)
    monkeypatch.setattr(
        reconstruction_impact,
        "leiden_labels",
        lambda adata, *, resolution, random_state: pd.Series(["0", "0", "1", "1"], index=adata.obs_names),
    )
    monkeypatch.setattr(
        reconstruction_impact,
        "select_shared_feature_names",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("must not reselect features")),
    )

    result = run_partition_analysis(
        raw,
        recon,
        level1_col="Level1",
        feature_names=["g0", "g1"],
        resolution_candidates=[0.3],
        raw_qc_min_genes=1,
        raw_qc_min_cells=1,
    )

    assert result.feature_names == ["g0", "g1"]
    assert result.audit["feature_selection"] == "preselected_shared_features"


def test_sp_svc_partition_uses_raw_qc_and_raw_derived_features(monkeypatch):
    from revise.analysis import reconstruction_impact

    raw = AnnData(
        X=np.ones((4, 2)),
        obs=pd.DataFrame({"Level1": ["A", "A", "B", "B"]}, index=["u0", "u1", "u2", "u3"]),
        var=pd.DataFrame(index=["g0", "g1"]),
    )
    recon = raw.copy()
    captured = {}

    def fake_filter(left, right, **kwargs):
        captured["qc"] = kwargs
        return left, right, {"n_units": left.n_obs, "excluded_raw_qc_units": 0}

    def fake_raw_hvg(adata, **kwargs):
        captured["raw_hvg"] = (adata, kwargs)
        return ["g0", "g1"]

    monkeypatch.setattr(reconstruction_impact, "filter_paired_sp_svc_inputs", fake_filter)
    monkeypatch.setattr(reconstruction_impact, "select_raw_hvg_feature_names", fake_raw_hvg)
    monkeypatch.setattr(reconstruction_impact, "prepare_leiden_graph", lambda adata, **_: adata)
    monkeypatch.setattr(
        reconstruction_impact,
        "leiden_labels",
        lambda adata, *, resolution, random_state: pd.Series(["0", "0", "1", "1"], index=adata.obs_names),
    )

    result = run_partition_analysis(
        raw,
        recon,
        level1_col="Level1",
        route_kind="sp_svc",
        resolution_candidates=[0.3],
        raw_qc_min_genes=2,
        raw_qc_min_cells=2,
    )

    assert captured["qc"] == {"min_genes": 2, "min_cells": 2}
    assert captured["raw_hvg"][0] is raw
    assert result.audit["feature_selection"] == "raw_canonical_seurat_v3_shared"


def test_sc_svc_partition_uses_fixed_resolution_and_omits_identity_expression_edge(monkeypatch):
    from revise.analysis import reconstruction_impact

    raw = AnnData(
        X=np.ones((4, 2)),
        obs=pd.DataFrame({"Level1": ["Fibroblast"] * 4}, index=["u0", "u1", "u2", "u3"]),
        var=pd.DataFrame(index=["g0", "g1"]),
    )
    recon = raw.copy()
    recon.obs["SVC_cluster"] = ["a", "a", "b", "b"]
    monkeypatch.setattr(reconstruction_impact, "prepare_leiden_graph", lambda adata, **_: adata)
    monkeypatch.setattr(
        reconstruction_impact,
        "leiden_labels",
        lambda adata, *, resolution, random_state: pd.Series(["0", "0", "1", "1"], index=adata.obs_names),
    )

    result = run_partition_analysis(
        raw,
        recon,
        level1_col="Level1",
        final_cluster_key="SVC_cluster",
        route_kind="sc_svc",
        resolution_mode="fixed_within_level1",
        within_level1_resolution=0.5,
    )

    assert result.complexity_resolution == 0.5
    assert result.resolution_source == "matched_cluster_count"
    assert list(result.comparisons) == ["raw_to_final_svc"]
    assert result.representation_audit["spatial_expression_identical"] is True


def test_spatial_impact_keeps_window_coordinates_after_anatomy_context_join():
    ids = pd.Index(["u0", "u1", "u2", "u3"])
    full_coordinates = pd.DataFrame({"x": [0.0, 1.0, 40.0, 41.0], "y": [0.0, 1.0, 0.0, 1.0]}, index=ids)
    impact = compute_spatial_impact(
        full_coordinates=full_coordinates,
        full_level1_labels=pd.Series(["Tumor", "Intestinal Epithelial", "Other", "Other"], index=ids),
        paired_coordinates=full_coordinates,
        raw_labels=pd.Series(["0", "0", "1", "1"], index=ids),
        reconstructed_labels=pd.Series(["0", "1", "1", "1"], index=ids),
        unit_changed=pd.Series([False, True, False, False], index=ids),
        microns_per_coordinate=1.0,
        candidate_window_sides_um=[40.0],
        min_parent_units=2,
        rarefaction_draws=5,
        threshold_bootstraps=5,
    )

    assert {"window_x", "window_y", "level1_region"} <= set(impact.window_metrics.columns)
    assert set(impact.anatomy_context_summary["level1_region"]) == {
        "Interface",
        "Other",
    }
    assert "Overall" in set(impact.cluster_change_by_anatomy["level1_region"])
    assert "Overall" in set(impact.diversity_by_anatomy["level1_region"])
    assert "Overall" in set(impact.region_extent_by_anatomy["level1_region"])
    assert "in_state_region" in impact.window_metrics
    assert "in_gain_region" in impact.window_metrics


def test_spatial_impact_carries_raw_level2_through_paired_rarefaction():
    ids = pd.Index(["u0", "u1", "u2", "u3"])
    coordinates = pd.DataFrame(
        {"x": [0.0, 1.0, 2.0, 3.0], "y": [0.0, 0.0, 0.0, 0.0]}, index=ids
    )
    impact = compute_spatial_impact(
        full_coordinates=coordinates,
        full_level1_labels=pd.Series(["Tumor"] * 4, index=ids),
        paired_coordinates=coordinates,
        raw_labels=pd.Series(["0", "0", "1", "1"], index=ids),
        raw_level2_labels=pd.Series(["L2a", "L2b", "L2a", "L2b"], index=ids),
        reconstructed_labels=pd.Series(["R0", "R0", "R1", "R1"], index=ids),
        unit_changed=pd.Series([False, False, False, False], index=ids),
        microns_per_coordinate=1.0,
        candidate_window_sides_um=[16.0],
        min_parent_units=4,
        rarefaction_draws=5,
        threshold_bootstraps=5,
    )

    assert {
        "k_obs_level2",
        "neff_level2",
        "evenness_level2",
        "delta_neff_vs_raw_level2",
    } <= set(impact.window_metrics.columns)


def test_route_level_anatomy_is_independent_of_parent_specific_window_selection():
    ids = pd.Index([f"u{i}" for i in range(8)])
    full_coordinates = pd.DataFrame(
        {"x": [0, 1, 2, 3, 40, 41, 42, 43], "y": [0] * 8}, index=ids
    )
    labels = pd.Series(
        ["Tumor", "Tumor", "Intestinal Epithelial", "Intestinal Epithelial"] + ["Other"] * 4,
        index=ids,
    )
    anatomy = compute_anatomy_regions(
        full_coordinates=full_coordinates,
        full_level1_labels=labels,
        microns_per_coordinate=1.0,
        candidate_window_sides_um=[4.0, 16.0],
        min_parent_units=2,
        cell_equivalent_um=8.0,
    )
    parent_ids = ids[:4]
    impact = compute_spatial_impact(
        full_coordinates=full_coordinates,
        full_level1_labels=labels,
        paired_coordinates=full_coordinates.loc[parent_ids],
        raw_labels=pd.Series(["a", "a", "b", "b"], index=parent_ids),
        reconstructed_labels=pd.Series(["a", "b", "b", "b"], index=parent_ids),
        unit_changed=pd.Series([False, True, False, False], index=parent_ids),
        microns_per_coordinate=1.0,
        candidate_window_sides_um=[16.0],
        min_parent_units=2,
        rarefaction_draws=3,
        threshold_bootstraps=5,
        anatomy_analysis=anatomy,
    )

    assert anatomy.scale_audit["main_window_side_um"] == 4.0
    assert impact.scale_audit["main_window_side_um"] == 16.0
    assert impact.window_metrics["level1_region"].eq("Interface").all()


def test_parent_window_uses_its_observed_anatomy_units_not_an_empty_geometric_center():
    ids = pd.Index(["u0", "u1", "u2", "u3"])
    coordinates = pd.DataFrame(
        {"x": [0.0, 0.1, 15.0, 15.1], "y": [0.0, 0.1, 0.0, 0.1]}, index=ids
    )
    anatomy = compute_anatomy_regions(
        full_coordinates=coordinates,
        full_level1_labels=pd.Series(["Tumor", "Tumor", "Other", "Other"], index=ids),
        microns_per_coordinate=1.0,
        candidate_window_sides_um=[4.0],
        min_parent_units=2,
    )

    impact = compute_spatial_impact(
        full_coordinates=coordinates,
        full_level1_labels=pd.Series(["Tumor", "Tumor", "Other", "Other"], index=ids),
        paired_coordinates=coordinates,
        raw_labels=pd.Series(["a", "a", "b", "b"], index=ids),
        reconstructed_labels=pd.Series(["a", "b", "b", "b"], index=ids),
        unit_changed=pd.Series([False, True, False, False], index=ids),
        microns_per_coordinate=1.0,
        candidate_window_sides_um=[16.0],
        min_parent_units=2,
        rarefaction_draws=3,
        threshold_bootstraps=5,
        anatomy_analysis=anatomy,
    )

    assert impact.window_metrics["level1_region"].tolist() == ["Tumor"]


def test_route_level_anatomy_artifacts_are_written_once(tmp_path: Path):
    ids = pd.Index(["u0", "u1", "u2", "u3"])
    anatomy = compute_anatomy_regions(
        full_coordinates=pd.DataFrame({"x": [0, 1, 2, 3], "y": [0, 0, 0, 0]}, index=ids),
        full_level1_labels=pd.Series(["Tumor", "Tumor", "Other", "Other"], index=ids),
        microns_per_coordinate=1.0,
        candidate_window_sides_um=[4.0],
        min_parent_units=2,
    )

    destination = write_anatomy_artifacts(tmp_path, anatomy)

    assert destination == tmp_path / "anatomy"
    assert (destination / "scale_audit.json").is_file()
    assert (destination / "scale_decision.csv").is_file()
    assert (destination / "window_assignments.csv.gz").is_file()
    assert (destination / "region_summary.csv").is_file()


def test_raw_level2_mapping_artifacts_preserve_source_audit(tmp_path: Path):
    result = type("Mapping", (), {})()
    result.labels = pd.Series(["L2a", "L2b"], index=["u0", "u1"], name="raw_level2")
    result.posterior = pd.DataFrame(
        [[0.8, 0.2], [0.1, 0.9]], index=["u0", "u1"], columns=["L2a", "L2b"]
    )
    result.assignments = pd.DataFrame(
        {
            "unit_id": ["u0", "u1"],
            "raw_level2": ["L2a", "L2b"],
            "confidence": [0.8, 0.9],
            "mapping_method": ["pot", "pot"],
            "source": ["raw_expression_reference_level2"] * 2,
        },
        index=["u0", "u1"],
    )
    result.audit = {"source": "raw_expression_reference_level2", "method": "pot"}

    destination = write_raw_level2_artifacts(tmp_path, result)

    assert destination == tmp_path / "raw_level2"
    assert (destination / "assignments.csv.gz").is_file()
    assert (destination / "posterior.csv.gz").is_file()
    assert '"source": "raw_expression_reference_level2"' in (
        destination / "audit.json"
    ).read_text(encoding="utf-8")


@pytest.mark.filterwarnings("error:Mean of empty slice:RuntimeWarning")
def test_gain_region_threshold_only_uses_positive_matched_delta(monkeypatch):
    from revise.analysis import reconstruction_impact

    captured = []

    def record_threshold(values, **_kwargs):
        captured.append(np.asarray(values))
        return {"status": "no_stable_threshold", "threshold": None, "n_windows": len(values), "n_valid_bootstrap": 0}, pd.DataFrame()

    monkeypatch.setattr(reconstruction_impact, "select_region_threshold", record_threshold)
    ids = pd.Index([f"u{i}" for i in range(8)])
    coordinates = pd.DataFrame({"x": np.arange(8), "y": np.zeros(8)}, index=ids)
    compute_spatial_impact(
        full_coordinates=coordinates,
        full_level1_labels=pd.Series(["Tumor"] * 8, index=ids),
        paired_coordinates=coordinates,
        raw_labels=pd.Series(["a"] * 8, index=ids),
        reconstructed_labels=pd.Series(["a", "a", "b", "b", "a", "a", "b", "b"], index=ids),
        unit_changed=pd.Series([False] * 8, index=ids),
        microns_per_coordinate=1.0,
        candidate_window_sides_um=[16.0],
        min_parent_units=4,
        rarefaction_draws=5,
        threshold_bootstraps=5,
    )

    assert len(captured) == 2
    assert (captured[1] > 0).all()


def test_spatial_impact_propagates_and_audits_configured_spatial_seed(monkeypatch):
    from revise.analysis import reconstruction_impact
    from revise.analysis.basic import spatial_region

    captured = {"rarefaction": [], "threshold": []}
    original_rarefaction = spatial_region.compute_rarefied_window_diversity

    def record_rarefaction(*args, **kwargs):
        captured["rarefaction"].append(kwargs["random_state"])
        return original_rarefaction(*args, **kwargs)

    def record_threshold(values, **kwargs):
        captured["threshold"].append(kwargs["random_state"])
        return (
            {
                "status": "no_stable_threshold",
                "threshold": None,
                "n_windows": len(values),
                "n_valid_bootstrap": 0,
            },
            pd.DataFrame(),
        )

    monkeypatch.setattr(reconstruction_impact, "compute_rarefied_window_diversity", record_rarefaction)
    monkeypatch.setattr(reconstruction_impact, "select_region_threshold", record_threshold)
    ids = pd.Index([f"u{i}" for i in range(8)])
    coordinates = pd.DataFrame({"x": np.arange(8), "y": np.zeros(8)}, index=ids)

    impact = compute_spatial_impact(
        full_coordinates=coordinates,
        full_level1_labels=pd.Series(["Tumor"] * 8, index=ids),
        paired_coordinates=coordinates,
        raw_labels=pd.Series(["a"] * 8, index=ids),
        raw_level2_labels=pd.Series(["l2"] * 8, index=ids),
        reconstructed_labels=pd.Series(["a", "a", "b", "b", "a", "a", "b", "b"], index=ids),
        unit_changed=pd.Series([False] * 8, index=ids),
        microns_per_coordinate=1.0,
        candidate_window_sides_um=[16.0],
        min_parent_units=4,
        rarefaction_draws=5,
        threshold_bootstraps=5,
        random_state=17,
    )

    assert captured["rarefaction"] == [17, 17]
    assert captured["threshold"] == [17, 17]
    assert impact.scale_audit["rarefaction_random_state"] == 17
    assert impact.scale_audit["threshold_random_state"] == 17


def test_scale_sensitivity_uses_rarefied_three_assignment_metrics():
    ids = pd.Index([f"u{i}" for i in range(8)])
    coordinates = pd.DataFrame({"x": np.arange(8), "y": np.zeros(8)}, index=ids)

    impact = compute_spatial_impact(
        full_coordinates=coordinates,
        full_level1_labels=pd.Series(["Tumor"] * 8, index=ids),
        paired_coordinates=coordinates,
        raw_labels=pd.Series(["a"] * 8, index=ids),
        raw_level2_labels=pd.Series(["l2"] * 8, index=ids),
        reconstructed_labels=pd.Series(["a", "a", "b", "b", "a", "a", "b", "b"], index=ids),
        unit_changed=pd.Series([False] * 8, index=ids),
        microns_per_coordinate=1.0,
        candidate_window_sides_um=[16.0],
        min_parent_units=4,
        rarefaction_draws=5,
        threshold_bootstraps=5,
        random_state=17,
    )

    sensitivity = impact.scale_sensitivity
    assert sensitivity.loc[0, "rarefaction_draws"] == 5
    assert sensitivity.loc[0, "random_state"] == 17
    assert "median_delta_neff_vs_raw_leiden" in sensitivity
    assert "median_delta_neff_vs_raw_level2" in sensitivity
    main = impact.window_metrics.loc[impact.window_metrics["valid_window"]]
    for column in ("neff_recon", "delta_neff_vs_raw_leiden", "delta_neff_vs_raw_level2"):
        assert sensitivity.loc[0, f"median_{column}"] == pytest.approx(main[column].median())
