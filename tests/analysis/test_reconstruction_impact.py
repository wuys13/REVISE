from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise.analysis.reconstruction_impact import (
    compute_spatial_impact,
    file_sha256,
    load_reconstruction_impact_config,
    run_partition_analysis,
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
    assert config["spatial_region"]["cell_equivalent_um"] == 8.0
    assert config["spatial_region"]["main_window_multiplier"] == 5
    assert config["spatial_region"]["anatomy_region"]["tumor_label"] == "Tumor"
    assert config["output"]["dir"] == "output/reconstruction_impact/example"


def test_partition_analysis_uses_raw_level1_ari_resolution_then_shared_resolution(monkeypatch):
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

    result = run_partition_analysis(raw, recon, level1_col="Level1", resolution_candidates=[0.3, 0.5])

    assert result.resolution == 0.5
    assert result.resolution_source == "raw_level1_ari"
    assert result.sweep["ARI"].tolist() == [pytest.approx(-0.5), pytest.approx(1.0)]
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
    )

    assert result.feature_names == ["g0", "g1"]


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

    assert result.resolution == 0.5
    assert result.resolution_source == "fixed_within_level1"
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
        cell_equivalent_um=8.0,
        main_window_multiplier=5,
    )

    assert {"window_x", "window_y", "level1_region"} <= set(impact.window_metrics.columns)
