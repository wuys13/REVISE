from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise.analysis.reconstruction_impact import (
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
    assert config["spatial_region"]["window_side_length_um"] == 40.0
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


def test_partition_analysis_uses_fixed_within_level1_resolution_and_final_svc(monkeypatch):
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
        resolution_mode="fixed_within_level1",
        within_level1_resolution=0.5,
    )

    assert result.resolution == 0.5
    assert result.resolution_source == "fixed_within_level1"
    assert set(result.comparisons) == {
        "raw_to_recon_expression",
        "recon_expression_to_final_svc",
        "raw_to_final_svc",
    }
