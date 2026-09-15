"""Execute notebook spatial cells against a small, deterministic spatial case."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
from scipy import sparse

matplotlib.use("Agg")


ROOT = Path(__file__).resolve().parents[2]
BUILDER_PATH = ROOT / "reproduce/case/reconstruction_impact/build_notebooks.py"


def _builder():
    spec = importlib.util.spec_from_file_location("impact_notebook_builder", BUILDER_PATH)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _fixture():
    from revise.analysis.basic.partition_change import compare_partitions
    from revise.analysis.reconstruction_impact import RawLevel2Mapping

    ids, coordinates, raw_labels, recon_labels = [], [], [], []
    for window in range(201):
        x0, y0 = (window % 67) * 10.0, (window // 67) * 10.0
        labels = (("a", "a", "a", "a"), ("a", "b", "a", "b"), ("a", "b", "c", "a"))[window % 3]
        for cell, label in enumerate(labels):
            ids.append(f"u{window}_{cell}")
            coordinates.append((x0 + (cell % 2) * 0.5, y0 + (cell // 2) * 0.5))
            raw_labels.append(("a", "b", "a", "b")[cell])
            recon_labels.append(label)
    # Add one invalid (<4 units) parent window to check its State mask stays NA.
    for cell in range(3):
        ids.append(f"invalid_{cell}")
        coordinates.append((1000 + cell * .2, 1000.0))
        raw_labels.append("a")
        recon_labels.append("a")

    obs = pd.DataFrame({"Level1": "Fibroblast"}, index=ids)
    var = pd.DataFrame(index=["gene_a", "gene_b", "gene_c"])
    raw = ad.AnnData(sparse.csr_matrix(np.tile([[1, 2, 3]], (len(ids), 1))), obs=obs, var=var)
    reconstruction = ad.AnnData(sparse.csr_matrix(np.tile([[3, 2, 1]], (len(ids), 1))), obs=obs.copy(), var=var.copy())
    raw.obsm["spatial"] = np.asarray(coordinates)
    reconstruction.obsm["spatial"] = np.asarray(coordinates)
    comparison = compare_partitions(
        pd.Series(raw_labels, index=ids), pd.Series(recon_labels, index=ids)
    )
    level2 = pd.Series(["l2a" if value == "a" else "l2b" for value in raw_labels], index=ids, name="raw_level2")
    mapping = RawLevel2Mapping(
        labels=level2,
        posterior=pd.DataFrame({"l2a": level2.eq("l2a"), "l2b": level2.eq("l2b")}, index=ids, dtype=float),
        assignments=pd.DataFrame({"unit_id": ids, "raw_level2": level2, "confidence": 1.0}, index=ids),
        audit={
            "source": "synthetic",
            "method": "stub",
            "n_raw_units": len(ids),
            "n_reference_cells": len(ids),
            "n_reference_level2": 2,
            "n_mapped_level2": 2,
            "n_missing": 0,
        },
    )
    return raw, reconstruction, comparison, mapping


def test_spatial_and_region_cells_match_reference_and_keep_invalid_masks_na(tmp_path, monkeypatch):
    """The notebook's explicit stages reproduce the core reference semantics."""
    from revise.analysis.reconstruction_impact import compute_spatial_impact

    builder = _builder()
    raw, reconstruction, comparison, level2 = _fixture()
    ids = comparison.assignments.index
    notebook_helpers = importlib.import_module(
        "reproduce.case.reconstruction_impact.notebook_helpers"
    )
    monkeypatch.setattr(
        notebook_helpers,
        "raw_level2_mapping",
        lambda raw_parent, parent, **kwargs: level2,
    )

    namespace: dict[str, object] = {"__name__": "notebook_cell_test"}
    exec(builder.HELPERS, namespace)
    namespace.update(
        OUTPUT_DIR=tmp_path,
        ROOT=ROOT,
        CONFIG={
            "route_kind": "sc_svc",
            "partition_change": {"random_state": 42},
            "spatial_region": {
                "microns_per_coordinate": 1.0,
                "candidate_window_sides_um": [2.0],
                "threshold_bootstraps": 500,
                "min_parent_units": 4,
                "rarefaction_draws": 200,
                "anatomy_region": {"tumor_label": "Tumor", "normal_source_label": "Normal"},
            },
        },
        CELL_EQUIVALENT_UM=1.0,
        SEED=42,
        PARENT_SOURCE={parent: parent for parent in ("Fibroblast", "Mono_Macro", "T")},
        EMT_NAME="EMT",
        SAMPLE_ID="fixture",
        GMT_PATH=ROOT / "reproduce/case/pathway/h.all.v2025.1.Hs.symbols.gmt",
        SOURCE_FILES={"fixture": Path(__file__)},
        COMPARISONS=[],
        full_coordinates=namespace["coordinates"](raw),
        full_level1=raw.obs["Level1"].copy(),
        TISSUE_XLIM=(0.0, 1002.0),
        TISSUE_YLIM=(0.0, 1002.0),
        ORIGIN_UM=(0.0, 0.0),
        RUNS={
            parent: {
                "raw_expression": raw,
                "recon_expression": reconstruction,
                "comparison": comparison,
                "partition_coordinates": namespace["coordinates"](raw).loc[ids],
                "available": raw.n_obs,
                "sampled": raw.n_obs,
                "used": raw.n_obs,
                "partition": SimpleNamespace(matched_cluster_status="matched_cluster_count"),
            }
            for parent in ("Fibroblast", "Mono_Macro", "T")
        },
        AUCELL=pd.DataFrame(
            {
                "scope": np.repeat(["Fibroblast", "Mono_Macro", "T"], len(ids)),
                "unit_id": np.tile(ids, 3),
                "raw_score": 0.1,
                "reconstruction_score": 0.2,
                "raw_status": "computed",
                "reconstruction_status": "computed",
            }
        ),
    )
    namespace["figures"] = namespace["NotebookFigures"](tmp_path, namespace["TISSUE_XLIM"], namespace["TISSUE_YLIM"], namespace["ORIGIN_UM"], 1.0)
    namespace["figures"].save_figure = lambda fig, *args: namespace["plt"].close(fig)
    namespace["AUCELL"]["delta"] = 0.1
    namespace["display"] = lambda *args, **kwargs: None
    namespace["save_figure"] = lambda figure, *args, **kwargs: namespace["plt"].close(figure)
    namespace["plt"].show = lambda *args, **kwargs: None
    namespace["change_table"] = pd.DataFrame(
        {
            "scope": ["Fibroblast", "Mono_Macro", "T"],
            "status": "matched_cluster_count",
            "headline_eligible": True,
        }
    )

    for cell_type, source in builder.spatial_cells() + builder.region_cells():
        if cell_type == "code":
            exec(source, namespace)

    spatial = namespace["RUNS"]["Fibroblast"]["spatial"]
    expected = compute_spatial_impact(
        full_coordinates=namespace["full_coordinates"],
        full_level1_labels=namespace["full_level1"],
        paired_coordinates=namespace["coordinates"](raw).loc[ids],
        raw_labels=comparison.assignments.set_index("unit_id")["raw_cluster"],
        raw_level2_labels=level2.labels,
        reconstructed_labels=comparison.assignments.set_index("unit_id")["recon_cluster"],
        unit_changed=comparison.assignments.set_index("unit_id")["unit_changed"],
        microns_per_coordinate=1.0,
        candidate_window_sides_um=[2.0],
        min_parent_units=4,
        rarefaction_draws=200,
        threshold_bootstraps=500,
        cell_equivalent_um=1.0,
        anatomy_analysis=namespace["ANATOMY"],
        tumor_label="Tumor",
        normal_source_label="Normal",
    )
    actual = spatial["metrics"].sort_values("window_id").reset_index(drop=True)
    reference = expected.window_metrics.sort_values("window_id").reset_index(drop=True)
    for column in ("k_obs_raw", "k_obs_recon", "k_obs_level2", "neff_raw", "neff_recon", "evenness_recon", "delta_neff_vs_raw_leiden", "delta_neff_vs_raw_level2"):
        np.testing.assert_allclose(actual[column], reference[column], equal_nan=True)
    assert actual.loc[~actual["valid_window"], "in_state_region"].isna().all()
    assert actual.loc[actual["valid_window"], "in_state_region"].notna().all() == (spatial["state_threshold"]["status"] == "ok")
    assert (tmp_path / "analysis/local_state/Fibroblast/window_metrics_with_state_region.csv").is_file()
    saved = pd.read_csv(tmp_path / "analysis/local_state/Fibroblast/unit_window_assignments.csv.gz")
    assert saved["unit_id"].tolist() == ids.tolist()
    level2_labels = pd.read_csv(tmp_path / "analysis/local_state/Fibroblast/raw_level2_labels.csv")
    level2_posterior = pd.read_csv(tmp_path / "analysis/local_state/Fibroblast/raw_level2_posterior.csv.gz")
    assert level2_labels.columns[0] == level2_posterior.columns[0] == "unit_id"
    assert level2_labels["unit_id"].tolist() == level2_posterior["unit_id"].tolist() == ids.tolist()
    anatomy_decision = json.loads((tmp_path / "analysis/anatomy/window_decision.json").read_text())
    assert {"main_window_side_um", "origin_x_um", "origin_y_um"} <= set(anatomy_decision["scale"])
    diversity_by_anatomy = pd.read_csv(tmp_path / "analysis/local_state/Fibroblast/diversity_by_anatomy.csv")
    change_by_anatomy = pd.read_csv(tmp_path / "analysis/local_state/Fibroblast/change_by_anatomy.csv")
    assert "median_delta_neff_vs_raw_leiden" in diversity_by_anatomy
    assert {"scope", "scale_um", "matched_cluster_status", "headline_eligible"} <= set(change_by_anatomy)

    namespace.update(
        change_table=namespace["change_table"].assign(partition_value=1.0),
        MORAN_SUMMARY=pd.DataFrame({
            "scope": ["Fibroblast", "Mono_Macro", "T"],
            "comparison_id": ["Fibroblast", "Mono_Macro", "T"],
            "moran_value": 2.0,
        }),
        aucell_summary=pd.DataFrame({"scope": ["Fibroblast", "Mono_Macro", "T"], "aucell_value": 3.0}),
    )
    final_summary_cell = builder.final_cells()[1]
    assert final_summary_cell[0] == "code"
    exec(final_summary_cell[1], namespace)

    local_summary = pd.read_csv(tmp_path / "analysis/cross_parent_local_diversity.csv")
    cross_parent = pd.read_csv(tmp_path / "analysis/cross_parent_summary.csv")
    region_extent = pd.read_csv(tmp_path / "analysis/cross_parent_state_region_extent.csv")
    assert len(local_summary) == 18
    assert set(local_summary["baseline"]) == {"raw_leiden", "raw_level2"}
    assert {"state_threshold_status", "state_region_overall_area_fraction"} <= set(cross_parent)
    assert set(region_extent["scope"]) == {"Fibroblast", "Mono_Macro", "T"}
