from __future__ import annotations

import json
from pathlib import Path
import runpy


ROOT = Path(__file__).resolve().parents[2]
NOTEBOOK_DIR = ROOT / "reproduce" / "case" / "reconstruction_impact"
NOTEBOOKS = (
    "VisiumHD_sp_SVC_Reconstruction_Impact.ipynb",
    "Xenium_sc_SVC_Fibroblast_Reconstruction_Impact.ipynb",
)


def _source(path: Path) -> str:
    notebook = json.loads(path.read_text(encoding="utf-8"))
    return "\n".join("".join(cell.get("source", [])) for cell in notebook["cells"])


def _notebook(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_reconstruction_impact_notebooks_are_route_scoped_and_use_internal_authority():
    for name in NOTEBOOKS:
        path = NOTEBOOK_DIR / name
        assert path.is_file()
        source = _source(path)
        assert "REVISE_ANALYSIS_OUTPUT_ROOT" in source
        assert (
            "revise.analysis.basic" in source
            or "revise.analysis.reconstruction_impact" in source
        )
        assert "Evidence boundary" in source
        headings = [
            "## 1. Route semantics and carrier audit",
            "## 2. Partition complexity diagnostic",
            "## 3. Matched-complexity change and Level1 localization",
            "## 4. Level1 anatomy Regions",
            "## 5. Fibroblast internal diversity and Regions",
            "## 6. Mono_Macro internal diversity and Regions",
            "## 7. T internal diversity and Regions",
            "## 8. Cross-parent summary and evidence boundary",
        ]
        positions = [source.index(heading) for heading in headings]
        assert positions == sorted(positions)
        assert "linear_sum_assignment" not in source
        assert "shannon_entropy" not in source
        assert "Raw Level2" in source
        assert "Neff≥2" not in source
        assert "main_window_multiplier" not in source
        assert "matched-K" in source
        assert "Raw Level1" in source
        assert "Raw Leiden" in source
        assert "High-diversity Region" in source
        assert "Diversity-gain Region" not in source
        assert "IPython.display import Markdown" in source

        notebook = _notebook(path)
        assert all(cell.get("id") for cell in notebook["cells"])
        assert "invert_yaxis" not in source
        assert "pcolormesh" in source
        assert "2, 2" in source
        assert "cells per side" in source
        assert "anatomy regions" in source.lower()
        assert "Region" in source
        for cell in notebook["cells"]:
            if cell["cell_type"] != "code":
                continue
            cell_source = "".join(cell.get("source", []))
            assert cell_source.count("save_figure(") <= 1

    visium_source = _source(NOTEBOOK_DIR / NOTEBOOKS[0])
    xenium_source = _source(NOTEBOOK_DIR / NOTEBOOKS[1])
    assert 'route_kind="sp_svc"' in visium_source
    assert "deterministic_same_id_sample" in visium_source
    assert "USE_FULL_VISIUMHD_COHORT = False" in visium_source
    assert "VISIUMHD_SAMPLE_N_UNITS = 30_000" in visium_source
    assert "global_ids" in visium_source
    assert "Raw-QC retained units" in visium_source
    assert "raw_qc_min_genes" in visium_source
    assert "raw_qc_min_cells" in visium_source
    assert "Raw-derived Seurat-v3 HVGs" in visium_source
    assert "parent-internal reassignment" in visium_source
    assert 'ax.barh(plot_table["level1"]' in visium_source
    assert "focus_level1" in visium_source
    assert "raw_context[reconstructed_ids, partition.feature_names]" not in visium_source
    assert "representation_audit" in xenium_source
    assert all(parent in xenium_source for parent in ("Fibroblast", "Mono_Macro", "T"))
    assert "expression_h5ad" in xenium_source


def test_reconstruction_impact_notebooks_expose_window_evidence_before_each_parent_metric():
    for name in NOTEBOOKS:
        source = _source(NOTEBOOK_DIR / name)
        for parent in ("Fibroblast", "Mono_Macro", "T"):
            heading = f"## {5 + (parent != 'Fibroblast') + (parent == 'T')}. {parent} internal diversity and Regions"
            start = source.index(heading)
            end = source.find("\n## ", start + 1)
            block = source[start:] if end == -1 else source[start:end]
            assert block.index("Cohort and window decision") < block.index("Internal-state baselines")
            assert "Kobs" in block
            assert "Neff" in block
            assert "evenness" in block
            assert 'plot_metric_comparison(RUN, "k_obs"' in block
            assert 'plot_metric_comparison(RUN, "neff"' in block
            assert 'plot_metric_comparison(RUN, "evenness"' in block
            assert "High-diversity Region" in block

    helper_source = _source(NOTEBOOK_DIR / NOTEBOOKS[0])
    assert "plt.subplots(2, 3" in helper_source
    assert "plt.subplots(1, 2" in helper_source
    assert "region_overlay" not in helper_source

    assert not (NOTEBOOK_DIR / "01_partition_change.ipynb").exists()
    assert not (NOTEBOOK_DIR / "02_spatial_diversity.ipynb").exists()
    assert not (NOTEBOOK_DIR / "03_region_and_sensitivity.ipynb").exists()


def test_reconstruction_impact_docs_and_configs_keep_the_post_analysis_boundary():
    design_dir = ROOT / "docs" / "design" / "reconstruction-impact"
    assert (design_dir / "README.md").is_file()
    assert (design_dir / "outputs-and-test-plan.md").is_file()
    for name in (
        "reconstruction_impact_visiumhd_p1crc.yaml",
        "reconstruction_impact_xenium_p2crc_fibroblast.yaml",
    ):
        path = ROOT / "configs" / "analysis" / name
        assert path.is_file()
        text = path.read_text(encoding="utf-8")
        assert "schema_version: 1" in text
        assert "output/reconstruction_impact" in text
        assert "final_notebook:" in text
        assert "candidate_window_sides_um" in text
        assert "min_parent_units: 4" in text
        assert "raw_level2_mapping:" in text
        assert "adata_sc_all_reanno.h5ad" in text

    xenium = (ROOT / "configs" / "analysis" / "reconstruction_impact_xenium_p2crc_fibroblast.yaml").read_text(encoding="utf-8")
    assert all(parent in xenium for parent in ("Fibroblast", "Mono_Macro", "T"))
    assert "expr.h5ad" in xenium
    assert "method: tacco" in xenium
    assert "value: P2CRC" in xenium
    visium = (ROOT / "configs" / "analysis" / "reconstruction_impact_visiumhd_p1crc.yaml").read_text(encoding="utf-8")
    assert "method: pot" in visium
    assert "raw_qc_min_genes: 50" in visium
    assert "raw_qc_min_cells: 3" in visium
    assert "feature_selection: raw_canonical_seurat_v3_shared" in visium


def test_final_reconstruction_impact_notebook_is_written_only_under_output():
    runner = NOTEBOOK_DIR / "run_route_notebook.py"
    assert runner.is_file()
    source = runner.read_text(encoding="utf-8")
    assert "output/reconstruction_impact/<sample_id>/notebook" in source
    assert "final_notebook" in source
    assert "notebook" in source
    assert "reproduce/case/reconstruction_impact" in source

    route_paths = runpy.run_path(str(runner))["route_paths"]
    _, output_root, final = route_paths("visiumhd")
    assert output_root == (ROOT / "output/reconstruction_impact/P1CRC_VisiumHD").resolve()
    assert final.parent == output_root / "notebook"


def test_reconstruction_impact_output_contract_names_metric_layers():
    source = (ROOT / "revise" / "analysis" / "reconstruction_impact.py").read_text(
        encoding="utf-8"
    )
    for filename in (
        "anatomy_context_summary.csv",
        "cluster_change_by_anatomy.csv",
        "diversity_by_anatomy.csv",
        "region_extent_by_anatomy.csv",
        "gain_region_extent_by_anatomy.csv",
        "state_threshold.json",
        "gain_threshold.json",
        "matched_k_resolution_sweep.csv",
        "change_by_level1.csv",
    ):
        assert filename in source
