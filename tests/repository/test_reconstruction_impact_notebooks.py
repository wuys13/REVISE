from __future__ import annotations

import json
from pathlib import Path


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
            "## 4. Anatomy overview",
            "## 5. Fibroblast internal diversity and Regions",
            "## 6. Mono_Macro internal diversity and Regions",
            "## 7. T internal diversity and Regions",
            "## 8. Cross-parent summary and evidence boundary",
        ]
        positions = [source.index(heading) for heading in headings]
        assert positions == sorted(positions)
        assert "linear_sum_assignment" not in source
        assert "shannon_entropy" not in source
        assert "Level2" not in source
        assert "Neff≥2" not in source
        assert "main_window_multiplier" not in source
        assert "matched-K" in source
        assert "Raw-Level1 / Recon subtype / Recon−1" in source
        assert "Raw-Leiden / Recon / Delta" in source
        assert "State Region" in source
        assert "Gain Region" in source
        assert "IPython.display import Markdown" in source

        notebook = _notebook(path)
        assert all(cell.get("id") for cell in notebook["cells"])
        assert "invert_yaxis" in source
        assert "anatomy context" in source.lower()
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
    assert "raw_context[reconstructed_ids, partition.feature_names]" not in visium_source
    assert "representation_audit" in xenium_source
    assert all(parent in xenium_source for parent in ("Fibroblast", "Mono_Macro", "T"))
    assert "expression_h5ad" in xenium_source

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
        assert "candidate_window_sides_um" in text
        assert "min_parent_units: 4" in text

    xenium = (ROOT / "configs" / "analysis" / "reconstruction_impact_xenium_p2crc_fibroblast.yaml").read_text(encoding="utf-8")
    assert all(parent in xenium for parent in ("Fibroblast", "Mono_Macro", "T"))
    assert "expr.h5ad" in xenium


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
        "matched_k_resolution_sweep.csv",
        "change_by_level1.csv",
    ):
        assert filename in source
