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
        assert "## 1. Question and route semantics" in source
        assert "## 2. Input and spatial-scale audit" in source
        assert "## 3. Expression-state partition impact" in source
        assert "## 4. Level1 anatomy candidates" in source
        assert "## 5. Spatial partition change" in source
        assert "## 6. Local diversity and Region" in source
        assert "## 7. Sensitivity and take-home results" in source
        assert "linear_sum_assignment" not in source
        assert "shannon_entropy" not in source
        assert "Level2" not in source
        assert "REVISE_ANALYSIS_OUTPUT_ROOT" in source

    visium_source = _source(NOTEBOOK_DIR / NOTEBOOKS[0])
    xenium_source = _source(NOTEBOOK_DIR / NOTEBOOKS[1])
    assert "route_kind=CONFIG[\"route_kind\"]" in visium_source
    assert "deterministic_same_id_sample" in visium_source
    assert "representation_audit" in xenium_source
    assert "final SVC" in xenium_source

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
        assert "expr.h5ad" not in text
        assert "cell_equivalent_um: 8.0" in text
        assert "main_window_multiplier: 5" in text
