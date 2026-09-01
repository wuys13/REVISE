from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
NOTEBOOK_DIR = ROOT / "reproduce" / "case" / "reconstruction_impact"
NOTEBOOKS = (
    "01_partition_change.ipynb",
    "02_spatial_diversity.ipynb",
    "03_region_and_sensitivity.ipynb",
)


def _source(path: Path) -> str:
    notebook = json.loads(path.read_text(encoding="utf-8"))
    return "\n".join("".join(cell.get("source", [])) for cell in notebook["cells"])


def test_reconstruction_impact_notebooks_are_new_and_use_internal_authority():
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
        assert "linear_sum_assignment" not in source
        assert "shannon_entropy" not in source

    region_source = _source(NOTEBOOK_DIR / "03_region_and_sensitivity.ipynb")
    assert "anatomy_region_map.csv.gz" in region_source
    assert "High-diversity Region" in region_source

    diversity_source = _source(NOTEBOOK_DIR / "02_spatial_diversity.ipynb")
    assert "anatomy_coords" in diversity_source
    assert "anatomy_windows" in diversity_source
    assert "anatomy_region_map.csv.gz" in diversity_source
    assert "index_label='unit_id'" in diversity_source
    assert "grid_origin" in diversity_source


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
