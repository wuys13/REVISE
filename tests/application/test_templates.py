from __future__ import annotations

from pathlib import Path

import pytest
import yaml


ROOT = Path(__file__).resolve().parents[2]
PACKAGE_TEMPLATES = ROOT / "revise" / "application" / "templates"
SOURCE_TEMPLATES = ROOT / "configs" / "application"
TOP_LEVEL_KEYS = [
    "schema_version",
    "application",
    "paths",
    "algorithm",
    "inputs",
    "preprocessing",
    "global_anchoring",
    "local_refinement",
    "output",
    "execution",
]


EXPECTED = {
    "VisiumHD.yaml": {
        "svc_type": "sp-SVC",
        "mode": None,
        "output_name": "sp_SVC",
        "output_dir": "output/sp_SVC_case/P1CRC",
        "ot_method": "pot",
        "st_path": "raw_data/Real_application/P1CRC_HD.h5ad",
        "reference_path": "raw_data/Real_application/adata_sc_all_reanno.h5ad",
        "preprocessing": {
            "spatial": {
                "min_transcript_counts": None,
                "min_counts": 20,
                "min_cell_counts": 30,
            },
            "reference": {
                "min_transcript_counts": None,
                "min_genes": 20,
                "min_cell_counts": 50,
            },
        },
        "local_refinement": {"strength": 0.2},
    },
    "Xenium.yaml": {
        "svc_type": "sc-SVC",
        "mode": "cluster",
        "output_name": None,
        "output_dir": "output/P2CRC_Xenium",
        "ot_method": "tacco",
        "st_path": "raw_data/Real_application/P2CRC_Xenium.h5ad",
        "reference_path": "raw_data/Real_application/adata_sc_all_reanno.h5ad",
        "reference_filter": {"filter_column": "Patient", "filter_value": "P2CRC"},
        "preprocessing": {
            "spatial": {"min_transcript_counts": 60, "min_cell_counts": 100},
            "reference": {"min_transcript_counts": None, "min_cell_counts": 100},
        },
        "local_refinement": {
            "subtype_column": "Level2",
            "select_cell_type": "T",
            "alpha": 0.2,
            "resolutions": [0.6, 0.7, 0.8],
        },
    },
    "Visium.yaml": {
        "svc_type": "sc-SVC",
        "mode": "sr",
        "output_name": "REVISEVisiumMouseBrain_sc-SVC",
        "output_dir": "output/visium_mouse_brain_revise",
        "ot_method": "pot",
        "st_path": "raw_data/visium_mouse_brain/ST_mouse_brain_prepared.h5ad",
        "reference_path": "raw_data/visium_mouse_brain/scRNA_mouse_brain_prepared.h5ad",
        "pm_on_cell_path": "raw_data/visium_mouse_brain/PM_on_cell.csv",
        "preprocessing": {
            "spatial": {
                "min_transcript_counts": None,
                "min_counts": 20,
                "min_cell_counts": 20,
            },
            "reference": {
                "min_transcript_counts": None,
                "min_genes": 20,
                "min_cell_counts": 20,
            },
        },
        "local_refinement": {
            "strength": 0.0,
            "graph": {
                "method": "pca",
                "alpha": 0.2,
                "n_neighbors": 10,
                "exp_neighbors": 10,
                "spatial_neighbors": 10,
            },
            "match_spot_sum": True,
        },
    },
}


def test_official_application_template_sets_are_exact_and_mirrored():
    expected = set(EXPECTED)

    assert {path.name for path in PACKAGE_TEMPLATES.glob("*.yaml")} == expected
    assert {path.name for path in SOURCE_TEMPLATES.glob("*.yaml")} == expected


@pytest.mark.parametrize("filename", EXPECTED)
def test_packaged_template_is_canonical_and_source_mirror_is_byte_exact(filename):
    package_path = PACKAGE_TEMPLATES / filename
    source_path = SOURCE_TEMPLATES / filename

    package_bytes = package_path.read_bytes()
    assert source_path.read_bytes() == package_bytes

    document = yaml.safe_load(package_bytes)
    expected = EXPECTED[filename]
    assert list(document) == TOP_LEVEL_KEYS
    assert document["schema_version"] == 1
    application = {"svc_type": expected["svc_type"]}
    if expected["mode"] is not None:
        application["mode"] = expected["mode"]
    assert document["application"] == application
    assert document["paths"] == {"root_dir": "."}
    expected_inputs = {
        "st": {"path": expected["st_path"], "format": "h5ad"},
        "reference": {
            "path": expected["reference_path"],
            "format": "h5ad",
            **expected.get("reference_filter", {}),
        },
    }
    if "pm_on_cell_path" in expected:
        expected_inputs["pm_on_cell"] = {"path": expected["pm_on_cell_path"]}
    assert document["algorithm"] == {"ot_method": expected["ot_method"]}
    assert document["inputs"] == expected_inputs
    assert document["preprocessing"] == expected["preprocessing"]
    assert document["global_anchoring"] == {"broad_column": "Level1"}
    assert document["local_refinement"] == expected["local_refinement"]
    assert document["output"] == {
        "dir": expected.get("output_dir", "output"),
        **({} if expected["output_name"] is None else {"name": expected["output_name"]}),
    }
    assert document["execution"] == {"seed": 42}
    assert "base_config" not in package_bytes.decode("utf-8")


def test_xenium_template_explains_that_the_selected_broad_type_must_exist():
    text = (PACKAGE_TEMPLATES / "Xenium.yaml").read_text(encoding="utf-8")

    assert "broad label exists" in text


def test_visium_template_explains_that_zero_strength_runs_local_refinement():
    text = (PACKAGE_TEMPLATES / "Visium.yaml").read_text(encoding="utf-8")

    assert "LR still runs" in text
