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
    "global_anchoring",
    "local_refinement",
    "output",
    "execution",
]


EXPECTED = {
    "Xenium_T.yaml": {
        "svc_type": "sc-SVC",
        "sample_name": "P1CRC",
        "ot_method": "tacco",
        "st_path": "raw_data/Real_application/P1CRC_Xenium.h5ad",
        "reference_path": "raw_data/Real_application/adata_sc_all_reanno.h5ad",
        "local_refinement": {"subtype_column": "Level2", "select_cell_type": "T"},
    },
    "Xenium_Fib.yaml": {
        "svc_type": "sc-SVC",
        "sample_name": "P1CRC",
        "ot_method": "tacco",
        "st_path": "raw_data/Real_application/P1CRC_Xenium.h5ad",
        "reference_path": "raw_data/Real_application/adata_sc_all_reanno.h5ad",
        "local_refinement": {
            "subtype_column": "Level2",
            "select_cell_type": "Fibroblast",
        },
    },
    "Xenium_Mono.yaml": {
        "svc_type": "sc-SVC",
        "sample_name": "P1CRC",
        "ot_method": "tacco",
        "st_path": "raw_data/Real_application/P1CRC_Xenium.h5ad",
        "reference_path": "raw_data/Real_application/adata_sc_all_reanno.h5ad",
        "local_refinement": {
            "subtype_column": "Level2",
            "select_cell_type": "Mono_Macro",
        },
    },
    "VisiumHD.yaml": {
        "svc_type": "sp-SVC",
        "sample_name": "P1CRC",
        "ot_method": "pot",
        "st_path": "raw_data/Real_application/P1CRC_HD.h5ad",
        "reference_path": "raw_data/Real_application/adata_sc_all_reanno.h5ad",
        "local_refinement": {"strength": 0.2},
    },
    "Visium.yaml": {
        "svc_type": "sc-SVC-sr",
        "sample_name": "REVISEVisiumMouseBrain",
        "ot_method": "pot",
        "st_path": "raw_data/visium_mouse_brain/ST_mouse_brain_prepared.h5ad",
        "reference_path": "raw_data/visium_mouse_brain/scRNA_mouse_brain_prepared.h5ad",
        "local_refinement": {"strength": 0.0},
    },
}


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
    assert document["application"] == {
        "svc_type": expected["svc_type"],
        "sample_name": expected["sample_name"],
    }
    assert document["paths"] == {"root_dir": "."}
    assert document["algorithm"] == {"ot_method": expected["ot_method"]}
    assert document["inputs"] == {
        "mode": "direct",
        "st": {"path": expected["st_path"], "format": "h5ad"},
        "reference": {
            "path": expected["reference_path"],
            "format": "h5ad",
            "patient_key": "Patient",
        },
    }
    assert document["global_anchoring"] == {"broad_column": "Level1"}
    assert document["local_refinement"] == expected["local_refinement"]
    assert document["output"] == {"path": "output"}
    assert document["execution"] == {"action": "run", "seed": 42}
    assert "base_config" not in package_bytes.decode("utf-8")


@pytest.mark.parametrize("filename", ("Xenium_T.yaml", "Xenium_Fib.yaml", "Xenium_Mono.yaml"))
def test_xenium_templates_explain_that_the_selected_broad_type_must_exist(filename):
    text = (PACKAGE_TEMPLATES / filename).read_text(encoding="utf-8")

    assert "broad label exists" in text


def test_visium_template_explains_that_zero_strength_runs_local_refinement():
    text = (PACKAGE_TEMPLATES / "Visium.yaml").read_text(encoding="utf-8")

    assert "LR still runs" in text
