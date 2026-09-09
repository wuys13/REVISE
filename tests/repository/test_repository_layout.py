from __future__ import annotations

from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize(
    "relative",
    [
        "scripts",
        "tools",
        "release",
        "MIGRATION.md",
        "legacy-assets.json",
        "application_sp_SVC_recon.py",
        "application_sp_SVC_recon.sh",
        "application_sc_SVC_recon.py",
        "application_sc_SVC_recon.sh",
    ],
)
def test_main_repository_excludes_migration_and_custom_release_surfaces(relative):
    assert not (ROOT / relative).exists()


def test_build_and_ci_do_not_reference_removed_maintenance_surfaces():
    manifest = (ROOT / "MANIFEST.in").read_text(encoding="utf-8")
    ci = (ROOT / ".github/workflows/ci.yml").read_text(encoding="utf-8")

    for removed in ("scripts/", "tools/", "release/", "legacy-assets.json"):
        assert removed not in manifest
        assert removed not in ci


def test_repository_root_exposes_only_the_reconstruction_script():
    assert {path.name for path in ROOT.glob("*.py")} == {"reconstruct.py"}
    assert not (ROOT / "benchmark_main.sh").exists()


def test_biological_metrics_remain_owned_by_the_analysis_package():
    module = ROOT / "revise/analysis/biological_metrics.py"

    assert module.is_file()
    source = module.read_text(encoding="utf-8")
    for function in (
        "compute_conditional_moran_i",
        "compute_local_label_entropy",
        "compute_identity_metrics",
        "compute_tmp_mer",
    ):
        assert f"def {function}(" in source


def test_case_notebooks_are_the_canonical_application_gallery():
    case_dir = ROOT / "reproduce" / "case"
    canonical = {
        "CosMx_SMI_267T_not_sp_SVC.ipynb",
        "MERFISH_Allen_VISp_sc_SVC_cluster.ipynb",
        "SlideSeq_mouse_colon_sp_SVC.ipynb",
        "SlideSeq_mouse_olfactory_bulb_sp_SVC.ipynb",
        "StereoSeq_zebrafish_5hpf_sp_SVC.ipynb",
        "osmFISH_sc_SVC_cluster.ipynb",
        "Xenium_sc_SVC_T.ipynb",
        "Xenium_sc_SVC_Fibroblast.ipynb",
        "Xenium_sc_SVC_Monocyte.ipynb",
        "VisiumHD_sp_SVC.ipynb",
        "Visium_sc_SVC_mouse_brain.ipynb",
    }
    retired = {
        "sc_SVC_case_T_analysis.ipynb",
        "sc_SVC_case_Fibroblast_analysis.ipynb",
        "sc_SVC_case_Monocyte_analysis.ipynb",
        "sp_SVC_case.ipynb",
        "sc_SVC_case_T_recon.ipynb",
        "sc_SVC_case_Fibroblast_recon.ipynb",
        "sc_SVC_case_Monocyte_recon.ipynb",
        "application_sc_SVC_analysis_case.ipynb",
    }

    assert {path.name for path in case_dir.glob("*.ipynb")} == canonical
    assert all(not (case_dir / name).exists() for name in retired)


def test_public_notebook_indexes_use_only_canonical_application_names():
    indexes = [
        ROOT / "reproduce/README.md",
        ROOT / "docs/source/gallery.rst",
        ROOT / "docs/source/quickstart.rst",
    ]
    combined = "\n".join(path.read_text(encoding="utf-8") for path in indexes)

    for canonical in (
        "Xenium_sc_SVC_T.ipynb",
        "Xenium_sc_SVC_Fibroblast.ipynb",
        "Xenium_sc_SVC_Monocyte.ipynb",
        "VisiumHD_sp_SVC.ipynb",
        "Visium_sc_SVC_mouse_brain.ipynb",
    ):
        assert canonical in combined
    for retired in (
        "sc_SVC_case_T_recon.ipynb",
        "sc_SVC_case_Fibroblast_recon.ipynb",
        "sc_SVC_case_Monocyte_recon.ipynb",
        "application_sc_SVC_analysis_case.ipynb",
        "sc_SVC_sr_case_Visium_mouse_brain.ipynb",
    ):
        assert retired not in combined
