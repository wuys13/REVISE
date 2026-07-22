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
        ".github/workflows/release.yml",
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


def test_moved_application_notebook_uses_repository_relative_paths():
    notebook = (
        ROOT / "reproduce/case/application_sc_SVC_analysis_case.ipynb"
    ).read_text(encoding="utf-8")

    assert "../../raw_data/Real_application" in notebook
    assert "../../output/sc_SVC_case" in notebook
    assert "./reproduce/case" not in notebook
