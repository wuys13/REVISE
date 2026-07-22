from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def test_biological_metrics_is_an_installed_console_entry():
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")

    assert (
        'revise-compute-biological-metrics = "revise.analysis.cli:main"'
        in pyproject
    )


def test_biological_metrics_cli_delegates_to_package_implementation():
    module = ROOT / "revise/analysis/cli.py"

    assert module.is_file()
    source = module.read_text(encoding="utf-8")
    assert "def main(" in source
    assert "compute_conditional_moran_i(" in source
    assert "compute_local_label_entropy(" in source
    assert "compute_identity_metrics(" in source
    assert "compute_tmp_mer(" in source
