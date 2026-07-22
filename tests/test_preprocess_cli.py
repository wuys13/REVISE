from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_histology_preprocessor_is_an_installed_console_entry():
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")

    assert (
        'revise-build-histology-priors = "revise.preprocess.cli:main"'
        in pyproject
    )


def test_histology_preprocessor_package_cli_exposes_help():
    module = ROOT / "revise/preprocess/cli.py"

    assert module.is_file()
    source = module.read_text(encoding="utf-8")
    assert "def main(" in source
    assert '"--st-h5ad"' in source
    assert '"--mask"' in source
    assert '"--out-h5ad"' in source
