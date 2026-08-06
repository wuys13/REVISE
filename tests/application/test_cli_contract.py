from pathlib import Path
import subprocess
import sys

import pytest

from revise import __version__


ROOT = Path(__file__).resolve().parents[2]


def test_source_entrypoint_exposes_version_without_a_config():
    result = subprocess.run(
        [sys.executable, "reconstruct.py", "--version"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0
    assert result.stdout.strip() == f"revise-reconstruct {__version__}"


def test_help_shows_application_flow_groups():
    import reconstruct

    help_text = reconstruct.build_parser().format_help()
    for group in (
        "Application:",
        "Inputs:",
        "Shared OT:",
        "Global Anchoring:",
        "Local Refinement:",
        "Output:",
        "Execution:",
    ):
        assert group in help_text


def test_execution_error_is_one_actionable_line_without_traceback(monkeypatch, capsys):
    import reconstruct

    def fail(*args, **kwargs):
        raise RuntimeError("missing input data/sample.h5ad")

    monkeypatch.setattr(reconstruct, "run_application", fail)
    with pytest.raises(SystemExit) as exc_info:
        reconstruct.main(["--config", "VisiumHD.yaml"])

    assert exc_info.value.code == 1
    stderr = capsys.readouterr().err
    assert stderr.strip() == "revise-reconstruct: error: missing input data/sample.h5ad"
    assert "Traceback" not in stderr
