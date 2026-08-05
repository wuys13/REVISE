from pathlib import Path
from types import SimpleNamespace
import subprocess
import sys

import pytest

from revise import __version__


ROOT = Path(__file__).resolve().parents[2]


def test_package_cli_exposes_version_without_a_config():
    result = subprocess.run(
        [sys.executable, "-m", "revise.application.cli", "--version"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0
    assert result.stdout.strip() == f"revise-reconstruct {__version__}"


def test_application_cli_rejects_benchmark_only_flags():
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "revise.application.cli",
            "--config",
            "run.yaml",
            "--spot-size",
            "100",
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 2
    assert "unrecognized arguments: --spot-size 100" in result.stderr


def test_source_wrapper_delegates_to_the_package_cli():
    import reconstruct
    from revise.application import cli

    assert reconstruct.main is cli.main
    assert not hasattr(reconstruct, "build_parser")


def test_execution_error_is_one_actionable_line_without_traceback(
    monkeypatch,
    capsys,
):
    from revise.application import cli, service

    monkeypatch.setattr(
        cli,
        "_resolve_application_source",
        lambda path: SimpleNamespace(path=None, payload=b"", label="test"),
    )
    monkeypatch.setattr(cli, "load_application_request", lambda *a, **k: object())

    def fail(_request):
        raise RuntimeError("missing input data/sample.h5ad")

    monkeypatch.setattr(service, "execute_application", fail)

    with pytest.raises(SystemExit) as exc_info:
        cli.main(["--config", "VisiumHD.yaml"])

    assert exc_info.value.code == 1
    stderr = capsys.readouterr().err
    assert stderr.strip() == (
        "revise-reconstruct: error: missing input data/sample.h5ad"
    )
    assert "Traceback" not in stderr
