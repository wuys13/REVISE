from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def test_benchmark_entrypoints_delegate_to_package_owned_implementations():
    cli = ROOT / "revise" / "benchmark" / "cli.py"
    launcher = ROOT / "revise" / "benchmark" / "launcher.py"
    wrapper = ROOT / "reproduce" / "benchmark_main.py"
    shell = ROOT / "reproduce" / "benchmark_main.sh"
    wrapper_source = wrapper.read_text(encoding="utf-8")
    shell_source = shell.read_text(encoding="utf-8")

    assert cli.is_file()
    assert launcher.is_file()
    assert "from revise.benchmark.cli import main" in wrapper_source
    assert 'parents[1]' in wrapper_source
    assert "revise/benchmark/launcher.py" in shell_source
    assert not (ROOT / "benchmark_launcher.py").exists()
    assert not (ROOT / "benchmark_main.py").exists()
    assert not (ROOT / "benchmark_main.sh").exists()

    launcher_source = launcher.read_text(encoding="utf-8")
    assert '"-m",' in launcher_source
    assert '"revise.benchmark.cli",' in launcher_source
    assert '"benchmark_main.py"' not in launcher_source


def test_benchmark_python_wrapper_finds_checkout_package_outside_repository(
    tmp_path,
):
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "reproduce" / "benchmark_main.py"),
            "--help",
        ],
        cwd=tmp_path,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "--confounding" in result.stdout
