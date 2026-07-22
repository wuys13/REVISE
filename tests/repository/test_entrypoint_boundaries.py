from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def test_benchmark_entrypoints_delegate_to_package_owned_implementations():
    cli = ROOT / "revise" / "benchmark" / "cli.py"
    launcher = ROOT / "revise" / "benchmark" / "launcher.py"
    root_cli = (ROOT / "benchmark_main.py").read_text(encoding="utf-8")
    root_shell = (ROOT / "benchmark_main.sh").read_text(encoding="utf-8")

    assert cli.is_file()
    assert launcher.is_file()
    assert "from revise.benchmark.cli import main" in root_cli
    assert "revise/benchmark/launcher.py" in root_shell
    assert not (ROOT / "benchmark_launcher.py").exists()

    launcher_source = launcher.read_text(encoding="utf-8")
    assert '"-m",' in launcher_source
    assert '"revise.benchmark.cli",' in launcher_source
    assert '"benchmark_main.py"' not in launcher_source
