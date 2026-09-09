from pathlib import Path
import subprocess
import sys

import tomli


def test_installed_entrypoints_declared():
    root = Path(__file__).resolve().parents[2]
    project = tomli.loads((root / 'pyproject.toml').read_text())['project']
    assert project['scripts']['revise-batch-reconstruct'] == 'revise.batch.cli:main'
    assert project['scripts']['revise-prepare-sample'] == 'revise.batch.cli:prepare_main'


def test_source_help_from_unrelated_directory(tmp_path):
    script = Path(__file__).resolve().parents[2] / 'batch_reconstruct.py'
    result = subprocess.run([sys.executable, str(script), '--help'], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert '--config' in result.stdout
