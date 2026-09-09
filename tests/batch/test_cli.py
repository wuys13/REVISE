from pathlib import Path
import subprocess
import sys

import tomli
import pytest


def test_installed_entrypoints_declared():
    root = Path(__file__).resolve().parents[2]
    project = tomli.loads((root / 'pyproject.toml').read_text())['project']
    assert project['scripts']['revise-batch-reconstruct'] == 'revise.batch.cli:main'
    assert project['scripts']['revise-batch-analyze'] == 'revise.batch.cli:analysis_main'
    assert project['scripts']['revise-prepare-sample'] == 'revise.batch.cli:prepare_main'


@pytest.mark.parametrize('name', ['batch_reconstruct.py', 'batch_analyze.py'])
def test_source_help_from_unrelated_directory(tmp_path, name):
    script = Path(__file__).resolve().parents[2] / name
    result = subprocess.run([sys.executable, str(script), '--help'], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert '--config' in result.stdout
