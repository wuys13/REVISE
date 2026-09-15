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
    assert 'revise-prepare-sample' not in project['scripts']


@pytest.mark.parametrize('name', ['batch_reconstruct.py', 'batch_analyze.py'])
def test_source_help_from_unrelated_directory(tmp_path, name):
    script = Path(__file__).resolve().parents[2] / name
    result = subprocess.run([sys.executable, str(script), '--help'], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert '--config' in result.stdout


def test_legacy_batch_yaml_has_migration_error(tmp_path):
    from contextlib import redirect_stderr
    from io import StringIO
    from revise.batch.cli import main
    config = tmp_path / 'batch.yaml'
    config.write_text('schema_version: 1\ninput_root: data\noutput_root: results\n')
    stderr = StringIO()
    with redirect_stderr(stderr), pytest.raises(SystemExit) as error:
        main(['--config', str(config)])
    assert error.value.code == 2
    assert 'migrate' in stderr.getvalue()
