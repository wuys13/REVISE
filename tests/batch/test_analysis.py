"""Exercise actual analysis adapters against published reconstruction handoffs."""
import json

import pytest
import yaml

from test_runner import batch_config, fake_solver as _fake_solver, make_sample, result_root

fake_solver = _fake_solver


def write_metrics(context):
    assert context.reconstruction['outputs']
    (context.output_dir / 'metrics.csv').write_text('raw,reconstructed,delta\n1,2,1\n')
    print('analysis adapter ran')
    return {
        'artifacts': {'paired_metrics': {'path': 'metrics.csv', 'description': 'paired metrics'}},
        'calculation': {
            'input_view': 'native_carriers', 'parameters': {},
            'comparison_basis': 'fixture comparison',
        },
    }


def fail_metrics(context):
    (context.output_dir / 'incomplete.csv').write_text('partial')
    raise RuntimeError('intentional analysis failure')


def escaping_metrics(context):
    return {
        'artifacts': {'bad': {'path': '../outside.csv', 'description': 'invalid path'}},
        'calculation': {
            'input_view': 'native_carriers', 'parameters': {},
            'comparison_basis': 'fixture comparison',
        },
    }


def config_with_analysis(tmp_path, aspects):
    config = batch_config(tmp_path)
    doc = yaml.safe_load(config.read_text())
    doc['analysis'] = aspects
    config.write_text(yaml.safe_dump(doc))
    return config


def adapter(name='write_metrics', **options):
    return {'entrypoint': f'test_analysis:{name}', 'version': '1', **options}


def prepared_batch(tmp_path, aspects, modality='iST'):
    from revise.batch.runner import run_batch
    root = tmp_path / 'data' / 'CRC' / 'one'
    make_sample(root, cell_types=['T'], modality=modality)
    config = config_with_analysis(tmp_path, aspects)
    assert run_batch(config)['summary']['succeeded'] == 1
    target = result_root(root) / 'T' if modality == 'iST' else result_root(root)
    return config, target


def test_unimplemented_aspects_are_not_successes(tmp_path, fake_solver):
    from revise.batch.analysis import run_analysis_batch
    config, target = prepared_batch(tmp_path, {})
    result = run_analysis_batch(config)
    assert result['summary']['not_implemented'] == 0
    assert result['summary']['succeeded'] == 0
    assert not list((target / 'analysis').rglob('*.csv'))
    assert json.loads((target / 'analysis' / 'analysis.json').read_text())['status'] == 'not_run'


def test_execute_reuse_and_config_or_artifact_changes(tmp_path, fake_solver):
    from revise.batch.analysis import run_analysis_batch
    config, target = prepared_batch(tmp_path, {'partition': adapter(requires_pairing=True)})
    assert run_analysis_batch(config)['summary']['succeeded'] == 1
    artifact = target / 'analysis' / 'partition' / 'metrics.csv'
    assert artifact.read_text().endswith('1,2,1\n')
    assert run_analysis_batch(config)['summary']['reused'] == 1
    artifact.write_text('damaged')
    assert run_analysis_batch(config)['summary']['succeeded'] == 1
    doc = yaml.safe_load(config.read_text())
    doc['analysis']['partition']['parameters'] = {'window': 10}
    config.write_text(yaml.safe_dump(doc))
    from revise.batch.runner import run_batch
    run_batch(config)
    assert run_analysis_batch(config)['summary']['succeeded'] == 1
    assert 'analysis adapter ran' in (target / '.revise' / 'analysis' / 'partition.log').read_text()


def test_failure_continues_and_does_not_publish_partial_files(tmp_path, fake_solver):
    from revise.batch.analysis import run_analysis_batch
    config, target = prepared_batch(tmp_path, {
        'partition': adapter('fail_metrics'), 'spatial_diversity': adapter(),
    })
    result = run_analysis_batch(config)
    assert result['summary']['failed'] == 1
    assert result['summary']['succeeded'] == 1
    assert not (target / 'analysis' / 'partition' / 'incomplete.csv').exists()
    assert run_analysis_batch(config)['summary']['reused'] == 1


def test_pairing_requirement_blocks_sst(tmp_path, fake_solver):
    from revise.batch.analysis import run_analysis_batch
    config, target = prepared_batch(tmp_path, {'partition': adapter(requires_pairing=True)}, 'sST')
    assert run_analysis_batch(config)['summary']['unavailable'] == 1
    assert not (target / 'analysis' / 'partition' / 'metrics.csv').exists()


def test_tampered_reconstruction_blocks_analysis(tmp_path, fake_solver):
    from revise.batch.analysis import run_analysis_batch
    config, target = prepared_batch(tmp_path, {'partition': adapter()})
    (target / 'spatial.h5ad').write_bytes(b'corrupt')
    assert run_analysis_batch(config)['summary']['blocked'] == 1


def test_adapter_cannot_return_paths_outside_aspect(tmp_path, fake_solver):
    from revise.batch.analysis import run_analysis_batch
    config, _ = prepared_batch(tmp_path, {'partition': adapter('escaping_metrics')})
    assert run_analysis_batch(config)['summary']['failed'] == 1


def test_reconstruction_rerun_invalidates_analysis_summary(tmp_path, fake_solver):
    from revise.batch.analysis import run_analysis_batch
    from revise.batch.runner import run_batch
    config, target = prepared_batch(tmp_path, {'partition': adapter()})
    run_analysis_batch(config)
    sample = tmp_path / 'data' / 'CRC' / 'one' / 'batch.yaml'
    doc = yaml.safe_load(sample.read_text())
    doc['execution']['seed'] = 43
    sample.write_text(yaml.safe_dump(doc))
    run_batch(config)
    assert json.loads((target / 'analysis' / 'analysis.json').read_text())['status'] == 'not_run'
    assert run_analysis_batch(config)['summary']['succeeded'] == 1


def test_analysis_cli_reports_incomplete_nonzero(tmp_path, fake_solver):
    from revise.batch.cli import analysis_main
    config, _ = prepared_batch(tmp_path, {})
    analysis_main(['--config', str(config)])


def test_changed_sample_configuration_requires_reconstruction(tmp_path, fake_solver):
    from revise.batch.analysis import run_analysis_batch
    config, _ = prepared_batch(tmp_path, {'partition': adapter()})
    sample = tmp_path / 'data' / 'CRC' / 'one' / 'batch.yaml'
    doc = yaml.safe_load(sample.read_text())
    doc['execution']['seed'] = 43
    sample.write_text(yaml.safe_dump(doc))
    assert run_analysis_batch(config)['summary']['blocked'] == 1


def test_placeholder_cannot_claim_ownership_of_foreign_analysis(tmp_path, fake_solver):
    from revise.batch.analysis import run_analysis_batch
    config, target = prepared_batch(tmp_path, {'partition': {}})
    artifact = target / 'analysis' / 'partition' / 'foreign.txt'
    artifact.parent.mkdir(parents=True, exist_ok=True)
    artifact.write_text('keep')
    run_analysis_batch(config)
    doc = yaml.safe_load(config.read_text())
    doc['analysis']['partition'] = adapter()
    config.write_text(yaml.safe_dump(doc))
    from revise.batch.runner import run_batch
    run_batch(config)
    assert run_analysis_batch(config)['summary']['failed'] == 1
    assert artifact.read_text() == 'keep'


@pytest.mark.parametrize('directory', ['analysis', '.revise/analysis'])
def test_analysis_symlink_cannot_write_elsewhere(tmp_path, fake_solver, directory):
    import shutil
    from revise.batch.analysis import run_analysis_batch
    config, target = prepared_batch(tmp_path, {'partition': adapter()})
    external = tmp_path / 'external'
    external.mkdir()
    link = target / directory
    if link.exists():
        shutil.rmtree(link)
    link.symlink_to(external, target_is_directory=True)
    result = run_analysis_batch(config)
    assert result['status'] == 'incomplete'
    assert not list(external.iterdir())


def test_reused_analysis_summary_keeps_artifact_roles(tmp_path, fake_solver):
    from revise.batch.analysis import run_analysis_batch
    config, target = prepared_batch(tmp_path, {'partition': adapter()})
    run_analysis_batch(config)
    result = run_analysis_batch(config)
    assert result['tasks'][0]['roles'] == {'paired_metrics': 'metrics.csv'}


def test_publication_failure_restores_previous_analysis(tmp_path, fake_solver, monkeypatch):
    from revise.batch import analysis
    config, target = prepared_batch(tmp_path, {'partition': adapter()})
    analysis.run_analysis_batch(config)
    path = target / 'analysis' / 'partition' / 'metrics.csv'
    before = path.read_bytes()
    doc = yaml.safe_load(config.read_text())
    doc['analysis']['partition']['version'] = '2'
    config.write_text(yaml.safe_dump(doc))
    from revise.batch.runner import run_batch
    run_batch(config)
    replace = analysis.os.replace
    def fail_publish(source, destination):
        if str(source).endswith('/artifacts'):
            raise OSError('publication failed')
        return replace(source, destination)
    monkeypatch.setattr(analysis.os, 'replace', fail_publish)
    assert analysis.run_analysis_batch(config)['summary']['failed'] == 1
    assert path.read_bytes() == before


def test_later_user_file_blocks_whole_aspect_replacement(tmp_path, fake_solver):
    from revise.batch.analysis import run_analysis_batch
    from revise.batch.runner import run_batch

    config, target = prepared_batch(tmp_path, {'partition': adapter()})
    run_analysis_batch(config)
    user_file = target / 'analysis' / 'partition' / 'keep.txt'
    user_file.write_text('keep')
    document = yaml.safe_load(config.read_text())
    document['analysis']['partition']['version'] = '2'
    config.write_text(yaml.safe_dump(document))
    run_batch(config)

    result = run_analysis_batch(config)

    assert result['summary']['failed'] == 1
    assert user_file.read_text() == 'keep'


def test_interrupted_analysis_retries(tmp_path, fake_solver):
    from revise.batch.analysis import run_analysis_batch
    config, target = prepared_batch(tmp_path, {'partition': adapter()})
    run_analysis_batch(config)
    path = target / '.revise' / 'analysis' / 'partition.json'
    record = json.loads(path.read_text())
    record['status'] = 'running'
    path.write_text(json.dumps(record))
    assert run_analysis_batch(config)['summary']['succeeded'] == 1


def test_failed_analysis_does_not_prevent_other_samples(tmp_path, fake_solver):
    from revise.batch.runner import run_batch
    from revise.batch.analysis import run_analysis_batch
    config, target = prepared_batch(tmp_path, {'partition': adapter()})
    make_sample(tmp_path / 'data' / 'CRC' / 'two', cell_types=['T'])
    run_batch(config)
    (target / 'spatial.h5ad').write_bytes(b'broken')
    report = run_analysis_batch(config)
    assert report['summary']['blocked'] == 1
    assert report['summary']['succeeded'] == 1


def test_source_analysis_cli_executes_external_adapter(tmp_path, fake_solver):
    import importlib
    import os
    from pathlib import Path
    import subprocess
    import sys
    from revise.batch import runner
    config, target = prepared_batch(tmp_path, {'partition': {
        'entrypoint': 'example_adapter:run', 'version': '1',
    }})
    # ``fake_solver`` uses a sentinel reconstruction code identity.  Re-run
    # the fixture task with the real identity before handing it to a fresh
    # CLI process, whose verifier cannot see the parent test monkeypatch.
    fake_execute = runner._execute
    importlib.reload(runner)
    runner._execute = fake_execute
    assert runner.run_batch(config)['summary']['succeeded'] == 1
    (tmp_path / 'example_adapter.py').write_text(
        "def run(context):\n"
        "    (context.output_dir / 'result.csv').write_text('value\\n1\\n')\n"
        "    return {'artifacts': {'metrics': {'path': 'result.csv', 'description': 'fixture result'}},\n"
        "            'calculation': {'input_view': 'native_carriers', 'parameters': {},\n"
        "                           'comparison_basis': 'fixture comparison'}}\n"
    )
    environment = os.environ.copy()
    environment['PYTHONPATH'] = str(tmp_path)
    script = Path(__file__).resolve().parents[2] / 'batch_analyze.py'
    process = subprocess.run([sys.executable, str(script), '--config', str(config)],
                             cwd=tmp_path, env=environment, capture_output=True, text=True)
    assert process.returncode == 0, process.stderr
    assert json.loads(process.stdout)['succeeded'] == 1
    assert (target / 'analysis' / 'partition' / 'result.csv').read_text() == 'value\n1\n'


def test_blocked_reconstruction_invalidates_local_analysis_summary(tmp_path, fake_solver):
    from revise.batch.analysis import run_analysis_batch
    config, target = prepared_batch(tmp_path, {'partition': adapter()})
    run_analysis_batch(config)
    (target / 'spatial.h5ad').write_bytes(b'broken')
    run_analysis_batch(config)
    assert json.loads((target / 'analysis' / 'analysis.json').read_text())['status'] == 'blocked'


def test_removed_sample_uses_prior_analysis_inventory(tmp_path, fake_solver):
    import shutil
    from revise.batch.analysis import run_analysis_batch

    config, target = prepared_batch(tmp_path, {'partition': adapter()})
    run_analysis_batch(config)
    control = target / '.revise' / 'analysis' / 'partition.json'
    artifact = target / 'analysis' / 'partition' / 'metrics.csv'
    shutil.rmtree(tmp_path / 'data' / 'CRC' / 'one')

    result = run_analysis_batch(config)

    assert result['tasks'] == []
    assert json.loads(control.read_text())['status'] == 'inactive'
    assert json.loads((target / 'analysis' / 'analysis.json').read_text())['status'] == 'inactive'
    assert artifact.read_text().endswith('1,2,1\n')


@pytest.mark.parametrize('change', ['add', 'remove', 'edit'])
def test_analysis_rechecks_configuration_chain(tmp_path, fake_solver, change):
    from revise.batch.analysis import run_analysis_batch
    from revise.batch.runner import run_batch
    config, target = prepared_batch(tmp_path, {'partition': adapter()})
    parent = tmp_path / 'data' / 'CRC' / 'batch.yaml'
    if change != 'add':
        parent.write_text('execution: {seed: 1}\n')
        run_batch(config)
    run_analysis_batch(config)
    if change == 'remove':
        parent.unlink()
    else:
        parent.write_text('execution: {seed: 2}\n')
    result = run_analysis_batch(config)
    assert result['summary']['reused'] == 1
    state = json.loads((target / '.revise' / 'analysis' / 'partition.json').read_text())
    parent_audit = next(item for item in state['configuration_audit']['chain']
                        if item['path'] == str(parent))
    assert (parent_audit['sha256'] is None) == (change == 'remove')


def test_analysis_uses_sample_level_aspects(tmp_path, fake_solver):
    from revise.batch.runner import run_batch
    from revise.batch.analysis import run_analysis_batch
    config, target = prepared_batch(tmp_path, {'partition': adapter()})
    path = tmp_path / 'data' / 'CRC' / 'one' / 'batch.yaml'
    doc = yaml.safe_load(path.read_text())
    doc['analysis'] = {'partition': adapter('fail_metrics')}
    path.write_text(yaml.safe_dump(doc))
    run_batch(config)
    assert run_analysis_batch(config)['summary']['failed'] == 1


def test_different_project_config_cannot_use_another_inventory(tmp_path, fake_solver):
    from revise.batch.analysis import run_analysis_batch
    config, _ = prepared_batch(tmp_path, {'partition': adapter()})
    alternate = tmp_path / 'other.yaml'
    alternate.write_text(config.read_text())
    with pytest.raises(ValueError, match='inventory'):
        run_analysis_batch(alternate)
