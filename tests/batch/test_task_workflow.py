"""Single and batch execution share scoped results on real project files."""
import json

import yaml

from test_runner import batch_config, fake_solver as _fake_solver, make_sample

fake_solver = _fake_solver


def test_resource_change_reruns_only_its_aspect(tmp_path, fake_solver, monkeypatch):
    from revise.batch import run_analysis_task, run_reconstruction_task
    from revise.batch.analysis import run_analysis_batch

    make_sample(tmp_path / 'data' / 'CRC' / 'one', cell_types=['T'])
    config = batch_config(tmp_path)
    module = tmp_path / 'resource_adapter.py'
    module.write_text('''
def run(context):
    value = context.resources['marker'].read_text()
    (context.output_dir / 'value.txt').write_text(value)
    return {
        'artifacts': {'value': {'path': 'value.txt', 'description': 'resource content'}},
        'calculation': {'input_view': 'native_carriers', 'parameters': {},
                        'comparison_basis': 'resource_fixture_only'},
    }
''')
    monkeypatch.syspath_prepend(str(tmp_path))
    for name in ('first', 'second'):
        (tmp_path / f'{name}.txt').write_text(name)
    first = run_reconstruction_task(config, 'CRC/one', cell_type='T')
    assert first['status'] == 'succeeded'
    document = yaml.safe_load(config.read_text())
    document['analysis'] = {
        name: {'entrypoint': 'resource_adapter:run', 'version': '1',
               'resources': {'marker': f'{name}.txt'}}
        for name in ('first', 'second')
    }
    config.write_text(yaml.safe_dump(document))
    for aspect in ('first', 'second'):
        assert run_analysis_task(config, 'CRC/one', aspect, cell_type='T')['status'] == 'succeeded'
    assert not (tmp_path / 'results' / 'batch_status.json').exists()
    assert run_analysis_batch(config)['summary']['reused'] == 2
    (tmp_path / 'first.txt').write_text('changed')
    report = run_analysis_batch(config)
    assert report['summary']['succeeded'] == 1
    assert report['summary']['reused'] == 1
    assert run_reconstruction_task(config, 'CRC/one', cell_type='T')['status'] == 'reused'
    assert fake_solver == ['T']
    output = tmp_path / 'results' / 'CRC' / 'one' / 'T'
    assert (output / 'analysis/first/value.txt').read_text() == 'changed'
    assert (output / 'analysis/second/value.txt').read_text() == 'second'
    assert run_analysis_batch(config)['summary']['reused'] == 2
    assert json.loads((output / '.revise/analysis/second.json').read_text())['status'] == 'succeeded'


def test_comment_only_config_edit_keeps_analysis_reusable(tmp_path, fake_solver, monkeypatch):
    from pathlib import Path
    from revise.batch import run_analysis_task, run_reconstruction_task

    make_sample(tmp_path / 'data' / 'one', cell_types=['T'])
    config = batch_config(tmp_path)
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[2] / 'configs/batch'))
    document = yaml.safe_load(config.read_text())
    document['analysis'] = {'inventory': {'entrypoint': 'analysis_example:run', 'version': '1'}}
    config.write_text(yaml.safe_dump(document))
    assert run_reconstruction_task(config, 'one', cell_type='T')['status'] == 'succeeded'
    assert run_analysis_task(config, 'one', 'inventory', cell_type='T')['status'] == 'succeeded'
    config.write_text(config.read_text() + '\n# Clarification without a semantic change\n')
    result = run_analysis_task(config, 'one', 'inventory', cell_type='T')
    assert result['status'] == 'reused'
    assert fake_solver == ['T']
