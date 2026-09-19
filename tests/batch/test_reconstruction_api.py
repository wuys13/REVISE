"""Public reconstruction task API tests."""
import json

import pytest
import yaml

from test_runner import batch_config, fake_solver as _fake_solver, make_sample, result_root


fake_solver = _fake_solver


def test_single_task_runs_only_the_selected_sample_and_type(tmp_path, fake_solver):
    from revise.batch.runner import run_reconstruction_task

    make_sample(tmp_path / 'data' / 'CRC' / 'one', cell_types=['T', 'Macro'])
    make_sample(tmp_path / 'data' / 'CRC' / 'two', cell_types=['T'])
    config = batch_config(tmp_path)

    result = run_reconstruction_task(config, 'CRC/one', cell_type='T')

    assert result['status'] == 'succeeded'
    assert result['sample_id'] == 'CRC/one'
    assert result['cell_type'] == 'T'
    assert result['fingerprint']
    assert result['outputs']['spatial']['path'].endswith('/T/spatial.h5ad')
    assert result['artifacts']
    assert result['log_path'].endswith('/T/.revise/reconstruction.log')
    assert fake_solver == ['T']
    assert (result_root(tmp_path / 'data' / 'CRC' / 'one') / 'T' / 'reconstruction.json').exists()
    assert not (result_root(tmp_path / 'data' / 'CRC' / 'one') / 'Macro').exists()
    assert not (result_root(tmp_path / 'data' / 'CRC' / 'two')).exists()
    assert not (tmp_path / 'results' / 'batch_status.json').exists()


def test_explicit_single_type_remains_available_without_legacy_cell_types_list(
        tmp_path, fake_solver):
    from revise.batch.runner import run_reconstruction_task

    root = tmp_path / 'data' / 'one'
    document = make_sample(root, cell_types=None)
    document['output'] = {}
    (root / 'batch.yaml').write_text(yaml.safe_dump(document))

    result = run_reconstruction_task(batch_config(tmp_path), 'one', cell_type='T')

    assert result['status'] == 'succeeded'
    assert result['cell_type'] == 'T'
    assert result['outputs']['spatial']['path'].endswith('/T/spatial.h5ad')
    assert fake_solver == ['T']


def test_single_task_reuses_without_batch_inventory(tmp_path, fake_solver):
    from revise.batch.runner import run_reconstruction_task

    root = tmp_path / 'data' / 'one'
    make_sample(root, cell_types=['T'])
    config = batch_config(tmp_path)
    first = run_reconstruction_task(config, 'one', cell_type='T')
    second = run_reconstruction_task(config, 'one', cell_type='T')

    assert first['status'] == 'succeeded'
    assert second['status'] == 'reused'
    assert fake_solver == ['T']


@pytest.mark.parametrize('sample_id, cell_type', [
    ('missing', 'T'), ('one', 'Missing'), ('one', None),
])
def test_single_task_rejects_out_of_scope_selection(tmp_path, fake_solver, sample_id, cell_type):
    from revise.batch.runner import run_reconstruction_task

    make_sample(tmp_path / 'data' / 'one', cell_types=['T'])
    config = batch_config(tmp_path)

    with pytest.raises(ValueError):
        run_reconstruction_task(config, sample_id, cell_type=cell_type)
    assert not (tmp_path / 'results').exists()
    assert not fake_solver


def test_single_task_reports_disabled_task_without_execution(tmp_path, fake_solver):
    from revise.batch.runner import run_reconstruction_task

    root = tmp_path / 'data' / 'one'
    document = make_sample(root, cell_types=['T'])
    config = batch_config(tmp_path)
    run_reconstruction_task(config, 'one', cell_type='T')
    document['enabled'] = False
    (root / 'batch.yaml').write_text(yaml.safe_dump(document))

    result = run_reconstruction_task(config, 'one', cell_type='T')

    assert result['status'] == 'inactive'
    assert fake_solver == ['T']
    assert json.loads((result_root(root) / 'T' / 'reconstruction.json').read_text())['status'] == 'inactive'


def test_verification_helper_checks_current_task_without_batch_inventory(tmp_path, fake_solver):
    from revise.batch.config import resolve_sample
    from revise.batch.runner import run_reconstruction_task, verify_reconstruction_task

    root = tmp_path / 'data' / 'one'
    make_sample(root, cell_types=['T'])
    config = batch_config(tmp_path)
    run_reconstruction_task(config, 'one', cell_type='T')
    resolved = resolve_sample(config, root)

    record = verify_reconstruction_task(resolved, 'T')

    assert record['status'] == 'succeeded'
    assert record['sample_id'] == 'one'
    assert record['cell_type'] == 'T'
    (result_root(root) / 'T' / 'spatial.h5ad').write_bytes(b'tampered')
    with pytest.raises(ValueError):
        verify_reconstruction_task(resolved, 'T')
