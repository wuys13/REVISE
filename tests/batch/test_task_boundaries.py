"""Boundary regressions found while integrating the task executors."""
import pytest
import yaml

from test_runner import batch_config, fake_solver as _fake_solver, make_sample

fake_solver = _fake_solver


def test_input_root_is_not_a_single_sample(tmp_path, fake_solver):
    from revise.batch import run_reconstruction_task

    make_sample(tmp_path / 'data', cell_types=['T'])
    config = batch_config(tmp_path)
    with pytest.raises(ValueError, match='sample_id'):
        run_reconstruction_task(config, '.', cell_type='T')
    assert fake_solver == []


def test_output_root_symlink_is_rejected(tmp_path, fake_solver):
    from revise.batch import run_reconstruction_task

    make_sample(tmp_path / 'data/one', cell_types=['T'])
    config = batch_config(tmp_path)
    outside = tmp_path / 'outside'
    outside.mkdir()
    (tmp_path / 'results').symlink_to(outside, target_is_directory=True)
    with pytest.raises(ValueError, match='symlink'):
        run_reconstruction_task(config, 'one', cell_type='T')
    assert list(outside.iterdir()) == []


def test_output_root_change_during_execution_is_not_success(tmp_path, fake_solver, monkeypatch):
    from revise.batch import runner

    make_sample(tmp_path / 'data/one', cell_types=['T'])
    config = batch_config(tmp_path)
    execute = runner._execute

    def changed_root(*args):
        result = execute(*args)
        doc = yaml.safe_load(config.read_text())
        doc['output_root'] = 'new-results'
        config.write_text(yaml.safe_dump(doc))
        return result

    monkeypatch.setattr(runner, '_execute', changed_root)
    result = runner.run_reconstruction_task(config, 'one', cell_type='T')
    assert result['status'] == 'failed'
    assert 'root' in result['error'].lower()


def test_reference_change_during_handoff_is_not_success(tmp_path, fake_solver, monkeypatch):
    import anndata as ad
    from revise.batch import runner

    make_sample(tmp_path / 'data/one', cell_types=['T'])
    config = batch_config(tmp_path)
    handoff = runner._handoff

    def changed_reference(sample, *args):
        result = handoff(sample, *args)
        data = ad.read_h5ad(sample.reference_path)
        data.uns['changed_during_handoff'] = True
        data.write_h5ad(sample.reference_path)
        return result

    monkeypatch.setattr(runner, '_handoff', changed_reference)
    result = runner.run_reconstruction_task(config, 'one', cell_type='T')
    assert result['status'] == 'failed'
    assert 'input' in result['error'].lower()
