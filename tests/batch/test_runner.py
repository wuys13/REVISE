"""Batch lifecycle tests use real H5AD publication and a small solver stand-in."""
import json
from contextlib import redirect_stdout
from io import StringIO

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml


def make_sample(root, *, sample_id='sample', modality='iST', cell_types=None):
    root.mkdir(parents=True)
    obs = pd.DataFrame({'Level1': ['T', 'Macro', 'Fibroblast'], 'Level2': ['a', 'b', 'c']}, index=['u1', 'u2', 'u3'])
    data = ad.AnnData(np.ones((3, 3)), obs=obs)
    data.obsm['spatial'] = np.array([[0, 0], [1, 1], [2, 2]], dtype=float)
    data.write_h5ad(root / 'source.h5ad')
    data.write_h5ad(root / 'reference.h5ad')
    local = {'subtype_column': 'Level2', 'alpha': 0.2, 'resolutions': [0.5]}
    if cell_types is not None:
        local['cell_types'] = cell_types
    if modality == 'hST':
        local = {'strength': 0.2}
    if modality == 'sST':
        local = {'strength': 0.2, 'match_spot_sum': True}
        data.obs['n_cells'] = 2
        data.write_h5ad(root / 'source.h5ad')
    document = {
        'schema_version': 1, 'sample': {'id': sample_id, 'modality': modality},
        'inputs': {'st': {'path': 'source.h5ad', 'format': 'h5ad'},
                   'reference': {'path': 'reference.h5ad', 'format': 'h5ad'}},
        'preparation': {'spatial': {'matrix': 'X', 'spatial_key': 'spatial', 'coordinate_unit': 'um'},
                        'reference': {'matrix': 'X'}},
        'algorithm': {'ot_method': 'pot'},
        'preprocessing': {'spatial': {'min_transcript_counts': None, 'min_cell_counts': 1},
                          'reference': {'min_transcript_counts': None, 'min_cell_counts': 1}},
        'global_anchoring': {'broad_column': 'Level1'}, 'local_refinement': local,
        'execution': {'seed': 42}, 'output': {'ist_mapping': 'paired'} if modality == 'iST' else {},
    }
    (root / 'sample.yaml').write_text(yaml.safe_dump(document))
    return document


def result_root(root):
    parts = list(root.parts)
    parts[parts.index('data')] = 'results'
    return type(root)(*parts)


def batch_config(tmp_path):
    path = tmp_path / 'batch.yaml'
    path.write_text(yaml.safe_dump({'schema_version': 1, 'input_root': 'data', 'output_root': 'results'}))
    return path


@pytest.fixture
def fake_solver(monkeypatch):
    from revise.batch import runner
    from revise.application.config import compile_application_config, load_application_yaml
    from revise.application.publication import output_paths
    calls = []
    def execute(config_path, log_path):
        source, document = load_application_yaml(config_path)
        config = compile_application_config(document, source=source)
        calls.append(config.select_cell_type)
        if config.select_cell_type == 'Missing':
            log_path.write_text('missing selected cell type')
            return 1
        source_data = ad.read_h5ad(config.st_path)
        if config.select_cell_type:
            source_data = source_data[source_data.obs.Level1 == config.select_cell_type].copy()
        manifest = config.output_dir / 'engine' / 'provenance.json'
        manifest.parent.mkdir(parents=True, exist_ok=True)
        manifest.write_text(json.dumps({'status': 'succeeded'}))
        for role, path in output_paths(config).items():
            data = source_data.copy()
            data.obs['SVC_cluster'] = '0'
            if role == 'expression':
                data.obs_names = ['donor-' + x for x in data.obs_names]
                del data.obsm['spatial']
            if config.mode == 'sr':
                data.obs_names = ['generated-' + x for x in data.obs_names]
            data.uns['revise_reconstruction'] = {'run_manifest': str(manifest), 'output_role': role}
            path.parent.mkdir(parents=True, exist_ok=True)
            data.write_h5ad(path)
        log_path.write_text('Finished\n')
        return 0
    monkeypatch.setattr(runner, '_execute', execute)
    monkeypatch.setattr(runner, '_code_identity', lambda: {'sha256': 'fixed-test-engine'})
    return calls


def test_default_types_nested_discovery_and_safe_reuse(tmp_path, fake_solver, monkeypatch):
    from revise.batch.runner import run_batch
    root = tmp_path / 'data' / 'CRC' / 'one'
    make_sample(root)
    config = batch_config(tmp_path)
    monkeypatch.chdir(tmp_path.parent)
    result = run_batch(config)
    assert result['summary'] == {'succeeded': 3, 'failed': 0, 'reused': 0}
    assert fake_solver == ['T', 'Macro', 'Fibroblast']
    handoff = json.loads((result_root(root) / 'T' / 'reconstruction.json').read_text())
    assert handoff['status'] == 'succeeded'
    assert handoff['ist_mapping'] == 'paired'
    assert handoff['pairing']['status'] == 'available'
    assert handoff['outputs']['expression']['observation_role'] == 'reference_expression'
    assert handoff['analysis']['status'] == 'not_run'
    assert (result_root(root) / 'T' / 'analysis' / 'partition').is_dir()
    assert not (result_root(root) / 'T' / 'analysis' / 'partition' / 'raw').exists()
    assert run_batch(config)['summary'] == {'succeeded': 0, 'failed': 0, 'reused': 3}
    assert len(fake_solver) == 3
    (result_root(root) / 'T' / 'spatial.h5ad').write_bytes(b'corrupt')
    assert run_batch(config)['summary'] == {'succeeded': 1, 'failed': 0, 'reused': 2}


def test_failure_isolated_and_stale_success_invalidated(tmp_path, fake_solver):
    from revise.batch import runner
    root = tmp_path / 'data' / 'one'
    doc = make_sample(root, cell_types=['T', 'Missing', 'Fibroblast'])
    config = batch_config(tmp_path)
    assert runner.run_batch(config)['summary'] == {'succeeded': 2, 'failed': 1, 'reused': 0}
    assert runner.run_batch(config)['summary'] == {'succeeded': 0, 'failed': 1, 'reused': 2}
    doc['execution']['seed'] = 43
    (root / 'sample.yaml').write_text(yaml.safe_dump(doc))
    original = runner._execute
    runner._execute = lambda *args: 1
    try:
        assert runner.run_batch(config)['summary']['failed'] == 3
    finally:
        runner._execute = original
    assert json.loads((result_root(root) / 'T' / 'reconstruction.json').read_text())['status'] == 'failed'
    assert runner.run_batch(config)['summary'] == {'succeeded': 2, 'failed': 1, 'reused': 0}


def test_duplicate_ids_rejected_without_running(tmp_path, fake_solver):
    from revise.batch.runner import run_batch
    make_sample(tmp_path / 'data' / 'a')
    make_sample(tmp_path / 'data' / 'b')
    result = run_batch(batch_config(tmp_path))
    assert result['summary']['failed'] == 2
    assert fake_solver == []


@pytest.mark.parametrize('modality', ['hST', 'sST'])
def test_non_ist_output_and_pairing_roles(tmp_path, fake_solver, modality):
    from revise.batch.runner import run_batch
    root = tmp_path / 'data' / 'one'
    make_sample(root, modality=modality)
    assert run_batch(batch_config(tmp_path))['summary']['succeeded'] == 1
    assert (result_root(root) / 'SVC.h5ad').exists()
    handoff = json.loads((result_root(root) / 'reconstruction.json').read_text())
    assert handoff['pairing']['status'] == ('unavailable' if modality == 'sST' else 'available')


def test_interrupted_task_and_changed_code_rerun(tmp_path, fake_solver, monkeypatch):
    from revise.batch import runner
    root = tmp_path / 'data' / 'one'
    make_sample(root, cell_types=['T'])
    config = batch_config(tmp_path)
    runner.run_batch(config)
    state_path = result_root(root) / 'T' / '.revise' / 'task.json'
    state = json.loads(state_path.read_text())
    state['status'] = 'running'
    state_path.write_text(json.dumps(state))
    assert runner.run_batch(config)['summary']['succeeded'] == 1
    monkeypatch.setattr(runner, '_code_identity', lambda: {'sha256': 'changed'})
    assert runner.run_batch(config)['summary']['succeeded'] == 1


@pytest.mark.parametrize('labels', [['../escape'], ['inputs'], ['T', 't'], ['T', 'T']])
def test_unsafe_or_colliding_type_directories_rejected(tmp_path, fake_solver, labels):
    from revise.batch.runner import run_batch
    make_sample(tmp_path / 'data' / 'one', cell_types=labels)
    result = run_batch(batch_config(tmp_path))
    assert result['summary']['failed'] == 1
    assert not fake_solver


def test_preparation_failure_invalidates_old_handoff(tmp_path, fake_solver):
    from revise.batch.runner import run_batch
    root = tmp_path / 'data' / 'one'
    make_sample(root, cell_types=['T'])
    config = batch_config(tmp_path)
    run_batch(config)
    (root / 'source.h5ad').unlink()
    assert run_batch(config)['summary']['failed'] == 1
    assert json.loads((result_root(root) / 'T' / 'reconstruction.json').read_text())['status'] == 'failed'


def test_renamed_sample_and_removed_types_do_not_reuse_old_handoffs(tmp_path, fake_solver):
    from revise.batch.runner import run_batch
    root = tmp_path / 'data' / 'one'
    doc = make_sample(root, cell_types=['T', 'Macro'])
    config = batch_config(tmp_path)
    run_batch(config)
    doc['sample']['id'] = 'renamed'
    doc['local_refinement']['cell_types'] = ['T']
    (root / 'sample.yaml').write_text(yaml.safe_dump(doc))
    assert run_batch(config)['summary']['succeeded'] == 1
    assert json.loads((result_root(root) / 'T' / 'reconstruction.json').read_text())['sample_id'] == 'renamed'
    assert json.loads((result_root(root) / 'Macro' / 'reconstruction.json').read_text())['status'] == 'inactive'
    assert (result_root(root) / 'Macro' / 'spatial.h5ad').exists()


@pytest.mark.parametrize('obstruction', ['file', 'symlink'])
def test_bad_task_directory_does_not_stop_other_tasks(tmp_path, fake_solver, obstruction):
    from revise.batch.runner import run_batch
    root = tmp_path / 'data' / 'one'
    make_sample(root, cell_types=['T', 'Macro'])
    external = tmp_path / 'external'
    external.mkdir()
    result_root(root).mkdir(parents=True)
    if obstruction == 'file':
        (result_root(root) / 'T').write_text('keep me')
    else:
        (result_root(root) / 'T').symlink_to(external, target_is_directory=True)
    assert run_batch(batch_config(tmp_path))['summary'] == {'succeeded': 1, 'failed': 1, 'reused': 0}
    assert list(external.iterdir()) == []


def test_publication_cannot_overwrite_original_input(tmp_path, fake_solver):
    from revise.batch.runner import run_batch
    root = tmp_path / 'data' / 'one'
    doc = make_sample(root, modality='hST')
    (root / 'source.h5ad').rename(root / 'SVC.h5ad')
    doc['inputs']['st']['path'] = 'SVC.h5ad'
    (root / 'sample.yaml').write_text(yaml.safe_dump(doc))
    original = (root / 'SVC.h5ad').read_bytes()
    assert run_batch(batch_config(tmp_path))['summary']['succeeded'] == 1
    assert (root / 'SVC.h5ad').read_bytes() == original
    assert (result_root(root) / 'SVC.h5ad').exists()


def test_cli_nonzero_for_partial_failure(tmp_path, fake_solver):
    from revise.batch.cli import main
    make_sample(tmp_path / 'data' / 'one', cell_types=['Missing'])
    output = StringIO()
    with redirect_stdout(output), pytest.raises(SystemExit) as exc:
        main(['--config', str(batch_config(tmp_path))])
    assert exc.value.code == 1
    assert json.loads(output.getvalue())['failed'] == 1


def test_sample_status_symlink_does_not_write_outside_package(tmp_path, fake_solver):
    from revise.batch.runner import run_batch
    root = tmp_path / 'data' / 'one'
    make_sample(root, cell_types=['T'])
    external = tmp_path / 'external'
    external.mkdir()
    result_root(root).mkdir(parents=True)
    (result_root(root) / '.revise').symlink_to(external, target_is_directory=True)
    result = run_batch(batch_config(tmp_path))
    assert result['summary']['failed'] == 1
    assert not list(external.iterdir())


def test_changed_input_during_reconstruction_is_not_published_as_success(tmp_path, fake_solver, monkeypatch):
    from revise.batch import runner
    root = tmp_path / 'data' / 'one'
    make_sample(root, cell_types=['T'])
    execute = runner._execute
    def changing_source(config_path, log_path):
        result = execute(config_path, log_path)
        source = ad.read_h5ad(root / 'reference.h5ad')
        source.X *= 2
        source.write_h5ad(root / 'reference.h5ad')
        return result
    monkeypatch.setattr(runner, '_execute', changing_source)
    assert runner.run_batch(batch_config(tmp_path))['summary']['failed'] == 1
    assert json.loads((result_root(root) / 'T' / 'reconstruction.json').read_text())['status'] == 'failed'


def test_task_config_failure_invalidates_previous_success(tmp_path, fake_solver):
    from revise.batch.runner import run_batch
    root = tmp_path / 'data' / 'one'
    document = make_sample(root, cell_types=['T'])
    config = batch_config(tmp_path)
    run_batch(config)
    document['output']['name'] = 'forbidden'
    (root / 'sample.yaml').write_text(yaml.safe_dump(document))
    assert run_batch(config)['summary']['failed'] == 1
    assert json.loads((result_root(root) / 'T' / 'reconstruction.json').read_text())['status'] == 'failed'


def test_symlink_sample_yaml_cannot_relocate_package(tmp_path, fake_solver):
    from revise.batch.runner import run_batch
    external = tmp_path / 'external'
    make_sample(external, cell_types=['T'])
    root = tmp_path / 'data' / 'one'
    root.mkdir(parents=True)
    (root / 'sample.yaml').symlink_to(external / 'sample.yaml')
    assert run_batch(batch_config(tmp_path))['summary']['failed'] == 1
    assert not (external / 'inputs').exists()
    assert not fake_solver


def test_reconstruction_writes_only_prepared_inputs_to_input_tree(tmp_path, fake_solver):
    from revise.batch.runner import run_batch
    root = tmp_path / 'data' / 'CRC' / 'one'
    make_sample(root)
    run_batch(batch_config(tmp_path))
    assert sorted(p.name for p in root.iterdir()) == ['inputs', 'reference.h5ad', 'sample.yaml', 'source.h5ad']
    assert not (tmp_path / 'data' / 'batch_status.json').exists()
    assert (tmp_path / 'results' / 'batch_status.json').is_file()


@pytest.mark.parametrize('output', ['data', 'data/results', '.'])
def test_input_output_roots_must_not_overlap(tmp_path, fake_solver, output):
    from revise.batch.runner import run_batch
    make_sample(tmp_path / 'data' / 'one')
    config = batch_config(tmp_path)
    config.write_text(yaml.safe_dump({'schema_version': 1, 'input_root': 'data', 'output_root': output}))
    with pytest.raises(ValueError, match='overlap'):
        run_batch(config)
    assert not fake_solver


def test_blocked_sample_destination_does_not_stop_other_samples(tmp_path, fake_solver):
    from revise.batch.runner import run_batch
    make_sample(tmp_path / 'data' / 'one', sample_id='one', cell_types=['T'])
    make_sample(tmp_path / 'data' / 'two', sample_id='two', cell_types=['T'])
    (tmp_path / 'results').mkdir()
    (tmp_path / 'results' / 'one').write_text('keep')
    report = run_batch(batch_config(tmp_path))
    assert report['summary'] == {'succeeded': 1, 'failed': 1, 'reused': 0}
