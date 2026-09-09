from pathlib import Path

import pytest
import yaml


def write(path, document):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(document))
    return path


def project(tmp_path):
    (tmp_path / 'ST/CRC/S01').mkdir(parents=True)
    return write(tmp_path / 'batch.yaml', {'schema_version': 2, 'input_root': 'ST',
                 'output_root': 'output', 'modality': 'iST', 'coordinates': {'unit': 'pixel'},
                 'inputs': {'reference': {'path': 'SC/shared.h5ad', 'filter_column': 'donor', 'filter_value': 'A'}},
                 'local_refinement': {'cell_types': ['T', 'Macro'], 'subtype_column': 'Level2'}})


def test_inheritance_paths_lists_and_reference_replacement(tmp_path, monkeypatch):
    from revise.batch.config import resolve_sample
    path = project(tmp_path)
    write(tmp_path / 'ST/CRC/batch.yaml', {'inputs': {'reference': {'path': '../../SC/crc.h5ad'}},
                                         'local_refinement': {'cell_types': ['T']}})
    write(tmp_path / 'ST/CRC/S01/batch.yaml', {'execution': {'seed': 7}})
    monkeypatch.chdir(tmp_path.parent)
    result = resolve_sample(path, tmp_path / 'ST/CRC/S01')
    assert result.sample_id == 'CRC/S01'
    assert result.document['inputs']['reference'] == {'path': str(tmp_path / 'SC/crc.h5ad')}
    assert result.document['local_refinement'] == {'cell_types': ['T'], 'subtype_column': 'Level2'}
    assert result.document['execution']['seed'] == 7
    assert result.document['enabled'] is True
    assert len(result.config_chain) == 4
    assert result.config_chain[1] == {'path': str(tmp_path / 'ST/batch.yaml'), 'sha256': None}


def test_chain_detects_added_removed_and_changed_configuration(tmp_path):
    from revise.batch.config import resolve_sample
    path = project(tmp_path)
    sample = tmp_path / 'ST/CRC/S01'
    first = resolve_sample(path, sample)
    override = write(sample / 'batch.yaml', {'enabled': False})
    second = resolve_sample(path, sample)
    assert second.document['enabled'] is False
    assert second.config_chain != first.config_chain
    override.unlink()
    assert resolve_sample(path, sample).config_chain == first.config_chain


@pytest.mark.parametrize('override', [{'schema_version': 2}, {'input_root': 'x'}, {'output_root': 'x'},
                                     {'inputs': {'st': {'path': 'x'}}}, {'enabled': 'false'},
                                     {'preparation': {}}, {'unknown': 1}, {'coordinates': None}])
def test_invalid_child_fields(tmp_path, override):
    from revise.batch.config import resolve_sample
    path = project(tmp_path)
    write(tmp_path / 'ST/CRC/batch.yaml', override)
    with pytest.raises(ValueError):
        resolve_sample(path, tmp_path / 'ST/CRC/S01')


@pytest.mark.parametrize('version', [1, None, True, '2'])
def test_old_protocol_has_migration_message(tmp_path, version):
    from revise.batch.config import load_batch_config
    path = write(tmp_path / 'batch.yaml', {'schema_version': version, 'input_root': 'ST', 'output_root': 'output'})
    with pytest.raises(ValueError, match='schema_version.*2.*migrat'):
        load_batch_config(path)


def test_discovery_does_not_follow_directory_links(tmp_path):
    from revise.batch.config import discover_samples
    root = tmp_path / 'ST'
    for folder in ['CRC/S01', 'BRCA/S01']:
        target = root / folder
        target.mkdir(parents=True)
        (target / 'spatial.h5ad').touch()
    (root / 'alias').symlink_to(root / 'CRC', target_is_directory=True)
    assert discover_samples(root) == [root / 'BRCA/S01', root / 'CRC/S01']


@pytest.mark.parametrize('output', ['ST', 'ST/output', '.'])
def test_roots_must_be_disjoint(tmp_path, output):
    from revise.batch.config import load_batch_config
    path = project(tmp_path)
    doc = yaml.safe_load(path.read_text())
    doc['output_root'] = output
    write(path, doc)
    with pytest.raises(ValueError, match='overlap'):
        load_batch_config(path)


def test_sample_must_be_inside_input_root(tmp_path):
    from revise.batch.config import resolve_sample
    path = project(tmp_path)
    with pytest.raises(ValueError):
        resolve_sample(path, tmp_path)


def test_prior_path_is_resolved_where_declared(tmp_path):
    from revise.batch.config import resolve_sample
    path = project(tmp_path)
    write(tmp_path / 'ST/CRC/batch.yaml', {'inputs': {'pm_on_cell': {'path': '../../prior.csv'}}})
    result = resolve_sample(path, tmp_path / 'ST/CRC/S01')
    assert result.document['inputs']['pm_on_cell']['path'] == str(tmp_path / 'prior.csv')


def test_null_input_spec_is_rejected(tmp_path):
    from revise.batch.config import resolve_sample
    path = project(tmp_path)
    write(tmp_path / 'ST/CRC/batch.yaml', {'inputs': {'pm_on_cell': None}})
    with pytest.raises(ValueError, match='pm_on_cell'):
        resolve_sample(path, tmp_path / 'ST/CRC/S01')


def test_symlink_configuration_is_rejected(tmp_path):
    from revise.batch.config import resolve_sample
    path = project(tmp_path)
    external = write(tmp_path / 'other.yaml', {'enabled': False})
    (tmp_path / 'ST/CRC/batch.yaml').symlink_to(external)
    with pytest.raises(ValueError, match='symlink'):
        resolve_sample(path, tmp_path / 'ST/CRC/S01')


def test_project_path_with_parent_segments_is_not_loaded_twice(tmp_path):
    from revise.batch.config import resolve_sample
    root = tmp_path / 'ST'
    sample = root / 'S01'
    sample.mkdir(parents=True)
    path = write(root / 'batch.yaml', {'schema_version': 2, 'input_root': '.', 'output_root': '../output'})
    result = resolve_sample(root / 'S01/../batch.yaml', sample)
    assert result.project_config == path
    assert len(result.config_chain) == 2


def test_configuration_hash_uses_the_parsed_snapshot(tmp_path, monkeypatch):
    from revise.batch.config import resolve_sample
    path = project(tmp_path)
    original = Path.read_bytes
    reads = []

    def track(candidate):
        if candidate.name == 'batch.yaml':
            reads.append(candidate)
        return original(candidate)

    monkeypatch.setattr(Path, 'read_bytes', track)
    monkeypatch.setattr(Path, 'read_text', lambda *args, **kwargs: pytest.fail('separate text read'))
    resolve_sample(path, tmp_path / 'ST/CRC/S01')
    assert reads.count(path) == 1


def test_malformed_yaml_is_a_configuration_error(tmp_path):
    from revise.batch.config import load_batch_config
    path = tmp_path / 'batch.yaml'
    path.write_text('schema_version: [')
    with pytest.raises(ValueError, match='batch.yaml'):
        load_batch_config(path)
