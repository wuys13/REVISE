from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml
from scipy import sparse


def fixture_sample(tmp_path):
    root = tmp_path / 'cancer' / 'sample-a'
    root.mkdir(parents=True)
    raw = tmp_path / 'raw.h5ad'
    ref = tmp_path / 'shared.h5ad'
    a = ad.AnnData(sparse.csr_matrix([[1., 2.], [3., 4.]]),
                   obs=pd.DataFrame(index=['c1', 'c2']), var=pd.DataFrame(index=['g1', 'g2']))
    a.obsm['xy'] = np.array([[1., 2.], [3., 4.]])
    a.layers['counts'] = a.X * 2
    a.uns['original_provenance'] = 'retain me'
    a.write_h5ad(raw)
    a.obs['Level1'] = ['T', 'Mono_Macro']
    a.obs['Level2'] = ['T1', 'M1']
    a.write_h5ad(ref)
    doc = {'schema_version': 1, 'sample': {'id': 'sample-a', 'modality': 'iST'},
           'inputs': {'st': {'path': '../../raw.h5ad', 'format': 'h5ad'},
                      'reference': {'path': str(ref), 'format': 'h5ad'}},
           'global_anchoring': {'broad_column': 'Level1'},
           'local_refinement': {'subtype_column': 'Level2'},
           'preparation': {'spatial': {'matrix': 'layers/counts', 'spatial_key': 'xy', 'coordinate_unit': 'pixel'},
                           'reference': {'matrix': 'X'}}}
    path = root / 'sample.yaml'
    path.write_text(yaml.safe_dump(doc))
    return path, doc, raw, ref


def prepare(path):
    from revise.batch.sample import prepare_sample
    return prepare_sample(path)


def test_nested_package_preserves_sources_and_shared_reference(tmp_path, monkeypatch):
    path, doc, raw, ref = fixture_sample(tmp_path)
    raw_before = raw.read_bytes()
    monkeypatch.chdir(tmp_path.parent)
    result = prepare(path)
    assert result.sample_id == 'sample-a'
    assert result.modality == 'iST'
    assert result.root == path.parent
    assert result.document == doc
    assert result.reference_path == ref
    a = ad.read_h5ad(result.st_path)
    np.testing.assert_array_equal(a.X.toarray(), [[2, 4], [6, 8]])
    np.testing.assert_array_equal(a.obsm['spatial'], a.obsm['xy'])
    assert a.uns['original_provenance'] == 'retain me'
    assert raw.read_bytes() == raw_before
    assert result.metadata['coordinates']['physical_distance_available'] is False
    assert result.metadata['sources']['spatial']['sha256']


def test_mapping_preserves_original_labels_and_source(tmp_path):
    path, doc, raw, ref = fixture_sample(tmp_path)
    doc['preparation']['reference']['label_mapping'] = {'Level1': {'Mono_Macro': 'Macro'}}
    path.write_text(yaml.safe_dump(doc))
    result = prepare(path)
    a = ad.read_h5ad(result.reference_path)
    assert a.obs.Level1.tolist() == ['T', 'Macro']
    assert a.obs['Level1_original'].tolist() == ['T', 'Mono_Macro']
    assert ad.read_h5ad(ref).obs.Level1.tolist() == ['T', 'Mono_Macro']


@pytest.mark.parametrize('field', ['matrix', 'spatial_key', 'coordinate_unit'])
def test_requires_explicit_spatial_metadata(tmp_path, field):
    path, doc, _, _ = fixture_sample(tmp_path)
    del doc['preparation']['spatial'][field]
    path.write_text(yaml.safe_dump(doc))
    with pytest.raises(ValueError, match=field):
        prepare(path)


@pytest.mark.parametrize('problem', ['duplicate_obs', 'duplicate_var', 'nan_matrix', 'negative_matrix', 'nan_coordinates', 'missing_label'])
def test_invalid_inputs_fail_before_publication(tmp_path, problem):
    path, doc, raw, ref = fixture_sample(tmp_path)
    source = ref if problem == 'missing_label' else raw
    a = ad.read_h5ad(source)
    if problem == 'duplicate_obs':
        a.obs_names = ['x', 'x']
    elif problem == 'duplicate_var':
        a.var_names = ['x', 'x']
    elif problem == 'nan_matrix':
        a.layers['counts'].data[0] = np.nan
    elif problem == 'negative_matrix':
        a.layers['counts'].data[0] = -1
    elif problem == 'nan_coordinates':
        a.obsm['xy'][0, 0] = np.nan
    else:
        del a.obs['Level1']
    a.write_h5ad(source)
    with pytest.raises(ValueError):
        prepare(path)
    assert not (path.parent / 'inputs' / 'spatial.h5ad').exists()


def test_cache_checks_output_and_source_content(tmp_path):
    path, _, raw, _ = fixture_sample(tmp_path)
    first = prepare(path)
    timestamp = first.st_path.stat().st_mtime_ns
    assert prepare(path).st_path.stat().st_mtime_ns == timestamp
    first.st_path.write_bytes(b'broken cache')
    repaired = prepare(path)
    assert ad.read_h5ad(repaired.st_path).n_obs == 2
    source = ad.read_h5ad(raw)
    source.layers['counts'] *= 3
    source.write_h5ad(raw)
    changed = prepare(path)
    assert changed.metadata['sources'] != first.metadata['sources']
    np.testing.assert_array_equal(ad.read_h5ad(changed.st_path).X.toarray(), [[6, 12], [18, 24]])


def test_never_overwrites_source_at_canonical_destination(tmp_path):
    path, doc, raw, _ = fixture_sample(tmp_path)
    destination = path.parent / 'inputs' / 'spatial.h5ad'
    destination.parent.mkdir()
    destination.write_bytes(raw.read_bytes())
    doc['inputs']['st']['path'] = 'inputs/spatial.h5ad'
    path.write_text(yaml.safe_dump(doc))
    before = destination.read_bytes()
    with pytest.raises(ValueError, match='source'):
        prepare(path)
    assert destination.read_bytes() == before


def test_directory_identity_tracks_contents(tmp_path):
    from revise.batch.sample import file_identity
    store = tmp_path / 'store.zarr'
    store.mkdir()
    (store / 'a').write_bytes(b'a')
    first = file_identity(store)
    (store / 'b').write_bytes(b'b')
    assert file_identity(store)['sha256'] != first['sha256']


def test_failed_publication_restores_previous_package(tmp_path, monkeypatch):
    from revise.batch import sample
    path, doc, _, _ = fixture_sample(tmp_path)
    first = prepare(path)
    before = first.st_path.read_bytes()
    manifest = path.parent / 'inputs' / 'preparation.json'
    previous_manifest = manifest.read_bytes()
    doc['preparation']['spatial']['matrix'] = 'X'
    path.write_text(yaml.safe_dump(doc))
    replace = sample.os.replace

    def fail_manifest(source, destination):
        if Path(source).name == 'preparation.json' and '.prepare-' in str(source):
            raise OSError('injected publication failure')
        return replace(source, destination)

    monkeypatch.setattr(sample.os, 'replace', fail_manifest)
    with pytest.raises(OSError, match='injected'):
        prepare(path)
    assert first.st_path.read_bytes() == before
    assert manifest.read_bytes() == previous_manifest


def test_spatialdata_requires_explicit_table_and_element(tmp_path):
    path, doc, _, _ = fixture_sample(tmp_path)
    doc['inputs']['st']['format'] = 'spatialdata'
    path.write_text(yaml.safe_dump(doc))
    with pytest.raises(ValueError, match='table'):
        prepare(path)
    doc['inputs']['st']['spatialdata'] = {'table': 'counts'}
    path.write_text(yaml.safe_dump(doc))
    with pytest.raises(ValueError, match='element'):
        prepare(path)


def test_reference_filter_is_preserved_and_validated(tmp_path):
    path, doc, _, _ = fixture_sample(tmp_path)
    doc['inputs']['reference'].update(filter_column='Level1', filter_value='T')
    path.write_text(yaml.safe_dump(doc))
    result = prepare(path)
    assert ad.read_h5ad(result.reference_path).n_obs == 2
    assert result.document['inputs']['reference']['filter_value'] == 'T'
    doc['inputs']['reference']['filter_value'] = 'absent'
    path.write_text(yaml.safe_dump(doc))
    with pytest.raises(ValueError, match='filter'):
        prepare(path)


def test_cache_invalidation_includes_preparation_implementation(tmp_path, monkeypatch):
    from revise.batch import sample
    path, _, _, _ = fixture_sample(tmp_path)
    first = prepare(path)
    identity = sample.file_identity

    def changed_implementation(candidate):
        result = identity(candidate)
        if Path(candidate).resolve() == Path(sample.__file__).resolve():
            result['sha256'] = 'new-code-identity'
        return result

    monkeypatch.setattr(sample, 'file_identity', changed_implementation)
    second = prepare(path)
    assert first.metadata['preparation_code'] != second.metadata['preparation_code']


@pytest.mark.parametrize('value', [None, 2, True, '1'])
def test_requires_supported_schema_version(tmp_path, value):
    path, doc, _, _ = fixture_sample(tmp_path)
    if value is None:
        doc.pop('schema_version', None)
    else:
        doc['schema_version'] = value
    path.write_text(yaml.safe_dump(doc))
    with pytest.raises(ValueError, match='schema_version'):
        prepare(path)


@pytest.mark.parametrize('block', ['sample', 'preparation', 'inputs', 'global_anchoring', 'local_refinement'])
def test_document_blocks_must_be_mappings(tmp_path, block):
    path, doc, _, _ = fixture_sample(tmp_path)
    doc[block] = []
    path.write_text(yaml.safe_dump(doc))
    with pytest.raises(ValueError, match=block):
        prepare(path)


@pytest.mark.parametrize('role', ['spatial', 'reference'])
def test_never_overwrites_unowned_prepared_destination(tmp_path, role):
    path, doc, _, _ = fixture_sample(tmp_path)
    doc['preparation']['reference']['matrix'] = 'layers/counts'
    path.write_text(yaml.safe_dump(doc))
    destination = path.parent / 'inputs' / f'{role}.h5ad'
    destination.parent.mkdir()
    destination.write_bytes(b'original user data')
    with pytest.raises(ValueError, match='unowned'):
        prepare(path)
    assert destination.read_bytes() == b'original user data'


def test_source_cannot_alias_preparation_manifest(tmp_path):
    path, doc, raw, _ = fixture_sample(tmp_path)
    destination = path.parent / 'inputs' / 'preparation.json'
    destination.parent.mkdir()
    destination.write_bytes(raw.read_bytes())
    doc['inputs']['st']['path'] = 'inputs/preparation.json'
    path.write_text(yaml.safe_dump(doc))
    before = destination.read_bytes()
    with pytest.raises(ValueError, match='source'):
        prepare(path)
    assert destination.read_bytes() == before


@pytest.mark.parametrize('inside_package', [False, True])
def test_preparation_rejects_symlink_inputs_directory(tmp_path, inside_package):
    path, _, _, _ = fixture_sample(tmp_path)
    destination = (path.parent if inside_package else tmp_path) / 'redirected-inputs'
    destination.mkdir()
    (path.parent / 'inputs').symlink_to(destination, target_is_directory=True)
    with pytest.raises(ValueError, match='symlink'):
        prepare(path)
    assert list(destination.iterdir()) == []


@pytest.mark.parametrize('name', ['spatial.h5ad', 'reference.h5ad', 'preparation.json'])
@pytest.mark.parametrize('dangling', [False, True])
def test_preparation_rejects_symlink_outputs_before_cache_reuse(tmp_path, name, dangling):
    path, doc, _, _ = fixture_sample(tmp_path)
    doc['preparation']['reference']['matrix'] = 'layers/counts'
    path.write_text(yaml.safe_dump(doc))
    prepare(path)
    destination = path.parent / 'inputs' / name
    external = tmp_path / name
    before = destination.read_bytes()
    if dangling:
        destination.unlink()
    else:
        destination.rename(external)
    destination.symlink_to(external)
    with pytest.raises(ValueError, match='symlink'):
        prepare(path)
    assert destination.is_symlink()
    assert not external.exists() if dangling else external.read_bytes() == before
