import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml
from scipy import sparse


def fixture_sample(tmp_path):
    root = tmp_path / 'ST/CRC/S01'
    root.mkdir(parents=True)
    st = root / 'spatial.h5ad'
    ref = tmp_path / 'shared.h5ad'
    a = ad.AnnData(sparse.csr_matrix([[1., 2.], [3., 4.]]),
                   obs=pd.DataFrame(index=['c1', 'c2']), var=pd.DataFrame(index=['g1', 'g2']))
    a.obsm['spatial'] = np.array([[1., 2.], [3., 4.]])
    a.write_h5ad(st)
    a.obs['Level1'] = ['T', 'Macro']
    a.obs['Level2'] = ['T1', 'M1']
    a.write_h5ad(ref)
    doc = {'schema_version': 2, 'input_root': 'ST', 'output_root': 'output',
           'modality': 'iST', 'coordinates': {'unit': 'pixel'},
           'inputs': {'reference': {'path': 'shared.h5ad'}},
           'global_anchoring': {'broad_column': 'Level1'},
           'local_refinement': {'subtype_column': 'Level2'}}
    path = tmp_path / 'batch.yaml'
    path.write_text(yaml.safe_dump(doc))
    return path, doc, st, ref


def read(path, st):
    from revise.batch.config import resolve_sample
    from revise.batch.sample import read_sample
    return read_sample(resolve_sample(path, st.parent))


def test_shared_inputs_are_read_only_without_new_files(tmp_path):
    path, _, st, ref = fixture_sample(tmp_path)
    before = {str(p): p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
    sample = read(path, st)
    assert sample.sample_id == 'CRC/S01'
    assert sample.st_path == st
    assert sample.reference_path == ref
    assert sample.metadata['outputs'] == sample.metadata['sources']
    assert sample.metadata['coordinates']['physical_distance_available'] is False
    assert sample.metadata['configuration']['project_config'] == str(path)
    assert before == {str(p): p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}


def test_partially_missing_reference_subtypes_are_left_for_per_type_eligibility(tmp_path):
    path, _, st, ref = fixture_sample(tmp_path)
    data = ad.read_h5ad(ref)
    data.obs['Level2'] = data.obs['Level2'].astype(object)
    data.obs.loc['c2', 'Level2'] = None
    data.write_h5ad(ref)

    sample = read(path, st)

    assert sample.reference_path == ref


@pytest.mark.parametrize('modality', ['hST', 'sST'])
def test_non_ist_routes_still_reject_missing_reference_broad_labels(tmp_path, modality):
    path, doc, st, ref = fixture_sample(tmp_path)
    doc['modality'] = modality
    doc['local_refinement'] = {'strength': 0.2}
    path.write_text(yaml.safe_dump(doc))
    data = ad.read_h5ad(ref)
    data.obs['Level1'] = data.obs['Level1'].astype(object)
    data.obs.loc['c2', 'Level1'] = None
    data.write_h5ad(ref)

    with pytest.raises(ValueError, match='non-null'):
        read(path, st)


@pytest.mark.parametrize('problem', ['duplicate_obs', 'duplicate_var', 'nan_matrix', 'negative_matrix',
                                      'nan_coordinates', 'missing_coordinates', 'missing_label', 'no_shared_genes'])
def test_invalid_inputs(tmp_path, problem):
    path, _, st, ref = fixture_sample(tmp_path)
    source = ref if problem in {'missing_label', 'no_shared_genes'} else st
    a = ad.read_h5ad(source)
    if problem == 'duplicate_obs':
        a.obs_names = ['x', 'x']
    elif problem == 'duplicate_var':
        a.var_names = ['x', 'x']
    elif problem == 'nan_matrix':
        a.X.data[0] = np.nan
    elif problem == 'negative_matrix':
        a.X.data[0] = -1
    elif problem == 'nan_coordinates':
        a.obsm['spatial'][0, 0] = np.nan
    elif problem == 'missing_coordinates':
        del a.obsm['spatial']
    elif problem == 'no_shared_genes':
        a.var_names = ['other1', 'other2']
    else:
        del a.obs['Level1']
    a.write_h5ad(source)
    with pytest.raises(ValueError):
        read(path, st)


def test_reference_filter_is_validated_without_filtering_source(tmp_path):
    path, doc, st, ref = fixture_sample(tmp_path)
    doc['inputs']['reference'].update(filter_column='Level1', filter_value='T')
    path.write_text(yaml.safe_dump(doc))
    assert read(path, st).reference_path == ref
    assert ad.read_h5ad(ref).n_obs == 2
    doc['inputs']['reference']['filter_value'] = 'absent'
    path.write_text(yaml.safe_dump(doc))
    with pytest.raises(ValueError, match='filter'):
        read(path, st)


@pytest.mark.parametrize('coordinates', [{}, {'unit': 'guess'}, {'unit': 'pixel', 'microns_per_coordinate': -1},
                                         {'unit': 'um', 'microns_per_coordinate': 2}])
def test_coordinates_must_be_explicit_and_consistent(tmp_path, coordinates):
    path, doc, st, _ = fixture_sample(tmp_path)
    doc['coordinates'] = coordinates
    path.write_text(yaml.safe_dump(doc))
    with pytest.raises(ValueError):
        read(path, st)
