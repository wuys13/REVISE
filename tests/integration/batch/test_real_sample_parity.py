"""Opt-in real-data route parity smoke tests; not full scientific protocols.

Set REVISE_REAL_DATA_ROOT to REVISE/raw_data. Tests read only backed metadata
and at most 128 spatial observations x 128 genes plus 160 reference cells.
All derived inputs and solver outputs are confined to pytest's temporary path.
"""
from copy import deepcopy
import json
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
from scipy import sparse
import yaml


DATA_ROOT = os.environ.get('REVISE_REAL_DATA_ROOT')
pytestmark = pytest.mark.skipif(not DATA_ROOT, reason='set REVISE_REAL_DATA_ROOT for real-data parity smoke tests')


def _subsets(modality):
    root = Path(DATA_ROOT)
    spatial_paths = {
        'hST': 'Real_application/P1CRC_HD.h5ad',
        'iST': 'Real_application/P2CRC_Xenium.h5ad',
        'sST': 'Sim2Real-ST-P5CRC/spot/part1/spot_50/xenium_spot.h5ad',
    }
    spatial_path = root / spatial_paths[modality]
    reference_path = root / 'Real_application/adata_sc_all_reanno.h5ad'
    assert spatial_path.is_file(), f'Missing required real fixture: {spatial_path}'
    assert reference_path.is_file(), f'Missing required real fixture: {reference_path}'
    spatial = ad.read_h5ad(spatial_path, backed='r')
    reference = ad.read_h5ad(reference_path, backed='r')
    try:
        genes = spatial.var_names.intersection(reference.var_names, sort=False)[:128]
        assert len(genes) >= 32, 'Real fixtures must share at least 32 genes'
        # Stratification prevents accidental loss of all T cells or reference subtypes.
        if modality == 'iST':
            rows = np.flatnonzero(spatial.obs['Level1'].eq('T').to_numpy())[:128]
            reference_rows = []
            t_cells = reference.obs['Level1'].eq('T')
            for subtype in reference.obs.loc[t_cells, 'Level2'].unique():
                reference_rows.extend(np.flatnonzero((t_cells & reference.obs.Level2.eq(subtype)).to_numpy())[:32])
            reference_rows = np.sort(reference_rows[:160])
        else:
            rows = np.arange(min(8 if modality == 'sST' else 128, spatial.n_obs))
            reference_rows = []
            for label in ['T', 'Mono/Macro', 'Fibroblast', 'Tumor', 'B']:
                reference_rows.extend(np.flatnonzero(reference.obs.Level1.eq(label).to_numpy())[:32])
            reference_rows = np.sort(reference_rows)
        assert len(rows) >= 8 and len(reference_rows) >= 32
        st = spatial[rows, genes].to_memory()
        ref = reference[reference_rows, genes].to_memory()
    finally:
        spatial.file.close()
        reference.file.close()
    if modality == 'sST':
        # AnnData slicing does not subset mappings in uns; preserve only active spots.
        st.uns['all_cells_in_spot'] = {
            name: st.uns['all_cells_in_spot'][name] for name in st.obs_names
        }
    st.uns['parity_smoke_source'] = str(spatial_path)
    ref.uns['parity_smoke_source'] = str(reference_path)
    # Explicit minimal fixture thresholds: remove empty measured rows, no inference.
    st = st[np.asarray(st.X.sum(axis=1)).ravel() > 0].copy()
    ref = ref[np.asarray(ref.X.sum(axis=1)).ravel() > 0].copy()
    assert 8 <= st.n_obs <= 128 and 32 <= ref.n_obs <= 160
    return st, ref


def _sample(root, modality):
    root.mkdir(parents=True)
    st, ref = _subsets(modality)
    st.write_h5ad(root / 'source.h5ad')
    ref.write_h5ad(root / 'reference.h5ad')
    local = {'strength': 0.2}
    if modality == 'iST':
        local = {'subtype_column': 'Level2', 'cell_types': ['T'], 'alpha': 0.2, 'resolutions': [0.5]}
    elif modality == 'sST':
        local['match_spot_sum'] = True
    document = {
        'schema_version': 1, 'sample': {'id': f'real-{modality}', 'modality': modality},
        'inputs': {'st': {'path': 'source.h5ad', 'format': 'h5ad'},
                   'reference': {'path': 'reference.h5ad', 'format': 'h5ad'}},
        'preparation': {'spatial': {'matrix': 'X', 'spatial_key': 'spatial', 'coordinate_unit': 'pixel'},
                        'reference': {'matrix': 'X'}},
        'algorithm': {'ot_method': 'tacco' if modality == 'iST' else 'pot'},
        'preprocessing': {
            'spatial': {'min_transcript_counts': None, 'min_counts': 1, 'min_cell_counts': 1},
            'reference': {'min_transcript_counts': None, 'min_genes': 1, 'min_cell_counts': 1}},
        'global_anchoring': {'broad_column': 'Level1'}, 'local_refinement': local,
        'execution': {'seed': 42}, 'output': {'ist_mapping': 'paired'} if modality == 'iST' else {},
    }
    path = root / 'sample.yaml'
    path.write_text(yaml.safe_dump(document))
    return path


@pytest.mark.parametrize('modality', ['hST', 'iST', 'sST'])
def test_real_batch_matches_direct_single_run(tmp_path, modality):
    from reconstruct import run_application
    from revise.application.config import compile_application_config, load_application_yaml
    from revise.application.publication import output_paths
    from revise.batch.runner import application_document, run_batch
    from revise.batch.sample import prepare_sample

    sample_path = _sample(tmp_path / 'data' / 'CRC' / modality, modality)
    batch_path = tmp_path / 'batch.yaml'
    batch_path.write_text(yaml.safe_dump({'schema_version': 1, 'input_root': 'data', 'output_root': 'results'}))
    report = run_batch(batch_path)
    if report['summary']['failed']:
        logs = {str(path): path.read_text()[-12000:] for path in (tmp_path / 'results').rglob('reconstruction.log')}
        pytest.fail(f'Actual batch route failed: {json.dumps(report)}\nLogs: {logs}')
    assert report['summary'] == {'succeeded': 1, 'failed': 0, 'reused': 0}

    prepared = prepare_sample(sample_path)
    document = deepcopy(application_document(prepared, 'T' if modality == 'iST' else None, tmp_path / 'direct'))
    document['output']['dir'] = str(tmp_path / 'direct').lstrip('/')
    direct_yaml = tmp_path / 'direct.yaml'
    direct_yaml.write_text(yaml.safe_dump(document))
    run_application(direct_yaml)
    source, loaded = load_application_yaml(direct_yaml)
    direct_config = compile_application_config(loaded, source=source)
    batch_root = tmp_path / 'results' / 'CRC' / modality
    if modality == 'iST':
        batch_root /= 'T'
    handoff = json.loads((batch_root / 'reconstruction.json').read_text())
    assert handoff['analysis']['status'] == 'not_run'
    for role, direct_path in output_paths(direct_config).items():
        actual = ad.read_h5ad(handoff['outputs'][role]['path'])
        expected = ad.read_h5ad(direct_path)
        assert actual.obs_names.equals(expected.obs_names)
        assert actual.var_names.equals(expected.var_names)
        actual_x = actual.X.toarray() if sparse.issparse(actual.X) else actual.X
        expected_x = expected.X.toarray() if sparse.issparse(expected.X) else expected.X
        np.testing.assert_allclose(actual_x, expected_x, rtol=1e-6, atol=1e-7)
        if 'SVC_cluster' in actual.obs or 'SVC_cluster' in expected.obs:
            assert actual.obs.SVC_cluster.astype(str).tolist() == expected.obs.SVC_cluster.astype(str).tolist()
        assert ('spatial' in actual.obsm) == ('spatial' in expected.obsm)
        if 'spatial' in actual.obsm:
            np.testing.assert_array_equal(actual.obsm['spatial'], expected.obsm['spatial'])
