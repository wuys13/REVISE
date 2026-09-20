from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import unquote

import numpy as np
import pandas as pd
import pytest
import yaml
from anndata import AnnData, read_h5ad

from revise.application.delivery import sample_document
from revise.application.publication import application_metadata, output_paths, publish_outputs
from tests.application.test_ist_publication import config
from tests.application.test_publication import _ctx


def delivery_fixture(tmp_path):
    cfg = config(tmp_path, 'random')
    cfg = replace(cfg, select_cell_type=None, output_dir=cfg.output_root,
                  delivery_sample_id='CRC/S01',
                  delivery_coordinates={'unit': 'um', 'microns_per_coordinate': 1.0})
    raw = AnnData(np.array([[1., 2., 3.], [4., 5., 6.], [0., 0., 0.]]),
                  obs=pd.DataFrame({'Level1': ['old', 'B', 'excluded'],
                                    'Level2': ['old-subtype', None, 'original']},
                                   index=['s0', 's1', 'qc-excluded']),
                  var=pd.DataFrame(index=['g1', 'g2', 'raw-only']))
    raw.obsm['spatial'] = np.array([[1., 2.], [3., 4.], [5., 6.]])
    spatial = raw[:2, :2].copy()
    spatial.obs['Level1'] = ['T', 'B']
    spatial.obs['Level2'] = ['t1', 'b1']
    spatial.obs['SVC_cluster'] = ['T:0', 'B:0']
    expression = spatial.copy()
    expression.obs_names = ['d0', 'd1']
    outputs = {'sc_svc_spatial': spatial, 'sc_svc_expr': expression}
    ctx = _ctx(tmp_path, outputs, application_config_metadata=application_metadata(cfg, paths=output_paths(cfg)))
    ctx.runner = SimpleNamespace(st_adata=spatial)
    ctx.artifacts = {'processed_cell_types': ['T', 'B'], 'skipped_cell_types': [
        {'cell_type': 'Mono', 'reason': 'insufficient_valid_reference_subtypes', 'valid_subtype_count': 1}]}
    return cfg, raw, ctx


def test_full_raw_inferences_and_consumer_document(tmp_path):
    cfg, raw, ctx = delivery_fixture(tmp_path)
    result = publish_outputs(cfg, output_paths(cfg), ctx, raw=raw)
    ctx.pending_publication[0]()
    paths = output_paths(cfg)
    stored = read_h5ad(paths['raw'])
    np.testing.assert_array_equal(stored.X, raw.X)
    np.testing.assert_array_equal(stored.obsm['spatial'], raw.obsm['spatial'])
    assert stored.obs_names.equals(raw.obs_names)
    assert stored.var_names.equals(raw.var_names)
    assert list(stored.obs.Level1) == list(raw.obs.Level1)
    assert stored.obs.loc['s0', 'revise_Level1'] == 'T'
    assert stored.obs.loc['s0', 'revise_Level2'] == 't1'
    assert pd.isna(stored.obs.loc['qc-excluded', 'revise_Level1'])
    assert pd.isna(stored.obs.loc['qc-excluded', 'revise_Level2'])
    assert 'revise_Level1' not in raw.obs
    assert result.shape == (2, 2)
    assert result.uns['revise_reconstruction']['skipped_cell_types']['0']['cell_type'] == 'Mono'
    assert list(result.obs.revise_Level1) == ['T', 'B']
    assert stored.uns['revise_delivery']['original_label_conflicts']['Level1'] == 1
    document = yaml.safe_load(paths['sample_config'].read_text())
    assert unquote(document['sample_id']) == 'CRC/S01'
    assert document['files'] == {'raw': 'raw.h5ad', 'svc': 'SVC.h5ad'}
    assert document['spatial']['unit'] == 'micron'
    assert document['expression']['svc']['scale'] == 'unknown'


def test_delivery_normalizes_inferred_labels_before_comparing_original_annotations(tmp_path):
    cfg, raw, ctx = delivery_fixture(tmp_path)
    raw.obs["Level1"] = ["T", "Mono_Macro", "excluded"]
    ctx.runner.st_adata.obs["Level1"] = ["T", "Mono/Macro"]
    ctx.svc.artifacts["outputs"]["sc_svc_spatial"].obs["Level1"] = ["T", "Mono/Macro"]

    publish_outputs(cfg, output_paths(cfg), ctx, raw=raw)
    ctx.pending_publication[0]()

    stored = read_h5ad(output_paths(cfg)["raw"])
    assert stored.obs.loc["s1", "revise_Level1"] == "Mono_Macro"
    assert stored.uns["revise_delivery"]["original_label_conflicts"]["Level1"] == 0


@pytest.mark.parametrize('failure', ['write', 'install', 'record'])
@pytest.mark.parametrize('mapping', ['random', 'within_cluster', 'outside_cluster'])
def test_delivery_failure_preserves_all_previous_bytes(monkeypatch, tmp_path, failure, mapping):
    import revise.application.publication as publication

    cfg, raw, ctx = delivery_fixture(tmp_path)
    cfg = replace(cfg, ist_mapping=mapping)
    paths = output_paths(cfg)
    cfg.output_dir.mkdir(parents=True)
    old = {role: f'old-{role}'.encode() for role in paths}
    for role, path in paths.items():
        path.write_bytes(old[role])
    if failure == 'write':
        original = AnnData.write_h5ad
        def write(data, filename, *args, **kwargs):
            if '.raw.h5ad.' in str(filename):
                raise OSError('injected Raw write failure')
            return original(data, filename, *args, **kwargs)
        monkeypatch.setattr(AnnData, 'write_h5ad', write)
    elif failure == 'install':
        original = publication.os.replace
        def install(source, target):
            if Path(target) == paths['sample_config'] and str(source).endswith('.tmp.yaml'):
                raise OSError('injected YAML install failure')
            return original(source, target)
        monkeypatch.setattr(publication.os, 'replace', install)
    else:
        def record(_):
            raise OSError('injected provenance failure')
        ctx.record_artifact = record
    with pytest.raises(OSError, match='injected'):
        publish_outputs(cfg, paths, ctx, raw=raw)
    assert {role: path.read_bytes() for role, path in paths.items()} == old
    assert not list(cfg.output_dir.glob('*.backup'))
    assert not list(cfg.output_dir.glob('*.tmp.*'))


@pytest.mark.parametrize('problem', ['duplicate', 'foreign', 'reserved'])
def test_invalid_raw_or_join_never_replaces_outputs(tmp_path, problem):
    cfg, raw, ctx = delivery_fixture(tmp_path)
    if problem == 'duplicate':
        raw.obs_names = ['s0', 's0', 'other']
    elif problem == 'foreign':
        ctx.runner.st_adata.obs_names = ['foreign', 's1']
    else:
        raw.obs['revise_Level2'] = 'user-owned'
    with pytest.raises(ValueError):
        publish_outputs(cfg, output_paths(cfg), ctx, raw=raw)
    assert not cfg.output_dir.exists()


def test_sample_identity_encoding_is_reversible_and_distinct(tmp_path):
    cfg, raw, ctx = delivery_fixture(tmp_path)
    ids = ['CRC/S01', 'CRC%2FS01', 'CRC_S01', '.hidden', 'a\\b']
    encoded = [sample_document(cfg, output_paths(cfg), raw, raw, sample_id=value)['sample_id'] for value in ids]
    assert len(set(encoded)) == len(ids)
    assert [unquote(value) for value in encoded] == ids
    assert all('/' not in value and '\\' not in value and not value.startswith('.') for value in encoded)


def test_source_changed_during_reconstruction_rejects_delivery(monkeypatch, tmp_path):
    import reconstruct

    cfg, raw, ctx = delivery_fixture(tmp_path)
    cfg.st_path.parent.mkdir(parents=True)
    raw.write_h5ad(cfg.st_path)
    def run(self, **kwargs):
        changed = raw.copy()
        changed.X[0, 0] = 99
        changed.write_h5ad(cfg.st_path)
        kwargs['finalize_callback'](ctx)
    monkeypatch.setattr(reconstruct.REVISEPipeline, 'run', run)
    with pytest.raises(ValueError, match='changed during reconstruction'):
        reconstruct.reconstruct(raw[:2].copy(), raw.copy(), cfg)
    assert not cfg.output_dir.exists()


def test_explicit_original_object_supports_in_memory_delivery(monkeypatch, tmp_path):
    import reconstruct

    cfg, raw, ctx = delivery_fixture(tmp_path)
    def run(self, **kwargs):
        kwargs['finalize_callback'](ctx)
        ctx.pending_publication[0]()
    monkeypatch.setattr(reconstruct.REVISEPipeline, 'run', run)
    result = reconstruct.reconstruct(raw[:2].copy(), raw.copy(), cfg, raw_adata=raw)
    assert result.n_obs == 2
    np.testing.assert_array_equal(read_h5ad(output_paths(cfg)['raw']).X, raw.X)


def test_explicit_raw_is_frozen_before_aliased_working_input_mutates(monkeypatch, tmp_path):
    import reconstruct

    cfg, raw, ctx = delivery_fixture(tmp_path)
    original = raw.X.copy()
    def run(self, **kwargs):
        kwargs['st_adata'].X[:] = 10000
        kwargs['finalize_callback'](ctx)
        ctx.pending_publication[0]()
    monkeypatch.setattr(reconstruct.REVISEPipeline, 'run', run)
    reconstruct.reconstruct(raw, raw.copy(), cfg, raw_adata=raw)
    np.testing.assert_array_equal(read_h5ad(output_paths(cfg)['raw']).X, original)


def test_missing_original_source_cannot_be_presented_as_complete_raw(tmp_path):
    import reconstruct

    cfg, raw, ctx = delivery_fixture(tmp_path)
    with pytest.raises(ValueError, match='requires an original spatial object'):
        reconstruct.reconstruct(raw[:2].copy(), raw.copy(), cfg)


def test_publication_cannot_replace_original_input(tmp_path):
    cfg, raw, ctx = delivery_fixture(tmp_path)
    cfg = replace(cfg, st_path=cfg.output_dir / 'raw.h5ad')
    cfg.output_dir.mkdir(parents=True)
    raw.write_h5ad(cfg.st_path)
    previous = cfg.st_path.read_bytes()
    with pytest.raises(ValueError, match='aliases original input'):
        publish_outputs(cfg, output_paths(cfg), ctx, raw=raw)
    assert cfg.st_path.read_bytes() == previous


def test_original_raw_reconstruction_metadata_is_not_replaced(tmp_path):
    cfg, raw, ctx = delivery_fixture(tmp_path)
    original = {'run_manifest': '/original/history.json', 'output_paths': {'svc': '/original/svc.h5ad'}}
    raw.uns['revise_reconstruction'] = original
    publish_outputs(cfg, output_paths(cfg), ctx, raw=raw)
    ctx.pending_publication[0]()
    stored = read_h5ad(output_paths(cfg)['raw'])
    assert stored.uns['revise_reconstruction'] == original
    assert stored.uns['revise_delivery']['publication']['output_role'] == 'raw'
    assert stored.uns['revise_delivery']['publication']['run_manifest'] == str(ctx.run_dir / 'provenance.json')


@pytest.mark.parametrize('route, mode, output_key', [
    ('sp-SVC', None, 'sp_svc'), ('sc-SVC', 'sr', 'sc_svc_dec')])
def test_other_routes_publish_raw_without_inventing_cell_correspondence(tmp_path, route, mode, output_key):
    cfg, raw, ctx = delivery_fixture(tmp_path)
    cfg = replace(cfg, svc_type=route, mode=mode, subtype_column=None)
    svc = ctx.runner.st_adata.copy()
    svc.obs = svc.obs.drop(columns=['Level2', 'SVC_cluster'])
    if mode == 'sr':
        svc.obs['cell_type'] = svc.obs.pop('Level1')
        svc.obs['spot_name'] = svc.obs_names
        svc.obs_names = ['generated-0', 'generated-1']
    ctx.svc.artifacts['outputs'] = {output_key: svc}
    ctx.svc.provenance = {'primary_output_key': output_key}
    paths = output_paths(cfg)
    publish_outputs(cfg, paths, ctx, raw=raw)
    ctx.pending_publication[0]()
    published = read_h5ad(paths['raw'])
    assert published.obs_names.equals(raw.obs_names)
    assert published.obs.loc['s0', 'revise_Level1'] == 'T'
    assert 'revise_Level2' not in published.obs
    declaration = yaml.safe_load(paths['sample_config'].read_text())
    assert 'reconstruction' not in declaration['columns']
    assert 'subtype' not in declaration['columns']
    assert 'subtype' not in published.uns['revise_delivery']['annotation_sources']
    assert declaration['expression']['svc']['scale'] == 'unknown'
    if mode == 'sr':
        delivered_svc = read_h5ad(paths['svc'])
        assert declaration['columns']['broad'] == 'revise_Level1'
        assert delivered_svc.obs['revise_Level1'].astype(str).tolist() == svc.obs['cell_type'].astype(str).tolist()
        assert 'cell_type' in delivered_svc.obs
        assert declaration['expression']['svc']['identity'] == 'sst_parent_spot_corrected_expression'
        assert 'not raw counts' in declaration['provenance']['svc_expression_processing']


def test_declared_missing_spatial_key_fails_before_publishing(tmp_path):
    cfg, raw, ctx = delivery_fixture(tmp_path)
    cfg = replace(cfg, delivery_coordinates={'key': 'nonexistent', 'unit': 'micron'})
    with pytest.raises(ValueError, match='missing declared spatial coordinates'):
        publish_outputs(cfg, output_paths(cfg), ctx, raw=raw)
    assert not output_paths(cfg)['sample_config'].exists()
