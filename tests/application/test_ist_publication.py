from dataclasses import replace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, read_h5ad
from scipy import sparse

from revise.application.config import ApplicationConfigError, compile_application_config, load_application_yaml
from revise.application.publication import application_metadata, output_paths, publish_outputs
from tests.application.test_request import _document, _write_config
from tests.application.test_publication import _ctx


def config(tmp_path, mapping='paired'):
    document = _document('sc-SVC', 'cluster')
    document['output'] = {'dir': 'out', 'ist_mapping': mapping}
    source, document = load_application_yaml(_write_config(tmp_path, document))
    return compile_application_config(document, source=source, cwd=tmp_path)


def carriers():
    spatial = AnnData(sparse.csr_matrix([[0], [0], [0]]), obs=pd.DataFrame({'SVC_cluster': ['b', 'a', 'a']}, index=['s0', 's1', 's2']))
    spatial.obsm['spatial'] = np.array([[1, 2], [3, 4], [5, 6]])
    expression = AnnData(sparse.csr_matrix([[2., 0], [4., 6], [8., 10]]), obs=pd.DataFrame({'SVC_cluster': ['a', 'a', 'b']}, index=['d2', 'd1', 'd3']), var=pd.DataFrame(index=['g1', 'g2']))
    return {'sc_svc_spatial': spatial, 'sc_svc_expr': expression}


@pytest.mark.parametrize('mapping', ['paired', 'mean', 'random'])
def test_modes_compile(tmp_path, mapping):
    assert config(tmp_path, mapping).ist_mapping == mapping


@pytest.mark.parametrize('route, mode', [('sp-SVC', None), ('sc-SVC', 'sr')])
def test_non_ist_rejects_explicit_mapping(tmp_path, route, mode):
    document = _document(route, mode)
    document['output']['ist_mapping'] = 'paired'
    source, document = load_application_yaml(_write_config(tmp_path, document))
    with pytest.raises(ApplicationConfigError, match='only valid for sc-SVC cluster'):
        compile_application_config(document, source=source, cwd=tmp_path)


def test_mean_is_sparse_and_uses_spatial_axis(tmp_path):
    cfg = config(tmp_path, 'mean')
    inputs = carriers()
    result = publish_outputs(cfg, output_paths(cfg), _ctx(tmp_path, inputs))
    assert output_paths(cfg)['svc'].name == 'SVC.h5ad'
    assert sparse.issparse(result.X)
    np.testing.assert_allclose(result.X.toarray(), [[8, 10], [3, 3], [3, 3]])
    assert list(result.obs_names) == ['s0', 's1', 's2']
    assert list(result.var_names) == ['g1', 'g2']
    np.testing.assert_array_equal(result.obsm['spatial'], inputs['sc_svc_spatial'].obsm['spatial'])
    assert result.uns['revise_reconstruction']['ist_mapping'] == 'mean'


def test_random_donors_are_seeded_and_reference_order_independent(tmp_path):
    cfg = config(tmp_path, 'random')
    inputs = carriers()
    first = publish_outputs(cfg, output_paths(cfg), _ctx(tmp_path, inputs))
    inputs['sc_svc_expr'] = inputs['sc_svc_expr'][[2, 1, 0]].copy()
    second = publish_outputs(cfg, output_paths(cfg), _ctx(tmp_path, inputs))
    assert list(first.obs.revise_ist_donor_id) == list(second.obs.revise_ist_donor_id)
    np.testing.assert_array_equal(first.X.toarray(), second.X.toarray())
    for i, donor in enumerate(first.obs.revise_ist_donor_id):
        np.testing.assert_array_equal(first.X[i].toarray(), inputs['sc_svc_expr'][donor].X.toarray())
    assert first.uns['revise_reconstruction']['effective_seed'] == cfg.seed


def test_switch_mode_cleans_owned_outputs_only_after_commit_and_can_rollback(tmp_path):
    paired = config(tmp_path)
    previous = _ctx(tmp_path, carriers(), application_config_metadata=application_metadata(paired, paths=output_paths(paired)))
    publish_outputs(paired, output_paths(paired), previous)
    previous.pending_publication[0]()
    original = {p: p.read_bytes() for p in output_paths(paired).values()}
    mean = replace(paired, ist_mapping='mean')
    ctx = _ctx(tmp_path, carriers())
    publish_outputs(mean, output_paths(mean), ctx)
    ctx.pending_publication[1]()
    assert not output_paths(mean)['svc'].exists()
    assert all(p.read_bytes() == content for p, content in original.items())
    ctx = _ctx(tmp_path, carriers())
    publish_outputs(mean, output_paths(mean), ctx)
    ctx.pending_publication[0]()
    assert all(not p.exists() for p in original)
    assert read_h5ad(output_paths(mean)['svc']).shape == (3, 2)


def test_unowned_alternative_is_preserved(tmp_path):
    cfg = config(tmp_path, 'mean')
    cfg.output_dir.mkdir(parents=True)
    alien = cfg.output_dir / 'spatial.h5ad'
    alien.write_bytes(b'user file')
    ctx = _ctx(tmp_path, carriers())
    publish_outputs(cfg, output_paths(cfg), ctx)
    ctx.pending_publication[0]()
    assert alien.read_bytes() == b'user file'


def test_overflow_does_not_publish_nonfinite_mean(tmp_path):
    cfg = config(tmp_path, 'mean')
    inputs = carriers()
    inputs['sc_svc_expr'].X = np.array([[1e308, 1], [1e308, 2], [1, 3]])
    with pytest.raises(ValueError, match='output X must contain only finite'):
        publish_outputs(cfg, output_paths(cfg), _ctx(tmp_path, inputs))
    assert not output_paths(cfg)['svc'].exists()


def test_bad_cluster_sets_do_not_replace_prior_outputs(tmp_path):
    cfg = config(tmp_path, 'mean')
    inputs = carriers()
    ctx = _ctx(tmp_path, inputs)
    publish_outputs(cfg, output_paths(cfg), ctx)
    ctx.pending_publication[0]()
    before = output_paths(cfg)['svc'].read_bytes()
    inputs['sc_svc_expr'].obs['SVC_cluster'] = ['a', 'a', 'c']
    with pytest.raises(ValueError, match='cluster sets must match'):
        publish_outputs(cfg, output_paths(cfg), _ctx(tmp_path, inputs))
    assert output_paths(cfg)['svc'].read_bytes() == before


def test_failed_switch_restores_old_pair(monkeypatch, tmp_path):
    import revise.application.publication as publication

    paired = config(tmp_path)
    ctx = _ctx(tmp_path, carriers(), application_config_metadata=application_metadata(paired, paths=output_paths(paired)))
    publish_outputs(paired, output_paths(paired), ctx)
    ctx.pending_publication[0]()
    before = {p: p.read_bytes() for p in output_paths(paired).values()}
    mean = replace(paired, ist_mapping='mean')
    target = output_paths(mean)['svc']
    original_replace = publication.os.replace

    def fail_install(source, destination):
        if destination == target and str(source).endswith('.tmp.h5ad'):
            raise OSError('injected installation failure')
        return original_replace(source, destination)

    monkeypatch.setattr(publication.os, 'replace', fail_install)
    with pytest.raises(OSError, match='injected installation'):
        publish_outputs(mean, output_paths(mean), _ctx(tmp_path, carriers()))
    assert not target.exists()
    assert all(p.read_bytes() == contents for p, contents in before.items())


@pytest.mark.parametrize('mapping', ['bogus', None, 2])
def test_invalid_mapping_rejected(tmp_path, mapping):
    with pytest.raises(ApplicationConfigError, match='output.ist_mapping'):
        config(tmp_path, mapping)


@pytest.mark.parametrize('foreign_field', ['output_dir', 'output_paths', 'resolved_inputs'])
def test_foreign_same_type_alternative_is_preserved(tmp_path, foreign_field):
    from revise.application.publication import application_metadata

    cfg = config(tmp_path, 'paired')
    paths = output_paths(cfg)
    metadata = application_metadata(cfg, paths=paths)
    if foreign_field == 'output_dir':
        metadata[foreign_field] = str(tmp_path / 'another_sample' / 'T')
    elif foreign_field == 'output_paths':
        metadata[foreign_field] = {role: str(tmp_path / 'another_sample' / p.name) for role, p in paths.items()}
    else:
        metadata[foreign_field]['st'] = str(tmp_path / 'another_input.h5ad')
    ctx = _ctx(tmp_path, carriers(), application_config_metadata=metadata)
    publish_outputs(cfg, paths, ctx)
    ctx.pending_publication[0]()
    before = {p: p.read_bytes() for p in paths.values()}
    mean = replace(cfg, ist_mapping='mean')
    ctx = _ctx(tmp_path, carriers())
    publish_outputs(mean, output_paths(mean), ctx)
    ctx.pending_publication[0]()
    assert all(p.exists() and p.read_bytes() == content for p, content in before.items())
