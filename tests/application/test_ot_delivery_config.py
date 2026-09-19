"""New assembly settings and explicit expression-source declarations."""
from dataclasses import replace

import numpy as np
import pytest
from anndata import AnnData

from revise.application.config import ApplicationConfigError, compile_application_config, load_application_yaml
from revise.application.expression import bind_expression_sources, consumer_declaration
from revise.application.publication import application_metadata, output_paths
from tests.application.test_request import _document, _write_config


def compile_doc(tmp_path, doc, **kwargs):
    source, document = load_application_yaml(_write_config(tmp_path, doc))
    return compile_application_config(document, source=source, cwd=tmp_path, **kwargs)


@pytest.mark.parametrize('mode', ['within_cluster', 'outside_cluster'])
def test_new_modes_effective_settings(tmp_path, mode):
    doc = _document('sc-SVC', 'cluster')
    doc['output']['ist_mapping'] = mode
    cfg = compile_doc(tmp_path, doc)
    assert cfg.ist_ot == {'spatial_weight': .2, 'max_cost_entries': 2_000_000, 'gene_block_size': 256}
    first = application_metadata(cfg, paths=output_paths(cfg))
    changed = replace(cfg, ist_ot={**cfg.ist_ot, 'spatial_weight': .4})
    assert first['effective_request_hash'] != application_metadata(changed, paths=output_paths(changed))['effective_request_hash']
    assert first['effective_request']['output']['ist_ot']['method'] == 'tacco'


@pytest.mark.parametrize('options', [{'spatial_weight': 1.1}, {'max_cost_entries': 0},
                                   {'gene_block_size': True}, {'unknown': 1}])
def test_invalid_ot_options(tmp_path, options):
    doc = _document('sc-SVC', 'cluster')
    doc['output'].update(ist_mapping='within_cluster', ist_ot=options)
    with pytest.raises(ApplicationConfigError):
        compile_doc(tmp_path, doc)


def test_expression_override_cannot_inherit_old_reference_declaration(tmp_path):
    doc = _document('sc-SVC', 'cluster')
    declaration = {'identity': 'measured_expression', 'scale': 'untransformed_nonnegative'}
    doc['inputs']['st']['expression'] = declaration
    doc['inputs']['reference']['expression'] = declaration
    cfg = compile_doc(tmp_path, doc)
    assert cfg.reference_expression == declaration
    replacement = tmp_path / 'new_reference.h5ad'
    AnnData(np.ones((2, 2))).write_h5ad(replacement)
    cfg = compile_doc(tmp_path, doc, reference_override=replacement)
    assert cfg.reference_expression is None
    assert cfg.st_expression == declaration


def test_source_declaration_and_unknown_behavior(tmp_path):
    cfg = compile_doc(tmp_path, _document('sc-SVC', 'cluster'))
    raw = AnnData(np.array([[.1, .2]]))
    ref = raw.copy()
    assert consumer_declaration(None, raw)['scale'] == 'unknown'
    ref.uns['revise_expression'] = {'identity': 'measured_reference', 'scale': 'untransformed_nonnegative'}
    bound = bind_expression_sources(cfg, raw, ref)
    assert bound.st_expression is None
    assert consumer_declaration(bound.reference_expression, ref) == {'matrix': 'X', 'identity': 'measured_reference'}
    ref.X[0, 0] = -1
    with pytest.raises(ValueError, match='finite and nonnegative'):
        bind_expression_sources(cfg, raw, ref)


def test_log_declaration_rejected(tmp_path):
    doc = _document('sc-SVC', 'cluster')
    doc['inputs']['st']['expression'] = {'identity': 'measured', 'scale': 'log1p'}
    with pytest.raises(ApplicationConfigError, match='log expression'):
        compile_doc(tmp_path, doc)


@pytest.mark.parametrize('failure', ['cost_limit', 'solver'])
def test_ot_failure_preserves_complete_previous_delivery(tmp_path, monkeypatch, failure):
    from revise.application.publication import publish_outputs
    from revise.utils.provenance import sha256_file
    from tests.application.test_ist_ot import _carriers
    from tests.application.test_ist_publication import config
    from tests.application.test_publication import _ctx

    cfg = config(tmp_path, 'random')
    cfg = replace(cfg, select_cell_type=None, output_dir=cfg.output_root, broad_column='parent')
    spatial, reference = _carriers()
    raw = spatial.copy()
    outputs = {'sc_svc_spatial': spatial, 'sc_svc_expr': reference}
    paths = output_paths(cfg)
    context = _ctx(tmp_path, outputs, application_config_metadata=application_metadata(cfg, paths=paths))
    publish_outputs(cfg, paths, context, raw=raw)
    context.pending_publication[0]()
    before = {role: sha256_file(path) for role, path in paths.items()}
    cfg = replace(cfg, ist_mapping='within_cluster', ist_ot={'max_cost_entries': 1} if failure == 'cost_limit' else None)
    if failure == 'solver':
        def fail(*args, **kwargs):
            raise RuntimeError('injected TACCO failure')
        monkeypatch.setattr('revise.application.ist_ot.OTKernel.couple', fail)
    context = _ctx(tmp_path, outputs, application_config_metadata=application_metadata(cfg, paths=paths))
    with pytest.raises((ValueError, RuntimeError), match='max_cost_entries|injected TACCO'):
        publish_outputs(cfg, paths, context, raw=raw)
    assert before == {role: sha256_file(path) for role, path in paths.items()}
