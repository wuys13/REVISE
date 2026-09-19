"""Optional joint check against the independently developed consumer checkout."""
from dataclasses import replace
import json
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise.application.publication import application_metadata, output_paths, publish_outputs
from revise.utils.provenance import sha256_file
from tests.application.test_ist_publication import config
from tests.application.test_publication import _ctx


def test_published_sample_is_consumed_by_real_analysis_entrypoint(tmp_path, monkeypatch):
    consumer = Path(os.environ.get('REVISE_ANALYSIS_ROOT',
                    str(Path(__file__).resolve().parents[3] / 'REVISE_Analysis_Agent')))
    if not (consumer / 'revise_analysis/io.py').is_file():
        pytest.skip('Set REVISE_ANALYSIS_ROOT to the independent analysis checkout')
    monkeypatch.syspath_prepend(str(consumer))
    from revise_analysis.io import load_sample
    from revise_analysis.runner import run_analysis

    cfg = config(tmp_path, 'random')
    cfg = replace(cfg, select_cell_type=None, output_dir=cfg.output_root,
                  delivery_sample_id='synthetic/whole-sample',
                  delivery_coordinates={'unit': 'um', 'microns_per_coordinate': 1.0})
    n = 36
    raw = AnnData(np.ones((n + 1, 3)),
                  obs=pd.DataFrame({'Level1': ['original'] * (n + 1)},
                                   index=[f's{i}' for i in range(n + 1)]),
                  var=pd.DataFrame(index=['g1', 'g2', 'raw-only']))
    raw.obsm['spatial'] = np.array([(i % 6, i // 6) for i in range(n + 1)], dtype=float)
    spatial = raw[:n, :2].copy()
    spatial.obs['Level1'] = ['T'] * n
    spatial.obs['Level2'] = ['t1' if i % 2 else 't2' for i in range(n)]
    spatial.obs['SVC_cluster'] = [f'T:{i % 3}' for i in range(n)]
    donor = spatial.copy()
    donor.obs_names = [f'd{i}' for i in range(n)]
    paths = output_paths(cfg)
    ctx = _ctx(tmp_path, {'sc_svc_spatial': spatial, 'sc_svc_expr': donor},
               application_config_metadata=application_metadata(cfg, paths=paths))
    ctx.runner = SimpleNamespace(st_adata=spatial)
    publish_outputs(cfg, paths, ctx, raw=raw)
    ctx.pending_publication[0]()
    hashes = {role: sha256_file(path) for role, path in paths.items()}
    elsewhere = tmp_path / 'elsewhere'
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    sample = load_sample(paths['sample_config'])
    assert sample.raw.shape == (37, 3)
    assert sample.svc.shape == (36, 2)
    assert sample.broad_key == 'revise_Level1'
    result = run_analysis(paths['sample_config'], 'reconstruction_impact', tmp_path / 'analysis', {
        'scopes': ['T'], 'min_window_units': 2, 'n_window_draws': 8,
        'region_n_bootstrap': 5, 'window_side_microns': 3.0,
        'anatomy_window_side_microns': 3.0, 'window_scale_candidates_microns': [2.0, 3.0],
    })
    assert not result.get('stage_errors'), result.get('stage_errors')
    assert result['status'] in {'partial', 'succeeded'}, result
    stages = {item['stage']: item for item in result['stages']}
    assert stages['diversity']['status'] == 'completed', stages['diversity']
    destination = tmp_path / 'analysis' / sample.sample_id / 'reconstruction_impact'
    assert (destination / 'tables/window_diversity_svc_T.csv').is_file()
    saved = json.loads((destination / 'result.json').read_text())
    assert saved['status'] == result['status']
    assert (destination / 'report.html').is_file()
    assert hashes == {role: sha256_file(path) for role, path in paths.items()}


def test_declared_linear_delivery_runs_real_expression_analysis(tmp_path, monkeypatch):
    consumer = Path(os.environ.get('REVISE_ANALYSIS_ROOT',
                    str(Path(__file__).resolve().parents[3] / 'REVISE_Analysis_Agent')))
    if not (consumer / 'revise_analysis/io.py').is_file():
        pytest.skip('Set REVISE_ANALYSIS_ROOT to the independent analysis checkout')
    monkeypatch.syspath_prepend(str(consumer))
    from revise_analysis.io import load_sample
    from revise_analysis.runner import run_analysis

    declaration = {'identity': 'measured_expression', 'scale': 'untransformed_nonnegative'}
    cfg = config(tmp_path, 'mean')
    cfg = replace(cfg, select_cell_type=None, output_dir=cfg.output_root,
                  st_expression=declaration, reference_expression=declaration)
    rng = np.random.default_rng(42)
    raw = AnnData(rng.uniform(.01, .9, size=(16, 4)),
                  obs=pd.DataFrame(index=[f's{i}' for i in range(16)]),
                  var=pd.DataFrame(index=[f'g{i}' for i in range(4)]))
    raw.obsm['spatial'] = np.array([(i % 4, i // 4) for i in range(16)], dtype=float)
    spatial = raw.copy()
    spatial.obs['Level1'] = 'T'
    spatial.obs['Level2'] = [f't{i % 2}' for i in range(16)]
    spatial.obs['SVC_cluster'] = [f'c{i % 4}' for i in range(16)]
    donor = spatial.copy()
    donor.obs_names = [f'd{i}' for i in range(16)]
    paths = output_paths(cfg)
    ctx = _ctx(tmp_path, {'sc_svc_spatial': spatial, 'sc_svc_expr': donor},
               application_config_metadata=application_metadata(cfg, paths=paths))
    ctx.runner = SimpleNamespace(st_adata=spatial)
    publish_outputs(cfg, paths, ctx, raw=raw)
    ctx.pending_publication[0]()
    hashes = {role: sha256_file(path) for role, path in paths.items()}
    sample = load_sample(paths['sample_config'])
    assert sample.expression_unavailable('raw') is None
    assert sample.expression_unavailable('svc') is None
    result = run_analysis(paths['sample_config'], 'spatial_autocorrelation', tmp_path / 'expression',
                          {'moran_n_neighbors': 3, 'random_state': 42})
    assert result['status'] == 'succeeded', result
    destination = tmp_path / 'expression' / sample.sample_id / 'spatial_autocorrelation'
    for side in ('raw', 'svc'):
        table = pd.read_csv(destination / f'tables/moran_{side}.csv')
        assert not table.empty
        if 'status' in table:
            assert table['status'].isin(['computed', 'ok']).all(), table
    assert hashes == {role: sha256_file(path) for role, path in paths.items()}
