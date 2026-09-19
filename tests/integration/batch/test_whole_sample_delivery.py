"""Real application/batch lifecycle with deterministic, inexpensive kernels."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import yaml
from anndata import AnnData, read_h5ad


@pytest.mark.parametrize("prepared_reference", [False, True])
@pytest.mark.parametrize("mapping", ["random", "within_cluster", "outside_cluster"])
def test_direct_and_batch_whole_sample_delivery_match(monkeypatch, tmp_path, prepared_reference, mapping):
    from reconstruct import run_application
    from revise.backend import kernels
    from revise.batch import runner
    from revise.batch.config import resolve_sample
    from revise.batch.sample import read_sample

    sample_root = tmp_path / 'data' / 'CRC' / 'S01'
    sample_root.mkdir(parents=True)
    raw = AnnData(np.ones((7, 3)),
                  obs=pd.DataFrame({'transcript_counts': [10] * 6 + [0]},
                                   index=[f's{i}' for i in range(7)]),
                  var=pd.DataFrame(index=['g1', 'g2', 'g3']))
    raw.obsm['spatial'] = np.arange(14).reshape(7, 2).astype(float)
    raw.write_h5ad(sample_root / 'spatial.h5ad')
    reference = AnnData(np.arange(21).reshape(7, 3).astype(float) + 1,
                        obs=pd.DataFrame({'Level1': ['A', 'A', 'A', 'B', 'B', 'C', 'C'],
                                          'Level2': ['a1', 'a2', None, 'b1', 'b2', 'c1', '']},
                                         index=[f'd{i}' for i in range(7)]), var=raw.var.copy())
    reference.write_h5ad(sample_root / 'reference.h5ad')
    reference_config = None
    if prepared_reference:
        from revise.reference_preparation import prepare_reference

        reference.obs['pair'] = 'one'
        reference.write_h5ad(sample_root / 'reference.h5ad')
        preparation = tmp_path / 'prepare.yaml'
        preparation.write_text(yaml.safe_dump({
            'schema_version': 1, 'mode': 'paired',
            'source': str(sample_root / 'reference.h5ad'),
            'pair_column': 'pair', 'pair_key': 'one', 'output_dir': 'prepared',
        }))
        reference_config = prepare_reference(preparation).reference_config_path
    project = tmp_path / 'batch.yaml'
    project.write_text(yaml.safe_dump({'schema_version': 2, 'input_root': 'data', 'output_root': 'results'}))
    document = {
        'modality': 'iST', 'inputs': {'reference': {'path': 'reference.h5ad', 'format': 'h5ad'}},
        'coordinates': {'unit': 'um'}, 'algorithm': {'ot_method': 'pot'},
        'preprocessing': {'spatial': {'min_transcript_counts': 1, 'min_cell_counts': 1},
                          'reference': {'min_transcript_counts': None, 'min_cell_counts': 1}},
        'global_anchoring': {'broad_column': 'Level1'},
        'local_refinement': {'subtype_column': 'Level2', 'alpha': .2, 'resolutions': [.5]},
        'execution': {'seed': 42},
        'output': {'ist_mapping': mapping},
    }
    if reference_config is not None:
        document['inputs']['reference'] = {'config': str(reference_config)}
    (sample_root / 'batch.yaml').write_text(yaml.safe_dump(document))
    calls = []
    def ga(spatial, ref, **kwargs):
        calls.append('GA')
        assert ref.n_obs == 7  # missing subtypes remain available to broad GA
        result = spatial.copy()
        result.obs['Level1'] = ['A', 'A', 'B', 'B', 'C', 'C']
        result.obs['confidence'] = 1.0
        return result
    monkeypatch.setattr(kernels, 'build_kernel', lambda *args, **kwargs: SimpleNamespace(run=ga))
    def anchor(self, spatial, ref, **kwargs):
        label = kwargs['cell_type_col']
        assert ref.obs[label].notna().all()
        result = spatial.copy()
        labels = ref.obs[label].unique()
        result.obs[label] = pd.Categorical([labels[i % len(labels)] for i in range(result.n_obs)])
        return result
    def graph(self, spatial, resolutions, subtype):
        calls.append(str(spatial.obs.Level1.iloc[0]))
        result = spatial.copy()
        result.obs['leiden_0.5'] = pd.Categorical(['0', '1'])
        return result, None, .5
    monkeypatch.setattr(kernels.LocalAnchoringKernel, 'run', anchor)
    monkeypatch.setattr(kernels.GraphClusterKernel, 'run', graph)

    sample = read_sample(resolve_sample(project, sample_root))
    direct_doc = runner.application_document(sample, None, tmp_path / 'direct')
    direct_path = tmp_path / 'direct.yaml'
    direct_path.write_text(yaml.safe_dump(direct_doc))
    run_application(direct_path, reference_config=reference_config)
    assert calls == ['GA', 'A', 'B']
    calls.clear()
    def execute(path, log, *, reference_config=None):
        run_application(path, reference_config=reference_config)
        log.write_text('in-process application integration fixture\n')
        return 0
    monkeypatch.setattr(runner, '_execute', execute)
    result = runner.run_batch(project)
    assert result['summary']['succeeded'] == 1, result
    assert calls == ['GA', 'A', 'B']
    batch_root = tmp_path / 'results' / 'CRC' / 'S01'
    for name in ('raw.h5ad', 'SVC.h5ad'):
        direct = read_h5ad(tmp_path / 'direct' / name)
        batch = read_h5ad(batch_root / name)
        np.testing.assert_array_equal(direct.X, batch.X)
        np.testing.assert_array_equal(direct.obsm['spatial'], batch.obsm['spatial'])
        pd.testing.assert_frame_equal(direct.obs, batch.obs)
        assert direct.var_names.equals(batch.var_names)
    raw_result = read_h5ad(batch_root / 'raw.h5ad')
    assert raw_result.n_obs == 7
    assert pd.isna(raw_result.obs.loc['s6', 'revise_Level1'])
    assert pd.isna(raw_result.obs.loc['s4', 'revise_Level2'])
    assert read_h5ad(batch_root / 'SVC.h5ad').n_obs == 4
    if reference_config is not None:
        evidence = raw_result.uns['revise_reconstruction']['effective_request']['inputs']['reference_preparation']
        assert evidence['verification'] == 'verified_report'
        assert evidence['config_path'] == str(reference_config)
    assert runner.run_batch(project)['summary']['reused'] == 1
