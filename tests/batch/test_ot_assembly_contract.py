"""Batch reuse must distinguish assembly modes, options and source declarations."""
import yaml

from revise.batch.runner import run_batch
from tests.batch.test_runner import batch_config, fake_solver, make_sample


def test_mode_options_and_expression_declarations_invalidate_reuse(tmp_path, fake_solver):
    root = tmp_path / 'data' / 'sample'
    doc = make_sample(root, cell_types=None, ist_mapping='random')
    project = batch_config(tmp_path)
    assert run_batch(project)['summary']['succeeded'] == 1
    assert run_batch(project)['summary']['reused'] == 1
    for output in (
        {'ist_mapping': 'within_cluster'},
        {'ist_mapping': 'within_cluster', 'ist_ot': {'spatial_weight': .4}},
        {'ist_mapping': 'outside_cluster'},
    ):
        doc['output'] = output
        (root / 'batch.yaml').write_text(yaml.safe_dump(doc))
        assert run_batch(project)['summary']['succeeded'] == 1
        assert run_batch(project)['summary']['reused'] == 1
    doc['inputs']['st'] = {'expression': {'identity': 'measured_expression', 'scale': 'untransformed_nonnegative'}}
    doc['inputs']['reference']['expression'] = {'identity': 'measured_reference', 'scale': 'untransformed_nonnegative'}
    (root / 'batch.yaml').write_text(yaml.safe_dump(doc))
    assert run_batch(project)['summary']['succeeded'] == 1
    assert len(fake_solver) == 5
