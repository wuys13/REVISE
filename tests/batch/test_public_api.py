"""The package exposes the same task executors used by batch orchestration."""
import inspect

import pytest


@pytest.mark.parametrize('name,arguments', [
    ('run_reconstruction_task', ['config_path', 'sample_id', 'cell_type']),
    ('run_analysis_task', ['config_path', 'sample_id', 'aspect', 'cell_type']),
])
def test_public_task_signature(name, arguments):
    import revise.batch as batch

    assert hasattr(batch, name), f'Missing public task API: {name}'
    signature = inspect.signature(getattr(batch, name))
    assert list(signature.parameters) == arguments
    assert signature.parameters['cell_type'].kind is inspect.Parameter.KEYWORD_ONLY
    assert signature.parameters['cell_type'].default is None
