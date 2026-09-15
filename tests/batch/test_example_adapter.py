"""The copyable adapter example returns a reloadable, truthful inventory."""
import csv
import importlib.util
from pathlib import Path
from types import SimpleNamespace


def test_example_adapter_records_actual_carriers(tmp_path):
    path = Path(__file__).resolve().parents[2] / 'configs/batch/analysis_example.py'
    assert path.is_file(), 'Missing copyable analysis adapter example'
    spec = importlib.util.spec_from_file_location('batch_analysis_example', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    context = SimpleNamespace(
        output_dir=tmp_path,
        parameters={},
        reconstruction={'outputs': {
            'spatial': {'shape': [2, 3], 'path': '/input/spatial.h5ad',
                        'observation_role': 'spatial_units'},
            'expression': {'shape': [7, 5], 'path': '/input/expr.h5ad',
                           'observation_role': 'reference_expression'},
        }},
    )
    result = module.run(context)
    artifact = result['artifacts']['carrier_inventory']
    with (tmp_path / artifact['path']).open() as stream:
        rows = list(csv.DictReader(stream))
    assert rows[1]['observation_role'] == 'reference_expression'
    assert rows[1]['observations'] == '7'
    assert rows[0]['source_path'] == '/input/spatial.h5ad'
    assert result['calculation']['input_view'] == 'native_carriers'
    assert result['calculation']['comparison_basis'] == 'inventory_only_no_scientific_comparison'
