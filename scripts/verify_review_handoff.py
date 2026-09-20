"""Retain small engineering acceptance artifacts; never runs real reconstruction."""
from __future__ import annotations

import argparse
import base64
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
TESTS = [
    'tests/application/test_sample_delivery.py',
    'tests/backend/test_sc_sr_guidance.py',
    'tests/application/test_ist_ot.py',
    'tests/application/test_ot_delivery_config.py',
    'tests/reference_preparation/test_confidence_shared.py',
    'tests/integration/test_analysis_delivery.py',
    'tests/batch/test_ot_assembly_contract.py',
    'tests/analysis/test_assembly_comparison.py',
    'tests/analysis/test_prepare_assembly_baseline.py',
]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path, help='New evidence directory; existing paths are refused')
    args = parser.parse_args()
    destination = args.output.expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=False)
    os.chdir(ROOT)
    sys.path.insert(0, str(ROOT))
    sys.modules['readline'] = None
    os.environ.update({
        'PYTEST_DISABLE_PLUGIN_AUTOLOAD': '1', 'OMP_NUM_THREADS': '1',
        'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'NUMBA_DISABLE_JIT': '1',
        'MPLBACKEND': 'Agg', 'MPLCONFIGDIR': str(destination / 'mpl'),
        'NUMBA_CACHE_DIR': str(destination / 'numba'), 'PYTHONPATH': str(ROOT),
    })
    summary = {'kind': 'synthetic engineering acceptance', 'python': sys.executable,
               'real_reconstruction': False, 'scientific_acceptance': False,
               'output': str(destination), 'tests': TESTS}
    command = [sys.executable, '-c',
               'import sys; sys.modules["readline"]=None; import pytest; raise SystemExit(pytest.main(sys.argv[1:]))',
               *TESTS, '-q', '--tb=short', '--basetemp', str(destination / 'fixtures'),
               '--junitxml', str(destination / 'junit.xml')]
    summary['command'] = command
    with (destination / 'tests.log').open('w') as stream:
        result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
    summary['test_exit_code'] = result.returncode
    if (destination / 'junit.xml').exists():
        suites = ET.parse(destination / 'junit.xml').getroot().findall('testsuite')
        summary['counts'] = {key: sum(int(s.attrib.get(key, 0)) for s in suites)
                             for key in ['tests', 'failures', 'errors', 'skipped']}
    (destination / 'summary.json').write_text(json.dumps(summary, indent=2))
    if result.returncode or summary.get('counts', {}).get('skipped', 0):
        print(json.dumps(summary, indent=2))
        return result.returncode or 1

    import runpy
    import nbformat
    from nbclient import NotebookClient
    from nbconvert import HTMLExporter

    helpers = runpy.run_path(str(ROOT / 'tests/analysis/test_assembly_comparison.py'))
    inputs_dir = destination / 'notebook-inputs'
    inputs_dir.mkdir()
    methods, baseline = helpers['_write_fixture'](inputs_dir)
    config = {'method_paths': {k: str(v) for k, v in methods.items()},
              'baseline_path': str(baseline), 'cell_types': ['T'],
              'broad_column': 'Level1', 'baseline_subtype_column': 'SVC_cluster',
              'spatial_key': 'spatial', 'coordinate_unit': 'pixel',
              'microns_per_coordinate': 0.2125, 'seed': 23,
              'resolutions': [0.6, 0.7], 'plot_resolution': 0.6,
              'output_dir': str(destination / 'notebook-report')}
    config_path = destination / 'notebook-config.json'
    config_path.write_text(json.dumps(config, indent=2))
    os.environ['REVISE_ASSEMBLY_COMPARISON_CONFIG'] = str(config_path)
    source = ROOT / 'reproduce/case/assembly_comparison.ipynb'
    before = digest(source)
    notebook = nbformat.read(source, as_version=4)
    notebook.cells.insert(0, nbformat.v4.new_code_cell('%matplotlib inline'))
    notebook.cells.append(nbformat.v4.new_code_cell(
        'from pathlib import Path\n'
        f'evidence_output = Path({str(destination)!r})\n'
        "comparison.metrics.to_csv(evidence_output / 'metrics.csv', index=False)\n"
        "comparison.coverage.to_csv(evidence_output / 'coverage.csv', index=False)\n"
        "for broad_type, result in comparison.by_type.items():\n"
        "    for (method, resolution), table in result.contingencies.items():\n"
        "        table.to_csv(evidence_output / f'contingency-{broad_type}-{method}-{resolution}.csv')\n"))
    try:
        executed = NotebookClient(notebook, timeout=180, kernel_name='python3',
                                  resources={'metadata': {'path': str(ROOT)}}).execute()
        nbformat.write(executed, destination / 'comparison.executed.ipynb')
        html, _ = HTMLExporter().from_notebook_node(executed)
        (destination / 'comparison.html').write_text(html)
        images = 0
        for cell in executed.cells:
            for out in cell.get('outputs', []):
                png = out.get('data', {}).get('image/png')
                if png:
                    images += 1
                    (destination / f'comparison-{images}.png').write_bytes(base64.b64decode(png))
        assert images > 0, 'Notebook did not render any PNG figures'
        assert digest(source) == before, 'Source notebook changed'
        summary.update({'notebook': 'passed', 'resolutions': [0.6, 0.7],
                        'source_notebook_unchanged': True, 'png_count': images})
    except Exception as exc:
        summary.update({'notebook': 'failed', 'notebook_error': repr(exc)})
        raise
    finally:
        summary['artifacts'] = {str(p.relative_to(destination)): digest(p)
                                for p in sorted(destination.rglob('*'))
                                if p.is_file() and not any(x in p.relative_to(destination).parts
                                                         for x in ['mpl', 'numba'])
                                and p.name != 'summary.json'}
        (destination / 'summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps({k: v for k, v in summary.items() if k != 'artifacts'}, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
