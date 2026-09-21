"""Execute an unchanged source Notebook against an explicit acceptance project.

Run in the chosen repository's scientific Python environment. Output contains
the executed Notebook, HTML export and source-hash/execution record; scientific
artifacts remain owned by the Notebook's existing workflow.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

sys.modules['readline'] = None
import nbformat
from nbclient import NotebookClient
from nbconvert import HTMLExporter

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--repo', type=Path, required=True)
parser.add_argument('--source', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--project', type=Path)
parser.add_argument('--comparison-config', type=Path)
args = parser.parse_args()
source = (args.repo / args.source).resolve()
output = args.output.resolve()
if (output / (source.stem + '.executed.ipynb')).exists():
    raise FileExistsError(f'Choose a new output directory: {output}')
output.mkdir(parents=True, exist_ok=True)
source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
os.environ['OUTPUT_DIR'] = str(output)
if args.project:
    os.environ['PROJECT_YAML'] = str(args.project.resolve())
if args.comparison_config:
    os.environ['REVISE_ASSEMBLY_COMPARISON_CONFIG'] = str(args.comparison_config.resolve())
nb = nbformat.read(source, as_version=4)
nb.cells.insert(0, nbformat.v4.new_code_cell(
    '# Execution-only environment initialization; scientific source cells follow unchanged.\n'
    'import sys\nsys.modules["readline"] = None\n%matplotlib inline'
))
record = {'source': str(source), 'source_sha256': source_hash,
          'kind': 'acceptance notebook execution', 'kernel': 'python3',
          'scientific_acceptance': False}
try:
    NotebookClient(nb, timeout=1800, kernel_name='python3',
                   resources={'metadata': {'path': str(args.repo.resolve())}}).execute()
    record['status'] = 'passed'
except BaseException as exc:
    record.update(status='failed', error=repr(exc))
    raise
finally:
    executed = output / (source.stem + '.executed.ipynb')
    nbformat.write(nb, executed)
    html, _ = HTMLExporter().from_notebook_node(nb)
    (output / (source.stem + '.html')).write_text(html)
    record['source_unchanged'] = source_hash == hashlib.sha256(source.read_bytes()).hexdigest()
    record['code_cells'] = sum(cell.cell_type == 'code' for cell in nb.cells) - 1
    record['cell_errors'] = [out.get('ename') for cell in nb.cells for out in cell.get('outputs', []) if out.output_type == 'error']
    (output / 'notebook_execution.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(record, indent=2))
