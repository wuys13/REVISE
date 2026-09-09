"""Sequential reconstruction of standard ST samples, with verified reuse."""
from __future__ import annotations

from collections import Counter
from contextlib import contextmanager
from copy import deepcopy
from hashlib import sha256
from importlib.metadata import PackageNotFoundError, version
import json
import os
from pathlib import Path
import subprocess
import sys
from tempfile import NamedTemporaryFile
import traceback

import yaml

from .config import ResolvedSample, discover_samples, load_batch_config, resolve_sample
from .sample import file_identity, read_sample


DEFAULT_CELL_TYPES = ['T', 'Macro', 'Fibroblast']
ROUTES = {'hST': {'svc_type': 'sp-SVC'},
          'iST': {'svc_type': 'sc-SVC', 'mode': 'cluster'},
          'sST': {'svc_type': 'sc-SVC', 'mode': 'sr'}}


def _json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with NamedTemporaryFile(mode='w', dir=path.parent, delete=False, encoding='utf-8') as handle:
            temporary = Path(handle.name)
            json.dump(value, handle, indent=2, allow_nan=False)
            handle.write('\n')
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text())
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError):
        return {}


def _digest(value) -> str:
    return sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def _code_identity() -> dict:
    """Hash reconstruction runtime files, excluding analysis-only code."""
    package = Path(__file__).resolve().parents[1]
    excluded = {
        Path('batch') / '__init__.py',
        Path('batch') / 'analysis.py',
        Path('batch') / 'cli.py',
        Path('batch') / 'inputs.py',
    }
    files = [path for path in sorted(package.rglob('*.py'))
             if path.relative_to(package) not in excluded]
    files.extend(sorted(package.rglob('*.yaml')))
    files = [path for path in files if 'analysis' not in path.relative_to(package).parts]
    import reconstruct
    files.append(Path(reconstruct.__file__).resolve())
    digest = sha256()
    for path in files:
        digest.update(str(path.relative_to(package.parent)).encode())
        digest.update(path.read_bytes())
    dependencies = {}
    for name in ('numpy', 'scipy', 'anndata', 'scanpy', 'POT', 'tacco', 'igraph', 'leidenalg'):
        try:
            dependencies[name] = version(name)
        except PackageNotFoundError:
            dependencies[name] = None
    return {'sha256': digest.hexdigest(), 'python': sys.version, 'dependencies': dependencies}


def _reconstruction_document(document: dict) -> dict:
    """Return only settings that can affect the existing reconstruction engine."""
    result = deepcopy(document)
    result.pop('analysis', None)
    return result


def _reconstruction_inputs(sample, document: dict) -> dict[str, dict[str, str]]:
    identities = {
        'spatial': file_identity(sample.st_path),
        'reference': file_identity(sample.reference_path),
    }
    prior = document.get('inputs', {}).get('pm_on_cell')
    if prior is not None:
        identities['pm_on_cell'] = file_identity(Path('/') / prior['path'])
    return identities


def _reconstruction_fingerprint(sample, document: dict, identities: dict, code: dict,
                                *, cell_type: str | None = None) -> str:
    return _digest({
        'fingerprint_version': 2,
        'sample': sample.sample_id,
        'inputs': identities,
        'config': _reconstruction_document(document),
        'project_config': str(Path(sample.metadata['configuration']['project_config']).resolve()),
        'cell_type': cell_type,
        'code': code,
    })


def _cell_types(document: dict) -> list[str | None]:
    modality = document['modality']
    if modality not in ROUTES:
        raise ValueError('modality must be hST, iST, or sST')
    if modality != 'iST':
        return [None]
    labels = document.get('local_refinement', {}).get('cell_types', DEFAULT_CELL_TYPES)
    if not isinstance(labels, list) or not labels:
        raise ValueError('local_refinement.cell_types must be a non-empty list')
    reserved = {'inputs', 'analysis', '.', '..', '.revise'}
    for label in labels:
        if (not isinstance(label, str) or not label.strip() or label != label.strip()
                or '/' in label or '\\' in label or any(ord(c) < 32 for c in label)
                or label.casefold() in reserved or label.startswith('.')):
            raise ValueError(f'Unsafe cell type directory: {label!r}; use an explicit label mapping')
    if len({label.casefold() for label in labels}) != len(labels):
        raise ValueError('cell_types contain duplicate or colliding directory names')
    return labels


def application_document(sample, cell_type: str | None, output_root: Path) -> dict:
    """Translate package paths to the unchanged single-run root-path contract."""
    doc = deepcopy(sample.document)
    for field in ('modality', 'coordinates', 'enabled', 'analysis', 'input_root', 'output_root'):
        doc.pop(field, None)
    doc['schema_version'] = 1
    doc['application'] = dict(ROUTES[sample.modality])
    doc['paths'] = {'root_dir': '/'}
    doc['inputs']['st'] = {'path': str(sample.st_path).lstrip('/'), 'format': 'h5ad'}
    doc['inputs']['reference']['path'] = str(sample.reference_path).lstrip('/')
    doc['inputs']['reference']['format'] = 'h5ad'
    prior = doc['inputs'].get('pm_on_cell')
    if prior is not None:
        source = Path(prior['path'])
        if not source.is_absolute():
            source = sample.root / source
        prior['path'] = str(source.resolve()).lstrip('/')
    local = doc['local_refinement']
    local.pop('cell_types', None)
    output = doc.setdefault('output', {})
    if set(output) - {'ist_mapping'}:
        raise ValueError('sample output only accepts ist_mapping; output locations are batch-owned')
    output['dir'] = str(output_root).lstrip('/')
    if sample.modality == 'iST':
        local['select_cell_type'] = cell_type
        output.setdefault('ist_mapping', 'paired')
    else:
        output['name'] = 'SVC'
    return doc


def _execute(config_path: Path, log_path: Path) -> int:
    environment = os.environ.copy()
    package_root = str(Path(__file__).resolve().parents[2])
    environment['PYTHONPATH'] = os.pathsep.join(filter(None, (package_root, environment.get('PYTHONPATH'))))
    with log_path.open('w') as log:
        process = subprocess.run(
            [sys.executable, '-m', 'reconstruct', '--config', str(config_path)],
            stdout=log, stderr=subprocess.STDOUT, env=environment,
        )
    return process.returncode


def _reusable(state: dict, fingerprint: str) -> bool:
    if state.get('status') != 'succeeded' or state.get('fingerprint') != fingerprint:
        return False
    artifacts = state.get('artifacts', [])
    if not artifacts:
        return False
    try:
        return all(file_identity(Path(item['path'])) == item for item in artifacts)
    except (OSError, ValueError, KeyError):
        return False


def _handoff(sample, config, paths: dict, fingerprint: str) -> dict:
    import anndata as ad
    import numpy as np

    outputs = {}
    provenance = {}
    spatial_path = paths.get('spatial', paths.get('svc'))
    with_source = ad.read_h5ad(sample.st_path, backed='r')
    raw_ids = with_source.obs_names.copy()
    raw_coordinates = np.asarray(with_source.obsm['spatial'])
    with_source.file.close()
    pairing = {'status': 'unavailable', 'reason': 'generated_spatial_units' if sample.modality == 'sST' else 'unverified',
               'raw_observations': len(raw_ids), 'id_key': 'obs_names'}
    for role, path in paths.items():
        data = ad.read_h5ad(path, backed='r')
        try:
            record = file_identity(path)
            record.update(shape=list(data.shape), observation_role=(
                'reference_expression' if role == 'expression' else
                'generated_spatial_units' if sample.modality == 'sST' else 'spatial_units'))
            outputs[role] = record
            metadata = data.uns.get('revise_reconstruction', {})
            manifest = metadata.get('run_manifest')
            if manifest:
                provenance[role] = file_identity(Path(manifest))
            if path == spatial_path and sample.modality != 'sST':
                ids = data.obs_names
                indexer = raw_ids.get_indexer(ids)
                coordinates = data.obsm.get('spatial')
                paired = (ids.is_unique and len(ids) > 0 and (indexer >= 0).all()
                          and coordinates is not None
                          and np.asarray(coordinates).shape == raw_coordinates[indexer].shape
                          and np.array_equal(np.asarray(coordinates), raw_coordinates[indexer]))
                pairing.update(status='available' if paired else 'unavailable',
                               reason='same_ids_and_coordinates' if paired else 'ids_or_coordinates_not_preserved',
                               reconstructed_observations=len(ids),
                               observation_ids_sha256=_digest(ids.tolist()))
        finally:
            data.file.close()
    mapping = config.ist_mapping if sample.modality == 'iST' else None
    return {
        'schema_version': 1, 'status': 'succeeded', 'sample_id': sample.sample_id,
        'directory': str(config.output_dir),
        'modality': sample.modality, 'cell_type': config.select_cell_type,
        'ist_mapping': mapping, 'fingerprint': fingerprint,
        'inputs': sample.metadata, 'outputs': outputs, 'provenance': provenance,
        'coordinates': sample.metadata['coordinates'], 'pairing': pairing,
        'expression_semantics': ('reference_expression_by_cluster_' + mapping
                                 if mapping in {'mean', 'random'} else 'native_reconstruction_carriers'),
        'analysis': {'status': 'not_run', 'state_path': 'analysis/analysis.json',
                     'storage': 'paired observations/windows share tables; independent axes use separate files',
                     'current_ist_spatial_carrier_contract': mapping == 'paired' if sample.modality == 'iST' else None},
    }


def _configuration_audit(sample) -> dict:
    """Re-read the full source chain for the mutable audit record."""
    configuration = sample.metadata['configuration']
    resolved = resolve_sample(configuration['project_config'], configuration['sample_dir'])
    audit = deepcopy(configuration)
    audit.update(chain=resolved.config_chain, effective=resolved.document,
                 input_root=str(resolved.input_root), output_root=str(resolved.output_root))
    if (_digest(_reconstruction_document(sample.document))
            != _digest(_reconstruction_document(resolved.document))):
        raise ValueError('Sample configuration changed; rerun reconstruction first')
    return audit


def _task_result(result: dict, status: str, *, state: dict | None = None,
                 handoff: dict | None = None, error: str | None = None) -> dict:
    record = dict(result, status=status)
    if state is not None:
        if state.get('fingerprint') is not None:
            record['fingerprint'] = state['fingerprint']
        if state.get('artifacts'):
            record['artifacts'] = state['artifacts']
    if handoff is not None:
        if handoff.get('outputs'):
            record['outputs'] = handoff['outputs']
        if handoff.get('provenance'):
            record['provenance'] = handoff['provenance']
    if error is not None:
        record['error'] = error
    record['log_path'] = str(Path(result['directory']) / '.revise' / 'reconstruction.log')
    return record


def _run_task(sample, cell_type: str | None, code: dict, output_root: Path) -> dict:
    from revise.application.config import compile_application_config, load_application_yaml
    from revise.application.publication import output_paths

    root = output_root / cell_type if cell_type is not None else output_root
    control = root / '.revise'
    state_path = control / 'task.json'
    handoff_path = root / 'reconstruction.json'
    result = {'sample_id': sample.sample_id, 'cell_type': cell_type, 'directory': str(root)}
    state = dict(result, status='running')
    writable = False
    try:
        for directory in (root, control, root / 'analysis'):
            if directory.is_symlink() or not directory.resolve().is_relative_to(output_root):
                raise ValueError(f'Package output directory must not be a symlink: {directory}')
        document = application_document(sample, cell_type, output_root)
        audit = _configuration_audit(sample)
        if (Path(audit['output_root']) / sample.sample_id != output_root
                or Path(audit['input_root']) / sample.sample_id != sample.root):
            raise ValueError('Input or output root changed during execution')
        sample.metadata = deepcopy(sample.metadata)
        sample.metadata['configuration'] = audit
        source_paths = [Path(item['path']) for item in sample.metadata['sources'].values()]
        source_paths.extend((sample.st_path, sample.reference_path))
        prior = document['inputs'].get('pm_on_cell')
        if prior is not None:
            source_paths.append(Path('/') / prior['path'])
        targets = [root / name for name in ('SVC.h5ad', 'spatial.h5ad', 'expr.h5ad', 'reconstruction.json')]
        for original in source_paths:
            if original.is_relative_to(control):
                raise ValueError(f'Original input is inside task control directory: {original}')
            for target in targets:
                if (target.resolve() == original.resolve()
                        or (original.is_dir() and target.resolve().is_relative_to(original.resolve()))
                        or (target.exists() and original.exists() and os.path.samefile(target, original))):
                    raise ValueError(f'Publication destination aliases original input: {original}')
        control.mkdir(parents=True, exist_ok=True)
        writable = True
        config_path = control / 'application.yaml'
        config_path.write_text(yaml.safe_dump(document, sort_keys=False))
        source, loaded = load_application_yaml(config_path)
        config = compile_application_config(loaded, source=source)
        identities = _reconstruction_inputs(sample, document)
        fingerprint = _reconstruction_fingerprint(
            sample, sample.document, identities, code, cell_type=cell_type)
        previous = _read_json(state_path)
        if _reusable(previous, fingerprint):
            previous['configuration_audit'] = audit
            _json(state_path, previous)
            handoff = _read_json(handoff_path)
            return _task_result(result, 'reused', state=previous, handoff=handoff)
        state.update(fingerprint=fingerprint, code=code, inputs=identities,
                     configuration_audit=audit)
        _json(state_path, state)
        _json(root / 'analysis' / 'analysis.json', {'status': 'not_run', 'reconstruction_fingerprint': fingerprint})
        _json(handoff_path, dict(result, schema_version=1, status='running'))
        returncode = _execute(config_path, control / 'reconstruction.log')
        if returncode:
            raise RuntimeError(f'Reconstruction exited {returncode}; see {control / "reconstruction.log"}')
        paths = output_paths(config)
        handoff = _handoff(sample, config, paths, fingerprint)
        _json(handoff_path, handoff)
        artifacts = [file_identity(path) for path in paths.values()]
        artifacts.extend(handoff['provenance'].values())
        artifacts.append(file_identity(handoff_path))
        if _code_identity() != code:
            raise ValueError('Reconstruction code changed during execution; retry with stable code')
        verify_configuration(sample.metadata['configuration'])
        if any(file_identity(Path(item['path'])) != item for item in identities.values()):
            raise ValueError('Input changed during reconstruction; retry with stable inputs')
        _json(state_path, dict(state, status='succeeded', artifacts=artifacts))
        return _task_result(result, 'succeeded', state=dict(state, status='succeeded', artifacts=artifacts),
                            handoff=handoff)
    except Exception as exc:
        error = f'{type(exc).__name__}: {exc}'
        if writable:
            with (control / 'reconstruction.log').open('a') as log:
                traceback.print_exc(file=log)
            _json(state_path, dict(state, status='failed', error=error))
            _json(handoff_path, dict(result, schema_version=1, status='failed', error=error))
        else:
            _invalidate_handoffs(output_root, status='failed', only=root)
        return _task_result(result, 'failed', state=state, error=error)


def _invalidate_handoffs(root: Path, *, status: str, active: set[Path] | None = None,
                         only: Path | None = None) -> None:
    """Invalidate tracked results without deleting scientific data or foreign files."""
    for handoff in [root / 'reconstruction.json', *root.glob('*/reconstruction.json')]:
        if only is not None and handoff.parent != only:
            continue
        if active is not None and handoff.parent in active:
            continue
        if handoff.is_symlink() or handoff.parent.is_symlink():
            continue
        state_path = handoff.parent / '.revise' / 'task.json'
        if state_path.parent.is_symlink():
            continue
        record = _read_json(handoff)
        state = _read_json(state_path)
        if record.get('schema_version') == 1 and state.get('directory') == str(handoff.parent):
            record.update(status=status, reason='sample_validation_failed' if status == 'failed' else 'task_no_longer_requested')
            state['status'] = status
            _json(handoff, record)
            _json(state_path, state)


@contextmanager
def _batch_lock(root: Path):
    import fcntl
    with (root / '.revise-batch.lock').open('a') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError('Another batch is already using this output_root') from exc
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def verify_configuration(configuration: dict):
    """Re-resolve the complete hierarchy, including previously absent YAML files."""
    resolved = resolve_sample(configuration['project_config'], configuration['sample_dir'])
    for key in ('input_root', 'output_root'):
        if key in configuration and configuration[key] != str(getattr(resolved, key)):
            raise ValueError('Input or output root changed; rerun reconstruction first')
    effective = configuration.get('effective')
    if (not isinstance(effective, dict)
            or _digest(_reconstruction_document(resolved.document))
            != _digest(_reconstruction_document(effective))
            or not resolved.document.get('enabled', True)):
        raise ValueError('Sample configuration changed; rerun reconstruction first')
    return resolved


def _output_directory(output: Path, relative: Path) -> Path:
    directory = output
    for part in relative.parts:
        directory = directory / part
        if directory.is_symlink() or (directory.exists() and not directory.is_dir()):
            raise ValueError(f'Output directory must be an ordinary directory: {directory}')
    return directory


def _validate_task_selection(resolved: ResolvedSample, cell_type: str | None) -> None:
    labels = _cell_types(resolved.document)
    if resolved.document.get('modality') == 'iST':
        if cell_type is None or cell_type not in labels:
            raise ValueError(f'cell_type must be one of the configured iST types: {labels}')
    elif cell_type is not None:
        raise ValueError('cell_type is only valid for iST samples')


def _selected_sample(input_root: Path, sample_id: str) -> Path:
    if (not isinstance(sample_id, str) or not sample_id.strip() or sample_id != sample_id.strip()
            or '\\' in sample_id or Path(sample_id).is_absolute()):
        raise ValueError('sample_id must be a relative path under input_root')
    relative = Path(sample_id)
    if not relative.parts or any(part in {'', '.', '..'} for part in relative.parts):
        raise ValueError('sample_id must be a normalized relative path under input_root')
    candidate = input_root / relative
    resolved = candidate.resolve()
    if not resolved.is_relative_to(input_root) or resolved.relative_to(input_root).as_posix() != sample_id:
        raise ValueError('sample_id must be a relative path under input_root')
    current = input_root
    for part in relative.parts:
        current /= part
        if current.is_symlink():
            raise ValueError(f'Selected sample must not contain symlinked directories: {candidate}')
    if not candidate.is_dir() or not (candidate / 'spatial.h5ad').is_file():
        raise ValueError(f'Unknown sample_id: {sample_id}')
    samples = discover_samples(input_root)
    if any(other != candidate and (other in candidate.parents or candidate in other.parents)
           for other in samples):
        raise ValueError('A selected sample cannot contain or be contained by another sample')
    return candidate


def _task_root(output: Path, sample_id: str, cell_type: str | None) -> Path:
    destination = _output_directory(output, Path(sample_id))
    return destination / cell_type if cell_type is not None else destination


def run_reconstruction_task(config_path: str | Path, sample_id: str, *, cell_type: str | None = None) -> dict:
    """Run exactly one standard ST reconstruction task under the shared lock."""
    _, input_root, output = load_batch_config(config_path)
    config_path = Path(config_path).resolve()
    sample_dir = _selected_sample(input_root, sample_id)
    resolved = resolve_sample(config_path, sample_dir)
    _validate_task_selection(resolved, cell_type)
    destination = output / Path(sample_id)
    task_root = _task_root(output, sample_id, cell_type)
    output.mkdir(parents=True, exist_ok=True)
    with _batch_lock(output):
        if not resolved.document.get('enabled', True):
            _invalidate_handoffs(destination, status='inactive', only=task_root)
            return {'sample_id': sample_id, 'cell_type': cell_type,
                    'directory': str(task_root), 'status': 'inactive'}
        try:
            sample = read_sample(resolved)
        except Exception as exc:
            _invalidate_handoffs(destination, status='failed', only=task_root)
            return _task_result(
                {'sample_id': sample_id, 'cell_type': cell_type, 'directory': str(task_root)},
                'failed', error=f'{type(exc).__name__}: {exc}')
        return _run_task(sample, cell_type, _code_identity(), destination)


def _artifacts_belong_to(root: Path, state: dict) -> bool:
    for item in state.get('artifacts', []):
        try:
            path = Path(item['path'])
            if path.is_symlink() or not path.resolve().is_relative_to(root.resolve()):
                return False
        except (KeyError, OSError, ValueError):
            return False
    return True


def verify_reconstruction_task(resolved: ResolvedSample, cell_type: str | None = None) -> dict:
    """Return a current handoff for one task without acquiring the output lock."""
    if not isinstance(resolved, ResolvedSample):
        raise TypeError('resolved must be a ResolvedSample')
    _validate_task_selection(resolved, cell_type)
    if not resolved.document.get('enabled', True):
        raise ValueError('Sample is disabled')
    sample = read_sample(resolved)
    destination = _output_directory(resolved.output_root, Path(resolved.sample_id))
    root = _task_root(resolved.output_root, resolved.sample_id, cell_type)
    control = root / '.revise'
    for directory in (root, control, root / 'analysis'):
        if directory.is_symlink() or not directory.resolve().is_relative_to(destination.resolve()):
            raise ValueError(f'Package output directory must not be a symlink: {directory}')
    state = _read_json(control / 'task.json')
    if state.get('directory') != str(root) or not _artifacts_belong_to(root, state):
        raise ValueError('Reconstruction task ownership or artifacts are invalid')
    document = application_document(sample, cell_type, destination)
    identities = _reconstruction_inputs(sample, document)
    fingerprint = _reconstruction_fingerprint(
        sample, sample.document, identities, _code_identity(), cell_type=cell_type)
    if not _reusable(state, fingerprint) or state.get('inputs') != identities:
        raise ValueError('Reconstruction is not successful or its artifacts changed')
    handoff = _read_json(root / 'reconstruction.json')
    if (handoff.get('status') != 'succeeded' or handoff.get('fingerprint') != fingerprint
            or handoff.get('directory') != str(root)
            or handoff.get('sample_id') != resolved.sample_id
            or handoff.get('cell_type') != cell_type):
        raise ValueError('Reconstruction handoff is not current')
    configuration = handoff.get('inputs', {}).get('configuration')
    if not isinstance(configuration, dict):
        raise ValueError('Reconstruction handoff is missing configuration provenance')
    if (Path(configuration.get('project_config', '')).resolve() != resolved.project_config
            or Path(configuration.get('sample_dir', '')).resolve() != resolved.root):
        raise ValueError('Reconstruction handoff belongs to another project or sample')
    verify_configuration(configuration)
    if handoff.get('inputs', {}).get('sources') != sample.metadata.get('sources'):
        raise ValueError('Reconstruction input provenance is not current')
    return handoff


def run_batch(config_path: str | Path) -> dict:
    """Read standard inputs and reconstruct into a separate output tree."""
    _, root, output = load_batch_config(config_path)
    config_path = Path(config_path).resolve()
    output.mkdir(parents=True, exist_ok=True)
    with _batch_lock(output):
        return _run_samples(config_path, root, output)


def _run_samples(config_path: Path, root: Path, output: Path) -> dict:
    paths = discover_samples(root)
    report_path = output / 'batch_status.json'
    sample_ids = {path.relative_to(root).as_posix() for path in paths}
    # Removed samples must not retain a current result in the previous inventory.
    for previous in _read_json(report_path).get('tasks', []):
        if previous.get('sample_id') in sample_ids or not previous.get('directory'):
            continue
        directory = Path(previous['directory'])
        try:
            _output_directory(output, directory.relative_to(output))
        except ValueError:
            continue
        _invalidate_handoffs(directory, status='inactive', only=directory)
    if not paths:
        raise ValueError(f'No spatial.h5ad samples under {root}; migrate old sample.yaml packages to schema_version: 2')
    report = {'schema_version': 1, 'status': 'running', 'project_config': str(config_path), 'input_root': str(root),
              'output_root': str(output), 'tasks': [], 'summary': {}}
    _json(report_path, report)
    code = _code_identity()
    for path in paths:
        sample_id = path.relative_to(root).as_posix()
        destination = output / path.relative_to(root)
        sample_control = destination / '.revise'
        safe_destination = False
        try:
            _output_directory(output, path.relative_to(root))
            safe_destination = path != root
            if sample_control.is_symlink() or (sample_control.exists() and not sample_control.is_dir()):
                raise ValueError(f'Sample control directory must be an ordinary directory: {sample_control}')
            if path == root or any(other != path and (other in path.parents or path in other.parents) for other in paths):
                raise ValueError('A sample must be below input_root and cannot contain another sample')
            resolved = resolve_sample(config_path, path)
            if not resolved.document.get('enabled', True):
                _invalidate_handoffs(destination, status='inactive')
                if sample_control.exists():
                    _json(sample_control / 'sample.json', {'sample_id': sample_id, 'status': 'inactive'})
                continue
            labels = _cell_types(resolved.document)
            sample = read_sample(resolved)
        except Exception as exc:
            task = {'sample_id': sample_id, 'sample_dir': str(path), 'status': 'failed', 'error': f'{type(exc).__name__}: {exc}'}
            report['tasks'].append(task)
            if safe_destination:
                _invalidate_handoffs(destination, status='failed')
            if safe_destination and not sample_control.is_symlink() and (not sample_control.exists() or sample_control.is_dir()):
                _json(sample_control / 'sample.json', task)
        else:
            active = {destination / label if label is not None else destination for label in labels}
            _invalidate_handoffs(destination, status='inactive', active=active)
            _json(sample_control / 'sample.json', {'sample_id': sample_id, 'status': 'validated'})
            for label in labels:
                report['tasks'].append(_run_task(sample, label, code, destination))
                _json(report_path, report)
        _json(report_path, report)
    summary = Counter(task['status'] for task in report['tasks'])
    report['summary'] = {status: summary[status] for status in ('succeeded', 'failed', 'reused')}
    report['status'] = 'completed_with_failures' if summary['failed'] else 'completed'
    _json(report_path, report)
    return report
