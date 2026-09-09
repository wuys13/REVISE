"""Execute configured analysis aspects against verified reconstruction handoffs."""
from __future__ import annotations

from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import dataclass, field
from importlib import import_module
import inspect
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import traceback
from typing import Any

from .config import discover_samples, load_batch_config, resolve_sample
from .runner import (
    _batch_lock, _cell_types, _digest, _json, _output_directory, _read_json,
    _selected_sample, _task_root, _validate_task_selection, verify_reconstruction_task,
)
from .sample import file_identity


_ANALYSIS_FIELDS = {
    'entrypoint', 'version', 'parameters', 'resources', 'enabled', 'requires_pairing',
}
_RESERVED_ASPECT_NAMES = {'inputs', 'analysis', '.', '..', '.revise'}
_UNSET = object()
_ADAPTER_SOURCES: dict[str, dict] = {}


def _safe_aspect_name(name: Any) -> bool:
    """Return whether *name* is safe to use as one output directory."""
    return (
        isinstance(name, str)
        and bool(name.strip())
        and name == name.strip()
        and name.casefold() not in _RESERVED_ASPECT_NAMES
        and not name.startswith('.')
        and '/' not in name
        and '\\' not in name
        and not any(ord(char) < 32 for char in name)
    )


@dataclass(frozen=True)
class AnalysisContext:
    """Adapter context with lazy input views and resolved resource paths.

    Existing adapters can keep using ``reconstruction``, ``output_dir`` and
    ``parameters``.  Input views are materialized only when ``inputs`` is read;
    ``resources`` contains paths resolved by the batch configuration loader.
    """

    reconstruction: dict
    output_dir: Path
    parameters: dict
    resources: dict[str, Path] = field(default_factory=dict)
    _inputs_value: Any = field(default=_UNSET, init=False, repr=False, compare=False)

    @property
    def inputs(self) -> Any:
        """Load ``AnalysisInputs`` only when a module requests it."""
        value = object.__getattribute__(self, '_inputs_value')
        if value is _UNSET:
            from .inputs import AnalysisInputs
            value = AnalysisInputs.from_reconstruction(self.reconstruction)
            object.__setattr__(self, '_inputs_value', value)
        return value

    def close(self) -> None:
        """Release any input views materialized by an adapter."""
        value = object.__getattribute__(self, '_inputs_value')
        if value is not _UNSET:
            value.close()


def _specifications(batch: dict) -> dict:
    """Validate and return only explicitly configured analysis aspects."""
    specs = batch.get('analysis', {})
    if specs is None:
        specs = {}
    if not isinstance(specs, dict):
        raise ValueError('analysis must be a mapping of safe aspect names to settings')
    seen_names = {}
    for name, spec in specs.items():
        if not _safe_aspect_name(name):
            raise ValueError(f'Unsafe analysis aspect directory: {name!r}')
        folded_name = name.casefold()
        if folded_name in seen_names:
            raise ValueError(
                f'Analysis aspect names collide case-insensitively: '
                f'{seen_names[folded_name]!r} and {name!r}'
            )
        seen_names[folded_name] = name
        if not isinstance(spec, dict) or set(spec) - _ANALYSIS_FIELDS:
            raise ValueError(f'Invalid analysis settings for {name}')
        entrypoint = spec.get('entrypoint')
        if entrypoint is not None:
            if (not isinstance(entrypoint, str) or entrypoint.count(':') != 1
                    or any(not part.strip() for part in entrypoint.split(':', 1))):
                raise ValueError(f'{name}.entrypoint must be module:function')
            if not isinstance(spec.get('version'), str) or not spec['version'].strip():
                raise ValueError(f'{name}.version must identify the algorithm and its dependencies')
        elif 'version' in spec:
            raise ValueError(f'{name}.version requires an entrypoint')
        if not isinstance(spec.get('parameters', {}), dict):
            raise ValueError(f'{name}.parameters must be a mapping')
        if type(spec.get('enabled', True)) is not bool:
            raise ValueError(f'{name}.enabled must be true or false')
        if type(spec.get('requires_pairing', False)) is not bool:
            raise ValueError(f'{name}.requires_pairing must be true or false')
        resources = spec.get('resources', {})
        if not isinstance(resources, dict):
            raise ValueError(f'{name}.resources must be a mapping')
        for resource_name, resource in resources.items():
            if not isinstance(resource_name, str) or not resource_name.strip():
                raise ValueError(f'{name}.resources names must be nonempty strings')
            if not isinstance(resource, str) or not resource.strip():
                raise ValueError(f'{name}.resources.{resource_name} must identify a path')
    return specs


def _resource_paths(spec: dict) -> tuple[dict[str, Path], dict[str, dict]]:
    paths: dict[str, Path] = {}
    identities: dict[str, dict] = {}
    for name, value in spec.get('resources', {}).items():
        path = Path(value)
        path = path.resolve()
        if not path.is_file():
            raise FileNotFoundError(f'Analysis resource does not exist: {path}')
        paths[name] = path
        identities[name] = file_identity(path)
    return paths, identities


def _analysis_code_identity(source: Path) -> dict:
    files = [
        Path(__file__).resolve(), source.resolve(),
        Path(__file__).with_name('runner.py').resolve(),
    ]
    input_helper = Path(__file__).with_name('inputs.py')
    if input_helper.is_file():
        files.append(input_helper.resolve())
    assembly_helper = Path(__file__).parents[1] / 'application' / 'ist_assembly.py'
    if assembly_helper.is_file():
        files.append(assembly_helper.resolve())
    unique = []
    for path in files:
        if path not in unique:
            unique.append(path)
    return {str(path): file_identity(path) for path in unique}


def _analysis_fingerprint(record: dict, spec: dict, resources: dict, source: Path) -> str:
    reconstruction = {
        'fingerprint': record.get('fingerprint'),
        'sample_id': record.get('sample_id'),
        'cell_type': record.get('cell_type'),
        'modality': record.get('modality'),
        'inputs': record.get('inputs', {}).get('sources', {}),
        'outputs': record.get('outputs', {}),
        'provenance': record.get('provenance', {}),
        'pairing': record.get('pairing', {}),
        'coordinates': record.get('coordinates', {}),
        'expression_semantics': record.get('expression_semantics'),
    }
    return _digest({
        'schema_version': 2,
        'aspect': spec,
        'resources': resources,
        'reconstruction': reconstruction,
        'code': _analysis_code_identity(source),
    })


def _artifact_records(previous: dict) -> list[dict]:
    # The structured result schema deliberately does not reinterpret the old
    # role-to-path state as a current success.  A prior run without
    # ``artifact_records`` must execute again under the current contract.
    records = previous.get('artifact_records', [])
    if isinstance(records, dict):
        records = list(records.values())
    return records if isinstance(records, list) else []


def _analysis_reusable(previous: dict, fingerprint: str, destination: Path) -> bool:
    if (previous.get('status') != 'succeeded'
            or previous.get('fingerprint') != fingerprint
            or previous.get('published_directory') != str(destination)):
        return False
    records = _artifact_records(previous)
    if not records:
        return False
    try:
        destination = destination.resolve()
        for item in records:
            path = Path(item['path'])
            if (not path.is_absolute() or path.is_symlink()
                    or not path.resolve().is_relative_to(destination)):
                return False
            expected = {'path': str(path), 'sha256': item['sha256']}
            if file_identity(path) != expected:
                return False
        return True
    except (OSError, ValueError, KeyError, TypeError):
        return False


def _has_unowned_files(destination: Path, previous: dict) -> bool:
    """Refuse to replace files that the prior analysis did not publish."""
    if not destination.exists():
        return False
    try:
        expected = {
            Path(item['path']).resolve()
            for item in _artifact_records(previous)
        }
        for path in destination.rglob('*'):
            if path.is_symlink() or (path.is_file() and path.resolve() not in expected):
                return True
    except (OSError, KeyError, TypeError, ValueError):
        return True
    return False


def _normalise_adapter_result(payload: Any) -> tuple[dict, dict]:
    if not isinstance(payload, dict) or not payload or 'artifacts' not in payload:
        raise ValueError('Analysis adapter must return structured artifacts and calculation')
    declared = payload['artifacts']
    calculation = payload.get('calculation')
    if not isinstance(declared, dict) or not declared:
        raise ValueError('Analysis result artifacts must be a nonempty mapping')
    if not isinstance(calculation, dict):
        raise ValueError('Structured analysis results require a calculation mapping')
    if not isinstance(calculation.get('input_view'), str) or not calculation['input_view'].strip():
        raise ValueError('Analysis calculation.input_view must be a nonempty string')
    if not isinstance(calculation.get('parameters'), dict):
        raise ValueError('Analysis calculation.parameters must be a mapping')
    if not isinstance(calculation.get('comparison_basis'), str) or not calculation['comparison_basis'].strip():
        raise ValueError('Analysis calculation.comparison_basis must be a nonempty string')
    normalised = {}
    for role, value in declared.items():
        if not isinstance(role, str) or not role.strip():
            raise ValueError('Analysis artifact roles must be nonempty strings')
        if not isinstance(value, dict):
            raise ValueError('Analysis artifacts must map roles to path and description mappings')
        path, metadata = value.get('path'), dict(value)
        if not isinstance(path, str) or not path.strip():
            raise ValueError(f'Analysis artifact {role!r} must declare a path')
        description = metadata.get('description')
        if not isinstance(description, str) or not description.strip():
            raise ValueError(f'Analysis artifact {role!r}.description must be a nonempty string')
        metadata['path'] = path
        normalised[role] = metadata
    return normalised, calculation


def _stage_artifacts(stage: Path, destination: Path, declared: dict) -> tuple[list[dict], dict]:
    records = []
    roles = {}
    for role, metadata in declared.items():
        relative = Path(metadata['path'])
        path = stage / relative
        if (relative.is_absolute() or not path.resolve().is_relative_to(stage.resolve())
                or not path.is_file() or path.is_symlink()):
            raise ValueError(f'Analysis artifact must be a file inside its aspect: {relative}')
        published = destination / relative
        identity = file_identity(path)
        identity['path'] = str(published)
        identity['role'] = role
        if 'description' in metadata:
            identity['description'] = metadata['description']
        records.append(identity)
        roles[role] = metadata['path']
    if any(path.is_symlink() for path in stage.rglob('*')):
        raise ValueError('Analysis artifacts must not contain symlinks')
    return records, roles


def _publish(stage: Path, destination: Path, temporary: Path) -> Path | None:
    backup = temporary / 'previous'
    if destination.exists():
        os.replace(destination, backup)
    try:
        os.replace(stage, destination)
    except BaseException:
        if backup.exists():
            os.replace(backup, destination)
        raise
    return backup if backup.exists() else None


def _restore_publication(destination: Path, backup: Path | None, temporary: Path) -> None:
    if destination.exists():
        os.replace(destination, temporary / 'unpublished')
    if backup is not None and backup.exists():
        os.replace(backup, destination)


def _configuration_audit(resolved) -> dict:
    """Record the current source chain without making it part of reuse identity."""
    return {
        'project_config': str(resolved.project_config),
        'sample_dir': str(resolved.root),
        'input_root': str(resolved.input_root),
        'output_root': str(resolved.output_root),
        'chain': resolved.config_chain,
        'effective': resolved.document,
    }


def _run_aspect(root: Path, record: dict, name: str, spec: dict, resolved,
                cell_type: str | None, inactive_reason: str | None = None) -> dict:
    result = {
        'aspect': name,
        'directory': str(root),
        'reconstruction_fingerprint': record.get('fingerprint'),
        'configuration_audit': _configuration_audit(resolved),
    }
    control = root / '.revise' / 'analysis'
    destination = root / 'analysis' / name
    state_path = control / f'{name}.json'
    result['log_path'] = str(control / f'{name}.log')
    writable = False
    try:
        _output_directory(root, Path('.revise/analysis'))
        _output_directory(root, Path('analysis') / name)
        control.mkdir(parents=True, exist_ok=True)
        writable = True
        previous = _read_json(state_path)
        if previous.get('published_directory') == str(destination):
            result['published_directory'] = str(destination)
        for key in ('artifact_records', 'roles'):
            if key in previous:
                result[key] = previous[key]
        if inactive_reason is not None:
            result.update(status='inactive', reason=inactive_reason, specification=spec)
        elif spec.get('enabled', True) is False:
            result.update(status='inactive', reason='disabled', specification=spec)
        elif 'entrypoint' not in spec:
            result.update(status='not_implemented', reason='No analysis adapter configured', specification=spec)
        elif spec.get('requires_pairing', False) and record.get('pairing', {}).get('status') != 'available':
            result.update(status='unavailable', reason=record.get('pairing', {}).get('reason', 'pairing unavailable'),
                          specification=spec)
        else:
            module_name, function_name = spec['entrypoint'].split(':')
            module = import_module(module_name)
            adapter = getattr(module, function_name)
            source = inspect.getsourcefile(adapter)
            if not callable(adapter) or source is None:
                raise ValueError('Analysis adapter must be a Python function with identifiable source')
            source_path = Path(source).resolve()
            source_identity = file_identity(source_path)
            previous_source = _ADAPTER_SOURCES.get(str(source_path))
            if previous_source is not None and previous_source != source_identity:
                raise RuntimeError('Analysis adapter source changed; restart the interpreter before retrying')
            _ADAPTER_SOURCES[str(source_path)] = source_identity
            resources, resource_identities = _resource_paths(spec)
            fingerprint = _analysis_fingerprint(record, spec, resource_identities, source_path)
            result.update(fingerprint=fingerprint, specification=spec, resources=resource_identities)
            if _analysis_reusable(previous, fingerprint, destination):
                reused = dict(previous, status='reused')
                reused['configuration_audit'] = _configuration_audit(resolved)
                reused['log_path'] = str(control / f'{name}.log')
                persisted = dict(reused, status='succeeded')
                _json(state_path, persisted)
                return reused
            if destination.exists() and any(destination.iterdir()):
                if previous.get('published_directory') != str(destination):
                    raise ValueError(f'Refusing to replace unowned analysis directory: {destination}')
                if _has_unowned_files(destination, previous):
                    raise ValueError(f'Refusing to replace analysis directory containing unowned files: {destination}')
            _json(state_path, dict(result, status='running'))
            destination.parent.mkdir(parents=True, exist_ok=True)
            with TemporaryDirectory(prefix=f'.{name}-', dir=destination.parent) as temporary_name:
                temporary = Path(temporary_name)
                stage = temporary / 'artifacts'
                stage.mkdir()
                context = AnalysisContext(
                    record, stage, spec.get('parameters', {}),
                    resources,
                )
                with (control / f'{name}.log').open('w') as log, redirect_stdout(log), redirect_stderr(log):
                    try:
                        payload = adapter(context)
                    finally:
                        context.close()
                declared, calculation = _normalise_adapter_result(payload)
                records, roles = _stage_artifacts(stage, destination, declared)
                current = resolve_sample(resolved.project_config, resolved.root)
                current_spec = _specifications(current.document).get(name)
                if current_spec != spec:
                    raise ValueError('Analysis configuration changed during analysis')
                _, current_resources = _resource_paths(current_spec)
                if _analysis_fingerprint(record, current_spec, current_resources, source_path) != fingerprint:
                    raise ValueError('Analysis resources or code changed during analysis')
                if verify_reconstruction_task(current, cell_type) != record:
                    raise ValueError('Reconstruction changed during analysis')
                current_audit = _configuration_audit(current)
                candidate = dict(
                    result, status='succeeded', artifacts=declared,
                    artifact_records=records, roles=roles,
                    published_directory=str(destination), calculation=calculation,
                    configuration_audit=current_audit,
                )
                json.dumps(candidate, sort_keys=True, allow_nan=False)
                backup = _publish(stage, destination, temporary)
                try:
                    _json(state_path, candidate)
                except BaseException:
                    _restore_publication(destination, backup, temporary)
                    raise
                result = candidate
                return result
        _json(state_path, result)
    except (Exception, SystemExit) as exc:
        result.update(status='failed', error=f'{type(exc).__name__}: {exc}')
        if writable:
            with (control / f'{name}.log').open('a') as log:
                traceback.print_exc(file=log)
            _json(state_path, result)
    return result


def _task_result(result: dict, sample_id: str, cell_type: str | None) -> dict:
    return dict(result, sample_id=sample_id, cell_type=cell_type)


def _update_single_summary(root: Path, record: dict, result: dict) -> None:
    """Record one aspect without claiming that sibling aspects ran."""
    try:
        _output_directory(root, Path('analysis'))
        _output_directory(root, Path('.revise'))
    except (OSError, ValueError):
        return
    path = root / 'analysis' / 'analysis.json'
    if not root.is_dir():
        return
    status = 'partial' if result.get('status') in {'succeeded', 'reused', 'inactive'} else 'incomplete'
    _json(path, {
        'status': status,
        'reconstruction_fingerprint': record.get('fingerprint'),
        'aspects': [result],
    })


def run_analysis_task(config_path: str | Path, sample_id: str, aspect: str, *,
                      cell_type: str | None = None) -> dict:
    """Run exactly one configured analysis aspect for one reconstruction task."""
    config_path = Path(config_path).resolve()
    _, inputs, output = load_batch_config(config_path)
    sample_dir = _selected_sample(inputs, sample_id)
    resolved = resolve_sample(config_path, sample_dir)
    _validate_task_selection(resolved, cell_type)
    selected_cell_type = cell_type
    specs = _specifications(resolved.document)
    if not _safe_aspect_name(aspect) or aspect not in specs:
        raise ValueError(f'Analysis aspect {aspect!r} is not configured for {sample_id}')
    root = _task_root(output, resolved.sample_id, selected_cell_type)
    if not output.is_dir():
        if not resolved.document.get('enabled', True):
            return _task_result({'aspect': aspect, 'directory': str(root), 'status': 'inactive',
                                 'reason': 'sample disabled'},
                                resolved.sample_id, selected_cell_type)
        return _task_result({'aspect': aspect, 'directory': str(root), 'status': 'failed',
                             'error': 'Reconstruction output root does not exist'},
                            resolved.sample_id, selected_cell_type)
    with _batch_lock(output):
        try:
            _output_directory(output, root.relative_to(output))
            if not resolved.document.get('enabled', True):
                if not root.is_dir():
                    result = {
                        'aspect': aspect, 'directory': str(root), 'status': 'inactive',
                        'reason': 'sample disabled',
                        'configuration_audit': _configuration_audit(resolved),
                        'log_path': str(root / '.revise' / 'analysis' / f'{aspect}.log'),
                    }
                    return _task_result(result, resolved.sample_id, selected_cell_type)
                result = _run_aspect(
                    root, _read_json(root / 'reconstruction.json'), aspect,
                    specs[aspect], resolved, selected_cell_type,
                    inactive_reason='sample disabled',
                )
                _update_single_summary(root, _read_json(root / 'reconstruction.json'), result)
                return _task_result(result, resolved.sample_id, selected_cell_type)
            record = verify_reconstruction_task(resolved, selected_cell_type)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            result = {'aspect': aspect, 'directory': str(root), 'status': 'failed',
                      'error': f'{type(exc).__name__}: {exc}'}
            _update_single_summary(root, _read_json(root / 'reconstruction.json'), result)
            return _task_result(result, resolved.sample_id, selected_cell_type)
        result = _run_aspect(root, record, aspect, specs[aspect], resolved, selected_cell_type)
        _update_single_summary(root, record, result)
        return _task_result(result, resolved.sample_id, selected_cell_type)


def _blocked_task(task: dict, error: Exception) -> dict:
    return {
        'sample_id': task.get('sample_id'),
        'cell_type': task.get('cell_type'),
        'status': 'blocked',
        'error': str(error),
    }


def _analysis_tasks(config_path: Path, input_root: Path, output: Path) -> list[dict]:
    """Discover current tasks; an inventory is only an ownership audit."""
    inventory_path = output / 'batch_status.json'
    inventory = _read_json(inventory_path)
    if inventory:
        if (inventory.get('project_config') != str(config_path)
                or inventory.get('input_root') != str(input_root)
                or inventory.get('output_root') != str(output)):
            raise ValueError('A completed reconstruction batch inventory is required')
    prior = {
        (task.get('sample_id'), task.get('cell_type')): task
        for task in inventory.get('tasks', []) if isinstance(task, dict)
    }
    tasks = []
    for sample_dir in discover_samples(input_root):
        sample_id = sample_dir.relative_to(input_root).as_posix()
        try:
            resolved = resolve_sample(config_path, sample_dir)
            if not resolved.document.get('enabled', True):
                continue
            labels = _cell_types(resolved.document)
        except Exception as exc:
            tasks.append({'sample_id': sample_id, 'cell_type': None, 'status': 'failed',
                          'error': f'{type(exc).__name__}: {exc}'})
            continue
        for cell_type in labels:
            root = _task_root(output, sample_id, cell_type)
            task = dict(prior.get((sample_id, cell_type), {}))
            task.update(sample_id=sample_id, cell_type=cell_type,
                        status='succeeded', directory=str(root))
            tasks.append(task)
    return tasks


def _invalidate_removed_aspects(root: Path, current: dict[str, dict]) -> None:
    control = root / '.revise' / 'analysis'
    if control.is_symlink() or not control.is_dir():
        return
    for state_path in control.glob('*.json'):
        name = state_path.stem
        if state_path.is_symlink() or name == 'analysis' or name in current:
            continue
        state = _read_json(state_path)
        if state:
            state.update(status='inactive', reason='aspect no longer configured')
            _json(state_path, state)


def _invalidate_removed_tasks(output: Path, previous: dict, current: list[dict]) -> None:
    """Invalidate tasks present in the previous analysis report only."""
    active = {(task.get('sample_id'), task.get('cell_type')) for task in current}
    for task in previous.get('tasks', []):
        key = (task.get('sample_id'), task.get('cell_type'))
        if key in active or not task.get('directory'):
            continue
        root = Path(task['directory'])
        try:
            _output_directory(output, root.relative_to(output))
            _output_directory(root, Path('analysis'))
            _output_directory(root, Path('.revise'))
        except (OSError, ValueError):
            continue
        control = root / '.revise' / 'analysis'
        if not control.is_symlink() and control.is_dir():
            for state_path in control.glob('*.json'):
                if state_path.stem == 'analysis' or state_path.is_symlink():
                    continue
                state = _read_json(state_path)
                if state:
                    state.update(status='inactive', reason='task no longer configured')
                    _json(state_path, state)
        _json(root / 'analysis' / 'analysis.json', {
            'status': 'inactive', 'reason': 'task no longer configured',
        })


def run_analysis_batch(config_path: str | Path) -> dict:
    """Run explicitly configured aspects against currently discovered tasks."""
    _, inputs, output = load_batch_config(config_path)
    config_path = Path(config_path).resolve()
    if not output.is_dir():
        raise ValueError('Run batch reconstruction before batch analysis')
    with _batch_lock(output):
        tasks = _analysis_tasks(config_path, inputs, output)
        previous_analysis = _read_json(output / 'analysis_status.json')
        if previous_analysis and any(
                previous_analysis.get(key) not in (None, expected)
                for key, expected in (
                    ('project_config', str(config_path)),
                    ('input_root', str(inputs)),
                    ('output_root', str(output)),
                )):
            previous_analysis = {}
        _invalidate_removed_tasks(output, previous_analysis, tasks)
        report = {
            'schema_version': 1, 'status': 'running', 'project_config': str(config_path),
            'input_root': str(inputs),
            'output_root': str(output), 'tasks': [],
        }
        report_path = output / 'analysis_status.json'
        _json(report_path, report)
        for task in tasks:
            root = None
            summary_safe = False
            try:
                if task.get('status') == 'failed':
                    raise ValueError(task.get('error', 'Sample configuration could not be resolved'))
                if task.get('status') not in {'succeeded', 'reused'}:
                    raise ValueError('Reconstruction task did not succeed')
                sample_id = task['sample_id']
                sample_dir = _selected_sample(inputs, sample_id)
                resolved = resolve_sample(config_path, sample_dir)
                _validate_task_selection(resolved, task.get('cell_type'))
                selected_cell_type = task.get('cell_type')
                root = _task_root(output, resolved.sample_id, selected_cell_type)
                _output_directory(output, root.relative_to(output))
                _output_directory(root, Path('analysis'))
                _output_directory(root, Path('.revise'))
                summary_safe = True
                if Path(task.get('directory', root)).resolve() != root.resolve():
                    raise ValueError('Reconstruction inventory task belongs to another output directory')
                record = verify_reconstruction_task(resolved, selected_cell_type)
                specs = _specifications(resolved.document)
                _invalidate_removed_aspects(root, specs)
            except (OSError, ValueError, KeyError, TypeError) as exc:
                if summary_safe and root is not None:
                    _json(root / 'analysis' / 'analysis.json', {'status': 'blocked', 'error': str(exc)})
                report['tasks'].append(_blocked_task(task, exc))
            else:
                summary_path = root / 'analysis' / 'analysis.json'
                _json(summary_path, {'status': 'running', 'reconstruction_fingerprint': record.get('fingerprint')})
                results = [
                    _run_aspect(root, record, name, spec, resolved, selected_cell_type)
                    for name, spec in specs.items()
                ]
                complete = all(item['status'] in {'succeeded', 'reused', 'inactive'} for item in results)
                _json(summary_path, {
                    'status': 'completed' if complete else 'incomplete',
                    'reconstruction_fingerprint': record.get('fingerprint'),
                    'aspects': results,
                })
                report['tasks'].extend(
                    _task_result(item, sample_id, selected_cell_type) for item in results
                )
            _json(report_path, report)
        counts = Counter(item['status'] for item in report['tasks'])
        report['summary'] = {
            name: counts[name] for name in (
                'succeeded', 'reused', 'failed', 'blocked', 'unavailable',
                'not_implemented', 'inactive',
            )
        }
        report['status'] = 'incomplete' if any(
            counts[name] for name in ('failed', 'blocked', 'unavailable', 'not_implemented')
        ) else 'completed'
        _json(report_path, report)
        return report
