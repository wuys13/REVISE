"""Batch analysis adapters consume verified reconstruction handoffs."""
from __future__ import annotations

from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import dataclass
from importlib import import_module
import inspect
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import traceback

from .runner import (
    ANALYSIS_ASPECTS, _batch_lock, _code_identity, _digest, _json,
    _output_directory, _read_json, _reusable, load_batch_config,
)
from .sample import file_identity


@dataclass(frozen=True)
class AnalysisContext:
    """Adapter input. Write only beneath output_dir; return role -> relative file."""
    reconstruction: dict
    output_dir: Path
    parameters: dict


def _specifications(batch: dict) -> dict:
    specs = batch.get('analysis', {})
    if specs == {}:
        specs = {name: {} for name in ANALYSIS_ASPECTS}
    if not isinstance(specs, dict) or set(specs) - set(ANALYSIS_ASPECTS):
        raise ValueError(f'analysis must map aspects from {ANALYSIS_ASPECTS} to adapter settings')
    for name, spec in specs.items():
        if not isinstance(spec, dict) or set(spec) - {'entrypoint', 'version', 'parameters', 'requires_pairing'}:
            raise ValueError(f'Invalid analysis settings for {name}')
        if 'entrypoint' in spec:
            if not isinstance(spec['entrypoint'], str) or spec['entrypoint'].count(':') != 1:
                raise ValueError(f'{name}.entrypoint must be module:function')
            if not isinstance(spec.get('version'), str) or not spec['version'].strip():
                raise ValueError(f'{name}.version must identify the algorithm and its dependencies')
        if not isinstance(spec.get('parameters', {}), dict) or type(spec.get('requires_pairing', False)) is not bool:
            raise ValueError(f'Invalid parameters or requires_pairing for {name}')
    return specs


def _verified_handoff(root: Path) -> dict:
    state = _read_json(root / '.revise' / 'task.json')
    if not _reusable(state, state.get('fingerprint')):
        raise ValueError('Reconstruction is not successful or its artifacts changed')
    record = _read_json(root / 'reconstruction.json')
    if record.get('status') != 'succeeded' or record.get('fingerprint') != state.get('fingerprint'):
        raise ValueError('Reconstruction handoff is not current')
    identities = list(state.get('inputs', {}).values()) + list(record['inputs']['sources'].values())
    if any(file_identity(Path(item['path'])) != item for item in identities):
        raise ValueError('Reconstruction inputs changed; rerun reconstruction first')
    return record


def _run_aspect(root: Path, record: dict, name: str, spec: dict, code: dict) -> dict:
    result = {'aspect': name, 'directory': str(root), 'reconstruction_fingerprint': record['fingerprint']}
    control = root / '.revise' / 'analysis'
    destination = root / 'analysis' / name
    state_path = control / f'{name}.json'
    writable = False
    try:
        _output_directory(root, Path('.revise/analysis'))
        _output_directory(root, Path('analysis') / name)
        control.mkdir(parents=True, exist_ok=True)
        writable = True
        previous = _read_json(state_path)
        if previous.get('published_directory') == str(destination):
            result['published_directory'] = str(destination)
        if 'entrypoint' not in spec:
            result.update(status='not_implemented', reason='No analysis adapter configured')
        elif spec.get('requires_pairing', False) and record['pairing']['status'] != 'available':
            result.update(status='unavailable', reason=record['pairing']['reason'])
        else:
            module_name, function_name = spec['entrypoint'].split(':')
            adapter = getattr(import_module(module_name), function_name)
            source = inspect.getsourcefile(adapter)
            if not callable(adapter) or source is None:
                raise ValueError('Analysis adapter must be a Python function with identifiable source')
            fingerprint = _digest({'reconstruction': record, 'specification': spec,
                                   'code': code, 'adapter': file_identity(Path(source))})
            result.update(fingerprint=fingerprint, specification=spec)
            if _reusable(previous, fingerprint):
                return dict(previous, status='reused')
            # Nonempty aspect directories must have been published by this runner.
            if destination.exists() and any(destination.iterdir()):
                if previous.get('published_directory') != str(destination):
                    raise ValueError(f'Refusing to replace unowned analysis directory: {destination}')
            _json(state_path, dict(result, status='running'))
            destination.parent.mkdir(parents=True, exist_ok=True)
            with TemporaryDirectory(prefix=f'.{name}-', dir=destination.parent) as temporary:
                stage = Path(temporary) / 'artifacts'
                stage.mkdir()
                with (control / f'{name}.log').open('w') as log, redirect_stdout(log), redirect_stderr(log):
                    artifacts = adapter(AnalysisContext(record, stage, spec.get('parameters', {})))
                if not isinstance(artifacts, dict) or not artifacts:
                    raise ValueError('Analysis adapter must return nonempty role -> relative file mapping')
                records = []
                for role, relative in artifacts.items():
                    if not isinstance(role, str) or not isinstance(relative, str):
                        raise ValueError('Analysis artifact roles and paths must be strings')
                    path = stage / relative
                    if Path(relative).is_absolute() or not path.resolve().is_relative_to(stage.resolve()) or not path.is_file():
                        raise ValueError(f'Analysis artifact must be a file inside its aspect: {relative}')
                    records.append({**file_identity(path), 'path': str(destination / relative)})
                if any(path.is_symlink() for path in stage.rglob('*')):
                    raise ValueError('Analysis artifacts must not contain symlinks')
                if _verified_handoff(root) != record:
                    raise ValueError('Reconstruction changed during analysis')
                backup = Path(temporary) / 'previous'
                if destination.exists():
                    os.replace(destination, backup)
                try:
                    os.replace(stage, destination)
                except BaseException:
                    if backup.exists():
                        os.replace(backup, destination)
                    raise
                result.update(status='succeeded', artifacts=records, roles=artifacts, published_directory=str(destination))
        _json(state_path, result)
    except (Exception, SystemExit) as exc:
        result.update(status='failed', error=f'{type(exc).__name__}: {exc}')
        if writable:
            with (control / f'{name}.log').open('a') as log:
                traceback.print_exc(file=log)
            _json(state_path, result)
    return result


def run_analysis_batch(config_path: str | Path) -> dict:
    """Run configured aspects against the latest reconstruction batch inventory."""
    batch, inputs, output = load_batch_config(config_path)
    specs = _specifications(batch)
    if not output.is_dir():
        raise ValueError('Run batch reconstruction before batch analysis')
    with _batch_lock(output):
        inventory = _read_json(output / 'batch_status.json')
        if (inventory.get('input_root') != str(inputs) or inventory.get('output_root') != str(output)
                or inventory.get('status') not in {'completed', 'completed_with_failures'}):
            raise ValueError('A completed reconstruction batch inventory is required')
        report = {'schema_version': 1, 'status': 'running', 'input_root': str(inputs),
                  'output_root': str(output), 'tasks': []}
        report_path = output / 'analysis_status.json'
        _json(report_path, report)
        code = _code_identity()
        for task in inventory['tasks']:
            root = None
            summary_safe = False
            try:
                if task['status'] not in {'succeeded', 'reused'}:
                    raise ValueError('Reconstruction task did not succeed')
                root = Path(task['directory'])
                _output_directory(output, root.relative_to(output))
                _output_directory(root, Path('analysis'))
                _output_directory(root, Path('.revise'))
                summary_safe = True
                record = _verified_handoff(root)
            except (OSError, ValueError, KeyError) as exc:
                if summary_safe:
                    _json(root / 'analysis' / 'analysis.json', {'status': 'blocked', 'error': str(exc)})
                report['tasks'].append({'sample_id': task.get('sample_id'), 'cell_type': task.get('cell_type'),
                                        'status': 'blocked', 'error': str(exc)})
            else:
                summary_path = root / 'analysis' / 'analysis.json'
                _json(summary_path, {'status': 'running', 'reconstruction_fingerprint': record['fingerprint']})
                results = [_run_aspect(root, record, name, spec, code) for name, spec in specs.items()]
                complete = all(item['status'] in {'succeeded', 'reused'} for item in results)
                _json(summary_path, {'status': 'completed' if complete else 'incomplete',
                                     'reconstruction_fingerprint': record['fingerprint'], 'aspects': results})
                report['tasks'].extend(dict(item, sample_id=task['sample_id'], cell_type=task['cell_type']) for item in results)
            _json(report_path, report)
        counts = Counter(item['status'] for item in report['tasks'])
        report['summary'] = {name: counts[name] for name in (
            'succeeded', 'reused', 'failed', 'blocked', 'unavailable', 'not_implemented')}
        report['status'] = 'incomplete' if any(counts[name] for name in (
            'failed', 'blocked', 'unavailable', 'not_implemented')) else 'completed'
        _json(report_path, report)
        return report
