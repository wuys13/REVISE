"""Discover standard ST directories and resolve their inherited configuration."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from hashlib import sha256
import os
from pathlib import Path

import yaml


ROOT_FIELDS = {'schema_version', 'input_root', 'output_root'}
MAPPING_FIELDS = {'coordinates', 'inputs', 'algorithm', 'preprocessing', 'global_anchoring',
                  'local_refinement', 'output', 'execution', 'analysis'}
FIELDS = ROOT_FIELDS | MAPPING_FIELDS | {'modality', 'enabled'}


@dataclass
class ResolvedSample:
    root: Path
    sample_id: str
    document: dict
    config_chain: list[dict]
    project_config: Path
    input_root: Path
    output_root: Path


def _read(path: Path, *, project: bool) -> tuple[dict, dict]:
    if path.is_symlink():
        raise ValueError(f'Configuration must not be a symlink: {path}')
    content = path.read_bytes()
    try:
        document = yaml.safe_load(content)
    except yaml.YAMLError as error:
        raise ValueError(f'Invalid YAML in {path}: {error}') from error
    if not isinstance(document, dict):
        raise ValueError(f'{path} must contain a mapping')
    if project and (type(document.get('schema_version')) is not int or document['schema_version'] != 2):
        raise ValueError('schema_version must be 2; migrate legacy sample.yaml packages to standard ST directories')
    forbidden = set(document) - FIELDS if project else set(document) - (FIELDS - ROOT_FIELDS)
    if forbidden:
        raise ValueError(f'Unsupported configuration fields in {path}: {sorted(forbidden)}')
    for key in MAPPING_FIELDS & document.keys():
        if not isinstance(document[key], dict):
            raise ValueError(f'{key} must be a mapping in {path}')
    if 'enabled' in document and type(document['enabled']) is not bool:
        raise ValueError('enabled must be true or false')
    inputs = document.get('inputs', {})
    if set(inputs) - {'reference', 'pm_on_cell'}:
        raise ValueError('inputs accepts only reference and pm_on_cell; ST is discovered as spatial.h5ad')
    for role, spec in inputs.items():
        if not isinstance(spec, dict):
            raise ValueError(f'inputs.{role} must be a mapping')
        value = spec.get('path')
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f'inputs.{role}.path must be explicit')
        spec['path'] = str((path.parent / value).resolve())
    analysis = document.get('analysis', {})
    if not isinstance(analysis, dict):
        raise ValueError(f'analysis must be a mapping in {path}')
    for aspect, spec in analysis.items():
        if not isinstance(spec, dict):
            raise ValueError(f'analysis.{aspect} must be a mapping in {path}')
        resources = spec.get('resources', {})
        if not isinstance(resources, dict):
            raise ValueError(f'analysis.{aspect}.resources must be a mapping in {path}')
        for name, value in resources.items():
            if not isinstance(name, str) or not name.strip():
                raise ValueError(f'analysis.{aspect}.resources names must be non-empty strings')
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f'analysis.{aspect}.resources.{name}.path must be explicit')
            resources[name] = str((path.parent / value).resolve())
    return document, {'path': str(path), 'sha256': sha256(content).hexdigest()}


def _load_project(path: Path) -> tuple[dict, Path, Path, dict]:
    document, identity = _read(path, project=True)
    roots = []
    for key in ('input_root', 'output_root'):
        value = document.get(key)
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f'{key} must be explicit')
        root = path.parent / value
        if root.is_symlink():
            raise ValueError(f'{key} must not be a symlink: {root}')
        roots.append(root.resolve())
    input_root, output_root = roots
    if input_root == output_root or input_root in output_root.parents or output_root in input_root.parents:
        raise ValueError('input_root and output_root must not overlap')
    if not input_root.is_dir():
        raise ValueError(f'input_root must be an existing directory: {input_root}')
    return document, input_root, output_root, identity


def _project_path(path: Path) -> Path:
    path = Path(path)
    if path.is_symlink():
        raise ValueError(f'Configuration must not be a symlink: {path}')
    return path.resolve()


def load_batch_config(path: Path) -> tuple[dict, Path, Path]:
    document, input_root, output_root, _ = _load_project(_project_path(path))
    return document, input_root, output_root


def _merge(base: dict, override: dict, trail: tuple = ()) -> dict:
    result = deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict) and trail + (key,) != ('inputs', 'reference'):
            result[key] = _merge(result[key], value, trail + (key,))
        else:
            result[key] = deepcopy(value)
    return result


def resolve_sample(config_path: Path, sample_dir: Path) -> ResolvedSample:
    project_path = _project_path(config_path)
    document, input_root, output_root, project_identity = _load_project(project_path)
    root = Path(sample_dir).resolve()
    relative = root.relative_to(input_root)
    paths = [project_path]
    current = input_root
    for part in [None, *relative.parts]:
        if part is not None:
            current /= part
        candidate = current / 'batch.yaml'
        if candidate != project_path:
            paths.append(candidate)
    chain = [project_identity]
    for path in paths[1:]:
        if path.is_symlink():
            raise ValueError(f'Configuration must not be a symlink: {path}')
        if path.exists():
            override, identity = _read(path, project=False)
            document = _merge(document, override)
            chain.append(identity)
        else:
            chain.append({'path': str(path), 'sha256': None})
    for key in ROOT_FIELDS:
        document.pop(key, None)
    document.setdefault('enabled', True)
    return ResolvedSample(root, relative.as_posix(), document, chain, project_path, input_root, output_root)


def discover_samples(input_root: Path) -> list[Path]:
    root = Path(input_root).resolve()
    result = []
    for directory, children, files in os.walk(root, followlinks=False):
        children[:] = sorted(name for name in children if not (Path(directory) / name).is_symlink())
        if 'spatial.h5ad' in files:
            result.append(Path(directory))
    return sorted(result)
