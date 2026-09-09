"""Prepare explicitly described inputs without changing their source files."""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import os
from pathlib import Path
import tempfile
from typing import Any

import numpy as np
from scipy import sparse
import yaml

from revise.io.input_service import REVISEInputService


@dataclass
class PreparedSample:
    sample_id: str
    modality: str
    root: Path
    document: dict[str, Any]
    st_path: Path
    reference_path: Path
    metadata: dict[str, Any]


def file_identity(path: Path) -> dict[str, str]:
    """Hash a file or a directory store, including relative store filenames."""
    path = Path(path).resolve()
    if not path.exists():
        raise FileNotFoundError(path)
    digest = sha256()
    files = sorted(p for p in path.rglob('*') if p.is_file()) if path.is_dir() else [path]
    for item in files:
        if path.is_dir():
            name = item.relative_to(path).as_posix().encode()
            digest.update(len(name).to_bytes(8, 'big'))
            digest.update(name)
            digest.update(item.stat().st_size.to_bytes(8, 'big'))
        with item.open('rb') as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                digest.update(block)
    return {'path': str(path), 'sha256': digest.hexdigest()}


def _required(mapping: dict, key: str, context: str):
    value = mapping.get(key)
    if value is None or value == '':
        raise ValueError(f'{context}.{key} must be explicit')
    return value


def _validate(adata, matrix: str, role: str):
    for axis in ('obs_names', 'var_names'):
        names = getattr(adata, axis)
        if not names.is_unique or names.isna().any() or any(not str(x).strip() for x in names):
            raise ValueError(f'{role}.{axis} must be unique and non-null')
    if not adata.n_obs or not adata.n_vars:
        raise ValueError(f'{role} must have nonempty observations and genes')
    if matrix == 'X':
        selected = adata.X
    elif isinstance(matrix, str) and matrix.startswith('layers/') and matrix[7:] in adata.layers:
        selected = adata.layers[matrix[7:]]
    else:
        raise ValueError(f'{role}.matrix must select X or an existing layers/<name>')
    if selected is None:
        raise ValueError(f'{role}.matrix is missing')
    values = selected.data if sparse.issparse(selected) else np.asarray(selected)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError(f'{role}.matrix must be finite and nonnegative')
    if matrix != 'X':
        adata.X = selected.copy()


def prepare_sample(path: Path) -> PreparedSample:
    """Read sample.yaml and publish validated H5AD inputs with a reusable manifest.

    No expression normalization, coordinate scaling, or label inference is done.
    Reference filtering is validated here and applied by the application runner.
    """
    path = Path(path).resolve()
    root = path.parent
    if (root / 'inputs').is_symlink():
        raise ValueError('Prepared inputs directory must not be a symlink')
    document = yaml.safe_load(path.read_text())
    if not isinstance(document, dict):
        raise ValueError('sample.yaml must contain a mapping')
    if type(document.get('schema_version')) is not int or document['schema_version'] != 1:
        raise ValueError('schema_version must be 1')
    mapping_blocks = ('sample', 'preparation', 'preparation.spatial', 'preparation.reference',
                      'inputs', 'inputs.st', 'inputs.st.spatialdata', 'inputs.reference',
                      'algorithm', 'preprocessing', 'preprocessing.spatial', 'preprocessing.reference',
                      'global_anchoring', 'local_refinement', 'output', 'execution')
    for field in mapping_blocks:
        value = document
        for key in field.split('.'):
            value = value.get(key, {})
        if not isinstance(value, dict):
            raise ValueError(f'{field} must be a mapping')
    sample = document.get('sample', {})
    sample_id = _required(sample, 'id', 'sample')
    modality = _required(sample, 'modality', 'sample')
    if not isinstance(sample_id, str) or not sample_id.strip():
        raise ValueError('sample.id must be a nonempty string')
    if modality not in {'hST', 'iST', 'sST'}:
        raise ValueError('sample.modality must be hST, iST, or sST')
    preparation = document.get('preparation', {})
    spatial = preparation.get('spatial', {})
    reference = preparation.get('reference', {})
    st_matrix = _required(spatial, 'matrix', 'preparation.spatial')
    ref_matrix = _required(reference, 'matrix', 'preparation.reference')
    spatial_key = _required(spatial, 'spatial_key', 'preparation.spatial')
    unit = _required(spatial, 'coordinate_unit', 'preparation.spatial')
    if unit not in {'um', 'pixel'}:
        raise ValueError('coordinate_unit must be um or pixel')
    scale = spatial.get('microns_per_coordinate', 1.0 if unit == 'um' else None)
    if scale is not None and (not isinstance(scale, (int, float)) or not np.isfinite(scale) or scale <= 0):
        raise ValueError('microns_per_coordinate must be positive and finite')
    if unit == 'um' and scale != 1.0:
        raise ValueError('um coordinates require microns_per_coordinate=1')
    inputs = document.get('inputs', {})
    st_spec = inputs.get('st', {})
    ref_spec = inputs.get('reference', {})
    st_source = (root / _required(st_spec, 'path', 'inputs.st')).resolve()
    ref_source = (root / _required(ref_spec, 'path', 'inputs.reference')).resolve()
    st_format = _required(st_spec, 'format', 'inputs.st')
    if st_format not in {'h5ad', 'spatialdata'}:
        raise ValueError('inputs.st.format must explicitly be h5ad or spatialdata')
    if ref_spec.get('format', 'h5ad') != 'h5ad':
        raise ValueError('inputs.reference.format must be h5ad')
    selectors = st_spec.get('spatialdata', {})
    if st_format == 'spatialdata':
        _required(selectors, 'table', 'inputs.st.spatialdata')
        _required(selectors, 'element', 'inputs.st.spatialdata')
    mappings = reference.get('label_mapping', {})
    if not isinstance(mappings, dict) or any(not isinstance(v, dict) for v in mappings.values()):
        raise ValueError('reference.label_mapping must map columns to label mappings')
    broad = _required(document.get('global_anchoring', {}), 'broad_column', 'global_anchoring')
    subtype = document.get('local_refinement', {}).get('subtype_column')
    if modality == 'iST' and not subtype:
        raise ValueError('local_refinement.subtype_column must be explicit for iST')
    columns = list(dict.fromkeys([broad] + ([subtype] if subtype else [])))
    normalized_ref = ref_matrix != 'X' or bool(mappings)
    st_path = root / 'inputs' / 'spatial.h5ad'
    ref_path = root / 'inputs' / 'reference.h5ad' if normalized_ref else ref_source
    manifest = root / 'inputs' / 'preparation.json'
    targets = [st_path] + ([ref_path] if normalized_ref else [])
    for target in targets + [manifest]:
        if target.is_symlink():
            raise ValueError(f'Prepared destination must not be a symlink: {target}')
        if any(target.resolve() == source or (target.exists() and os.path.samefile(target, source))
               for source in (st_source, ref_source)):
            raise ValueError(f'Prepared destination aliases an original source: {target}')
    sources = {'spatial': file_identity(st_source), 'reference': file_identity(ref_source)}
    coordinates = {'source_key': spatial_key, 'unit': unit, 'microns_per_coordinate': scale,
                   'physical_distance_available': scale is not None}
    signature = {'cache_version': 1, 'preparation_code': file_identity(Path(__file__)),
                 'sources': sources, 'transforms': preparation,
                 'selectors': st_spec, 'reference_spec': ref_spec, 'required_labels': columns,
                 'coordinates': coordinates}
    cached = {}
    if manifest.exists():
        try:
            cached = json.loads(manifest.read_text())
            if not isinstance(cached, dict):
                cached = {}
            if all(cached.get(key) == value for key, value in signature.items()):
                actual = {'spatial': file_identity(st_path), 'reference': file_identity(ref_path)}
                if cached.get('outputs') == actual:
                    return PreparedSample(sample_id, modality, root, document, st_path, ref_path, cached)
        except (OSError, ValueError, TypeError):
            pass
    recorded_outputs = cached.get('outputs', {})
    for target in targets:
        if not target.exists():
            continue
        owned = (cached.get('cache_version') == 1
                 and isinstance(recorded_outputs, dict)
                 and any(isinstance(record, dict)
                         and record.get('path') == str(target)
                         and isinstance(record.get('sha256'), str)
                         and len(record['sha256']) == 64
                         for record in recorded_outputs.values()))
        if not owned:
            raise ValueError(f'Refusing to overwrite unowned prepared destination: {target}')
    service = REVISEInputService({'input_format': st_format,
                                 'spatialdata_table': selectors.get('table'),
                                 'spatialdata_spatial_element': selectors.get('element')})
    st = service.read_st_adata(st_source)
    ref = service.read_sc_ref_adata(ref_source)
    _validate(st, st_matrix, 'spatial')
    _validate(ref, ref_matrix, 'reference')
    if spatial_key not in st.obsm:
        raise ValueError(f'spatial_key {spatial_key!r} is missing from obsm')
    xy = np.asarray(st.obsm[spatial_key])
    if xy.ndim != 2 or xy.shape[0] != st.n_obs or xy.shape[1] < 2 or not np.isfinite(xy).all():
        raise ValueError('spatial coordinates must be finite with shape (n_obs, >=2)')
    st.obsm['spatial'] = xy.copy()
    for column in set(columns) | set(mappings):
        if column not in ref.obs or ref.obs[column].isna().any() or ref.obs[column].astype(str).str.strip().eq('').any():
            raise ValueError(f'reference label column {column!r} must exist and be non-null')
    for column, mapping in mappings.items():
        original = f'{column}_original'
        if original in ref.obs:
            raise ValueError(f'Cannot preserve labels: column {original!r} already exists')
        ref.obs[original] = ref.obs[column].copy()
        ref.obs[column] = ref.obs[column].astype(object).map(lambda value: mapping.get(value, value))
        if ref.obs[column].isna().any() or ref.obs[column].astype(str).str.strip().eq('').any():
            raise ValueError(f'label_mapping produces missing labels for {column!r}')
    filter_column = ref_spec.get('filter_column')
    filter_value = ref_spec.get('filter_value')
    if (filter_column is None) != (filter_value is None):
        raise ValueError('reference filter_column and filter_value must be supplied together')
    if filter_column is not None:
        if filter_column not in ref.obs or not ref.obs[filter_column].eq(filter_value).any():
            raise ValueError('reference filter must select at least one row')
    if st.var_names.intersection(ref.var_names).empty:
        raise ValueError('spatial and reference inputs must share genes')
    if sources != {'spatial': file_identity(st_source), 'reference': file_identity(ref_source)}:
        raise ValueError('Input source changed during preparation; retry with stable inputs')
    metadata = dict(signature)
    metadata['label_original_columns'] = {column: f'{column}_original' for column in mappings}
    metadata['coordinate_transform'] = 'copied_to_obsm_spatial_without_scaling'
    (root / 'inputs').mkdir(exist_ok=True)
    # Stage every artifact before replacement; rollback prior files on publication failure.
    with tempfile.TemporaryDirectory(prefix='.prepare-', dir=root / 'inputs') as temp:
        stage = Path(temp)
        staged = {st_path: stage / 'spatial.h5ad'}
        st.write_h5ad(staged[st_path])
        if normalized_ref:
            staged[ref_path] = stage / 'reference.h5ad'
            ref.write_h5ad(staged[ref_path])
        metadata['outputs'] = {
            'spatial': {**file_identity(staged[st_path]), 'path': str(st_path)},
            'reference': {**file_identity(staged[ref_path] if normalized_ref else ref_path), 'path': str(ref_path)},
        }
        staged[manifest] = stage / 'preparation.json'
        staged[manifest].write_text(json.dumps(metadata, indent=2, sort_keys=True) + '\n')
        backups = {}
        installed = []
        try:
            for destination, pending in staged.items():
                if destination.exists():
                    backup = stage / (destination.name + '.backup')
                    os.replace(destination, backup)
                    backups[destination] = backup
                os.replace(pending, destination)
                installed.append(destination)
        except BaseException:
            for destination in installed:
                destination.unlink(missing_ok=True)
            for destination, backup in backups.items():
                os.replace(backup, destination)
            raise
    return PreparedSample(sample_id, modality, root, document, st_path, ref_path, metadata)
