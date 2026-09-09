"""Read and validate standard H5AD inputs without creating or changing files."""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
from scipy import sparse

from .config import ResolvedSample


@dataclass
class BatchSample:
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


def _validate(adata, role: str) -> None:
    for axis in ('obs_names', 'var_names'):
        names = getattr(adata, axis)
        if not names.is_unique or names.isna().any() or any(not str(x).strip() for x in names):
            raise ValueError(f'{role}.{axis} must be unique and non-null')
    if not adata.n_obs or not adata.n_vars:
        raise ValueError(f'{role} must have nonempty observations and genes')
    if adata.X is None:
        raise ValueError(f'{role}.X is missing')
    for start in range(0, adata.n_obs, 4096):
        matrix = adata.X[start:start + 4096]
        values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(f'{role}.X must be finite and nonnegative')


def read_sample(resolved: ResolvedSample) -> BatchSample:
    """Validate input eligibility; QC, filtering and reconstruction stay in the engine."""
    document = resolved.document
    modality = document.get('modality')
    if modality not in {'hST', 'iST', 'sST'}:
        raise ValueError('modality must be hST, iST, or sST')
    coordinates = document.get('coordinates', {})
    if set(coordinates) - {'unit', 'microns_per_coordinate'}:
        raise ValueError('coordinates accepts unit and microns_per_coordinate only')
    unit = coordinates.get('unit')
    if unit not in {'um', 'pixel'}:
        raise ValueError('coordinates.unit must explicitly be um or pixel')
    scale = coordinates.get('microns_per_coordinate', 1.0 if unit == 'um' else None)
    if scale is not None and (type(scale) not in (int, float) or not np.isfinite(scale) or scale <= 0):
        raise ValueError('microns_per_coordinate must be positive and finite')
    if unit == 'um' and scale != 1.0:
        raise ValueError('um coordinates require microns_per_coordinate=1')
    ref_spec = document.get('inputs', {}).get('reference', {})
    if not ref_spec.get('path'):
        raise ValueError('inputs.reference.path must be explicit')
    if set(ref_spec) - {'path', 'format', 'filter_column', 'filter_value'} or ref_spec.get('format', 'h5ad') != 'h5ad':
        raise ValueError('reference must be standard H5AD with optional filter_column/filter_value')
    st_path = resolved.root / 'spatial.h5ad'
    ref_path = Path(ref_spec['path'])
    broad = document.get('global_anchoring', {}).get('broad_column')
    subtype = document.get('local_refinement', {}).get('subtype_column')
    if not isinstance(broad, str) or not broad.strip():
        raise ValueError('global_anchoring.broad_column must be explicit')
    if modality == 'iST' and (not isinstance(subtype, str) or not subtype.strip()):
        raise ValueError('local_refinement.subtype_column must be explicit for iST')
    sources = {'spatial': file_identity(st_path), 'reference': file_identity(ref_path)}
    st = ad.read_h5ad(st_path, backed='r')
    try:
        ref = ad.read_h5ad(ref_path, backed='r')
        try:
            _validate(st, 'spatial')
            _validate(ref, 'reference')
            xy = np.asarray(st.obsm.get('spatial'))
            if xy.ndim != 2 or xy.shape[0] != st.n_obs or xy.shape[1] < 2 or not np.isfinite(xy).all():
                raise ValueError('spatial coordinates must be finite with shape (n_obs, >=2)')
            for column in dict.fromkeys([broad] + ([subtype] if subtype else [])):
                if column not in ref.obs or ref.obs[column].isna().any() or ref.obs[column].astype(str).str.strip().eq('').any():
                    raise ValueError(f'reference label column {column!r} must exist and be non-null')
            column, value = ref_spec.get('filter_column'), ref_spec.get('filter_value')
            if (column is None) != (value is None):
                raise ValueError('reference filter_column and filter_value must be supplied together')
            if column is not None and (column not in ref.obs or not ref.obs[column].eq(value).any()):
                raise ValueError('reference filter must select at least one row')
            if st.var_names.intersection(ref.var_names).empty:
                raise ValueError('spatial and reference inputs must share genes')
        finally:
            ref.file.close()
    finally:
        st.file.close()
    if sources != {'spatial': file_identity(st_path), 'reference': file_identity(ref_path)}:
        raise ValueError('Input source changed during validation; retry with stable inputs')
    metadata = {
        'sources': sources, 'outputs': sources,
        'coordinates': {'source_key': 'spatial', 'unit': unit, 'microns_per_coordinate': scale,
                        'physical_distance_available': scale is not None},
        'configuration': {'project_config': str(resolved.project_config), 'sample_dir': str(resolved.root),
                          'chain': resolved.config_chain, 'effective': document},
    }
    return BatchSample(resolved.sample_id, modality, resolved.root, document, st_path, ref_path, metadata)
