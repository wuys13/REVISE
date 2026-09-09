"""Assemble paired iST carriers without changing reconstruction algorithms."""

from __future__ import annotations

import numpy as np
from anndata import AnnData
from scipy import sparse

from revise.utils.provenance import hash_jsonable


def _cluster_keys(adata, source):
    if 'SVC_cluster' not in adata.obs:
        raise ValueError(f'{source} is missing SVC_cluster')
    labels = adata.obs['SVC_cluster']
    if labels.isna().any():
        raise ValueError(f'{source} contains null SVC_cluster labels')
    return [(type(value), value) for value in labels.to_numpy(dtype=object)]


def assemble_ist(spatial, expression, *, mapping: str, seed: int):
    """Assign cluster means or sampled donors to the spatial observation axis."""
    spatial_keys = _cluster_keys(spatial, 'spatial carrier')
    expression_keys = _cluster_keys(expression, 'expression carrier')
    if not spatial_keys or set(spatial_keys) != set(expression_keys):
        raise ValueError('spatial and expression SVC cluster sets must match exactly and be nonempty')
    if not expression.var_names.is_unique:
        raise ValueError('expression var_names must be unique')
    values = expression.X.data if sparse.issparse(expression.X) else expression.X
    if not np.isfinite(values).all():
        raise ValueError('expression X must contain only finite values')
    donors = {}
    for index, key in enumerate(expression_keys):
        donors.setdefault(key, []).append(index)
    metadata = {'ist_mapping': mapping, 'expression_source': 'expression_carrier.X_as_is'}
    obs = spatial.obs.copy(deep=True)
    if mapping == 'mean':
        if sparse.issparse(expression.X):
            # A sparse averaging operator avoids densifying cluster-by-gene blocks.
            rows, columns, weights = [], [], []
            clusters = list(donors)
            for row, key in enumerate(clusters):
                indices = donors[key]
                rows.extend([row] * len(indices))
                columns.extend(indices)
                weights.extend([1.0 / len(indices)] * len(indices))
            averaging = sparse.csr_matrix((weights, (rows, columns)), shape=(len(clusters), expression.n_obs))
            means = averaging @ expression.X
            lookup = {key: i for i, key in enumerate(clusters)}
            output_x = means[[lookup[key] for key in spatial_keys]].tocsr()
        else:
            means = {key: expression.X[indices].mean(axis=0) for key, indices in donors.items()}
            output_x = np.vstack([means[key] for key in spatial_keys])
    elif mapping == 'random':
        if not expression.obs_names.is_unique or any(not str(name) for name in expression.obs_names):
            raise ValueError('expression obs_names must be unique and nonempty for random mapping')
        for indices in donors.values():
            indices.sort(key=lambda index: str(expression.obs_names[index]))
        rng = np.random.default_rng(seed)
        chosen = [donors[key][int(rng.integers(len(donors[key])))] for key in spatial_keys]
        donor_ids = expression.obs_names[chosen].tolist()
        output_x = expression.X[chosen].copy()
        obs['revise_ist_donor_id'] = donor_ids
        metadata.update(effective_seed=seed, donor_column='revise_ist_donor_id', donor_sha256=hash_jsonable(donor_ids))
    else:
        raise ValueError('iST assembly mapping must be mean or random')
    output_values = output_x.data if sparse.issparse(output_x) else output_x
    if not np.isfinite(output_values).all():
        raise ValueError("iST output X must contain only finite values")
    return AnnData(
        X=output_x, obs=obs, var=expression.var.copy(deep=True),
        uns={'revise_reconstruction': metadata},
        obsm={key: value.copy() for key, value in spatial.obsm.items()},
        obsp={key: value.copy() for key, value in spatial.obsp.items()},
        varm={key: value.copy() for key, value in expression.varm.items()},
        varp={key: value.copy() for key, value in expression.varp.items()},
    )
