"""Assemble paired or optimal-transport iST output carriers."""

from __future__ import annotations

import numpy as np
from anndata import AnnData
from scipy import sparse

from revise.utils.provenance import hash_jsonable

from .ist_ot import (
    compile_ot_options,
    display_label,
    first_seen,
    overlap_gene_positions,
    smooth_spatial_overlap,
    solve_group,
    typed_labels,
    validate_matrix,
)


def _cluster_keys(adata, source):
    return typed_labels(adata, "SVC_cluster", source)


def _cluster_means(expression):
    """Compute sparse-preserving donor means in first-seen cluster order."""
    expression_keys = _cluster_keys(expression, 'expression carrier')
    if not expression_keys:
        raise ValueError('expression carrier SVC cluster set must be nonempty')
    if not expression.var_names.is_unique:
        raise ValueError('expression var_names must be unique')
    values = expression.X.data if sparse.issparse(expression.X) else expression.X
    if not np.isfinite(values).all():
        raise ValueError('expression X must contain only finite values')
    donors = {}
    for index, key in enumerate(expression_keys):
        donors.setdefault(key, []).append(index)
    clusters = list(donors)
    if sparse.issparse(expression.X):
        # A sparse averaging operator avoids densifying cluster-by-gene blocks.
        rows, columns, weights = [], [], []
        for row, key in enumerate(clusters):
            indices = donors[key]
            rows.extend([row] * len(indices))
            columns.extend(indices)
            weights.extend([1.0 / len(indices)] * len(indices))
        averaging = sparse.csr_matrix((weights, (rows, columns)), shape=(len(clusters), expression.n_obs))
        means = (averaging @ expression.X).tocsr()
    else:
        means = np.vstack([expression.X[donors[key]].mean(axis=0) for key in clusters])
    return means, clusters


def cluster_means(expression):
    """Return donor expression means and their first-seen cluster labels."""
    means, clusters = _cluster_means(expression)
    output_values = means.data if sparse.issparse(means) else means
    if not np.isfinite(output_values).all():
        raise ValueError("cluster means must contain only finite values")
    return means, [key[1] for key in clusters]


def _cluster_broad_mapping(adata, cluster_keys, broad_column, source):
    broad_keys = typed_labels(adata, broad_column, source)
    result = {}
    for cluster, broad in zip(cluster_keys, broad_keys):
        previous = result.setdefault(cluster, broad)
        if previous != broad:
            raise ValueError(
                f"{source} SVC_cluster {cluster[1]!r} has inconsistent "
                f"{broad_column} labels"
            )
    return broad_keys, result


def _cluster_profiles(expression, expression_keys, overlap_positions):
    clusters = first_seen(expression_keys)
    donors = {cluster: [] for cluster in clusters}
    for index, cluster in enumerate(expression_keys):
        donors[cluster].append(index)
    overlap = expression.X[:, overlap_positions]
    overlap_profiles = []
    full_profiles = []
    masses = []
    for cluster in clusters:
        indices = donors[cluster]
        overlap_values = overlap[indices]
        masses.append(float(np.asarray(overlap_values.sum()).item()))
        if sparse.issparse(overlap_values):
            overlap_profiles.append(np.asarray(overlap_values.mean(axis=0)).ravel())
            full_profiles.append(expression.X[indices].mean(axis=0))
        else:
            overlap_profiles.append(np.asarray(overlap_values).mean(axis=0))
            full_profiles.append(np.asarray(expression.X[indices]).mean(axis=0))
    overlap_result = np.vstack(overlap_profiles)
    if sparse.issparse(expression.X):
        full_result = sparse.vstack(
            [sparse.csr_matrix(row) for row in full_profiles]
        ).tocsr()
    else:
        full_result = np.vstack(full_profiles)
    return clusters, overlap_result, full_result, np.asarray(masses, dtype=np.float64)


def _assemble_ot(spatial, expression, *, mapping, broad_column, ot_options):
    options = compile_ot_options(ot_options)
    validate_matrix(spatial.X, "spatial X")
    validate_matrix(expression.X, "expression X")
    spatial_keys = typed_labels(spatial, "SVC_cluster", "spatial carrier")
    expression_keys = typed_labels(expression, "SVC_cluster", "expression carrier")
    if not spatial_keys:
        raise ValueError("spatial carrier must contain at least one observation")
    if not expression_keys:
        raise ValueError("expression carrier must contain at least one observation")
    spatial_gene_positions, expression_gene_positions = overlap_gene_positions(
        spatial, expression
    )
    smoothed, smoothing_metadata = smooth_spatial_overlap(
        spatial,
        spatial_gene_positions,
        spatial_keys,
        spatial_weight=options.spatial_weight,
    )
    expression_overlap = expression.X[:, expression_gene_positions]
    group_keys: list[tuple[type, object]]
    expression_candidates: dict[tuple[type, object], list[int]] = {}
    profile_candidates = None
    if mapping == "within_cluster":
        group_keys = first_seen(spatial_keys)
        for index, cluster in enumerate(expression_keys):
            expression_candidates.setdefault(cluster, []).append(index)
        missing = [key for key in group_keys if key not in expression_candidates]
        if missing:
            labels = ", ".join(repr(key[1]) for key in missing)
            raise ValueError(
                "expression carrier has no within-cluster candidates for: "
                f"{labels}"
            )
    else:
        spatial_broad, _ = _cluster_broad_mapping(
            spatial, spatial_keys, broad_column, "spatial carrier"
        )
        _, expression_cluster_broad = _cluster_broad_mapping(
            expression, expression_keys, broad_column, "expression carrier"
        )
        group_keys = first_seen(spatial_broad)
        clusters, overlap_profiles, full_profiles, cluster_masses = _cluster_profiles(
            expression, expression_keys, expression_gene_positions
        )
        by_broad: dict[tuple[type, object], list[int]] = {}
        for index, cluster in enumerate(clusters):
            by_broad.setdefault(expression_cluster_broad[cluster], []).append(index)
        missing = [key for key in group_keys if key not in by_broad]
        if missing:
            labels = ", ".join(repr(key[1]) for key in missing)
            raise ValueError(
                "expression carrier has no outside-cluster candidates for "
                f"{broad_column}: {labels}"
            )
        profile_candidates = (
            spatial_broad,
            by_broad,
            overlap_profiles,
            full_profiles,
            cluster_masses,
        )

    output = np.empty((spatial.n_obs, expression.n_vars), dtype=np.float64)
    labels = []
    spatial_sizes = []
    candidate_sizes = []
    timings = []
    for group in group_keys:
        if mapping == "within_cluster":
            spatial_indices = np.asarray(
                [i for i, key in enumerate(spatial_keys) if key == group]
            )
            candidate_indices = np.asarray(expression_candidates[group])
            candidate_overlap = expression_overlap[candidate_indices]
            candidate_full = expression.X[candidate_indices]
            candidate_overlap_dense = (
                candidate_overlap.toarray()
                if sparse.issparse(candidate_overlap)
                else np.asarray(candidate_overlap, dtype=np.float64)
            )
            target_mass = np.asarray(
                candidate_overlap_dense.sum(axis=1), dtype=np.float64
            )
        else:
            (
                spatial_broad,
                by_broad,
                overlap_profiles,
                full_profiles,
                cluster_masses,
            ) = profile_candidates
            spatial_indices = np.asarray(
                [i for i, key in enumerate(spatial_broad) if key == group]
            )
            candidate_indices = np.asarray(by_broad[group])
            candidate_overlap_dense = overlap_profiles[candidate_indices]
            candidate_full = full_profiles[candidate_indices]
            target_mass = cluster_masses[candidate_indices]
        elapsed = solve_group(
            smoothed[spatial_indices],
            candidate_overlap_dense,
            candidate_full,
            target_mass,
            options=options,
            group_label=display_label(group),
            output=output,
            output_rows=spatial_indices,
        )
        labels.append(display_label(group))
        spatial_sizes.append(int(spatial_indices.size))
        candidate_sizes.append(int(candidate_indices.size))
        timings.append(float(elapsed))

    metadata = {
        "ist_mapping": mapping,
        "expression_source": "expression_carrier.X_as_is",
        "ot_effective_options": options.metadata(),
        "ot_group_labels": labels,
        "ot_group_spatial_sizes": spatial_sizes,
        "ot_group_candidate_sizes": candidate_sizes,
        "ot_group_timings_seconds": timings,
        "ot_smoothing": smoothing_metadata,
        "ot_output_bytes": int(output.nbytes),
        "ot_source_mass": "spatial_overlap_expression_total",
        "ot_target_mass": (
            "sum_of_donor_overlap_expression_per_cluster"
            if mapping == "outside_cluster" else "donor_overlap_expression_total"
        ),
    }
    if mapping == "outside_cluster":
        metadata["ot_broad_column"] = broad_column
    return output, metadata


def assemble_ist(
    spatial,
    expression,
    *,
    mapping: str,
    seed: int,
    broad_column: str = "Level1",
    ot_options: dict | None = None,
):
    """Assign cluster means or sampled donors to the spatial observation axis."""
    obs = spatial.obs.copy(deep=True)
    if mapping in {"within_cluster", "outside_cluster"}:
        output_x, metadata = _assemble_ot(
            spatial,
            expression,
            mapping=mapping,
            broad_column=broad_column,
            ot_options=ot_options,
        )
    else:
        spatial_keys = _cluster_keys(spatial, 'spatial carrier')
        expression_keys = _cluster_keys(expression, 'expression carrier')
        if not spatial_keys or set(spatial_keys) != set(expression_keys):
            raise ValueError(
                'spatial and expression SVC cluster sets must match exactly '
                'and be nonempty'
            )
        if not expression.var_names.is_unique:
            raise ValueError('expression var_names must be unique')
        validate_matrix(expression.X, "expression X")
        metadata = {
            'ist_mapping': mapping,
            'expression_source': 'expression_carrier.X_as_is',
        }
    if mapping == 'mean':
        means, clusters = _cluster_means(expression)
        lookup = {key: i for i, key in enumerate(clusters)}
        positions = [lookup[key] for key in spatial_keys]
        output_x = means[positions].tocsr() if sparse.issparse(means) else means[positions]
    elif mapping == 'random':
        donors = {}
        for index, key in enumerate(expression_keys):
            donors.setdefault(key, []).append(index)
        if not expression.obs_names.is_unique or any(not str(name) for name in expression.obs_names):
            raise ValueError('expression obs_names must be unique and nonempty for random mapping')
        for indices in donors.values():
            indices.sort(key=lambda index: str(expression.obs_names[index]))
        rng = np.random.default_rng(seed)
        chosen = [donors[key][int(rng.integers(len(donors[key])))] for key in spatial_keys]
        donor_ids = expression.obs_names[chosen].tolist()
        output_x = expression.X[chosen].copy()
        obs['revise_ist_donor_id'] = donor_ids
        metadata.update(
            effective_seed=seed,
            donor_column='revise_ist_donor_id',
            donor_sha256=hash_jsonable(donor_ids),
        )
    elif mapping not in {"within_cluster", "outside_cluster"}:
        raise ValueError(
            'iST assembly mapping must be mean, random, within_cluster, or outside_cluster'
        )
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
