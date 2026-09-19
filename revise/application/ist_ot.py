"""Optimal-transport assembly primitives for iST output carriers."""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter
from typing import Any, Iterable

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from revise.backend.kernels.ot import OTKernel
from revise.backend.ops.distance import bhattacharyya_distance


_DEFAULT_OPTIONS = {
    "spatial_weight": 0.2,
    "max_cost_entries": 2_000_000,
    "gene_block_size": 256,
}


@dataclass(frozen=True)
class OTOptions:
    spatial_weight: float
    max_cost_entries: int
    gene_block_size: int

    def metadata(self) -> dict[str, float | int | str]:
        return {
            "method": "tacco",
            "spatial_weight": self.spatial_weight,
            "max_cost_entries": self.max_cost_entries,
            "gene_block_size": self.gene_block_size,
        }


def compile_ot_options(options: dict[str, Any] | None) -> OTOptions:
    supplied = {} if options is None else options
    if not isinstance(supplied, dict):
        raise ValueError("ot_options must be a dictionary or None")
    unknown = sorted(set(supplied) - set(_DEFAULT_OPTIONS))
    if unknown:
        raise ValueError(f"ot_options contains unsupported key(s): {', '.join(unknown)}")
    values = {**_DEFAULT_OPTIONS, **supplied}

    spatial_weight = values["spatial_weight"]
    if isinstance(spatial_weight, bool) or not isinstance(
        spatial_weight,
        (int, float, np.integer, np.floating),
    ):
        raise ValueError("ot_options.spatial_weight must be a finite number in [0, 1]")
    spatial_weight = float(spatial_weight)
    if not np.isfinite(spatial_weight) or not 0 <= spatial_weight <= 1:
        raise ValueError("ot_options.spatial_weight must be a finite number in [0, 1]")

    integers: dict[str, int] = {}
    for key in ("max_cost_entries", "gene_block_size"):
        value = values[key]
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, np.integer))
            or value <= 0
        ):
            raise ValueError(f"ot_options.{key} must be a positive integer")
        integers[key] = int(value)
    return OTOptions(spatial_weight, integers["max_cost_entries"], integers["gene_block_size"])


def typed_labels(adata: AnnData, column: str, source: str) -> list[tuple[type, Any]]:
    if column not in adata.obs:
        raise ValueError(f"{source} is missing {column}")
    keys: list[tuple[type, Any]] = []
    for value in adata.obs[column].to_numpy(dtype=object):
        try:
            is_null = bool(pd.isna(value))
        except (TypeError, ValueError):
            is_null = True
        if is_null or (isinstance(value, str) and not value.strip()):
            raise ValueError(f"{source} contains null {column} labels")
        try:
            hash(value)
        except TypeError as exc:
            raise ValueError(
                f"{source} contains an invalid unhashable {column} label"
            ) from exc
        keys.append((type(value), value))
    return keys


def validate_matrix(matrix: Any, name: str) -> None:
    if getattr(matrix, "ndim", None) != 2:
        raise ValueError(f"{name} must be a two-dimensional matrix")
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    try:
        finite = np.isfinite(values)
    except TypeError as exc:
        raise ValueError(f"{name} must contain numeric values") from exc
    if not np.all(finite):
        raise ValueError(f"{name} must contain only finite values")
    try:
        negative = values < 0
    except TypeError as exc:
        raise ValueError(f"{name} must contain numeric values") from exc
    if np.any(negative):
        raise ValueError(f"{name} must contain only non-negative values")


def _dense(matrix: Any) -> np.ndarray:
    values = matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)
    return np.asarray(values, dtype=np.float64)


def overlap_gene_positions(
    spatial: AnnData,
    expression: AnnData,
) -> tuple[np.ndarray, np.ndarray]:
    if not spatial.var_names.is_unique:
        raise ValueError("spatial var_names must be unique")
    if not expression.var_names.is_unique:
        raise ValueError("expression var_names must be unique")
    expression_lookup = {name: index for index, name in enumerate(expression.var_names)}
    spatial_positions: list[int] = []
    expression_positions: list[int] = []
    for index, name in enumerate(spatial.var_names):
        if name in expression_lookup:
            spatial_positions.append(index)
            expression_positions.append(expression_lookup[name])
    if not spatial_positions:
        raise ValueError("spatial and expression carriers have no overlapping genes")
    return np.asarray(spatial_positions), np.asarray(expression_positions)


def smooth_spatial_overlap(
    spatial: AnnData,
    overlap_positions: np.ndarray,
    cluster_keys: list[tuple[type, Any]],
    *,
    spatial_weight: float,
) -> tuple[np.ndarray, dict[str, int]]:
    """Smooth ST expression over a freshly built same-cluster Squidpy graph."""
    if "spatial" not in spatial.obsm:
        raise ValueError("spatial carrier is missing obsm['spatial'] coordinates")
    coordinates = np.asarray(spatial.obsm["spatial"])
    if (
        coordinates.ndim != 2
        or coordinates.shape[0] != spatial.n_obs
        or coordinates.shape[1] < 2
    ):
        raise ValueError("spatial coordinates must have shape (n_obs, at least 2)")
    if not np.all(np.isfinite(coordinates)):
        raise ValueError("spatial coordinates must contain only finite values")

    values = _dense(spatial.X[:, overlap_positions])
    if spatial.n_obs == 0:
        raise ValueError("spatial carrier must contain at least one observation")
    if spatial.n_obs == 1:
        return values.copy(), {"same_cluster_edges": 0, "isolated_rows": 1}

    # Construct a minimal copy so no stored graph or Visium metadata can change
    # Squidpy's generic-coordinate convention or be overwritten on the input.
    work = AnnData(
        X=np.zeros((spatial.n_obs, 1), dtype=np.float32),
        obs=pd.DataFrame(index=spatial.obs_names.copy()),
        obsm={"spatial": coordinates.copy()},
    )
    import squidpy as sq

    connectivities, _ = sq.gr.spatial_neighbors(
        work,
        coord_type="generic",
        n_neighs=min(6, spatial.n_obs - 1),
        copy=True,
    )
    graph = sparse.coo_matrix(connectivities)
    keep = np.fromiter(
        (
            row != column and cluster_keys[row] == cluster_keys[column]
            for row, column in zip(graph.row, graph.col)
        ),
        dtype=bool,
        count=graph.nnz,
    )
    same_cluster = sparse.csr_matrix(
        (graph.data[keep], (graph.row[keep], graph.col[keep])),
        shape=graph.shape,
    )
    same_cluster.eliminate_zeros()
    row_sums = np.asarray(same_cluster.sum(axis=1)).ravel()
    isolated = row_sums <= 0
    if np.any(~isolated):
        same_cluster = sparse.diags(
            np.divide(1.0, row_sums, out=np.zeros_like(row_sums), where=~isolated)
        ) @ same_cluster
    neighbor_values = same_cluster @ values
    smoothed = (1.0 - spatial_weight) * values + spatial_weight * neighbor_values
    smoothed[isolated] = values[isolated]
    validate_matrix(smoothed, "smoothed spatial overlap expression")
    return smoothed, {
        "same_cluster_edges": int(same_cluster.nnz),
        "isolated_rows": int(isolated.sum()),
    }


def _positive_masses(values: np.ndarray, name: str) -> np.ndarray:
    masses = np.asarray(values.sum(axis=1), dtype=np.float64).ravel()
    if not np.all(np.isfinite(masses)):
        raise ValueError(f"{name} overlap expression totals must be finite")
    if np.any(masses <= 0):
        raise ValueError(f"{name} overlap expression totals must be strictly positive")
    return masses


def _normalized_cost(
    source_overlap: np.ndarray,
    candidate_overlap: np.ndarray,
    *,
    options: OTOptions,
    group_label: str,
) -> np.ndarray:
    _check_cost_allocation(
        source_overlap.shape[0],
        candidate_overlap.shape[0],
        options=options,
        group_label=group_label,
    )
    distance = bhattacharyya_distance(candidate_overlap, source_overlap)
    cost = np.asarray(distance.T, dtype=np.float64)
    expected = (source_overlap.shape[0], candidate_overlap.shape[0])
    if cost.shape != expected:
        raise ValueError(f"Bhattacharyya cost has shape {cost.shape}, expected {expected}")
    if not np.all(np.isfinite(cost)):
        raise ValueError("Bhattacharyya cost must contain only finite values")
    if np.any(cost < 0):
        raise ValueError("Bhattacharyya cost must contain only non-negative values")
    maximum = float(cost.max(initial=0.0))
    return cost / maximum if maximum > 0 else cost


def _check_cost_allocation(
    source_count: int,
    candidate_count: int,
    *,
    options: OTOptions,
    group_label: str,
) -> None:
    entries = int(source_count) * int(candidate_count)
    if entries > options.max_cost_entries:
        raise ValueError(
            f"OT cost allocation for {group_label} requires {entries} entries, "
            f"exceeding max_cost_entries={options.max_cost_entries}; use "
            "outside_cluster to reduce donor candidates when scientifically "
            "appropriate, or raise max_cost_entries intentionally"
        )


def _weighted_full_expression(
    row_weights: np.ndarray,
    candidate_full: Any,
    *,
    gene_block_size: int,
    output: np.ndarray,
    output_rows: np.ndarray,
) -> None:
    n_rows, n_candidates = row_weights.shape
    if candidate_full.shape[0] != n_candidates:
        raise ValueError("candidate expression rows do not match OT weights")
    if len(output_rows) != n_rows or output.shape[1] != candidate_full.shape[1]:
        raise ValueError("output positions do not match the expression projection")
    for start in range(0, candidate_full.shape[1], gene_block_size):
        stop = min(start + gene_block_size, candidate_full.shape[1])
        block = candidate_full[:, start:stop]
        product = row_weights @ block
        values = product.toarray() if sparse.issparse(product) else product
        validate_matrix(values, "iST OT output X")
        output[output_rows, start:stop] = values


def solve_group(
    source_overlap: np.ndarray,
    candidate_overlap: np.ndarray,
    candidate_full: Any,
    target_mass: np.ndarray,
    *,
    options: OTOptions,
    group_label: str,
    output: np.ndarray,
    output_rows: np.ndarray,
) -> float:
    """Solve one complete group OT problem and project all expression genes."""
    started = perf_counter()
    source_mass = _positive_masses(source_overlap, "spatial")
    target_mass = np.asarray(target_mass, dtype=np.float64).ravel()
    if target_mass.shape != (candidate_overlap.shape[0],):
        raise ValueError("candidate overlap totals do not match candidate profiles")
    if not np.all(np.isfinite(target_mass)) or np.any(target_mass <= 0):
        raise ValueError(
            "candidate overlap expression totals must be finite and strictly positive"
        )

    # Enforce the allocation cap before either distance construction or the
    # single-candidate shortcut, so the option retains one simple contract.
    _check_cost_allocation(
        source_overlap.shape[0],
        candidate_overlap.shape[0],
        options=options,
        group_label=group_label,
    )
    if candidate_overlap.shape[0] == 1:
        row_weights = np.ones((source_overlap.shape[0], 1), dtype=np.float64)
    else:
        cost = _normalized_cost(
            source_overlap, candidate_overlap, options=options, group_label=group_label
        )
        coupling = OTKernel.couple(source_mass, target_mass, cost, method="tacco")
        coupling = np.asarray(coupling, dtype=np.float64)
        if coupling.shape != cost.shape:
            raise ValueError(
                f"OT coupling has shape {coupling.shape}, expected {cost.shape}"
            )
        if not np.all(np.isfinite(coupling)) or np.any(coupling < 0):
            raise ValueError("OT coupling weights must be finite and non-negative")
        row_mass = coupling.sum(axis=1)
        if np.any(row_mass <= 0):
            raise ValueError("OT coupling must have positive mass in every spatial row")
        row_weights = coupling / row_mass[:, None]
    _weighted_full_expression(
        row_weights, candidate_full, gene_block_size=options.gene_block_size,
        output=output, output_rows=output_rows,
    )
    return perf_counter() - started


def first_seen(keys: Iterable[tuple[type, Any]]) -> list[tuple[type, Any]]:
    return list(dict.fromkeys(keys))


def display_label(key: tuple[type, Any]) -> str:
    return f"{key[0].__name__}:{key[1]}"
