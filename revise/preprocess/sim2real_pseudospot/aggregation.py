"""Real-cell aggregation into pseudo-spots."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse
from sklearn.neighbors import kneighbors_graph


def _nearest_neighbor_distance(coordinates: np.ndarray) -> float:
    if coordinates.shape[0] < 2:
        raise ValueError("At least two cells are required to construct pseudo-spots.")
    graph = kneighbors_graph(
        coordinates, n_neighbors=1, mode="distance", include_self=False
    )
    return float(graph.data.min())


def aggregate_real_cells(cells: AnnData, *, spot_size: int) -> tuple[AnnData, pd.DataFrame]:
    """Sum real-cell expression into the legacy square pseudo-spot grid."""
    if spot_size <= 0:
        raise ValueError("spot_size must be positive.")
    if "spatial" not in cells.obsm:
        raise ValueError("Input cells must provide obsm['spatial'] coordinates.")

    coordinates = np.asarray(cells.obsm["spatial"], dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 2:
        raise ValueError("obsm['spatial'] must have exactly two coordinate columns.")

    cell_distance = _nearest_neighbor_distance(coordinates)
    xmin, ymin = coordinates.min(axis=0) - 0.5 * cell_distance
    xmax, ymax = coordinates.max(axis=0) + 0.5 * cell_distance
    nx = int(math.ceil(((xmax - xmin) + cell_distance) / spot_size))
    ny = int(math.ceil(((ymax - ymin) + cell_distance) / spot_size))

    grid = ((coordinates - np.array([xmin, ymin])) / spot_size).astype(int)
    spot_indices = grid[:, 0] * ny + grid[:, 1]
    n_spots = nx * ny
    assignment = sparse.csr_matrix(
        (np.ones(cells.n_obs), (spot_indices, np.arange(cells.n_obs))),
        shape=(n_spots, cells.n_obs),
    )
    counts = sparse.csr_matrix(assignment @ cells.X)

    spot_names = [f"SPOT_{i}_{j}" for i in range(nx) for j in range(ny)]
    spot_coordinates = np.array(
        [
            (xmin + spot_size * (i + 0.5), ymin + spot_size * (j + 0.5))
            for i in range(nx)
            for j in range(ny)
        ]
    )
    cell_ids = (
        cells.obs["cell_id"].astype(str).tolist()
        if "cell_id" in cells.obs
        else cells.obs_names.astype(str).tolist()
    )
    mapping: dict[str, list[str] | None] = {}
    for index, name in enumerate(spot_names):
        member_indices = np.flatnonzero(spot_indices == index)
        mapping[name] = (
            [cell_ids[cell_index] for cell_index in member_indices]
            if member_indices.size
            else None
        )

    spots = AnnData(
        X=counts,
        obs=pd.DataFrame(index=spot_names),
        var=cells.var.copy(),
    )
    spots.obsm["spatial"] = spot_coordinates
    spots.uns["all_cells_in_spot"] = mapping
    nonempty = np.asarray(counts.getnnz(axis=1) > 0) | np.asarray(assignment.getnnz(axis=1) > 0)
    spots = spots[nonempty].copy()

    distribution = pd.Series(
        [len(members) if members is not None else 0 for members in mapping.values()]
    ).value_counts().sort_index().to_frame("count")
    distribution.index.name = None
    return spots, distribution
