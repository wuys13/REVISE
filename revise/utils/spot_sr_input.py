from __future__ import annotations

from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse
from scipy.spatial import cKDTree


ALL_CELLS_IN_SPOT_KEY = "all_cells_in_spot"
ESTIMATED_CELL_COUNT_COL = "estimated_cell_count"


def ensure_all_cells_in_spot(
    st_adata: AnnData,
    *,
    logger=None,
    real_adata: Optional[AnnData] = None,
    key: str = ALL_CELLS_IN_SPOT_KEY,
    target_median_cells: int = 4,
    min_cells_per_spot: int = 1,
    max_cells_per_spot: int = 12,
) -> AnnData:
    """Ensure a spot-to-virtual-cell mapping exists for sST inputs.

    The user-facing sST input contract stays small: users provide a
    spot-level AnnData plus a scRNA-seq reference. REVISE still needs a list of
    virtual cells inside each spot before spot expression can be redistributed.
    This helper therefore implements the fallback in one place:

    1. If ``st_adata.uns[key]`` already exists, validate that every active spot
       maps to a list-like collection of cell ids. Existing mappings are never
       overwritten because user-provided segmentation/deconvolution is more
       informative than any default heuristic.
    2. In benchmark mode, when ground-truth cell-level coordinates are
       available, assign real benchmark cell ids to the nearest spot. This keeps
       evaluation indices meaningful and avoids fabricating ids that cannot be
       matched back to the benchmark ground truth.
    3. In normal application mode, estimate the number of virtual cells from
       each spot's transcript count. This is intentionally a conservative
       allocation heuristic, not a claim that true cell boundaries are known.
    """

    if key in st_adata.uns and st_adata.uns[key] is not None:
        st_adata.uns[key] = _validate_all_cells_in_spot(
            st_adata.uns[key],
            spot_names=st_adata.obs_names.astype(str),
            key=key,
        )
        _write_estimated_cell_count(st_adata, st_adata.uns[key])
        return st_adata

    if real_adata is not None and _has_spatial_coordinates(real_adata):
        mapping, radius = _build_mapping_from_nearest_ground_truth_cells(st_adata, real_adata)
        st_adata.uns[key] = mapping
        _write_estimated_cell_count(st_adata, mapping)
        if logger is not None:
            logger.warning(
                "[spot-sr-input] st_adata.uns['%s'] is missing; inferred benchmark "
                "spot-cell mapping from nearest ground-truth cell coordinates "
                "(spots=%d, mapped_cells=%d, max_distance=%.4f).",
                key,
                st_adata.n_obs,
                sum(len(v) for v in mapping.values()),
                radius,
            )
        return st_adata

    mapping = _build_mapping_from_transcript_counts(
        st_adata,
        target_median_cells=target_median_cells,
        min_cells_per_spot=min_cells_per_spot,
        max_cells_per_spot=max_cells_per_spot,
    )
    st_adata.uns[key] = mapping
    _write_estimated_cell_count(st_adata, mapping)
    if logger is not None:
        counts = st_adata.obs[ESTIMATED_CELL_COUNT_COL]
        logger.warning(
            "[spot-sr-input] st_adata.uns['%s'] is missing; generated default "
            "virtual-cell ids from spot transcript counts "
            "(target_median_cells=%d, min=%d, median=%.2f, max=%d).",
            key,
            target_median_cells,
            int(counts.min()),
            float(np.median(counts)),
            int(counts.max()),
        )
    return st_adata


def _validate_all_cells_in_spot(
    raw_mapping: Mapping[Any, Any],
    *,
    spot_names: Iterable[str],
    key: str,
) -> Dict[str, List[str]]:
    if not isinstance(raw_mapping, Mapping):
        raise TypeError(
            f"st_adata.uns['{key}'] must be a mapping from spot id to a list of cell ids"
        )

    normalized: Dict[str, List[str]] = {}
    for spot, cells in raw_mapping.items():
        spot_id = str(spot)
        if isinstance(cells, str):
            # A bare string is iterable in Python and would otherwise be split
            # into characters by downstream code. Failing early makes malformed
            # user input obvious instead of silently creating bogus cells.
            raise TypeError(
                f"st_adata.uns['{key}'][{spot_id!r}] must be a list-like "
                "collection of cell ids, not a string"
            )
        if cells is None:
            cell_ids: List[str] = []
        else:
            cell_ids = [str(cell_id) for cell_id in list(cells)]
        normalized[spot_id] = cell_ids

    missing = [spot for spot in spot_names if spot not in normalized]
    if missing:
        preview = ", ".join(missing[:5])
        raise KeyError(
            f"st_adata.uns['{key}'] is missing mappings for {len(missing)} "
            f"active spots; examples: {preview}"
        )

    empty = [spot for spot in spot_names if len(normalized[spot]) == 0]
    if empty:
        preview = ", ".join(empty[:5])
        raise ValueError(
            f"st_adata.uns['{key}'] contains empty cell lists for {len(empty)} "
            f"active spots; examples: {preview}"
        )
    return normalized


def _build_mapping_from_transcript_counts(
    st_adata: AnnData,
    *,
    target_median_cells: int,
    min_cells_per_spot: int,
    max_cells_per_spot: int,
) -> Dict[str, List[str]]:
    if target_median_cells <= 0:
        raise ValueError("target_median_cells must be positive")
    if min_cells_per_spot <= 0:
        raise ValueError("min_cells_per_spot must be positive")
    if max_cells_per_spot < min_cells_per_spot:
        raise ValueError("max_cells_per_spot must be >= min_cells_per_spot")

    counts = _get_transcript_counts(st_adata)
    positive = counts[counts > 0]
    if positive.size == 0:
        raise ValueError(
            "Cannot infer all_cells_in_spot because every spot has zero transcript counts"
        )

    # Use the positive-count median rather than the mean so a few high-depth
    # spots do not inflate the scale and collapse most spots to one virtual
    # cell. The default maps the median expressed spot to four virtual cells,
    # then clamps the result to keep the fallback computationally bounded.
    scale = max(float(np.median(positive)) / float(target_median_cells), 1.0)
    estimated = np.rint(counts / scale).astype(int)
    estimated = np.clip(estimated, int(min_cells_per_spot), int(max_cells_per_spot))

    mapping: Dict[str, List[str]] = {}
    for spot_id, n_cells in zip(st_adata.obs_names.astype(str), estimated):
        mapping[str(spot_id)] = [
            f"{spot_id}_vc_{idx:02d}" for idx in range(int(n_cells))
        ]
    return mapping


def _build_mapping_from_nearest_ground_truth_cells(
    st_adata: AnnData,
    real_adata: AnnData,
) -> Tuple[Dict[str, List[str]], float]:
    spot_names = st_adata.obs_names.astype(str).to_numpy()
    spot_xy = _get_spatial_coordinates(st_adata)
    real_xy = _get_spatial_coordinates(real_adata)

    if spot_xy.shape[0] == 0:
        raise ValueError(
            "Cannot infer all_cells_in_spot from ground truth because st_adata has no spots"
        )
    if real_xy.shape[0] == 0:
        raise ValueError(
            "Cannot infer all_cells_in_spot from ground truth because real_adata has no cells"
        )

    tree = cKDTree(spot_xy)
    distances, indices = tree.query(real_xy, k=1)
    distances = np.asarray(distances, dtype=np.float64)
    indices = np.asarray(indices, dtype=np.int64)
    max_distance = _estimate_spot_assignment_radius(spot_xy)

    mapping: Dict[str, List[str]] = {str(spot): [] for spot in spot_names}
    if "cell_id" in real_adata.obs:
        raw_ids = real_adata.obs["cell_id"]
        real_ids = raw_ids.astype(str).to_numpy()
        if (
            raw_ids.isna().any()
            or pd.Series(real_ids).str.strip().eq("").any()
            or pd.Index(real_ids).duplicated().any()
        ):
            raise ValueError(
                "Ground-truth obs['cell_id'] must contain unique non-null values"
            )
    else:
        real_ids = real_adata.obs_names.astype(str).to_numpy()
    for cell_id, spot_idx, distance in zip(real_ids, indices, distances):
        if float(distance) <= max_distance:
            mapping[str(spot_names[int(spot_idx)])].append(str(cell_id))

    # A sampled benchmark can contain isolated spots whose nearest true cells
    # fall outside the inferred radius. Downstream code requires every active
    # spot to have at least one virtual cell, so fall back to one deterministic
    # synthetic id for those rare spots. These synthetic ids will not contribute
    # to benchmark metrics because they do not intersect the ground truth index.
    for spot in spot_names:
        spot_id = str(spot)
        if not mapping[spot_id]:
            mapping[spot_id] = [f"{spot_id}_vc_00"]

    return mapping, max_distance


def _estimate_spot_assignment_radius(spot_xy: np.ndarray) -> float:
    if spot_xy.shape[0] < 2:
        return float("inf")

    tree = cKDTree(spot_xy)
    distances, _indices = tree.query(spot_xy, k=2)
    nearest = np.asarray(distances[:, 1], dtype=np.float64)
    nearest = nearest[np.isfinite(nearest) & (nearest > 0)]
    if nearest.size == 0:
        return float("inf")

    # The Sim2Real spot generator uses square/circular spot neighborhoods, and
    # real Visium spots are circular. A 0.75 * center-spacing radius is a simple
    # conservative default: broad enough to recover edge cells in square spots,
    # but still tight enough to avoid assigning distant cells to a sampled spot.
    return float(np.median(nearest) * 0.75)


def _get_transcript_counts(st_adata: AnnData) -> np.ndarray:
    if "transcript_counts" in st_adata.obs:
        counts = np.asarray(st_adata.obs["transcript_counts"], dtype=np.float64)
    elif "total_counts" in st_adata.obs:
        counts = np.asarray(st_adata.obs["total_counts"], dtype=np.float64)
    else:
        row_sums = st_adata.X.sum(axis=1)
        if sparse.issparse(row_sums) and hasattr(row_sums, "A1"):
            counts = row_sums.A1
        else:
            counts = np.asarray(row_sums).ravel()
        st_adata.obs["transcript_counts"] = counts
    if counts.shape[0] != st_adata.n_obs:
        raise ValueError("transcript count vector length does not match st_adata.n_obs")
    return counts.astype(np.float64, copy=False)


def _write_estimated_cell_count(st_adata: AnnData, mapping: Mapping[str, List[str]]) -> None:
    counts = pd.Series(
        [len(mapping[str(spot)]) for spot in st_adata.obs_names.astype(str)],
        index=st_adata.obs_names,
        dtype="int64",
    )
    st_adata.obs[ESTIMATED_CELL_COUNT_COL] = counts


def _has_spatial_coordinates(adata: AnnData) -> bool:
    if "spatial" in adata.obsm:
        return True
    return {"x", "y"}.issubset(set(adata.obs.columns))


def _get_spatial_coordinates(adata: AnnData) -> np.ndarray:
    if "spatial" in adata.obsm:
        coords = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    elif {"x", "y"}.issubset(set(adata.obs.columns)):
        coords = adata.obs.loc[:, ["x", "y"]].to_numpy(dtype=np.float64)
    else:
        raise KeyError("AnnData must contain obsm['spatial'] or obs[['x', 'y']]")
    if coords.ndim != 2 or coords.shape[1] < 2:
        raise ValueError("spatial coordinates must have shape (n_obs, >=2)")
    if not np.all(np.isfinite(coords[:, :2])):
        raise ValueError("spatial coordinates must contain only finite values")
    return coords[:, :2]
