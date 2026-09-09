from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple, Union

import anndata as ad
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.spatial import cKDTree

from revise.utils.spot_sr_input import (
    ALL_CELLS_IN_SPOT_KEY,
    CELL_LOCATIONS_KEY,
    ensure_all_cells_in_spot,
    validate_cell_locations,
)


HISTOLOGY_PRIOR_KEY = "revise_histology_prior"

PathLike = Union[str, Path]

_SPOT_ID_COLUMNS = (
    "spot_id",
    "spot",
    "spot_name",
    "barcode",
    "barcodes",
    "obs_name",
    "id",
)
_COORDINATE_COLUMN_PAIRS = (
    ("x", "y"),
    ("spatial_x", "spatial_y"),
    ("center_x", "center_y"),
    ("image_x", "image_y"),
    ("pixel_x", "pixel_y"),
    ("pxl_col_in_fullres", "pxl_row_in_fullres"),
)


@dataclass
class HistologyPriorResult:
    """Outputs from histology preprocessing.

    ``mapping`` is written to ``st_adata.uns["all_cells_in_spot"]``. The cell
    and spot tables are returned for inspection or external reporting. The core
    REVISE engine consumes the mapping and, when available, the standardized
    segmented-cell centers in ``st_adata.uns["revise_cell_locations"]``.
    """

    mapping: Dict[str, List[str]]
    cell_table: pd.DataFrame
    spot_table: pd.DataFrame
    provenance: Dict[str, Any]


def read_histology_image(image_path: PathLike) -> np.ndarray:
    """Read a raw histology image from disk and validate its shape."""

    from skimage import io

    path = Path(image_path)
    if not path.exists():
        raise FileNotFoundError(f"Histology image does not exist: {path}")
    image = np.asarray(io.imread(path))
    if image.ndim not in (2, 3):
        raise ValueError(
            f"Histology image must be 2D grayscale or 3D multichannel; got shape {image.shape}"
        )
    if image.shape[0] == 0 or image.shape[1] == 0:
        raise ValueError(f"Histology image has an empty spatial dimension: {image.shape}")
    return image


def read_labeled_mask(mask_path: PathLike) -> np.ndarray:
    """Read a labeled segmentation mask where 0 is background and labels are cells."""

    from skimage import io

    path = Path(mask_path)
    if not path.exists():
        raise FileNotFoundError(f"Histology segmentation mask does not exist: {path}")
    return _coerce_labeled_mask(io.imread(path), source=str(path))


def extract_labeled_mask_cells(
    mask: np.ndarray,
    *,
    image: Optional[np.ndarray] = None,
    min_area: float = 1.0,
    cell_id_prefix: str = "histology_cell",
) -> pd.DataFrame:
    """Extract image-derived cell centroids and morphology summaries.

    Parameters
    ----------
    mask:
        2D integer label image. Label 0 is treated as background; every positive
        label is one segmented cell or nucleus.
    image:
        Optional raw histology image aligned to ``mask``. When provided, the
        output includes per-cell mean intensity from the grayscale image.
    min_area:
        Minimum labeled pixel area to retain.
    cell_id_prefix:
        Prefix used to create stable cell ids from mask labels.
    """

    from skimage import measure

    if min_area <= 0:
        raise ValueError("min_area must be positive")
    label_mask = _coerce_labeled_mask(mask, source="mask")
    intensity_image = None
    properties: Tuple[str, ...] = (
        "label",
        "area",
        "centroid",
        "bbox",
    )
    if image is not None:
        _validate_image_mask_shapes(image, label_mask)
        intensity_image = _as_grayscale_image(image)
        properties = properties + ("mean_intensity",)

    table = pd.DataFrame(
        measure.regionprops_table(
            label_mask,
            intensity_image=intensity_image,
            properties=properties,
        )
    )
    if table.empty:
        raise ValueError("Segmentation mask contains no positive cell labels")

    table = table.rename(
        columns={
            "label": "source_label",
            "centroid-1": "x",
            "centroid-0": "y",
            "bbox-0": "bbox_ymin",
            "bbox-1": "bbox_xmin",
            "bbox-2": "bbox_ymax",
            "bbox-3": "bbox_xmax",
        }
    )
    table = table.loc[table["area"].astype(float) >= float(min_area)].copy()
    if table.empty:
        raise ValueError(
            "No segmented cells remain after min_area filtering "
            f"(min_area={float(min_area)})"
        )

    table["source_label"] = table["source_label"].astype(int)
    table["cell_id"] = [
        f"{cell_id_prefix}_{label}" for label in table["source_label"].tolist()
    ]
    ordered_columns = [
        "cell_id",
        "source_label",
        "x",
        "y",
        "area",
        "bbox_xmin",
        "bbox_ymin",
        "bbox_xmax",
        "bbox_ymax",
    ]
    if "mean_intensity" in table.columns:
        ordered_columns.append("mean_intensity")
    return table.loc[:, ordered_columns].reset_index(drop=True)


def load_spot_table(
    spots_path: PathLike,
    *,
    spot_id_col: Optional[str] = None,
    x_col: Optional[str] = None,
    y_col: Optional[str] = None,
) -> pd.DataFrame:
    """Load spot image coordinates from a CSV table.

    Supported coordinate pairs include ``x/y``, ``spatial_x/spatial_y``,
    ``center_x/center_y``, and 10x Genomics Visium
    ``pxl_col_in_fullres/pxl_row_in_fullres``.
    """

    path = Path(spots_path)
    if not path.exists():
        raise FileNotFoundError(f"Spot table does not exist: {path}")
    table = pd.read_csv(path)
    if table.empty:
        raise ValueError(f"Spot table is empty: {path}")

    resolved_spot_id_col = spot_id_col or _find_column(table, _SPOT_ID_COLUMNS)
    if resolved_spot_id_col is None and table.columns[0].startswith("Unnamed:"):
        resolved_spot_id_col = table.columns[0]
    if resolved_spot_id_col is None:
        raise ValueError(
            "Spot table must contain a spot id column such as spot_id, barcode, or spot_name"
        )

    if x_col is None or y_col is None:
        resolved_x_col, resolved_y_col = _find_coordinate_pair(table)
    else:
        resolved_x_col, resolved_y_col = x_col, y_col
    _require_columns(table, [resolved_spot_id_col, resolved_x_col, resolved_y_col])

    spots = pd.DataFrame(
        {
            "spot_id": table[resolved_spot_id_col].astype(str),
            "x": pd.to_numeric(table[resolved_x_col], errors="raise"),
            "y": pd.to_numeric(table[resolved_y_col], errors="raise"),
        }
    )
    return _validate_spot_table(spots)


def spot_table_from_adata(st_adata: AnnData) -> pd.DataFrame:
    """Build a spot coordinate table from ``obsm["spatial"]`` or ``obs[["x", "y"]]``."""

    if "spatial" in st_adata.obsm:
        coords = np.asarray(st_adata.obsm["spatial"], dtype=np.float64)
    elif {"x", "y"}.issubset(set(st_adata.obs.columns)):
        coords = st_adata.obs.loc[:, ["x", "y"]].to_numpy(dtype=np.float64)
    else:
        raise KeyError(
            "st_adata must contain obsm['spatial'] or obs[['x', 'y']] when --spots is omitted"
        )
    if coords.ndim != 2 or coords.shape[1] < 2:
        raise ValueError("Spatial coordinates must have shape (n_spots, >=2)")
    spots = pd.DataFrame(
        {
            "spot_id": st_adata.obs_names.astype(str),
            "x": coords[:, 0],
            "y": coords[:, 1],
        }
    )
    return _validate_spot_table(spots)


def build_spot_cell_mapping(
    cell_table: pd.DataFrame,
    spot_table: pd.DataFrame,
    *,
    spot_radius: Optional[float] = None,
) -> Dict[str, List[str]]:
    """Assign each image-derived cell centroid to its nearest spot.

    If ``spot_radius`` is provided, cells farther than that radius are left
    unassigned. The returned mapping includes every spot in ``spot_table``;
    downstream completion fills any empty active spots with deterministic
    virtual ids.
    """

    _require_columns(cell_table, ["cell_id", "x", "y"])
    spots = _validate_spot_table(spot_table)
    cells = _validate_cell_table(cell_table)
    if spot_radius is not None and spot_radius <= 0:
        raise ValueError("spot_radius must be positive when provided")

    mapping: Dict[str, List[str]] = {spot_id: [] for spot_id in spots["spot_id"]}
    spot_xy = spots.loc[:, ["x", "y"]].to_numpy(dtype=np.float64)
    cell_xy = cells.loc[:, ["x", "y"]].to_numpy(dtype=np.float64)
    tree = cKDTree(spot_xy)
    distances, indices = tree.query(cell_xy, k=1)

    for cell_id, distance, spot_idx in zip(cells["cell_id"], distances, indices):
        if spot_radius is not None and float(distance) > float(spot_radius):
            continue
        spot_id = str(spots.iloc[int(spot_idx)]["spot_id"])
        mapping[spot_id].append(str(cell_id))

    return mapping


def build_histology_prior_from_tables(
    st_adata: AnnData,
    cell_table: pd.DataFrame,
    spot_table: pd.DataFrame,
    *,
    image: Optional[np.ndarray] = None,
    mask: Optional[np.ndarray] = None,
    image_path: Optional[PathLike] = None,
    mask_path: Optional[PathLike] = None,
    spots_path: Optional[PathLike] = None,
    spot_radius: Optional[float] = None,
    mapping_key: str = ALL_CELLS_IN_SPOT_KEY,
    target_median_cells: int = 4,
    min_cells_per_spot: int = 1,
    max_cells_per_spot: int = 12,
) -> HistologyPriorResult:
    """Write an image-derived spot-cell prior into an AnnData object.

    This mutates ``st_adata`` in place by setting ``uns[mapping_key]`` and
    ``uns["revise_histology_prior"]``. Spots not covered by the segmentation
    are filled with the same transcript-count virtual-cell fallback used by the
    sc-SVC sr-mode input service, and the number of fallback spots is recorded in
    provenance.
    """

    _validate_spot_coordinates_against_adata(st_adata, spot_table)

    image_mapping = build_spot_cell_mapping(
        cell_table,
        spot_table,
        spot_radius=spot_radius,
    )
    complete_mapping, fallback_spots = _complete_mapping_for_adata(
        st_adata,
        image_mapping,
        mapping_key=mapping_key,
        target_median_cells=target_median_cells,
        min_cells_per_spot=min_cells_per_spot,
        max_cells_per_spot=max_cells_per_spot,
    )

    st_adata.uns[mapping_key] = complete_mapping
    ensure_all_cells_in_spot(
        st_adata,
        key=mapping_key,
        target_median_cells=target_median_cells,
        min_cells_per_spot=min_cells_per_spot,
        max_cells_per_spot=max_cells_per_spot,
    )
    cell_lookup = _validate_cell_table(cell_table).set_index("cell_id")
    cell_to_spot = {
        str(cell_id): str(spot_id)
        for spot_id, cell_ids in complete_mapping.items()
        for cell_id in cell_ids
        if str(cell_id) in cell_lookup.index
    }
    locations = cell_lookup.loc[list(cell_to_spot), ["x", "y"]].copy()
    locations.insert(
        0,
        "spot_name",
        [cell_to_spot[cell_id] for cell_id in locations.index],
    )
    st_adata.uns[CELL_LOCATIONS_KEY] = validate_cell_locations(
        locations,
        all_cells_in_spot=complete_mapping,
    )

    active_spots = set(st_adata.obs_names.astype(str))
    histology_mapped_spots = [
        spot_id
        for spot_id, cells in image_mapping.items()
        if spot_id in active_spots and len(cells) > 0
    ]
    mapped_cell_count = sum(
        len(cells)
        for spot_id, cells in image_mapping.items()
        if spot_id in active_spots
    )
    provenance = {
        "mapping_key": mapping_key,
        "n_active_spots": int(st_adata.n_obs),
        "n_spots_in_coordinate_table": int(len(spot_table)),
        "n_segmented_cells": int(len(cell_table)),
        "n_image_mapped_spots": int(len(histology_mapped_spots)),
        "n_image_mapped_cells": int(mapped_cell_count),
        "n_fallback_spots": int(len(fallback_spots)),
        "cell_locations_key": CELL_LOCATIONS_KEY,
        "n_cell_locations": int(len(locations)),
        "n_spot_table_not_in_adata": int(
            len(set(spot_table["spot_id"].astype(str)) - active_spots)
        ),
        "fallback_policy": (
            "Spots without image-derived cells are filled with deterministic "
            "transcript-count virtual-cell ids from revise.utils.spot_sr_input."
        ),
    }
    if image_path is not None:
        provenance["image_path"] = str(image_path)
    if mask_path is not None:
        provenance["mask_path"] = str(mask_path)
    if spots_path is not None:
        provenance["spots_path"] = str(spots_path)
    if spot_radius is not None:
        provenance["spot_radius"] = float(spot_radius)
    if image is not None:
        provenance["image_summary"] = _image_summary(image)
    if mask is not None:
        provenance["mask_summary"] = _mask_summary(mask)
    st_adata.uns[HISTOLOGY_PRIOR_KEY] = provenance
    return HistologyPriorResult(
        mapping=complete_mapping,
        cell_table=cell_table.copy(),
        spot_table=spot_table.copy(),
        provenance=provenance,
    )


def build_histology_prior_h5ad(
    *,
    st_h5ad_path: PathLike,
    mask_path: PathLike,
    out_h5ad_path: PathLike,
    image_path: Optional[PathLike] = None,
    spots_path: Optional[PathLike] = None,
    spot_radius: Optional[float] = None,
    min_area: float = 1.0,
    cell_id_prefix: str = "histology_cell",
    mapping_key: str = ALL_CELLS_IN_SPOT_KEY,
    target_median_cells: int = 4,
    min_cells_per_spot: int = 1,
    max_cells_per_spot: int = 12,
) -> HistologyPriorResult:
    """Build and write histology-derived priors for an ST H5AD file."""

    st_h5ad = Path(st_h5ad_path)
    if not st_h5ad.exists():
        raise FileNotFoundError(f"ST H5AD does not exist: {st_h5ad}")

    st_adata = ad.read_h5ad(st_h5ad)
    mask = read_labeled_mask(mask_path)
    image = read_histology_image(image_path) if image_path is not None else None
    cell_table = extract_labeled_mask_cells(
        mask,
        image=image,
        min_area=min_area,
        cell_id_prefix=cell_id_prefix,
    )
    if spots_path is not None:
        spot_table = load_spot_table(spots_path)
    else:
        spot_table = spot_table_from_adata(st_adata)

    result = build_histology_prior_from_tables(
        st_adata,
        cell_table,
        spot_table,
        image=image,
        mask=mask,
        image_path=image_path,
        mask_path=mask_path,
        spots_path=spots_path,
        spot_radius=spot_radius,
        mapping_key=mapping_key,
        target_median_cells=target_median_cells,
        min_cells_per_spot=min_cells_per_spot,
        max_cells_per_spot=max_cells_per_spot,
    )
    out_h5ad = Path(out_h5ad_path)
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    st_adata.write_h5ad(out_h5ad)
    return result


def _complete_mapping_for_adata(
    st_adata: AnnData,
    image_mapping: Mapping[str, List[str]],
    *,
    mapping_key: str,
    target_median_cells: int,
    min_cells_per_spot: int,
    max_cells_per_spot: int,
) -> Tuple[Dict[str, List[str]], List[str]]:
    fallback_adata = st_adata.copy()
    fallback_adata.uns.pop(mapping_key, None)
    ensure_all_cells_in_spot(
        fallback_adata,
        key=mapping_key,
        target_median_cells=target_median_cells,
        min_cells_per_spot=min_cells_per_spot,
        max_cells_per_spot=max_cells_per_spot,
    )
    fallback_mapping = fallback_adata.uns[mapping_key]

    complete_mapping: Dict[str, List[str]] = {}
    fallback_spots: List[str] = []
    for spot_id in st_adata.obs_names.astype(str):
        image_cells = [str(cell_id) for cell_id in image_mapping.get(str(spot_id), [])]
        if image_cells:
            complete_mapping[str(spot_id)] = image_cells
        else:
            complete_mapping[str(spot_id)] = [
                str(cell_id) for cell_id in fallback_mapping[str(spot_id)]
            ]
            fallback_spots.append(str(spot_id))
    return complete_mapping, fallback_spots


def _coerce_labeled_mask(mask: np.ndarray, *, source: str) -> np.ndarray:
    array = np.asarray(mask)
    if array.ndim == 3 and array.shape[-1] == 1:
        array = array[:, :, 0]
    if array.ndim != 2:
        raise ValueError(
            f"{source} must be a 2D labeled mask with 0 as background; got shape {array.shape}"
        )
    if not np.issubdtype(array.dtype, np.integer):
        if not np.all(np.isfinite(array)):
            raise ValueError(f"{source} contains non-finite mask labels")
        rounded = np.rint(array)
        if not np.allclose(array, rounded):
            raise ValueError(f"{source} must contain integer labels")
        array = rounded
    label_mask = array.astype(np.int64, copy=False)
    if np.any(label_mask < 0):
        raise ValueError(f"{source} contains negative labels")
    if not np.any(label_mask > 0):
        raise ValueError(f"{source} contains no positive cell labels")
    return label_mask


def _as_grayscale_image(image: np.ndarray) -> np.ndarray:
    array = np.asarray(image)
    if array.ndim == 2:
        return array.astype(np.float64, copy=False)
    if array.ndim != 3:
        raise ValueError(f"Histology image must be 2D or 3D; got shape {array.shape}")
    if array.shape[-1] == 1:
        return array[:, :, 0].astype(np.float64, copy=False)
    if array.shape[-1] >= 3:
        rgb = array[:, :, :3].astype(np.float64, copy=False)
        return rgb[:, :, 0] * 0.299 + rgb[:, :, 1] * 0.587 + rgb[:, :, 2] * 0.114
    raise ValueError(f"Unsupported histology image channel shape: {array.shape}")


def _validate_image_mask_shapes(image: np.ndarray, mask: np.ndarray) -> None:
    if image.shape[0] != mask.shape[0] or image.shape[1] != mask.shape[1]:
        raise ValueError(
            "Histology image and segmentation mask must share spatial dimensions; "
            f"got image={image.shape[:2]} mask={mask.shape[:2]}"
        )


def _validate_spot_table(spot_table: pd.DataFrame) -> pd.DataFrame:
    _require_columns(spot_table, ["spot_id", "x", "y"])
    spots = spot_table.loc[:, ["spot_id", "x", "y"]].copy()
    spots["spot_id"] = spots["spot_id"].astype(str)
    spots["x"] = pd.to_numeric(spots["x"], errors="raise")
    spots["y"] = pd.to_numeric(spots["y"], errors="raise")
    if spots.empty:
        raise ValueError("Spot table contains no rows")
    if spots["spot_id"].duplicated().any():
        duplicated = spots.loc[spots["spot_id"].duplicated(), "spot_id"].iloc[0]
        raise ValueError(f"Spot table contains duplicated spot id: {duplicated}")
    if not np.isfinite(spots.loc[:, ["x", "y"]].to_numpy(dtype=np.float64)).all():
        raise ValueError("Spot table contains non-finite coordinates")
    return spots.reset_index(drop=True)


def _validate_cell_table(cell_table: pd.DataFrame) -> pd.DataFrame:
    _require_columns(cell_table, ["cell_id", "x", "y"])
    cells = cell_table.loc[:, ["cell_id", "x", "y"]].copy()
    cells["cell_id"] = cells["cell_id"].astype(str)
    cells["x"] = pd.to_numeric(cells["x"], errors="raise")
    cells["y"] = pd.to_numeric(cells["y"], errors="raise")
    if cells.empty:
        raise ValueError("Cell table contains no rows")
    if cells["cell_id"].duplicated().any():
        duplicated = cells.loc[cells["cell_id"].duplicated(), "cell_id"].iloc[0]
        raise ValueError(f"Cell table contains duplicated cell id: {duplicated}")
    if not np.isfinite(cells.loc[:, ["x", "y"]].to_numpy(dtype=np.float64)).all():
        raise ValueError("Cell table contains non-finite coordinates")
    return cells.reset_index(drop=True)


def _validate_spot_coordinates_against_adata(
    st_adata: AnnData,
    spot_table: pd.DataFrame,
) -> None:
    """Reject mixed coordinate frames before persisting segmented centers."""

    spots = _validate_spot_table(spot_table).set_index("spot_id")
    expected_spots = spot_table_from_adata(st_adata).set_index("spot_id")
    shared = [spot_id for spot_id in expected_spots.index if spot_id in spots.index]
    if not shared:
        return
    supplied = spots.loc[shared, ["x", "y"]].to_numpy(dtype=np.float64)
    expected = expected_spots.loc[shared, ["x", "y"]].to_numpy(dtype=np.float64)
    if not np.allclose(supplied, expected, rtol=1e-5, atol=1e-8):
        raise ValueError(
            "Spot-table x/y must match st_adata.obsm['spatial'] for shared spots "
            "so segmented and fallback cell centers use one coordinate frame"
        )


def _find_column(table: pd.DataFrame, candidates: Tuple[str, ...]) -> Optional[str]:
    columns_by_lower = {str(column).lower(): column for column in table.columns}
    for candidate in candidates:
        column = columns_by_lower.get(candidate.lower())
        if column is not None:
            return str(column)
    return None


def _find_coordinate_pair(table: pd.DataFrame) -> Tuple[str, str]:
    columns_by_lower = {str(column).lower(): str(column) for column in table.columns}
    for x_candidate, y_candidate in _COORDINATE_COLUMN_PAIRS:
        x_col = columns_by_lower.get(x_candidate.lower())
        y_col = columns_by_lower.get(y_candidate.lower())
        if x_col is not None and y_col is not None:
            return x_col, y_col
    raise ValueError(
        "Spot table must contain coordinate columns such as x/y or "
        "pxl_col_in_fullres/pxl_row_in_fullres"
    )


def _require_columns(table: pd.DataFrame, columns: List[str]) -> None:
    missing = [column for column in columns if column not in table.columns]
    if missing:
        raise KeyError(f"Missing required columns: {', '.join(missing)}")


def _image_summary(image: Optional[np.ndarray]) -> Optional[Dict[str, Any]]:
    if image is None:
        return None
    array = np.asarray(image)
    return {
        "shape": [int(dim) for dim in array.shape],
        "ndim": int(array.ndim),
        "dtype": str(array.dtype),
        "min": float(np.nanmin(array)),
        "max": float(np.nanmax(array)),
    }


def _mask_summary(mask: Optional[np.ndarray]) -> Optional[Dict[str, Any]]:
    if mask is None:
        return None
    label_mask = _coerce_labeled_mask(mask, source="mask")
    return {
        "shape": [int(dim) for dim in label_mask.shape],
        "dtype": str(label_mask.dtype),
        "n_labels": int(np.unique(label_mask[label_mask > 0]).size),
    }
