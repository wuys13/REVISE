from __future__ import annotations

import importlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional, Tuple, Union

import numpy as np
import pandas as pd
from anndata import AnnData


class SpatialDataInputError(ValueError):
    """Raised when a SpatialData object cannot be normalized for REVISE."""


@dataclass
class SpatialDataTableResult:
    adata: AnnData
    sdata: Any
    metadata: Dict[str, Any] = field(default_factory=dict)


class SpatialDataService:
    """Read SpatialData stores and expose REVISE's AnnData contract.

    SpatialData is used here as an internal spatial service layer. The selected
    table is copied into AnnData because the current reconstruction stack
    expects matrix-like AnnData inputs, while the original SpatialData object is
    kept as optional provenance for future write-back and visualization.
    """

    OBS_COORDINATE_PAIRS: Tuple[Tuple[str, str], ...] = (
        ("x", "y"),
        ("X", "Y"),
        ("spatial_x", "spatial_y"),
        ("global_x", "global_y"),
        ("center_x", "center_y"),
        ("x_centroid", "y_centroid"),
        ("array_col", "array_row"),
        ("pxl_col_in_fullres", "pxl_row_in_fullres"),
    )

    def __init__(self, logger=None):
        self.logger = logger

    def read_table(
        self,
        path: Union[str, Path],
        *,
        table_key: Optional[str] = None,
        spatial_element_key: Optional[str] = None,
        coordinate_system: str = "global",
    ) -> SpatialDataTableResult:
        spatialdata = self._import_spatialdata()
        read_zarr = getattr(spatialdata, "read_zarr", None)
        if read_zarr is None:
            raise ImportError("Installed spatialdata package does not expose read_zarr().")

        sdata = read_zarr(str(path))
        selected_key, table, selection_reason = self._select_table(sdata, table_key)
        adata = table.copy()
        source_report = {
            "format": "spatialdata",
            "path": str(path),
            "table_key": selected_key,
            "table_selection": selection_reason,
            "coordinate_system": coordinate_system,
        }
        coord_source = self._ensure_spatial_coordinates(
            adata,
            sdata=sdata,
            spatial_element_key=spatial_element_key,
        )
        source_report["spatial_coordinates"] = coord_source
        if spatial_element_key:
            source_report["spatial_element_key"] = spatial_element_key
        adata.uns.setdefault("revise_input", {})["spatialdata"] = source_report
        if self.logger is not None:
            self.logger.info(
                "[input] loaded SpatialData table=%s path=%s coordinates=%s",
                selected_key,
                path,
                coord_source,
            )
        return SpatialDataTableResult(adata=adata, sdata=sdata, metadata=source_report)

    @staticmethod
    def _import_spatialdata():
        try:
            return importlib.import_module("spatialdata")
        except ImportError as exc:  # pragma: no cover - depends on optional env
            raise ImportError(
                "SpatialData input requires optional dependencies. Install with "
                "`python -m pip install \"revise-svc[spatialdata]\"`."
            ) from exc

    def _select_table(self, sdata: Any, table_key: Optional[str]) -> Tuple[str, AnnData, str]:
        tables = getattr(sdata, "tables", None)
        if tables is None:
            raise SpatialDataInputError("SpatialData object has no tables attribute.")
        keys = list(tables.keys())
        if not keys:
            raise SpatialDataInputError("SpatialData object contains no AnnData tables.")

        if table_key:
            if table_key not in tables:
                raise SpatialDataInputError(
                    f"SpatialData table {table_key!r} was requested but available tables are {keys}."
                )
            return table_key, tables[table_key], "explicit"

        if "table" in tables:
            return "table", tables["table"], "default_table_key"

        if len(keys) == 1:
            key = keys[0]
            return key, tables[key], "single_table"

        raise SpatialDataInputError(
            "SpatialData object contains multiple tables. Set io.spatialdata_table "
            f"to one of: {keys}."
        )

    def _ensure_spatial_coordinates(
        self,
        adata: AnnData,
        *,
        sdata: Any,
        spatial_element_key: Optional[str],
    ) -> str:
        if "spatial" in adata.obsm:
            coords = np.asarray(adata.obsm["spatial"])
            if coords.ndim == 2 and coords.shape[0] == adata.n_obs and coords.shape[1] >= 2:
                return "table.obsm['spatial']"
            raise SpatialDataInputError(
                "Selected SpatialData table has obsm['spatial'], but it is not shaped "
                "(n_obs, >=2)."
            )

        obs_source = self._write_spatial_from_obs_columns(adata)
        if obs_source is not None:
            return obs_source

        element_source = self._write_spatial_from_linked_element(
            adata,
            sdata=sdata,
            spatial_element_key=spatial_element_key,
        )
        if element_source is not None:
            return element_source

        raise SpatialDataInputError(
            "Could not infer spatial coordinates for the selected SpatialData table. "
            "Provide table.obsm['spatial'], numeric x/y obs columns, or a linked "
            "Shapes element with matching instance ids."
        )

    def _write_spatial_from_obs_columns(self, adata: AnnData) -> Optional[str]:
        for x_key, y_key in self.OBS_COORDINATE_PAIRS:
            if x_key not in adata.obs or y_key not in adata.obs:
                continue
            try:
                coords = adata.obs[[x_key, y_key]].to_numpy(dtype=np.float64)
            except (TypeError, ValueError):
                continue
            if coords.shape == (adata.n_obs, 2) and np.isfinite(coords).all():
                adata.obsm["spatial"] = coords
                return f"table.obs[{x_key!r}, {y_key!r}]"
        return None

    def _write_spatial_from_linked_element(
        self,
        adata: AnnData,
        *,
        sdata: Any,
        spatial_element_key: Optional[str],
    ) -> Optional[str]:
        element_key, instance_key = self._resolve_region_and_instance_keys(
            adata,
            sdata=sdata,
            spatial_element_key=spatial_element_key,
        )
        if element_key is None:
            return None

        element = self._get_spatial_element(sdata, element_key)
        if element is None or not hasattr(element, "geometry"):
            return None

        coord_df = self._shape_centroid_coordinates(element)
        if coord_df.empty:
            return None

        if instance_key and instance_key in adata.obs:
            lookup_ids = pd.Index(adata.obs[instance_key].astype(str))
        else:
            lookup_ids = pd.Index(adata.obs_names.astype(str))

        missing = lookup_ids.difference(coord_df.index)
        if len(missing) > 0:
            preview = missing[:5].tolist()
            raise SpatialDataInputError(
                f"Linked SpatialData element {element_key!r} is missing {len(missing)} "
                f"table instances, for example {preview}."
            )

        adata.obsm["spatial"] = coord_df.loc[lookup_ids].to_numpy(dtype=np.float64)
        if instance_key:
            return f"shapes[{element_key!r}].centroid via obs[{instance_key!r}]"
        return f"shapes[{element_key!r}].centroid via obs_names"

    def _resolve_region_and_instance_keys(
        self,
        adata: AnnData,
        *,
        sdata: Any,
        spatial_element_key: Optional[str],
    ) -> Tuple[Optional[str], Optional[str]]:
        if spatial_element_key:
            return spatial_element_key, self._get_instance_key_from_table(adata)

        table_keys = self._get_table_annotation_keys(adata)
        if table_keys is None:
            return None, None

        regions, region_key, instance_key = table_keys
        region = self._select_region(regions, adata, region_key)
        if region is None or self._get_spatial_element(sdata, region) is None:
            return None, instance_key
        return region, instance_key

    @staticmethod
    def _get_instance_key_from_table(adata: AnnData) -> Optional[str]:
        for key in ("instance_id", "cell_id", "spot_id", "label"):
            if key in adata.obs:
                return key
        return None

    @staticmethod
    def _get_table_annotation_keys(adata: AnnData) -> Optional[Tuple[Iterable[str], str, str]]:
        try:
            models = importlib.import_module("spatialdata.models")
        except ImportError as exc:
            raise ImportError(
                "SpatialData table metadata requires a working spatialdata "
                "installation; run `python -m pip install "
                "\"revise-svc[spatialdata]\"`."
            ) from exc
        try:
            region, region_key, instance_key = models.get_table_keys(adata)
        except (AttributeError, KeyError, TypeError, ValueError):
            return None

        if isinstance(region, str):
            regions: Iterable[str] = [region]
        else:
            regions = region
        return regions, region_key, instance_key

    @staticmethod
    def _select_region(regions: Iterable[str], adata: AnnData, region_key: str) -> Optional[str]:
        region_list = [str(region) for region in regions]
        if len(region_list) == 1:
            return region_list[0]
        if region_key in adata.obs:
            values = pd.Index(adata.obs[region_key].astype(str).unique())
            usable = [region for region in region_list if region in values]
            if len(usable) == 1:
                return usable[0]
        return None

    @staticmethod
    def _get_spatial_element(sdata: Any, key: str) -> Any:
        if hasattr(sdata, "get"):
            try:
                element = sdata.get(key)
            except KeyError:
                element = None
            if element is not None:
                return element
        for attr in ("shapes", "labels", "points", "images"):
            container = getattr(sdata, attr, None)
            if container is not None and isinstance(container, Mapping) and key in container:
                return container[key]
        return None

    @staticmethod
    def _shape_centroid_coordinates(element: Any) -> pd.DataFrame:
        centroids = element.geometry.centroid
        return pd.DataFrame(
            {
                "x": np.asarray(centroids.x, dtype=np.float64),
                "y": np.asarray(centroids.y, dtype=np.float64),
            },
            index=pd.Index(element.index.astype(str)),
        )
