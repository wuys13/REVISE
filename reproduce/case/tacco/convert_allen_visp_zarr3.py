from __future__ import annotations

import importlib.metadata
import platform
import sys
from pathlib import Path

import anndata as ad
import numpy as np

EXPECTED_SHAPE = (2_389, 268)
EXPECTED_OBS_COLUMNS = ["cell_id", "region"]


def main(argv: list[str] | None = None) -> None:
    arguments = list(sys.argv[1:] if argv is None else argv)
    if len(arguments) != 2:
        raise SystemExit(
            "usage: convert_allen_visp_zarr3.py SOURCE_DATA_ZARR OUTPUT_H5AD"
        )

    source_zarr = Path(arguments[0]).resolve()
    output_h5ad = Path(arguments[1]).resolve()
    if tuple(map(int, platform.python_version_tuple()[:2])) < (3, 11):
        raise RuntimeError(
            "The isolated Zarr3 converter requires Python >=3.11; "
            f"got {platform.python_version()}"
        )
    zarr_version = importlib.metadata.version("zarr")
    anndata_version = importlib.metadata.version("anndata")
    if not zarr_version.startswith("3."):
        raise RuntimeError(f"The converter requires Zarr 3, got {zarr_version}")
    if not anndata_version.startswith("0.12."):
        raise RuntimeError(
            f"The converter requires AnnData 0.12.x, got {anndata_version}"
        )
    if output_h5ad.exists():
        raise FileExistsError(f"Refusing to overwrite converter output: {output_h5ad}")

    table_path = source_zarr / "tables" / "table"
    table = ad.read_zarr(table_path)
    if tuple(table.shape) != EXPECTED_SHAPE:
        raise ValueError(f"Unexpected Allen VISp table shape: {table.shape}")
    if list(table.obs.columns) != EXPECTED_OBS_COLUMNS:
        raise ValueError(
            f"Unexpected Allen VISp obs columns: {list(table.obs.columns)!r}"
        )
    if list(table.layers) or list(table.obsm) or table.raw is not None:
        raise ValueError("Allen VISp source table unexpectedly has layers, obsm, or raw")
    if not table.obs_names.is_unique or not table.var_names.is_unique:
        raise ValueError("Allen VISp source table axes are not unique")
    values = np.asarray(table.X)
    if values.shape != EXPECTED_SHAPE:
        raise ValueError(f"Unexpected dense expression shape: {values.shape}")
    if values.dtype.kind not in {"i", "u"}:
        raise ValueError(f"Expected integer MERFISH counts, got {values.dtype}")
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Allen VISp source table contains invalid counts")

    blank_genes = table.var_names[
        table.var_names.str.startswith("Blank-")
    ].astype(str).tolist()
    if set(blank_genes) != {f"Blank-{index}" for index in range(1, 11)}:
        raise ValueError(f"Unexpected MERFISH blank controls: {blank_genes}")
    cell_ids = table.obs["cell_id"].to_numpy(dtype=int)
    if not np.array_equal(np.sort(cell_ids), np.arange(EXPECTED_SHAPE[0])):
        raise ValueError("Allen VISp table cell_id must be a permutation of 0..2388")

    table.uns["zarr3_conversion"] = {
        "source_zarr": str(source_zarr),
        "source_table": "tables/table",
        "python": platform.python_version(),
        "anndata": anndata_version,
        "zarr": zarr_version,
        "expression_source": "tables/table/X",
        "expression_dtype": str(values.dtype),
        "source_shape": list(EXPECTED_SHAPE),
        "source_obs_columns": EXPECTED_OBS_COLUMNS,
        "blank_controls": blank_genes,
    }
    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    table.write_h5ad(output_h5ad, compression="gzip")


if __name__ == "__main__":
    main()
