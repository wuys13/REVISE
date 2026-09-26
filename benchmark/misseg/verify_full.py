#!/usr/bin/env python3
"""Check that all three P1CRC runs produced usable, spatially aligned outputs."""

from __future__ import annotations

import argparse
import gzip
import json
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def matrix_market_shape(path: Path) -> tuple[int, int, int]:
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.startswith("%"):
                return tuple(map(int, line.split()))
    raise ValueError(f"Missing Matrix Market shape: {path}")


def h5ad_summary(path: Path, spatial_key: str) -> dict:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(path)
    obj = ad.read_h5ad(path, backed="r")
    try:
        if spatial_key not in obj.obsm:
            raise ValueError(f"{path}: missing obsm[{spatial_key!r}]")
        xy = np.asarray(obj.obsm[spatial_key])
        if xy.shape != (obj.n_obs, 2) or not np.isfinite(xy).all():
            raise ValueError(f"{path}: invalid spatial coordinates {xy.shape}")
        return {
            "path": str(path),
            "cells": int(obj.n_obs),
            "genes": int(obj.n_vars),
            "x_um_range": [float(xy[:, 0].min()), float(xy[:, 0].max())],
            "y_um_range": [float(xy[:, 1].min()), float(xy[:, 1].max())],
        }
    finally:
        obj.file.close()


def write_proseg_correspondence(cells: pd.DataFrame, common_path: Path, output_path: Path) -> int:
    common_obj = ad.read_h5ad(common_path, backed="r")
    try:
        common_cells = pd.DataFrame({
            "stardist_cell": common_obj.obs_names.astype(str).to_numpy(),
            "mask_label": pd.to_numeric(common_obj.obs["mask_label"], errors="raise").to_numpy(dtype=np.int64),
        })
    finally:
        common_obj.file.close()
    proseg_cells = pd.DataFrame({
        "proseg_cell": pd.to_numeric(cells["cell"], errors="raise").to_numpy(dtype=np.int64),
        "mask_label": pd.to_numeric(cells["original_cell_id"], errors="raise").to_numpy(dtype=np.int64),
    })
    if not np.array_equal(proseg_cells["proseg_cell"].to_numpy(), np.arange(len(cells))):
        raise ValueError("Proseg cell IDs are not aligned to counts matrix rows")
    correspondence = common_cells.merge(proseg_cells, on="mask_label", validate="one_to_one")
    if len(correspondence) != len(common_cells):
        raise ValueError(
            f"Only {len(correspondence)} of {len(common_cells)} input cells map to Proseg cells"
        )
    correspondence.to_parquet(output_path, index=False)
    return len(correspondence)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root
    results = root / "results"
    for stage in ("proseg", "aggregate", "resolvi", "split_inputs", "split", "split_h5ad"):
        if not (results / f".{stage}.done").exists():
            raise RuntimeError(f"Stage {stage} has no completion marker")

    mask_meta = json.loads((results / "proseg/p1crc_stardist_masks.npy.json").read_text())
    x0, y0, x1, y1 = mask_meta["crop_xyxy_fullres"]
    microns = mask_meta["microns_per_fullres_pixel"]
    proseg = results / "proseg/full"
    manifest = json.loads((proseg / "run_manifest.json").read_text())
    if manifest.get("returncode") != 0:
        raise RuntimeError(f"Proseg did not exit successfully: {manifest.get('returncode')}")
    cells = pd.read_parquet(proseg / "cell_metadata.parquet")
    count_shape = matrix_market_shape(proseg / "counts.mtx.gz")
    if count_shape[0] != len(cells):
        raise ValueError(f"Proseg counts have {count_shape[0]} rows for {len(cells)} cells")
    xy = cells[["centroid_x", "centroid_y"]].to_numpy()
    inside = (
        (xy[:, 0] >= x0 * microns) & (xy[:, 0] <= x1 * microns)
        & (xy[:, 1] >= y0 * microns) & (xy[:, 1] <= y1 * microns)
    )
    if inside.mean() < 0.95:
        raise ValueError(f"Only {inside.mean():.1%} of Proseg cells lie in the mask crop")

    common = h5ad_summary(results / "common/star_dist_cells_5um.h5ad", "X_spatial")
    correspondence_path = results / "common/proseg_cell_correspondence.parquet"
    matched_cells = write_proseg_correspondence(
        cells, results / "common/star_dist_cells_5um.h5ad", correspondence_path
    )
    latent_path = results / "resolvi/full/resolvi_latent.h5ad"
    resolvi = h5ad_summary(latent_path, "X_spatial")
    latent = ad.read_h5ad(latent_path, backed="r")
    try:
        if "X_resolVI" not in latent.obsm:
            raise ValueError("ResolVI latent embedding missing")
        latent_dims = int(latent.obsm["X_resolVI"].shape[1])
    finally:
        latent.file.close()
    expression = h5ad_summary(results / "resolvi/full/resolvi_expression.h5ad", "X_spatial")
    if (resolvi["cells"], resolvi["genes"]) != (expression["cells"], expression["genes"]):
        raise ValueError("ResolVI latent and expression outputs have different dimensions")
    split = h5ad_summary(results / "split/full/split_purified.h5ad", "X_spatial")
    prepared = json.loads((results / "split/full_input/manifest.json").read_text())
    if split["cells"] > prepared["spatial_cells"]:
        raise ValueError("SPLIT output has more cells than its input")

    report = {
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "mask_nuclei": mask_meta["num_nuclei"],
        "proseg": {
            "cells": len(cells), "genes": count_shape[1], "nonzeros": count_shape[2],
            "fraction_centroids_inside_mask_crop": float(inside.mean()),
        },
        "common_input": common,
        "common_to_proseg": {
            "matched_cells": matched_cells,
            "correspondence_path": str(correspondence_path),
        },
        "resolvi": {**resolvi, "latent_dimensions": latent_dims, "expression_path": expression["path"]},
        "split": {**split, "input_cells": prepared["spatial_cells"]},
    }
    path = results / "full_verification.json"
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
