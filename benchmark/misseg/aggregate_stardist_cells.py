#!/usr/bin/env python3
"""Aggregate P1CRC 2 µm bins into cells from the same StarDist prior as Proseg.

This provides a common *pre-Proseg* cell count matrix for ResolVI and SPLIT.
Nucleus labels are extended to the closest pixel up to a recorded distance;
this is an explicit analysis choice rather than a 10x cell segmentation output.
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import ndimage, sparse


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binned-outputs", type=Path, required=True)
    parser.add_argument("--mask", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-expansion-um", type=float, default=5.0)
    args = parser.parse_args()

    mask_meta = json.loads(args.mask.with_suffix(".json").read_text())
    with gzip.open(args.mask, "rb") as handle:
        mask = np.load(handle)
    if mask.shape != tuple(mask_meta["mask_shape_yx"]):
        raise ValueError("Mask dimensions differ from metadata")

    positions_path = args.binned_outputs / "square_002um/spatial/tissue_positions.parquet"
    positions = pd.read_parquet(
        positions_path,
        columns=["barcode", "in_tissue", "pxl_col_in_fullres", "pxl_row_in_fullres"],
    )
    positions = positions.loc[positions["in_tissue"] == 1].set_index("barcode")
    h5_path = args.binned_outputs / "square_002um/filtered_feature_bc_matrix.h5"
    with h5py.File(h5_path, "r") as handle:
        matrix = handle["matrix"]
        barcodes = np.char.decode(matrix["barcodes"][:], "utf-8")
        if len(barcodes) != len(positions):
            raise ValueError("Filtered matrix and under-tissue positions have different sizes")
        loc = positions.index.get_indexer(barcodes)
        if (loc < 0).any():
            raise ValueError(f"Missing {(loc < 0).sum()} filtered barcodes in positions")
        x_px = positions["pxl_col_in_fullres"].to_numpy()[loc]
        y_px = positions["pxl_row_in_fullres"].to_numpy()[loc]
        del positions, barcodes, loc

        x0, y0, _, _ = mask_meta["crop_xyxy_fullres"]
        factor = mask_meta["downsample"]
        col = np.floor((x_px - x0) / factor).astype(np.int32)
        row = np.floor((y_px - y0) / factor).astype(np.int32)
        inside = (row >= 0) & (row < mask.shape[0]) & (col >= 0) & (col < mask.shape[1])
        print(f"Bins inside mask crop: {inside.sum():,} / {len(inside):,}", flush=True)

        # Pixel distance is converted to µm with the same scale used by Proseg.
        distance, nearest = ndimage.distance_transform_edt(mask == 0, return_indices=True)
        selected = np.zeros(len(row), dtype=np.uint32)
        valid_indices = np.flatnonzero(inside)
        r, c = row[valid_indices], col[valid_indices]
        max_distance_px = args.max_expansion_um / (
            factor * mask_meta["microns_per_fullres_pixel"]
        )
        near = distance[r, c] <= max_distance_px
        selected[valid_indices[near]] = mask[nearest[0, r[near], c[near]], nearest[1, r[near], c[near]]]
        del distance, nearest, mask, row, col, inside, r, c, valid_indices, near

        labels = np.unique(selected)
        labels = labels[labels != 0]
        label_to_cell = np.full(int(labels.max()) + 1, -1, dtype=np.int32)
        label_to_cell[labels] = np.arange(len(labels), dtype=np.int32)
        cell_for_bin = np.full(len(selected), -1, dtype=np.int32)
        assigned = selected > 0
        cell_for_bin[assigned] = label_to_cell[selected[assigned]]
        del label_to_cell, selected
        print(f"Assigned bins: {assigned.sum():,}; cells: {len(labels):,}", flush=True)

        indptr = matrix["indptr"][:]
        genes = np.char.decode(matrix["features/name"][:], "utf-8")
        gene_ids = np.char.decode(matrix["features/id"][:], "utf-8")
        nnz_per_bin = np.diff(indptr).astype(np.int32)
        row_index = np.repeat(cell_for_bin, nnz_per_bin)
        keep = row_index >= 0
        gene_index = matrix["indices"][:].astype(np.int32, copy=False)[keep]
        counts = matrix["data"][:][keep]
        cell_index = row_index[keep]
        del row_index, keep, indptr, nnz_per_bin

    cell_counts = sparse.coo_matrix(
        (counts, (cell_index, gene_index)),
        shape=(len(labels), len(genes)),
        dtype=np.int32,
    ).tocsr()
    cell_counts.sum_duplicates()
    del counts, gene_index, cell_index

    bin_counts = np.bincount(cell_for_bin[assigned], minlength=len(labels))
    x_mean = np.bincount(cell_for_bin[assigned], weights=x_px[assigned], minlength=len(labels)) / bin_counts
    y_mean = np.bincount(cell_for_bin[assigned], weights=y_px[assigned], minlength=len(labels)) / bin_counts
    um_per_px = mask_meta["microns_per_fullres_pixel"]
    obs = pd.DataFrame(
        {
            "mask_label": labels,
            "n_2um_bins": bin_counts,
            "x_um": x_mean * um_per_px,
            "y_um": y_mean * um_per_px,
        },
        index=[f"stardist_{label}" for label in labels],
    )
    var = pd.DataFrame({"gene_id": gene_ids}, index=pd.Index(genes, name="gene"))
    result = ad.AnnData(X=cell_counts, obs=obs, var=var)
    result.var_names_make_unique()
    result.obsm["spatial"] = result.obs[["x_um", "y_um"]].to_numpy()
    result.obsm["X_spatial"] = result.obsm["spatial"].copy()
    result.uns["misseg_input"] = {
        "source": "P1CRC 2 um filtered Space Ranger matrix",
        "segmentation": "StarDist2D 2D_versatile_he nuclei with nearest-nucleus expansion",
        "max_expansion_um": args.max_expansion_um,
        "num_input_bins": len(cell_for_bin),
        "num_assigned_bins": int(assigned.sum()),
        "mask_metadata": mask_meta,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.write_h5ad(args.output, compression="gzip")
    print(f"Wrote {result.shape} with {result.X.nnz:,} nonzero entries to {args.output}", flush=True)


if __name__ == "__main__":
    main()
