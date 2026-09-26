#!/usr/bin/env python3
"""Create a small Visium HD binned_outputs directory for a Proseg smoke test."""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import io, sparse


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binned-outputs", type=Path, required=True)
    parser.add_argument("--mask", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    mask_meta = json.loads(args.mask.with_suffix(".json").read_text())
    x0, y0, x1, y1 = mask_meta["crop_xyxy_fullres"]

    original = args.binned_outputs / "square_002um"
    out = args.output / "square_002um"
    matrix_out = out / "raw_feature_bc_matrix"
    spatial_out = out / "spatial"
    matrix_out.mkdir(parents=True, exist_ok=True)
    spatial_out.mkdir(parents=True, exist_ok=True)

    positions = pd.read_parquet(original / "spatial/tissue_positions.parquet")
    positions = positions.loc[
        positions.in_tissue.eq(1)
        & positions.pxl_col_in_fullres.between(x0, x1)
        & positions.pxl_row_in_fullres.between(y0, y1)
    ].set_index("barcode")
    if not positions.index.is_unique:
        raise ValueError("Duplicate 2 µm barcodes")

    with h5py.File(original / "filtered_feature_bc_matrix.h5", "r") as handle:
        matrix = handle["matrix"]
        barcodes = np.char.decode(matrix["barcodes"][:], "utf-8")
        col = pd.Index(barcodes).get_indexer(positions.index)
        if (col < 0).any():
            raise ValueError("Not all selected positions are in the filtered matrix")
        col.sort()
        genes, n_bins = matrix["shape"][:]
        csc = sparse.csc_matrix(
            (matrix["data"][:], matrix["indices"][:], matrix["indptr"][:]),
            shape=(genes, n_bins),
        )
        subset = csc[:, col]
        selected_barcodes = barcodes[col]

    with gzip.open(matrix_out / "matrix.mtx.gz", "wb", compresslevel=3) as handle:
        io.mmwrite(handle, subset, field="integer")
    with gzip.open(matrix_out / "barcodes.tsv.gz", "wt") as handle:
        handle.write("\n".join(selected_barcodes) + "\n")
    shutil.copy2(original / "filtered_feature_bc_matrix/features.tsv.gz", matrix_out / "features.tsv.gz")
    shutil.copy2(original / "spatial/scalefactors_json.json", spatial_out / "scalefactors_json.json")
    positions.loc[selected_barcodes].reset_index().to_parquet(spatial_out / "tissue_positions.parquet")
    print(f"Wrote {len(col):,} bins, {genes:,} genes, {subset.nnz:,} nonzeros to {args.output}")


if __name__ == "__main__":
    main()
