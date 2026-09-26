#!/usr/bin/env python3
"""Convert SPLIT's Matrix Market output to an H5AD with matched coordinates."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--split-output", type=Path, required=True)
    parser.add_argument("--split-input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    genes = pd.read_csv(args.split_output / "genes.tsv", header=None)[0].astype(str).to_numpy()
    cells = pd.read_csv(args.split_output / "cells.tsv", header=None)[0].astype(str).to_numpy()
    matrix = io.mmread(args.split_output / "purified_counts.mtx")
    if not sparse.issparse(matrix):
        matrix = sparse.csr_matrix(matrix)
    matrix = matrix.T.tocsr().astype(np.float32)
    if matrix.shape != (len(cells), len(genes)):
        raise ValueError(f"Matrix shape {matrix.shape} does not match {len(cells)} cells and {len(genes)} genes")
    metadata = pd.read_csv(args.split_output / "cell_metadata.csv", index_col=0)
    metadata.index = metadata.index.astype(str)
    if not pd.Index(cells).isin(metadata.index).all():
        raise ValueError("SPLIT cell IDs missing from its metadata")
    obs = metadata.loc[cells].copy()
    coords = pd.read_csv(args.split_input / "spatial_coordinates.csv", index_col="cell_id")
    coords.index = coords.index.astype(str)
    if not pd.Index(cells).isin(coords.index).all():
        raise ValueError("SPLIT cell IDs missing from spatial input coordinates")
    xy = coords.loc[cells, ["x_um", "y_um"]].to_numpy(dtype=np.float32)
    obs[["x_um", "y_um"]] = xy

    result = ad.AnnData(X=matrix, obs=obs, var=pd.DataFrame(index=genes))
    result.obsm["spatial"] = xy
    result.obsm["X_spatial"] = xy
    result.uns["split_conversion"] = {
        "split_output": str(args.split_output),
        "split_input": str(args.split_input),
        "n_cells": int(result.n_obs),
        "n_genes": int(result.n_vars),
        "n_nonzero": int(matrix.nnz),
        "sum_purified_counts": float(matrix.sum()),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.write_h5ad(args.output, compression="lzf")
    (args.output.with_suffix(".json")).write_text(
        json.dumps(result.uns["split_conversion"], indent=2) + "\n"
    )
    print(json.dumps(result.uns["split_conversion"], indent=2), flush=True)


if __name__ == "__main__":
    main()
