#!/usr/bin/env python3
"""Export matched P1CRC spatial counts and stratified scRNA reference for RCTD/SPLIT."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse


def write_mtx(path: Path, matrix: sparse.spmatrix) -> None:
    with gzip.open(path, "wb", compresslevel=3) as handle:
        io.mmwrite(handle, matrix.tocoo(), field="integer")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spatial", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-spatial-cells", type=int, default=10000)
    parser.add_argument("--spatial-sampling", choices=("central", "random"), default="central")
    parser.add_argument("--max-reference-per-type", type=int, default=500)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--spatial-scale", type=float, default=1.0,
        help="Multiply obsm['spatial'] by this factor when X_spatial is absent",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    spatial = ad.read_h5ad(args.spatial)
    if "X_spatial" not in spatial.obsm:
        spatial.obsm["X_spatial"] = np.asarray(spatial.obsm["spatial"]) * args.spatial_scale
    if args.max_spatial_cells and spatial.n_obs > args.max_spatial_cells:
        if args.spatial_sampling == "central":
            xy = np.asarray(spatial.obsm["X_spatial"])
            # A contiguous central field preserves spatial neighborhoods for SPLIT.
            center = np.median(xy, axis=0)
            distance = np.sum(((xy - center) / np.std(xy, axis=0)) ** 2, axis=1)
            selected = np.sort(np.argpartition(distance, args.max_spatial_cells)[: args.max_spatial_cells])
        else:
            selected = np.sort(rng.choice(spatial.n_obs, args.max_spatial_cells, replace=False))
        spatial = spatial[selected].copy()

    reference = ad.read_h5ad(args.reference, backed="r")
    candidates = np.flatnonzero(reference.obs["Patient"].astype(str).to_numpy() == "P1CRC")
    labels = reference.obs.iloc[candidates]["Level1"].astype(str).to_numpy()
    selected_ref = []
    for label in sorted(set(labels)):
        idx = candidates[labels == label]
        chosen = rng.choice(idx, size=min(len(idx), args.max_reference_per_type), replace=False)
        selected_ref.extend(chosen.tolist())
    selected_ref = np.array(sorted(selected_ref), dtype=np.int64)
    ref = reference[selected_ref].to_memory()
    if not sparse.issparse(spatial.X) or not sparse.issparse(ref.X):
        raise ValueError("Expected sparse raw count matrices")
    if np.any(ref.X.data % 1) or np.any(spatial.X.data % 1):
        raise ValueError("Input matrices must contain integer UMI counts")
    ref_genes = pd.Index(ref.var_names.astype(str))
    common = pd.Index(spatial.var_names.astype(str)).intersection(ref_genes)
    if len(common) < 1000:
        raise ValueError(f"Only {len(common)} common genes")
    spatial = spatial[:, common].copy()
    ref = ref[:, common].copy()

    write_mtx(args.output / "spatial_genes_by_cells.mtx.gz", spatial.X.T)
    write_mtx(args.output / "reference_genes_by_cells.mtx.gz", ref.X.T)
    pd.Series(common).to_csv(args.output / "genes.tsv", index=False, header=False)
    pd.Series(spatial.obs_names).to_csv(args.output / "spatial_cells.tsv", index=False, header=False)
    pd.Series(ref.obs_names).to_csv(args.output / "reference_cells.tsv", index=False, header=False)
    pd.DataFrame({
        "cell_id": spatial.obs_names,
        "x_um": spatial.obsm["X_spatial"][:, 0],
        "y_um": spatial.obsm["X_spatial"][:, 1],
    }).to_csv(args.output / "spatial_coordinates.csv", index=False)
    safe_labels = ref.obs["Level1"].astype(str).str.replace(r"[^A-Za-z0-9_]", "_", regex=True)
    if safe_labels.nunique() != ref.obs["Level1"].nunique():
        raise ValueError("Sanitized cell type labels are not unique")
    pd.DataFrame({"cell_id": ref.obs_names, "cell_type": safe_labels.to_numpy()}).to_csv(
        args.output / "reference_labels.csv", index=False
    )
    manifest = {
        "spatial_input": str(args.spatial),
        "reference_input": str(args.reference),
        "spatial_cells": spatial.n_obs,
        "reference_cells": ref.n_obs,
        "genes": len(common),
        "cell_type_counts": ref.obs["Level1"].astype(str).value_counts().to_dict(),
        "max_spatial_cells": args.max_spatial_cells,
        "spatial_sampling": args.spatial_sampling,
        "max_reference_per_type": args.max_reference_per_type,
        "seed": args.seed,
        "spatial_scale": args.spatial_scale,
    }
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == "__main__":
    main()
