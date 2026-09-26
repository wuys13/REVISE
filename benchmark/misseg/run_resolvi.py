#!/usr/bin/env python3
"""Train ResolVI on P1CRC's common, pre-Proseg StarDist cell matrix."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc
import torch

# The inherited qz base environment contains an incompatible TensorFlow 2.8.
# ResolVI is PyTorch-only; prevent an optional torchmetrics import from loading it.
sys.modules["tensorflow"] = None
import scvi


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--max-cells", type=int, default=None, help="Pilot only; omit for full run")
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--skip-expression", action="store_true", help="Pilot only")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    scvi.settings.seed = args.seed
    adata = ad.read_h5ad(args.input)
    if args.max_cells and args.max_cells < adata.n_obs:
        rng = np.random.default_rng(args.seed)
        selected = np.sort(rng.choice(adata.n_obs, args.max_cells, replace=False))
        adata = adata[selected].copy()
    if "X_spatial" not in adata.obsm:
        adata.obsm["X_spatial"] = adata.obsm["spatial"].copy()
    sc.pp.filter_genes(adata, min_cells=3)
    if adata.n_obs < 20 or adata.n_vars < 100:
        raise ValueError(f"Too few observations/features after filtering: {adata.shape}")

    scvi.external.RESOLVI.setup_anndata(
        adata,
        prepare_data_kwargs={"n_neighbors": 10, "spatial_rep": "X_spatial"},
    )
    model = scvi.external.RESOLVI(adata, semisupervised=False)
    model.train(
        max_epochs=args.epochs,
        batch_size=args.batch_size,
    )
    model.save(args.output / "model", overwrite=True, save_anndata=False)
    adata.obsm["X_resolVI"] = model.get_latent_representation(adata).astype(np.float32)
    adata.uns["resolvi_run"] = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "scvi_version": scvi.__version__,
        "torch_version": str(torch.__version__),
        "input": str(args.input),
        "epochs": args.epochs,
        "batch_size": args.batch_size,
        "seed": args.seed,
        "semisupervised": False,
        "shape": list(adata.shape),
        "expression_output": "normalized decoded gene expression" if not args.skip_expression else None,
    }
    (args.output / "run_manifest.json").write_text(json.dumps(adata.uns["resolvi_run"], indent=2) + "\n")

    latent = adata.copy()
    latent.write_h5ad(args.output / "resolvi_latent.h5ad", compression="lzf")
    if not args.skip_expression:
        normalized = model.get_normalized_expression(
            adata, batch_size=args.batch_size, return_numpy=True
        )
        adata.X = np.asarray(normalized, dtype=np.float32)
        adata.write_h5ad(args.output / "resolvi_expression.h5ad", compression="lzf")
        print(f"Wrote decoded gene expression for {adata.shape}", flush=True)


if __name__ == "__main__":
    main()
