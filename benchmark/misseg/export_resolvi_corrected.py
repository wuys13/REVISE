#!/usr/bin/env python3
"""Export ResolVI posterior corrected count rates without holding the full matrix in RAM."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import torch

# The inherited qz environment contains an incompatible optional TensorFlow.
sys.modules["tensorflow"] = None
import scvi


def write_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="ResolVI latent H5AD with raw X")
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--block-size", type=int, default=2048)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--num-samples", type=int, default=3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--max-cells", type=int, default=None, help="Pilot only")
    args = parser.parse_args()
    if min(args.block_size, args.batch_size, args.num_samples) < 1:
        raise ValueError("block-size, batch-size, and num-samples must be positive")
    args.output.parent.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input)
    n_cells = min(adata.n_obs, args.max_cells) if args.max_cells else adata.n_obs
    n_genes = adata.n_vars
    if n_cells < 1 or n_genes < 1:
        raise ValueError(f"Empty input: {adata.shape}")
    scvi.settings.seed = args.seed
    model = scvi.external.RESOLVI.load(args.model, adata=adata)

    scratch = args.output.with_suffix(".float32.mmap")
    progress_path = args.output.with_suffix(".progress.json")
    if args.output.exists() and not progress_path.exists():
        raise FileExistsError(f"Completed output already exists: {args.output}")
    config = {
        "input": str(args.input.resolve()),
        "model": str(args.model.resolve()),
        "output": str(args.output.resolve()),
        "cells": n_cells,
        "genes": n_genes,
        "block_size": args.block_size,
        "batch_size": args.batch_size,
        "num_samples": args.num_samples,
        "seed": args.seed,
        "method": "median posterior px_rate from model_corrected",
        "scvi_version": scvi.__version__,
    }
    if progress_path.exists():
        progress = json.loads(progress_path.read_text())
        if any(progress.get(key) != value for key, value in config.items()):
            raise ValueError(f"Existing progress has different settings: {progress_path}")
        if not scratch.exists() or scratch.stat().st_size != n_cells * n_genes * 4:
            raise ValueError("Existing progress has no matching float32 scratch matrix")
        next_start = int(progress["next_start"])
        matrix = np.memmap(scratch, mode="r+", dtype=np.float32, shape=(n_cells, n_genes))
    else:
        if scratch.exists():
            raise FileExistsError(f"Scratch exists without a progress file: {scratch}")
        next_start = 0
        matrix = np.memmap(scratch, mode="w+", dtype=np.float32, shape=(n_cells, n_genes))
        write_json(progress_path, {**config, "next_start": next_start})

    for start in range(next_start, n_cells, args.block_size):
        end = min(start + args.block_size, n_cells)
        np.random.seed(args.seed + start)
        torch.manual_seed(args.seed + start)
        samples = model.sample_posterior(
            adata=adata,
            indices=np.arange(start, end),
            model=model.module.model_corrected,
            return_sites=["px_rate"],
            summary_fun={"post_sample_q50": np.median},
            num_samples=args.num_samples,
            summary_frequency=max(1, args.block_size // args.batch_size),
            batch_size=args.batch_size,
        )
        values = np.asarray(samples["post_sample_q50"]["px_rate"], dtype=np.float32)
        if values.shape != (end - start, n_genes):
            raise ValueError(f"Block {start}:{end} has unexpected shape {values.shape}")
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(f"Block {start}:{end} has invalid corrected rates")
        matrix[start:end] = values
        matrix.flush()
        write_json(progress_path, {**config, "next_start": end})
        print(f"Corrected counts: {end:,}/{n_cells:,} cells", flush=True)

    corrected = adata[:n_cells].copy()
    # AnnData's writer accepts ndarray views but does not register numpy.memmap.
    corrected.X = np.asarray(matrix)
    corrected.uns["resolvi_corrected_run"] = {
        **config,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "output_values": "continuous expected counts, not integer UMI",
    }
    corrected.write_h5ad(args.output, compression="lzf")
    check = ad.read_h5ad(args.output, backed="r")
    try:
        if check.shape != (n_cells, n_genes):
            raise ValueError(f"Written output has unexpected shape {check.shape}")
    finally:
        check.file.close()
    del corrected, matrix
    scratch.unlink()
    progress_path.unlink()
    print(f"Wrote {args.output} with shape {(n_cells, n_genes)}", flush=True)


if __name__ == "__main__":
    main()
