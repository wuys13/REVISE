"""Sample baseline images at held-out Xenium cells and score like REVISE.

The cell coordinates are used only after model inference. They must never be
passed to TESLA or iStar for training or prediction.
"""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import pearsonr
from skimage.metrics import structural_similarity


def truth_and_positions(case: Path) -> tuple[ad.AnnData, np.ndarray]:
    manifest = json.loads((case / "case_manifest.json").read_text())
    truth = ad.read_h5ad(manifest["source_truth"])
    locs = pd.read_csv(case / "truth-cell-locs.tsv", sep="\t", index_col=0)
    if len(locs) != truth.n_obs:
        raise ValueError("Truth cell coordinate count mismatch")
    return truth, locs[["x", "y"]].to_numpy(dtype=np.float64)


def sample_tesla(case: Path, xy: np.ndarray) -> tuple[np.ndarray, list[str], dict]:
    prediction = ad.read_h5ad(case / "tesla_superpixels.h5ad")
    grid_xy = prediction.obs[["image_col", "image_row"]].to_numpy(dtype=np.float64)
    distances, indices = cKDTree(grid_xy).query(xy)
    matrix = prediction.X[indices]
    matrix = matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)
    meta = {
        "sampling": "nearest TESLA superpixel to each held-out cell center",
        "median_distance_pixels": float(np.median(distances)),
        "max_distance_pixels": float(np.max(distances)),
    }
    return matrix, prediction.var_names.astype(str).tolist(), meta


def sample_istar(case: Path, xy: np.ndarray) -> tuple[np.ndarray, list[str], dict]:
    genes = [s for s in (case / "gene-names.txt").read_text().splitlines() if s]
    # iStar's histology embedding and gene output have one sample per 16 px.
    # Its raw image is padded to a multiple of 16 by preprocess.py.
    first = np.asarray(pickle.loads((case / "cnts-super" / f"{genes[0]}.pickle").read_bytes()))
    if first.ndim != 2:
        raise ValueError(f"Expected 2D iStar gene grid; found {first.shape}")
    # Grid (0, 0) describes the image patch [0:16, 0:16], not a point at (0, 0).
    ij = np.floor(xy[:, ::-1] / 16).astype(int)
    ij[:, 0] = np.clip(ij[:, 0], 0, first.shape[0] - 1)
    ij[:, 1] = np.clip(ij[:, 1], 0, first.shape[1] - 1)
    matrix = np.empty((len(xy), len(genes)), dtype=np.float32)
    for j, gene in enumerate(genes):
        grid = np.asarray(pickle.loads((case / "cnts-super" / f"{gene}.pickle").read_bytes()))
        if grid.shape != first.shape:
            raise ValueError(f"Grid shape mismatch for {gene}: {grid.shape}")
        matrix[:, j] = grid[ij[:, 0], ij[:, 1]]
    meta = {"sampling": "iStar 16-pixel patch containing each held-out cell center",
            "grid_shape": list(first.shape)}
    return matrix, genes, meta


def score(truth: ad.AnnData, prediction: ad.AnnData) -> pd.DataFrame:
    # Matches revise.analysis.metrics.compute_metric(normalize=True): library
    # normalization is performed before selecting the overlapping gene panel.
    truth = truth.copy()
    prediction = prediction.copy()
    sc.pp.normalize_total(truth, target_sum=1e4)
    sc.pp.normalize_total(prediction, target_sum=1e4)
    genes = truth.var_names.intersection(prediction.var_names)
    if len(genes) == 0:
        raise ValueError("No common genes")
    truth = truth[:, genes]
    prediction = prediction[:, genes]
    actual = truth.X.toarray() if sparse.issparse(truth.X) else np.asarray(truth.X)
    estimated = prediction.X.toarray() if sparse.issparse(prediction.X) else np.asarray(prediction.X)
    actual, estimated = actual.astype(float), estimated.astype(float)
    for matrix in (actual, estimated):
        lo = matrix.min(axis=0)
        span = matrix.max(axis=0) - lo
        span[span == 0] = 1
        matrix -= lo
        matrix /= span
    rows = []
    for j, gene in enumerate(genes):
        a, p = actual[:, j], estimated[:, j]
        pcc = float(pearsonr(a, p).statistic) if np.std(a) and np.std(p) else float("nan")
        mse = float(np.mean((a - p) ** 2))
        rows.append((gene, pcc, float(structural_similarity(a, p, data_range=1)),
                     mse, float(np.sqrt(mse) / np.mean(a)) if np.mean(a) else float("nan")))
    return pd.DataFrame(rows, columns=["Gene", "PCC", "SSIM", "MSE", "NRMSE"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--method", choices=["tesla", "istar"], required=True)
    args = parser.parse_args()
    truth, xy = truth_and_positions(args.case)
    if args.method == "tesla":
        data, genes, sampling = sample_tesla(args.case, xy)
    else:
        data, genes, sampling = sample_istar(args.case, xy)
    if not np.all(np.isfinite(data)):
        raise ValueError("Prediction contains nonfinite values")
    prediction = ad.AnnData(np.maximum(data, 0).astype(np.float32))
    prediction.var_names = genes
    prediction.obs_names = truth.obs_names.copy()
    prediction.obsm["spatial"] = truth.obsm["spatial"].copy()
    prediction.write_h5ad(args.case / f"{args.method}_predicted_cells.h5ad", compression="gzip")
    metrics = score(truth, prediction)
    metrics.to_csv(args.case / f"{args.method}_metrics_normalized.csv", index=False)
    summary = {
        "method": args.method, "case": str(args.case), "n_cells": truth.n_obs,
        "n_genes": len(metrics), "metrics": {
            column: float(np.nanmean(metrics[column]))
            for column in ("PCC", "SSIM", "MSE", "NRMSE")
        }, **sampling,
    }
    (args.case / f"{args.method}_metrics_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
