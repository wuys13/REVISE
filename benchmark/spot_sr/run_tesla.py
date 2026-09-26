"""Run the upstream TESLA interpolation on one prepared REVISE spot case."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import anndata as ad
import cv2
import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--tesla-repo", type=Path, required=True)
    parser.add_argument("--grid-microns", type=float, default=10.0)
    parser.add_argument("--neighbors", type=int, default=10)
    args = parser.parse_args()
    sys.path.insert(0, str(args.tesla_repo / "TESLA_package"))
    import TESLA as tesla

    manifest = json.loads((args.case / "case_manifest.json").read_text())
    meta = json.loads((args.case.parent / "image_transform.json").read_text())
    image = cv2.imread(str(args.case / "he-raw.jpg"), cv2.IMREAD_COLOR)
    if image is None:
        raise FileNotFoundError(args.case / "he-raw.jpg")
    cnts = pd.read_csv(args.case / "cnts.tsv", sep="\t", index_col=0)
    locs = pd.read_csv(args.case / "locs-raw.tsv", sep="\t", index_col=0).loc[cnts.index]
    spots = ad.AnnData(np.log1p(cnts.to_numpy(dtype=np.float32)))
    spots.obs_names = cnts.index.astype(str)
    spots.var_names = cnts.columns.astype(str)
    # TESLA's pixel_x is image row and pixel_y is image column.
    spots.obs["pixel_x"] = locs["y"].round().astype(int).to_numpy()
    spots.obs["pixel_y"] = locs["x"].round().astype(int).to_numpy()
    h, w = image.shape[:2]
    contour = np.array([[[1, 1]], [[w - 2, 1]], [[w - 2, h - 2]], [[1, h - 2]]], dtype=np.int32)
    grid_px = max(2, int(round(args.grid_microns / meta["target_mpp"])))
    start = time.time()
    prediction = tesla.imputation(
        img=image, raw=spots, cnt=contour, genes=spots.var_names.tolist(),
        shape="None", res=grid_px, s=1, k=2,
        num_nbs=min(args.neighbors, spots.n_obs),
    )
    prediction.obs["image_row"] = prediction.obs["x"].to_numpy()
    prediction.obs["image_col"] = prediction.obs["y"].to_numpy()
    prediction.X = np.expm1(np.asarray(prediction.X, dtype=np.float32))
    prediction.X = np.maximum(prediction.X, 0)
    prediction.write_h5ad(args.case / "tesla_superpixels.h5ad", compression="gzip")
    status = {
        "method": "TESLA", "case": str(args.case), "n_superpixels": prediction.n_obs,
        "n_genes": prediction.n_vars, "grid_microns": args.grid_microns,
        "runtime_seconds": time.time() - start,
        "upstream_revision": subprocess.check_output(
            ["git", "-C", str(args.tesla_repo), "rev-parse", "HEAD"], text=True
        ).strip(),
    }
    (args.case / "tesla_run.json").write_text(json.dumps(status, indent=2) + "\n")
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
