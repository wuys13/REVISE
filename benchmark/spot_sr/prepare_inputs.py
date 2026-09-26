"""Convert REVISE Sim2Real-ST spot cases to TESLA and iStar inputs.

Run on qz with the shared baseline environment. H&E keypoints are fitted as a
similarity transform from Xenium morphology pixels to H&E image pixels.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


XENIUM_MPP = 0.2125


def fit_alignment(path: Path) -> tuple[np.ndarray, np.ndarray, float]:
    landmarks = pd.read_csv(path)
    fixed = landmarks[["fixedX", "fixedY"]].to_numpy(dtype=float)
    aligned = landmarks[["alignmentX", "alignmentY"]].to_numpy(dtype=float)
    x0, y0 = fixed.mean(axis=0), aligned.mean(axis=0)
    u, singular, vt = np.linalg.svd((fixed - x0).T @ (aligned - y0))
    rotation = u @ vt  # reflection is permitted, as in Xenium Explorer
    scale = singular.sum() / np.square(fixed - x0).sum()
    linear = scale * rotation
    shift = y0 - x0 @ linear
    residual = np.linalg.norm(fixed @ linear + shift - aligned, axis=1)
    return linear, shift, float(np.sqrt(np.mean(residual**2)))


def to_image(xy_microns: np.ndarray, meta: dict) -> np.ndarray:
    xy = np.asarray(xy_microns, dtype=float) / XENIUM_MPP
    full = xy @ np.asarray(meta["linear"]) + np.asarray(meta["shift"])
    return (full - np.asarray(meta["crop_origin"])) * meta["resize_scale"]


def prepare_part(args: argparse.Namespace) -> None:
    import pyvips

    part_dir = args.root / "spot" / args.part
    cells = ad.read_h5ad(part_dir / "selected_xenium.h5ad", backed="r")
    xy = np.asarray(cells.obsm["spatial"], dtype=float)
    linear, shift, rmse = fit_alignment(args.alignment)
    if rmse > 150:
        raise ValueError(f"H&E alignment keypoint RMSE is {rmse:.1f} pixels")
    he_xy = xy / XENIUM_MPP @ linear + shift
    image = pyvips.Image.new_from_file(str(args.image), access="random")
    if image.bands > 3:
        image = image[:3]
    he_mpp = XENIUM_MPP / np.linalg.norm(linear[0])
    margin = args.margin_microns / he_mpp
    origin = np.maximum(np.floor(he_xy.min(axis=0) - margin).astype(int), 0)
    end = np.minimum(
        np.ceil(he_xy.max(axis=0) + margin).astype(int),
        np.array([image.width, image.height]),
    )
    wh = end - origin
    if np.any(wh <= 0):
        raise ValueError("Selected region does not overlap the H&E image")
    resize_scale = he_mpp / args.target_mpp
    out = args.output / args.part
    out.mkdir(parents=True, exist_ok=True)
    cropped = image.crop(int(origin[0]), int(origin[1]), int(wh[0]), int(wh[1]))
    cropped = cropped.resize(resize_scale, kernel="lanczos3")
    cropped.jpegsave(str(out / "he-raw.jpg"), Q=92)
    meta = {
        "part": args.part,
        "image": str(args.image),
        "alignment": str(args.alignment),
        "linear": linear.tolist(),
        "shift": shift.tolist(),
        "alignment_rmse_pixels": rmse,
        "crop_origin": origin.tolist(),
        "crop_size_native_pixels": wh.tolist(),
        "resize_scale": resize_scale,
        "target_mpp": args.target_mpp,
        "image_size": [cropped.width, cropped.height],
        "n_truth_cells": cells.n_obs,
    }
    (out / "image_transform.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps({k: meta[k] for k in ("part", "alignment_rmse_pixels", "image_size", "n_truth_cells")}, indent=2))


def prepare_case(args: argparse.Namespace) -> None:
    part_dir = args.root / "spot" / args.part
    part_out = args.output / args.part
    meta = json.loads((part_out / "image_transform.json").read_text())
    source = part_dir / f"spot_{args.size}" / "xenium_spot.h5ad"
    spots = ad.read_h5ad(source)
    cells = ad.read_h5ad(part_dir / "selected_xenium.h5ad", backed="r")
    case = part_out / f"spot_{args.size}"
    case.mkdir(parents=True, exist_ok=True)
    raw_image = case / "he-raw.jpg"
    if not raw_image.exists():
        raw_image.symlink_to(Path("..") / "he-raw.jpg")

    counts = spots.X.toarray() if sparse.issparse(spots.X) else np.asarray(spots.X)
    if np.any(counts < 0) or not np.all(np.isfinite(counts)):
        raise ValueError("Expected nonnegative, finite raw spot counts")
    names = spots.var_names.astype(str)
    pd.DataFrame(counts, index=spots.obs_names, columns=names).to_csv(case / "cnts.tsv", sep="\t")
    spot_xy = to_image(spots.obsm["spatial"], meta)
    pd.DataFrame(spot_xy, index=spots.obs_names, columns=["x", "y"]).to_csv(
        case / "locs-raw.tsv", sep="\t"
    )
    cell_xy = to_image(cells.obsm["spatial"], meta)
    cell_ids = cells.obs["cell_id"].astype(str) if "cell_id" in cells.obs else cells.obs_names.astype(str)
    pd.DataFrame(cell_xy, index=cell_ids, columns=["x", "y"]).to_csv(case / "truth-cell-locs.tsv", sep="\t")
    (case / "gene-names.txt").write_text("\n".join(names) + "\n")
    (case / "pixel-size-raw.txt").write_text(f'{meta["target_mpp"]}\n')
    (case / "pixel-size.txt").write_text(f'{meta["target_mpp"]}\n')
    # iStar models circular capture disks; the released pseudo-spots are square.
    (case / "radius-raw.txt").write_text(f'{args.size / (2 * meta["target_mpp"]):.6f}\n')
    width, height = meta["image_size"]
    if not (np.all((spot_xy[:, 0] >= 0) & (spot_xy[:, 0] < width)) and
            np.all((spot_xy[:, 1] >= 0) & (spot_xy[:, 1] < height))):
        raise ValueError("Some pseudo-spots fall outside the cropped H&E image")
    manifest = {
        "part": args.part,
        "spot_size_microns": args.size,
        "n_spots": spots.n_obs,
        "n_genes": spots.n_vars,
        "n_truth_cells": cells.n_obs,
        "source_spot": str(source),
        "source_truth": str(part_dir / "selected_xenium.h5ad"),
        "image_transform": str(part_out / "image_transform.json"),
    }
    (case / "case_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["part", "case"])
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--part", choices=["part1", "part2", "part3"], required=True)
    parser.add_argument("--size", type=int, choices=[50, 100, 150, 200])
    parser.add_argument("--image", type=Path)
    parser.add_argument("--alignment", type=Path)
    parser.add_argument("--target-mpp", type=float, default=0.5)
    parser.add_argument("--margin-microns", type=float, default=250.0)
    args = parser.parse_args()
    if args.mode == "part":
        if args.image is None or args.alignment is None:
            parser.error("part mode requires --image and --alignment")
        prepare_part(args)
    else:
        if args.size is None:
            parser.error("case mode requires --size")
        prepare_case(args)


if __name__ == "__main__":
    main()
