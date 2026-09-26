#!/usr/bin/env python3
"""Run patched Proseg on the matched P1CRC Visium HD 2 µm outputs."""

from __future__ import annotations

import argparse
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--binned-outputs", type=Path, required=True)
    parser.add_argument("--mask", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--nthreads", type=int, default=48)
    args = parser.parse_args()

    scales_path = args.binned_outputs / "square_002um/spatial/scalefactors_json.json"
    matrix_path = args.binned_outputs / "square_002um/raw_feature_bc_matrix/matrix.mtx.gz"
    mask_meta_path = args.mask.with_suffix(".json")
    for path in (args.binary, scales_path, matrix_path, args.mask, mask_meta_path):
        if not path.exists():
            raise FileNotFoundError(path)

    scales = json.loads(scales_path.read_text())
    mask_meta = json.loads(mask_meta_path.read_text())
    if abs(scales["microns_per_pixel"] - mask_meta["microns_per_fullres_pixel"]) > 1e-7:
        raise ValueError("Mask and Space Ranger microns-per-pixel values differ")
    if mask_meta["num_nuclei"] < 10000:
        raise ValueError("Too few nuclei for full P1CRC; check the mask")

    args.output.mkdir(parents=True, exist_ok=True)
    x_transform = [str(value) for value in mask_meta["cellpose_x_transform"]]
    y_transform = [str(value) for value in mask_meta["cellpose_y_transform"]]
    command = [
        str(args.binary),
        "--visiumhd",
        "--cellpose-masks", str(args.mask),
        "--cellpose-x-transform", *x_transform,
        "--cellpose-y-transform", *y_transform,
        "--nthreads", str(args.nthreads),
        "--output-path", str(args.output),
        "--output-spatialdata", "proseg-output.zarr",
        "--exclude-spatialdata-transcripts",
        "--output-counts", "counts.mtx.gz",
        "--output-expected-counts", "expected_counts.mtx.gz",
        "--output-cell-metadata", "cell_metadata.parquet",
        "--output-gene-metadata", "gene_metadata.parquet",
        str(args.binned_outputs),
    ]
    manifest = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "binary_version": subprocess.check_output([str(args.binary), "--version"], text=True).strip(),
        "command": command,
        "mask_metadata": mask_meta,
        "space_ranger_scalefactors": scales,
    }
    (args.output / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    with (args.output / "proseg.log").open("w") as log:
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
    manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
    manifest["returncode"] = completed.returncode
    (args.output / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if completed.returncode:
        raise SystemExit(f"Proseg failed with status {completed.returncode}; see {args.output / 'proseg.log'}")


if __name__ == "__main__":
    main()
