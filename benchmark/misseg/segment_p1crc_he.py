#!/usr/bin/env python3
"""Create a P1CRC H&E nuclei mask for Proseg from the matched Visium HD image.

The output mask is cropped and downsampled. The metadata JSON records the
pixel-to-micron affine transform required by Proseg's --cellpose-*-transform
arguments. Coordinates in Space Ranger's 2 µm tissue_positions.parquet and the
original BTF image were checked to share orientation before using this script.
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import cv2
import numpy as np
import tifffile
import zarr
from csbdeep.utils import normalize
from stardist.models import StarDist2D


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="Output .npy.gz mask")
    parser.add_argument("--microns-per-pixel", type=float, required=True)
    parser.add_argument("--x0", type=int, default=0)
    parser.add_argument("--y0", type=int, default=11000)
    parser.add_argument("--x1", type=int, default=26000)
    parser.add_argument("--y1", type=int, default=38000)
    parser.add_argument("--downsample", type=int, default=2)
    parser.add_argument("--block-size", type=int, default=2048)
    parser.add_argument("--prob-thresh", type=float, default=0.01)
    parser.add_argument("--nms-thresh", type=float, default=0.001)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tifffile.TiffFile(args.image) as tif:
        source = zarr.open(tif.series[0].aszarr(), mode="r")
        height, width, channels = source.shape
        if channels != 3:
            raise ValueError(f"Expected RGB BTF, got {source.shape}")
        if not (0 <= args.x0 < args.x1 <= width and 0 <= args.y0 < args.y1 <= height):
            raise ValueError(f"Crop exceeds image dimensions {source.shape}")
        image = np.asarray(
            source[args.y0 : args.y1 : args.downsample,
                   args.x0 : args.x1 : args.downsample, :]
        )

    print(f"Segmenting RGB crop {image.shape}, original origin=({args.x0}, {args.y0})", flush=True)
    model = StarDist2D.from_pretrained("2D_versatile_he")
    image_normalized = normalize(image, 5, 95)
    labels, _ = model.predict_instances_big(
        image_normalized,
        axes="YXC",
        block_size=args.block_size,
        prob_thresh=args.prob_thresh,
        nms_thresh=args.nms_thresh,
        min_overlap=128,
        context=128,
        normalizer=None,
        n_tiles=(2, 2, 1),
    )
    del image_normalized
    labels = labels.astype(np.uint32, copy=False)
    with gzip.open(args.output, "wb", compresslevel=3) as handle:
        np.save(handle, labels)

    scale = args.microns_per_pixel * args.downsample
    metadata = {
        "image": str(args.image),
        "mask": str(args.output),
        "image_shape_yxc": [height, width, channels],
        "mask_shape_yx": list(labels.shape),
        "crop_xyxy_fullres": [args.x0, args.y0, args.x1, args.y1],
        "downsample": args.downsample,
        "microns_per_fullres_pixel": args.microns_per_pixel,
        "cellpose_x_transform": [scale, 0.0, args.x0 * args.microns_per_pixel],
        "cellpose_y_transform": [0.0, scale, args.y0 * args.microns_per_pixel],
        "segmentation_model": "StarDist2D/2D_versatile_he",
        "prob_thresh": args.prob_thresh,
        "nms_thresh": args.nms_thresh,
        "num_nuclei": int(labels.max()),
        "foreground_fraction": float(np.count_nonzero(labels) / labels.size),
    }
    metadata_path = args.output.with_suffix(".json")
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")

    mid_y, mid_x = labels.shape[0] // 2, labels.shape[1] // 2
    y0, y1 = max(0, mid_y - 512), min(labels.shape[0], mid_y + 512)
    x0, x1 = max(0, mid_x - 512), min(labels.shape[1], mid_x + 512)
    preview = image[y0:y1, x0:x1].copy()
    roi = labels[y0:y1, x0:x1]
    boundary = np.zeros(roi.shape, dtype=bool)
    boundary[1:, :] |= (roi[1:, :] != roi[:-1, :]) & (roi[1:, :] != 0)
    boundary[:, 1:] |= (roi[:, 1:] != roi[:, :-1]) & (roi[:, 1:] != 0)
    preview[boundary] = [255, 0, 0]
    cv2.imwrite(str(args.output.with_suffix(".qc.jpg")), cv2.cvtColor(preview, cv2.COLOR_RGB2BGR))
    print(json.dumps(metadata, indent=2), flush=True)


if __name__ == "__main__":
    main()
