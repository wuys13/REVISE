#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict

from revise.preprocess.histology_priors import build_histology_prior_h5ad


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build optional histology-derived REVISE priors from a raw image, "
            "a labeled segmentation mask, and spot coordinates."
        )
    )
    parser.add_argument(
        "--st-h5ad",
        required=True,
        help="Input spot-level ST AnnData/H5AD file.",
    )
    parser.add_argument(
        "--out-h5ad",
        required=True,
        help="Output H5AD with uns['all_cells_in_spot'] populated.",
    )
    parser.add_argument(
        "--mask",
        required=True,
        help="2D labeled segmentation mask. Label 0 is background; positive labels are cells.",
    )
    parser.add_argument(
        "--image",
        default=None,
        help=(
            "Optional raw histology image aligned to --mask. When provided, "
            "the preprocessor records image shape/range and per-cell mean intensity."
        ),
    )
    parser.add_argument(
        "--spots",
        default=None,
        help=(
            "Optional CSV with spot_id plus image x/y coordinates. If omitted, "
            "coordinates are read from st_h5ad.obsm['spatial'] or obs[['x', 'y']]."
        ),
    )
    parser.add_argument(
        "--spot-radius",
        type=float,
        default=None,
        help="Optional maximum image-coordinate distance for assigning a segmented cell to a spot.",
    )
    parser.add_argument(
        "--min-area",
        type=float,
        default=1.0,
        help="Minimum labeled mask area retained as a cell.",
    )
    parser.add_argument(
        "--cell-id-prefix",
        default="histology_cell",
        help="Prefix for cell ids generated from mask labels.",
    )
    parser.add_argument(
        "--report-json",
        default=None,
        help="Optional path for a compact JSON preprocessing report.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = build_histology_prior_h5ad(
        st_h5ad_path=args.st_h5ad,
        mask_path=args.mask,
        out_h5ad_path=args.out_h5ad,
        image_path=args.image,
        spots_path=args.spots,
        spot_radius=args.spot_radius,
        min_area=args.min_area,
        cell_id_prefix=args.cell_id_prefix,
    )
    report = _build_report(result.provenance, result.mapping)
    if args.report_json is not None:
        report_path = Path(args.report_json)
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


def _build_report(
    provenance: Dict[str, Any],
    mapping: Dict[str, Any],
    *,
    max_examples: int = 5,
) -> Dict[str, Any]:
    example_mapping: Dict[str, Any] = {}
    for spot_id in list(mapping.keys())[:max_examples]:
        example_mapping[spot_id] = list(mapping[spot_id])[:max_examples]
    report = dict(provenance)
    report["example_mapping"] = example_mapping
    return report


if __name__ == "__main__":
    main()
