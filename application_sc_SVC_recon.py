#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path
from typing import List

from revise.framework import REVISEPipeline


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser("Run sc-SVC reconstruction pipeline (unified interface wrapper)")
    parser.add_argument("--config", default="revise/revise.yaml", help="Path to unified config")
    parser.add_argument("--platform", default="iST", choices=["iST", "sST"], help="Application platform route")
    parser.add_argument("--confounding", default=None, help="Route confounding factor; defaults by platform")
    parser.add_argument("--profile", default=None, help="Unified config profile; inferred from route when omitted")
    parser.add_argument("--sample-name", required=True, help="Patient/sample ID")
    parser.add_argument("--st-file", required=True, help="Spatial h5ad filename")
    parser.add_argument("--data-root", required=True, help="Directory containing input h5ad files")
    parser.add_argument(
        "--output-root",
        default="output/sc_SVC_case",
        help="Unified output root directory",
    )
    parser.add_argument("--sc-ref-file", required=True, help="sc reference h5ad filename")
    parser.add_argument("--patient-key", default="Patient", help="Patient key in sc reference obs")
    parser.add_argument("--select-ct", required=True, help="Cell type for local refinement")
    parser.add_argument("--cell-type-col", default="Level1", help="Cell type column name")
    parser.add_argument("--sub-cell-type-col", default="Level2", help="Sub-cell type column name")
    parser.add_argument(
        "--compatibility-mode",
        action="store_true",
        help="Write notebook-compatible output filenames (sc_SVC_expr.h5ad/sc_SVC_spatial.h5ad)",
    )
    return parser.parse_args()


def _build_set_overrides(args: argparse.Namespace) -> List[str]:
    set_overrides: List[str] = [
        f"sc.select_ct={args.select_ct}",
        f"columns.cell_type_col={args.cell_type_col}",
        f"columns.sub_cell_type_col={args.sub_cell_type_col}",
    ]
    return set_overrides


def _first_existing(paths: List[Path], label: str) -> Path:
    for path in paths:
        if path.exists():
            return path
    searched = ", ".join(str(path) for path in paths)
    raise FileNotFoundError(f"{label} output was not found; searched: {searched}")


def _publish_notebook_outputs(args: argparse.Namespace, run_dir: str | None) -> dict[str, str]:
    if not run_dir:
        raise RuntimeError("pipeline finished without a run_dir")

    run_path = Path(run_dir)
    expr_src = _first_existing(
        [
            run_path / "sc_SVC_expr.h5ad",
            run_path / "artifacts" / "sc_svc_expr.h5ad",
        ],
        "sc-SVC expression",
    )
    spatial_src = _first_existing(
        [
            run_path / "sc_SVC_spatial.h5ad",
            run_path / "artifacts" / "sc_svc_spatial.h5ad",
        ],
        "sc-SVC spatial",
    )

    data_type = Path(args.st_file).stem
    dst_dir = Path(args.output_root) / f"{args.sample_name}_{data_type}" / args.select_ct
    dst_dir.mkdir(parents=True, exist_ok=True)
    expr_dst = dst_dir / "sc_SVC_expr.h5ad"
    spatial_dst = dst_dir / "sc_SVC_spatial.h5ad"
    shutil.copy2(expr_src, expr_dst)
    shutil.copy2(spatial_src, spatial_dst)
    return {
        "sc_SVC_expr": str(expr_dst),
        "sc_SVC_spatial": str(spatial_dst),
    }


def main(args: argparse.Namespace) -> None:
    runtime_overrides = {
        "platform": args.platform,
        "compatibility_mode": bool(args.compatibility_mode),
    }
    if args.confounding:
        runtime_overrides["confounding"] = args.confounding

    pipeline = REVISEPipeline(config_path=args.config)
    svc = pipeline.run(
        profile=args.profile,
        runtime_overrides=runtime_overrides,
        io_overrides={
            "data_root": args.data_root,
            "output_root": args.output_root,
            "sample_name": args.sample_name,
            "st_file": args.st_file,
            "sc_ref_file": args.sc_ref_file,
            "patient_key": args.patient_key,
        },
        set_overrides=_build_set_overrides(args),
        dry_run=False,
    )
    summary = svc.summary()
    summary["run_dir"] = svc.provenance.get("run_dir")
    summary["notebook_outputs"] = _publish_notebook_outputs(args, summary["run_dir"])
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    cli_args = get_args()
    main(cli_args)
    print(f"Finish sample_name: {cli_args.sample_name} .....")
