#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

from revise.framework import REVISEPipeline


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser("REVISE application sp-SVC (unified interface wrapper)")
    parser.add_argument("--config", default="revise/revise.yaml", help="Path to unified config")
    parser.add_argument("--platform", default="hST", choices=["hST"], help="Application platform route")
    parser.add_argument("--confounding", default="bin2cell", choices=["bin2cell"], help="Route confounding factor")
    parser.add_argument("--profile", default="application_sp", help="Unified config profile")
    parser.add_argument("--data-root", required=True, type=str, help="Data root path")
    parser.add_argument("--sample-name", required=True, type=str, help="Sample name of datasets")
    parser.add_argument("--st-file", required=True, type=str, help="Spatial file name")
    parser.add_argument("--sc-ref-file", required=True, type=str, help="Single-cell reference file name")
    parser.add_argument("--patient-key", default="Patient", type=str, help="Patient key in sc-ref-file")
    parser.add_argument("--sample-size", default=None, type=int, help="Optional ST subsample size")
    parser.add_argument(
        "--output-root",
        default="output/sp_SVC_case",
        type=str,
        help="Unified output root path",
    )
    parser.add_argument(
        "--compatibility-mode",
        action="store_true",
        help="Also emit notebook-compatible filenames inside the run directory",
    )
    return parser.parse_args()


def _publish_notebook_outputs(args: argparse.Namespace, run_dir: str | None) -> dict[str, str]:
    if not run_dir:
        raise RuntimeError("pipeline finished without a run_dir")

    run_path = Path(run_dir)
    candidates = [
        run_path / "sp_SVC.h5ad",
        run_path / "artifacts" / "sp_svc.h5ad",
    ]
    src = next((path for path in candidates if path.exists()), None)
    if src is None:
        searched = ", ".join(str(path) for path in candidates)
        raise FileNotFoundError(f"sp-SVC output was not found; searched: {searched}")

    dst_dir = Path(args.output_root) / args.sample_name
    dst_dir.mkdir(parents=True, exist_ok=True)
    dst = dst_dir / "sp_SVC.h5ad"
    shutil.copy2(src, dst)
    return {"sp_SVC": str(dst)}


def main(args: argparse.Namespace) -> None:
    pipeline = REVISEPipeline(config_path=args.config)
    svc = pipeline.run(
        profile=args.profile,
        runtime_overrides={
            "platform": args.platform,
            "confounding": args.confounding,
            "compatibility_mode": bool(args.compatibility_mode),
        },
        io_overrides={
            "data_root": args.data_root,
            "output_root": args.output_root,
            "sample_name": args.sample_name,
            "st_file": args.st_file,
            "sc_ref_file": args.sc_ref_file,
            "patient_key": args.patient_key,
            "sample_size": args.sample_size,
        },
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
