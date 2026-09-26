"""Collect completed baseline runs into a comparison table."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for case in sorted(args.prepared.glob("part*/spot_*")):
        manifest_file = case / "case_manifest.json"
        if not manifest_file.exists():
            continue
        manifest = json.loads(manifest_file.read_text())
        for method in ("tesla", "istar"):
            result = case / f"{method}_metrics_summary.json"
            run = case / f"{method}_run.json"
            if not (result.exists() and run.exists()):
                continue
            result_data = json.loads(result.read_text())
            run_data = json.loads(run.read_text())
            rows.append({
                "part": manifest["part"],
                "spot_size_microns": manifest["spot_size_microns"],
                "method": method,
                "n_spots": manifest["n_spots"],
                "n_cells": result_data["n_cells"],
                "n_genes": result_data["n_genes"],
                "runtime_seconds": run_data.get("runtime_seconds"),
                **result_data["metrics"],
            })
    table = pd.DataFrame(rows)
    if table.empty:
        raise RuntimeError("No completed runs found")
    table = table.sort_values(["part", "spot_size_microns", "method"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, index=False)
    grouped = table.groupby(["spot_size_microns", "method"])[
        ["PCC", "SSIM", "MSE", "NRMSE"]
    ].agg(["mean", "std"])
    grouped.columns = [f"{metric}_{stat}" for metric, stat in grouped.columns]
    grouped = grouped.reset_index()
    grouped["n_regions"] = table.groupby(["spot_size_microns", "method"]).size().to_numpy()
    grouped_path = args.output.with_name("baselines_by_size.csv")
    grouped.to_csv(grouped_path, index=False)
    print(table.to_string(index=False))
    print(f"Saved {len(table)} completed runs to {args.output} and {grouped_path}")


if __name__ == "__main__":
    main()
