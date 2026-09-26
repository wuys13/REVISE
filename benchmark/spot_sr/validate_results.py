"""Validate the complete 3-region x 4-size x 2-method result set."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    errors = []
    for part in ("part1", "part2", "part3"):
        for size in (50, 100, 150, 200):
            case = args.prepared / part / f"spot_{size}"
            manifest = json.loads((case / "case_manifest.json").read_text())
            expected_genes = pd.read_csv(case / "cnts.tsv", sep="\t", nrows=0).columns[1:]
            for method in ("tesla", "istar"):
                label = f"{part}/spot_{size}/{method}"
                try:
                    run = json.loads((case / f"{method}_run.json").read_text())
                    summary = json.loads((case / f"{method}_metrics_summary.json").read_text())
                    metric = pd.read_csv(case / f"{method}_metrics_normalized.csv")
                    pred = ad.read_h5ad(case / f"{method}_predicted_cells.h5ad", backed="r")
                    if pred.n_obs != manifest["n_truth_cells"]:
                        raise ValueError("Predicted cell count mismatch")
                    if pred.n_vars != manifest["n_genes"] or pred.n_vars != 324:
                        raise ValueError("Predicted gene count mismatch")
                    if set(pred.var_names) != set(expected_genes):
                        raise ValueError("Predicted gene panel mismatch")
                    if len(metric) != pred.n_vars or set(metric["Gene"]) != set(pred.var_names):
                        raise ValueError("Metric gene panel mismatch")
                    if not np.isfinite(metric[["PCC", "SSIM", "MSE", "NRMSE"]].to_numpy()).all():
                        raise ValueError("Nonfinite metric")
                    for name in ("PCC", "SSIM", "MSE", "NRMSE"):
                        if not np.isclose(metric[name].mean(), summary["metrics"][name], atol=1e-8):
                            raise ValueError(f"Summary and per-gene {name} disagree")
                    if method == "istar" and (run.get("epochs"), run.get("n_states")) != (400, 5):
                        raise ValueError("Unexpected iStar training configuration")
                    if method == "tesla" and run.get("grid_microns") != 10:
                        raise ValueError("Unexpected TESLA grid size")
                    rows.append({"case": label, "cells": pred.n_obs, "genes": pred.n_vars,
                                 "PCC": summary["metrics"]["PCC"], "run_seconds": run["runtime_seconds"]})
                    pred.file.close()
                except Exception as exc:
                    errors.append({"case": label, "error": str(exc)})
    report = {"expected_runs": 24, "validated_runs": len(rows), "errors": errors,
              "results": rows}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: report[key] for key in ("expected_runs", "validated_runs", "errors")}, indent=2))
    if len(rows) != 24 or errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
