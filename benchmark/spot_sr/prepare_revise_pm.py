"""Build the REVISE spot PM sidecar using its published Visium prior recipe.

This is a timing-compatible substitute when the released spot archive lacks
PM_on_cell.csv. It does not reproduce an unpublished original PM matrix.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spot-root", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--force", action="store_true", help="Replace an existing PM sidecar")
    args = parser.parse_args()

    spot_root = args.spot_root
    output = spot_root / "PM_on_cell.csv"
    if output.exists() and not args.force:
        raise FileExistsError(f"Refusing to replace {output}; use --force explicitly")
    reference = ad.read_h5ad(spot_root / "real_sc_ref_all.h5ad", backed="r")
    cell_types = sorted(reference.obs["Level1"].astype(str).unique())
    frequencies = (
        reference.obs["Level1"]
        .astype(str)
        .value_counts(normalize=True)
        .reindex(cell_types)
        .to_numpy()
    )
    cell_ids: set[str] = set()
    for part in ("part1", "part2", "part3"):
        selected = ad.read_h5ad(spot_root / part / "selected_xenium.h5ad", backed="r")
        cell_ids.update(map(str, selected.obs_names))
    sorted_ids = sorted(cell_ids)
    rng = np.random.default_rng(args.seed)
    probabilities = np.tile(frequencies, (len(sorted_ids), 1))
    probabilities += rng.uniform(0, 1e-6, size=probabilities.shape)
    probabilities /= probabilities.sum(axis=1, keepdims=True)
    pd.DataFrame(probabilities, index=sorted_ids, columns=cell_types).to_csv(output)
    (spot_root / "PM_on_cell.provenance.json").write_text(
        json.dumps(
            {
                "source": "REVISE Visium example prior recipe",
                "note": "Synthetic reference-frequency prior; not the original spot benchmark PM",
                "seed": args.seed,
                "n_cells": len(sorted_ids),
                "n_cell_types": len(cell_types),
            },
            indent=2,
        ) + "\n"
    )
    print(f"Wrote {output}: {len(sorted_ids)} cells x {len(cell_types)} cell types")


if __name__ == "__main__":
    main()
