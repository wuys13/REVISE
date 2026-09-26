"""Reconstruct a shared P1CRC reference from released per-region subsets.

The upstream REVISE spot config expects `real_sc_ref_all.h5ad` in every part,
but this copy of the data contains only `real_sc_ref_part.h5ad`. The union of
their cell IDs provides a usable reference for the three regions. It has not
been proven byte-identical to the original full reference.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spot-root", type=Path, required=True)
    args = parser.parse_args()
    refs = [ad.read_h5ad(args.spot_root / f"part{i}" / "real_sc_ref_part.h5ad")
            for i in (1, 2, 3)]
    for ref in refs[1:]:
        if not refs[0].var_names.equals(ref.var_names):
            raise ValueError("Reference gene panels differ")
    # The released part3 reference is a subset of part1; part2 adds one cell
    # type absent from part1. Verify overlapping records before deduplicating.
    for i in range(len(refs)):
        for j in range(i + 1, len(refs)):
            common = refs[i].obs_names.intersection(refs[j].obs_names)
            if not refs[i][common].obs.equals(refs[j][common].obs):
                raise ValueError(f"Reference annotations differ for parts {i+1}, {j+1}")
            # Compare a deterministic subset of expression rows; full matrices
            # are otherwise copied unchanged from the released parts.
            check = common[: min(100, len(common))]
            delta = refs[i][check].X - refs[j][check].X
            unequal = delta.nnz if sparse.issparse(delta) else np.count_nonzero(delta)
            if unequal:
                raise ValueError(f"Reference expression differs for parts {i+1}, {j+1}")
    seen: set[str] = set()
    pieces = []
    for ref in refs:
        keep = np.array([name not in seen for name in ref.obs_names], dtype=bool)
        pieces.append(ref[keep].copy())
        seen.update(ref.obs_names[keep])
    combined = ad.concat(pieces, join="inner", merge="same")
    if combined.n_obs != len(seen) or not combined.obs_names.is_unique:
        raise ValueError("Deduplication failed")
    output = args.spot_root / "real_sc_ref_all.h5ad"
    if not output.exists():
        combined.write_h5ad(output, compression="gzip")
    for i in (1, 2, 3):
        link = args.spot_root / f"part{i}" / "real_sc_ref_all.h5ad"
        if not link.exists():
            link.symlink_to(Path("..") / output.name)
    patient_dir = args.spot_root / "P1CRC"
    patient_dir.mkdir(exist_ok=True)
    for i in (1, 2, 3):
        link = patient_dir / f"cut_part{i}"
        if not link.exists():
            link.symlink_to(Path("..") / f"part{i}")
    print({"output": str(output), "n_cells": combined.n_obs,
           "n_genes": combined.n_vars,
           "cell_types": pd.Series(combined.obs["Level1"]).value_counts().to_dict()})


if __name__ == "__main__":
    main()
