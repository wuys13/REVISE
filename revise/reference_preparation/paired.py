"""Extract one exact paired subset without changing its expression content."""

from __future__ import annotations

import os
from pathlib import Path
import tempfile
from typing import Any

import anndata as ad


def extract_pair(
    source: Path,
    pair_column: str,
    pair_key: str,
    output_path: Path,
) -> dict[str, Any]:
    """Write rows whose ``.obs[pair_column]`` exactly equal ``pair_key``."""
    source = Path(source)
    output_path = Path(output_path)
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite {output_path}")

    source_adata = ad.read_h5ad(source)
    if pair_column not in source_adata.obs:
        raise KeyError(f"pair column is missing from source .obs: {pair_column}")
    selected_mask = source_adata.obs[pair_column] == pair_key
    if hasattr(selected_mask, "fillna"):
        selected_mask = selected_mask.fillna(False)
    selected_count = int(selected_mask.sum())
    if selected_count == 0:
        raise ValueError(f"no source cells exactly match {pair_column}={pair_key!r}")

    selected = source_adata[selected_mask].copy()
    descriptor, temporary_name = tempfile.mkstemp(
        dir=output_path.parent, prefix=f".{output_path.name}.", suffix=".h5ad",
    )
    os.close(descriptor)
    temporary_path = Path(temporary_name)
    try:
        selected.write_h5ad(temporary_path)
        try:
            os.link(temporary_path, output_path)
        except FileExistsError:
            raise FileExistsError(f"refusing to overwrite {output_path}") from None
    finally:
        temporary_path.unlink(missing_ok=True)

    return {
        "source_path": str(source),
        "pair_column": pair_column,
        "pair_key": pair_key,
        "source_cells": int(source_adata.n_obs),
        "selected_cells": selected_count,
        "genes": int(source_adata.n_vars),
    }
