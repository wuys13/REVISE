from __future__ import annotations

import os
from pathlib import Path

import anndata as ad
import numpy as np


def data_root() -> Path:
    return Path(os.environ["REVISE_TACCO_DATA_ROOT"]).expanduser().resolve()


def _counts_adata(source: ad.AnnData, *, spatial: bool = False) -> ad.AnnData:
    counts = source.layers["counts"].copy()
    prepared = ad.AnnData(
        X=counts,
        obs=source.obs.copy(),
        var=source.var.copy(),
    )
    prepared.obs["transcript_counts"] = np.asarray(counts.sum(axis=1)).ravel()
    if spatial:
        prepared.obsm["spatial"] = source.obs[["spatial_x", "spatial_y"]].to_numpy()
    return prepared


def prepare() -> Path:
    root = data_root()
    raw_dir = root / "raw" / "zesta_zf5"
    prepared_dir = root / "prepared" / "zesta_zf5"
    st_source = ad.read_h5ad(raw_dir / "zf5_stereoseq.h5ad")
    reference_source = ad.read_h5ad(raw_dir / "zf5_scRNA.h5ad")
    st = _counts_adata(st_source, spatial=True)
    reference = _counts_adata(reference_source)
    prepared_dir.mkdir(parents=True, exist_ok=True)
    st.write_h5ad(prepared_dir / "zf5_stereoseq_counts.h5ad", compression="gzip")
    reference.write_h5ad(prepared_dir / "zf5_scRNA_counts.h5ad", compression="gzip")
    return prepared_dir


if __name__ == "__main__":
    prepare()
