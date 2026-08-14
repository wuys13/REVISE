from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from revise.analysis.metrics import compute_clustering_metrics


def compute_leiden_sweep(
    adata: AnnData,
    *,
    resolutions: Sequence[float],
    reference_col: str,
    random_state: int = 0,
) -> tuple[AnnData, pd.DataFrame]:
    """Run the preprocessing and Leiden sweep used by the sp-SVC analysis."""
    if reference_col not in adata.obs:
        raise KeyError(f"Column {reference_col!r} not found in adata.obs")

    resolution_values = list(resolutions)
    if not resolution_values or any(
        not np.isfinite(resolution) or resolution <= 0
        for resolution in resolution_values
    ):
        raise ValueError("resolutions must contain positive finite values")

    work = adata.copy()
    sc.pp.filter_cells(work, min_genes=50)
    sc.pp.filter_genes(work, min_cells=3)
    work = work[:, ~work.var_names.str.startswith("MT-")].copy()
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    if work.n_vars > 2000:
        sc.pp.highly_variable_genes(
            work,
            n_top_genes=2000,
            flavor="seurat_v3",
            subset=True,
        )
        work = work[:, work.var["highly_variable"]].copy()
    sc.tl.pca(work, svd_solver="arpack")
    sc.pp.neighbors(work, n_neighbors=10, n_pcs=40)

    rows = []
    for resolution in resolution_values:
        key = f"leiden_res_{resolution}"
        sc.tl.leiden(
            work,
            resolution=resolution,
            key_added=key,
            random_state=random_state,
        )
        ari, nmi = compute_clustering_metrics(work, key, reference_col)
        rows.append(
            {
                "resolution": resolution,
                "ARI": ari,
                "NMI": nmi,
                "cluster_num": work.obs[key].nunique(),
            }
        )

    return work, pd.DataFrame(
        rows,
        columns=["resolution", "ARI", "NMI", "cluster_num"],
    )
