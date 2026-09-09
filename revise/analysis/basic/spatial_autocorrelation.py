from __future__ import annotations

import pandas as pd
import scanpy as sc
import squidpy as sq
from anndata import AnnData


def compute_groupwise_spatial_autocorrelation(
    adata: AnnData,
    *,
    groupby: str,
    mode: str = "moran",
    min_obs: int = 51,
) -> pd.DataFrame:
    """Compute Moran I or Geary C after independently graphing each group."""
    if mode not in {"moran", "geary"}:
        raise ValueError("mode must be 'moran' or 'geary'")
    if groupby not in adata.obs:
        raise KeyError(f"Column {groupby!r} not found in adata.obs")
    if "spatial" not in adata.obsm:
        raise KeyError("Missing spatial coordinates: adata.obsm['spatial']")
    if min_obs < 1:
        raise ValueError("min_obs must be at least 1")

    group_values = sorted(pd.unique(adata.obs[groupby]), key=lambda value: str(value))
    subsets = [("All", adata)]
    subsets.extend(
        (group, adata[adata.obs[groupby] == group])
        for group in group_values
    )

    uns_key, statistic = ("moranI", "I") if mode == "moran" else ("gearyC", "C")
    columns: list[pd.Series] = []
    for group, subset in subsets:
        if subset.n_obs < min_obs:
            continue
        work = subset.copy()
        sc.pp.normalize_total(work, target_sum=1e4)
        sc.pp.log1p(work)
        sq.gr.spatial_neighbors(work)
        sq.gr.spatial_autocorr(
            work,
            mode=mode,
            genes=work.var_names,
            n_perms=1,
            n_jobs=1,
        )
        values = work.uns[uns_key][statistic].reindex(adata.var_names)
        values.name = group
        columns.append(values)

    return pd.concat(columns, axis=1) if columns else pd.DataFrame(index=adata.var_names)
