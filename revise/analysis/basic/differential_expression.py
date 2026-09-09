from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc


def get_degs(
    adata,
    groupby,
    method="t-test",
    fc_threshold=None,
    reference="rest",
):
    """Run differential expression and return a tidy statistics table."""
    adata = adata.copy()
    adata.obs[groupby] = adata.obs[groupby].astype("category")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    print("Conducting differential expression analysis...")
    sc.tl.rank_genes_groups(
        adata,
        reference=reference,
        groupby=groupby,
        method=method,
        use_raw=False,
    )

    result = pd.DataFrame()
    categories = adata.obs[groupby].cat.categories.tolist()
    if reference in categories:
        categories.remove(reference)
    for group in categories:
        group_df = pd.DataFrame(
            {
                "group": group,
                "gene": adata.uns["rank_genes_groups"]["names"][group],
                "logfoldchanges": adata.uns["rank_genes_groups"][
                    "logfoldchanges"
                ][group],
                "pvals": adata.uns["rank_genes_groups"]["pvals"][group],
                "pvals_adj": adata.uns["rank_genes_groups"]["pvals_adj"][group],
            }
        )
        result = pd.concat([result, group_df])

    if fc_threshold is not None:
        if fc_threshold > 0:
            result = result[result["logfoldchanges"] > fc_threshold]
        else:
            result = result[result["logfoldchanges"] < fc_threshold]

    result["log_q"] = -np.log10(result["pvals_adj"] + 1e-100)
    result.sort_values(by="log_q", ascending=False, inplace=True)
    result.reset_index(drop=True, inplace=True)
    return result


__all__ = ["get_degs"]
