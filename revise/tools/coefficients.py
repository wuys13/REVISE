import numpy as np
import pandas as pd
import scipy
import squidpy as sq


def compute_autocorr_metrics(adata_raw, connectivity):
    """
    Calculate Moran's I and Geary's C
    """
    adata_raw.obsp["connectivity"] = connectivity

    sq.gr.spatial_autocorr(
        adata_raw,
        connectivity_key="connectivity",
        mode="moran",
        genes=list(adata_raw.var_names),
        copy=False,
    )
    moran_df = adata_raw.uns["moranI"]
    if "names" in moran_df.columns:
        moran_res = moran_df.set_index("names")["I"]
    else:
        moran_res = moran_df["I"]

    sq.gr.spatial_autocorr(
        adata_raw,
        connectivity_key="connectivity",
        mode="geary",
        genes=list(adata_raw.var_names),
        copy=False,
    )
    geary_df = adata_raw.uns["gearyC"]
    if "names" in geary_df.columns:
        geary_res = geary_df.set_index("names")["C"]
    else:
        geary_res = geary_df["C"]

    df = pd.concat(
        [moran_res.rename("moranI"), geary_res.rename("gearyC")],
        axis=1,
    )
    return df


def get_weighted_align_score(adata, res, label="Level2"):
    """
    Calculate weighted alignment score between clusters and cell type labels.
    
    This function measures how well Leiden clusters align with soft cell type
    labels by computing the maximum soft sum for each cluster and averaging
    across all cells.
    
    Args:
        adata: AnnData object with clustering and cell type labels
        res: Resolution value for Leiden clustering (used to find leiden_{res} column)
        label: Column name in adata.obsm containing soft cell type labels
            
    Returns:
        float: Alignment score (0-1, higher is better) rounded to 4 decimal places
    """
    leiden_class = adata.obs[f'leiden_{res}'].to_numpy()
    unique_class = np.unique(leiden_class)
    gene_scores = adata.obsm[label].to_numpy()

    align_score = 0.0

    for class_name in unique_class:
        idx = np.where(leiden_class == class_name)[0]

        if len(idx) == 0:
            continue

        cluster_scores = gene_scores[idx]
        sum_scores = cluster_scores.sum(axis=0)
        max_sum = sum_scores.max()

        align_score += max_sum

    align_score = align_score / len(leiden_class)
    align_score = round(align_score, 4)

    return align_score


def get_spatial_score(adata, res):
    """Calculate for each cell the number of spatial neighbors with the same leiden label"""
    if "spatial_connectivities" not in adata.obsp:
        sq.gr.spatial_neighbors(adata)
    nn_graph_space = adata.obsp["spatial_connectivities"]
    labels = adata.obs[f"leiden_{res}"].to_numpy()

    same_label_matrix = (labels[:, None] == labels[None, :]).astype(int)
    same_label_matrix = scipy.sparse.csr_matrix(same_label_matrix)
    spatial_count = (nn_graph_space.multiply(same_label_matrix)).sum(axis=0).A1

    adata.obs[f"spatial_score_{res}"] = spatial_count

    return adata
