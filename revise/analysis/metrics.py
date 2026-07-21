from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from scipy.stats import pearsonr
from skimage.metrics import structural_similarity as ssim
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score


def _to_numpy_matrix(x):
    if issparse(x):
        arr = x.toarray()
    elif isinstance(x, pd.DataFrame):
        arr = x.values
    else:
        arr = np.asarray(x)

    if hasattr(arr, "A"):
        arr = np.asarray(arr)

    arr = arr.astype(float, copy=False)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    return arr


def normalize_data(data):
    min_vals = np.min(data, axis=0, keepdims=True)
    max_vals = np.max(data, axis=0, keepdims=True)

    range_vals = max_vals - min_vals
    range_vals[range_vals == 0] = 1

    normalized_data = (data - min_vals) / range_vals
    return normalized_data


def compute_metric(
    adata_gt,
    adata_pred,
    logger,
    adata_process=False,
    sample_ratio=None,
    gene_list=None,
    normalize=False,
):
    """Compute per-gene expression reconstruction metrics.

    This preserves the compatibility metric contract (PCC/SSIM/MSE/NRMSE) so parity
    reports and benchmark evaluation outputs remain directly comparable.
    """
    adata_gt = adata_gt.copy()
    adata_pred = adata_pred.copy()

    sc.pp.normalize_total(adata_gt, target_sum=1e4)
    sc.pp.normalize_total(adata_pred, target_sum=1e4)

    if sample_ratio is not None:
        logger.info("Sampling %s%% of cells.", sample_ratio * 100)
        indices = np.random.choice(adata_gt.n_obs, int(adata_gt.n_obs * sample_ratio), replace=False)
        adata_gt = adata_gt[indices, :]
        adata_pred = adata_pred[indices, :]

    if gene_list is not None:
        logger.info("Using %s genes.", len(gene_list))
        adata_gt = adata_gt[:, gene_list]
        adata_pred = adata_pred[:, gene_list]
    else:
        overlap_genes = list(adata_gt.var_names.intersection(adata_pred.var_names))
        if len(overlap_genes) == 0:
            raise ValueError("No overlap genes between ground-truth and prediction adatas.")
        adata_gt = adata_gt[:, overlap_genes]
        adata_pred = adata_pred[:, overlap_genes]

    if adata_process:
        logger.info("Normalizing and log-transforming data.")
        sc.pp.log1p(adata_gt)
        sc.pp.log1p(adata_pred)

    x_gt = _to_numpy_matrix(adata_gt.X)
    x_pred = _to_numpy_matrix(adata_pred.X)

    if normalize:
        x_gt = normalize_data(x_gt)
        x_pred = normalize_data(x_pred)

    pcc_values = []
    ssim_values = []
    mse_values = []
    nrmse_values = []

    genes = adata_gt.var_names

    for i, _gene in enumerate(genes):
        expr_gt = x_gt[:, i]
        expr_pred = x_pred[:, i]

        pcc, _ = pearsonr(expr_gt, expr_pred)
        pcc_values.append(pcc)

        if normalize:
            data_range = 1
        else:
            data_range = expr_gt.max() - expr_gt.min()
        ssim_value = ssim(expr_gt, expr_pred, data_range=data_range)
        ssim_values.append(ssim_value)

        mse = np.mean((expr_gt - expr_pred) ** 2)
        mse_values.append(mse)

        nrmse = np.sqrt(mse) / np.mean(expr_gt)
        nrmse_values.append(nrmse)

    metrics_df = pd.DataFrame(
        {
            "Gene": genes,
            "PCC": pcc_values,
            "SSIM": ssim_values,
            "MSE": mse_values,
            "NRMSE": nrmse_values,
        }
    )

    logger.info("PCC: %.4f, %.4f, %.4f", np.nanmean(pcc_values), np.nanmax(pcc_values), np.nanmin(pcc_values))
    logger.info("SSIM: %.4f, %.4f, %.4f", np.nanmean(ssim_values), np.nanmax(ssim_values), np.nanmin(ssim_values))
    logger.info("MSE: %.4f, %.4f, %.4f", np.nanmean(mse_values), np.nanmax(mse_values), np.nanmin(mse_values))
    logger.info(
        "NRMSE: %.4f, %.4f, %.4f",
        np.nanmean(nrmse_values),
        np.nanmax(nrmse_values),
        np.nanmin(nrmse_values),
    )

    return metrics_df


def compute_clustering_metrics(adata, pred_label_key, true_label_key):
    """Compute ARI/NMI from predicted and reference labels in `adata.obs`."""
    pred_labels = adata.obs[pred_label_key].values
    true_labels = adata.obs[true_label_key].values
    true_labels = pd.Categorical(true_labels).codes
    pred_labels = pd.Categorical(pred_labels).codes
    ari = adjusted_rand_score(true_labels, pred_labels)
    nmi = normalized_mutual_info_score(true_labels, pred_labels)
    return ari, nmi
