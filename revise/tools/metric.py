import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from scipy.stats import pearsonr
from skimage.metrics import structural_similarity as ssim
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score


def _to_numpy_matrix(X):
    if issparse(X):
        arr = X.toarray()
    elif isinstance(X, pd.DataFrame):
        arr = X.values
    else:
        arr = np.asarray(X)

    if hasattr(arr, "A"):
        arr = np.asarray(arr)

    arr = arr.astype(float, copy=False)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    return arr


def compute_metric(
        adata_gt,
        adata_pred,
        logger,
        adata_process=False,
        sample_ratio=None,
        gene_list=None,
        normalize=False
):
    """
    Compute evaluation metrics between ground truth and predicted expression.
    
    This function calculates per-gene metrics including Pearson correlation
    coefficient (PCC), structural similarity index (SSIM), mean squared error
    (MSE), and normalized root mean squared error (NRMSE).
    
    Args:
        adata_gt: Ground truth AnnData object
        adata_pred: Predicted AnnData object
        logger: Logger instance for logging metrics
        adata_process: If True, apply log1p transformation before computing metrics
        sample_ratio: If provided, randomly sample this fraction of cells for evaluation
        gene_list: If provided, only evaluate metrics for these genes
        normalize: If True, normalize data to [0, 1] range before computing metrics
        
    Returns:
        pd.DataFrame: DataFrame with columns:
            - Gene: Gene names
            - PCC: Pearson correlation coefficient
            - SSIM: Structural similarity index
            - MSE: Mean squared error
            - NRMSE: Normalized root mean squared error
            
    The function logs summary statistics (mean, max, min) for each metric.
    """
    adata_gt = adata_gt.copy()
    adata_pred = adata_pred.copy()

    sc.pp.normalize_total(adata_gt, target_sum=1e4)
    sc.pp.normalize_total(adata_pred, target_sum=1e4)

    if sample_ratio is not None:
        logger.info(f"Sampling {sample_ratio * 100}% of cells.")
        indices = np.random.choice(adata_gt.n_obs, int(adata_gt.n_obs * sample_ratio), replace=False)
        adata_gt = adata_gt[indices, :]
        adata_pred = adata_pred[indices, :]

    if gene_list is not None:
        logger.info(f"Using {len(gene_list)} genes.")
        adata_gt = adata_gt[:, gene_list]
        adata_pred = adata_pred[:, gene_list]
    else:
        overlap_genes = list(adata_gt.var_names.intersection(adata_pred.var_names))
        if len(overlap_genes) == 0:
            import pdb; pdb.set_trace()
            print(overlap_genes)
        adata_gt = adata_gt[:, overlap_genes]
        adata_pred = adata_pred[:, overlap_genes]

    if adata_process:
        logger.info("Normalizing and log-transforming data.")
        sc.pp.log1p(adata_gt)
        sc.pp.log1p(adata_pred)

    X_gt = _to_numpy_matrix(adata_gt.X)
    X_pred = _to_numpy_matrix(adata_pred.X)

    if normalize:
        X_gt = normalize_data(X_gt)
        X_pred = normalize_data(X_pred)

    pcc_values = []
    ssim_values = []
    mse_values = []
    nrmse_values = []

    genes = adata_gt.var_names

    for i, gene in enumerate(genes):
        expr_gt = X_gt[:, i]
        expr_pred = X_pred[:, i]

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

    metrics_df = pd.DataFrame({'Gene': genes,
                               'PCC': pcc_values,
                               'SSIM': ssim_values,
                               "MSE": mse_values,
                               "NRMSE": nrmse_values})

    logger.info(f"PCC: {np.nanmean(pcc_values):.4f}, {np.nanmax(pcc_values):.4f}, {np.nanmin(pcc_values):.4f}")
    logger.info(f"SSIM: {np.nanmean(ssim_values):.4f}, {np.nanmax(ssim_values):.4f}, {np.nanmin(ssim_values):.4f}")
    logger.info(f"MSE: {np.nanmean(mse_values):.4f}, {np.nanmax(mse_values):.4f}, {np.nanmin(mse_values):.4f}")
    logger.info(
        f"NRMSE: {np.nanmean(nrmse_values):.4f}, {np.nanmax(nrmse_values):.4f}, {np.nanmin(nrmse_values):.4f}")

    return metrics_df


def normalize_data(data):
    min_vals = np.min(data, axis=0, keepdims=True)
    max_vals = np.max(data, axis=0, keepdims=True)

    range_vals = max_vals - min_vals
    range_vals[range_vals == 0] = 1

    normalized_data = (data - min_vals) / range_vals
    return normalized_data


def compute_clustering_metrics(adata, pred_label_key, true_label_key):
    """
    Compute clustering evaluation metrics (ARI and NMI).
    
    Args:
        adata: AnnData object with predicted and true labels in obs
        pred_label_key: Column name in adata.obs for predicted cluster labels
        true_label_key: Column name in adata.obs for true cell type labels
        
    Returns:
        tuple: (ari, nmi)
            - ari: Adjusted Rand Index (higher is better, range: -1 to 1)
            - nmi: Normalized Mutual Information (higher is better, range: 0 to 1)
    """
    pred_labels = adata.obs[pred_label_key].values
    true_labels = adata.obs[true_label_key].values
    true_labels = pd.Categorical(true_labels).codes
    pred_labels = pd.Categorical(pred_labels).codes
    ari = adjusted_rand_score(true_labels, pred_labels)
    nmi = normalized_mutual_info_score(true_labels, pred_labels)
    return ari, nmi
