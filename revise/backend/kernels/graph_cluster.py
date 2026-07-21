import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import squidpy as sq
from anndata import AnnData
from sklearn.metrics.cluster import adjusted_rand_score, normalized_mutual_info_score
from tqdm import tqdm

from revise.backend.kernels.base import BaseKernel
from revise.backend.ops.posterior_conditioning import (
    condition_sparse_graph,
    get_posterior_matrix,
    posterior_conditioning_enabled,
    posterior_conditioning_mode,
)


class GraphClusterKernel(BaseKernel):
    """
    Graph-based clustering with spatial and alignment score evaluation.
    
    This class performs Leiden clustering at multiple resolutions and evaluates
    clustering quality using spatial coherence and alignment with cell type labels.
    """
    def __init__(self, config, logger):
        super().__init__(config, logger)

    def run(self, adata: AnnData, resolution, label):
        """
        Perform graph-based clustering and evaluate at multiple resolutions.
        
        Args:
            adata: AnnData object to cluster
            resolution: List of resolution values for Leiden clustering
            label: Column name in adata.obsm containing soft cell type labels
                (used for computing alignment scores)
                
        Returns:
            tuple: (adata, merge_df, best_res)
                - adata: AnnData with clustering results in obs
                - merge_df: DataFrame with metrics for each resolution
                - best_res: Best resolution based on alignment score
                
        For each resolution, computes:
        - Spatial score: Number of spatial neighbors with same cluster label
        - Alignment score: Agreement between clusters and cell type labels
        """
        adata = adata.copy()
        adata_raw = adata.copy()
        # Keep graph construction and Leiden partitioning deterministic so
        # old/new wrappers can be compared with strict equality.
        random_state = int(getattr(self.config, "rec_random_state", 0))

        def _log_graph_stats(name, graph):
            nnz = getattr(graph, "nnz", None)
            shape = getattr(graph, "shape", None)
            try:
                total = float(graph.sum())
            except Exception:
                total = None
            print(f"{name}: shape={shape} nnz={nnz} sum={total}")
        def _log_graph_diff(name, left, right):
            try:
                diff = left - right
            except Exception:
                print(f"{name}: diff=unavailable")
                return
            if hasattr(diff, "eliminate_zeros"):
                diff.eliminate_zeros()
            nnz = getattr(diff, "nnz", None)
            if nnz:
                max_abs = float(np.abs(diff.data).max())
                abs_sum = float(np.abs(diff.data).sum())
            else:
                max_abs = 0.0
                abs_sum = 0.0
            print(f"{name}: nnz={nnz} abs_sum={abs_sum} max_abs={max_abs}")

        if adata.n_obs < 2:
            raise ValueError("Graph clustering requires at least two spatial cells")

        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=min(100, adata.n_vars))
        if not bool(adata.var["highly_variable"].any()):
            adata.var["highly_variable"] = True
        adata = adata[:, adata.var['highly_variable']]
        n_pcs = min(30, max(0, adata.n_obs - 1), max(0, adata.n_vars - 1))
        if n_pcs >= 1:
            sc.pp.pca(adata, n_comps=n_pcs, random_state=random_state)
            sc.pp.neighbors(adata, n_pcs=n_pcs, random_state=random_state)
        else:
            sc.pp.neighbors(adata, use_rep="X", random_state=random_state)
        print(
            "Preprocess: n_obs=%s n_vars=%s n_hvg=%s n_pcs=%s"
            % (adata.n_obs, adata.n_vars, int(adata.var['highly_variable'].sum()), n_pcs)
        )
        nn_graph_genes = adata.obsp["connectivities"]
        # spatial proximity graph
        sq.gr.spatial_neighbors(adata)
        nn_graph_space = adata.obsp["spatial_connectivities"]
        joint_graph = (1 - self.config.rec_graph_alpha) * nn_graph_genes + self.config.rec_graph_alpha * nn_graph_space

        if self.config.rec_graph_method == "pca":
            adjacency_graph = nn_graph_genes
        elif self.config.rec_graph_method == "spatial":
            adjacency_graph = nn_graph_space
        elif self.config.rec_graph_method == "joint":
            adjacency_graph = joint_graph
        else:
            raise ValueError("neighbors_method must be pca or spatial")

        pc_mode = posterior_conditioning_mode(self.config)
        if pc_mode != "off":
            posterior_key = getattr(self.config, "posterior_conditioning_key", label)
            q = get_posterior_matrix(adata, posterior_key)
            if q is None and posterior_key != label:
                q = get_posterior_matrix(adata, label)
            if q is None:
                self.logger.warning(
                    "Posterior conditioning requested for graph clustering but neither obsm[%r] nor obsm[%r] "
                    "is available; using the unconditioned graph.",
                    posterior_key,
                    label,
                )
            elif posterior_conditioning_enabled(self.config, "cost"):
                adjacency_graph = condition_sparse_graph(
                    adjacency_graph,
                    q,
                    beta=getattr(self.config, "posterior_conditioning_beta", 1.0),
                    min_affinity=getattr(self.config, "posterior_conditioning_min_affinity", 0.05),
                )
                adata.obsp["posterior_conditioned_connectivities"] = adjacency_graph
                self.logger.info(
                    "Applied posterior-conditioned graph edges with key=%s, mode=%s",
                    posterior_key if posterior_key in adata.obsm else label,
                    pc_mode,
                )
            elif posterior_conditioning_enabled(self.config, "reference"):
                self.logger.info(
                    "posterior_conditioning.mode=%s has no entropic-reference analogue for graph clustering; "
                    "leaving graph unchanged.",
                    pc_mode,
                )

        print(
            "Graph config: method=%s alpha=%s resolutions=%s random_state=%s"
            % (self.config.rec_graph_method, self.config.rec_graph_alpha, resolution, random_state)
        )
        _log_graph_stats("nn_graph_genes", nn_graph_genes)
        _log_graph_stats("nn_graph_space", nn_graph_space)
        _log_graph_stats("adjacency_graph", adjacency_graph)
        _log_graph_diff("diff_genes_space", nn_graph_genes, nn_graph_space)
        _log_graph_diff("diff_adj_genes", adjacency_graph, nn_graph_genes)
        _log_graph_diff("diff_adj_space", adjacency_graph, nn_graph_space)

        # adjacency_graph = get_adjacency_graph(
        #     adata,
        #     "sc_app",
        #     self.config.rec_graph_method,
        #     self.config.rec_graph_alpha,
        #     self.config.rec_graph_exp_neighbor_num,
        #     self.config.rec_graph_spatial_neighbor_num
        # )

        merge_df = pd.DataFrame()

        for res in tqdm(resolution, desc="leiden"):
            sc.tl.leiden(
                adata,
                adjacency=adjacency_graph,
                resolution=res,
                key_added=f"leiden_{res}",
                random_state=random_state,
            )
            adata = get_spatial_score(adata, res=res)
            from revise.backend.ops.coefficients import get_weighted_align_score
            align_score = get_weighted_align_score(adata, res=res, label=label)
            mean_score = np.mean(adata.obs[f"spatial_score_{res}"])
            cluster_num = len(np.unique(adata.obs[f"leiden_{res}"]))
            # weighted_score = mean_score * np.log(cluster_num)
            df = pd.DataFrame({
                "resolution": res,
                "cluster_num": cluster_num,
                "mean_score": mean_score,
                "align_score": align_score,
            }, index=[0])
            self.logger.info(f"Resolution {res}: {cluster_num} clusters mean spatial score: {mean_score:.4f} {align_score}...")
            # Plotting is optional and should never block clustering.
            # Newer scanpy wrappers can raise on list-valued `color`, so we
            # render each panel separately and guard with best-effort fallback.
            if bool(getattr(self.config, "plot_flag", False)):
                try:
                    sc.pl.scatter(adata, x="x", y="y", color=f"leiden_{res}", show=False)
                    sc.pl.scatter(adata, x="x", y="y", color=f"spatial_score_{res}", show=False)
                except Exception as exc:
                    self.logger.warning("Skip cluster scatter plotting at resolution %s: %s", res, exc)
            merge_df = pd.concat([merge_df, df], axis=0)
        
        best_res = merge_df[merge_df["align_score"] == merge_df["align_score"].max()]["resolution"].values[-1]
        print(f"Best resolution: {best_res}")
        merge_df.reset_index(drop=True, inplace=True)
        adata_raw.obs = adata.obs.copy()
        return adata_raw, merge_df, best_res


def get_spatial_score(adata, res):

    nn_graph_space = adata.obsp["spatial_connectivities"]

    labels = adata.obs[f"leiden_{res}"].to_numpy()

    edges = nn_graph_space.tocoo(copy=True)
    edges.sum_duplicates()
    same_label = labels[edges.row] == labels[edges.col]
    spatial_count = np.bincount(
        edges.row,
        weights=edges.data * same_label,
        minlength=nn_graph_space.shape[0],
    )

    adata.obs[f"spatial_score_{res}"] = spatial_count

    return adata


def get_weighted_align_score(adata, res, label="Level2"):
    # leiden cluster
    leiden_class = adata.obs[f'leiden_{res}'].to_numpy()
    unique_class = np.unique(leiden_class)

    # soft label
    gene_scores_raw = adata.obsm[label]
    gene_scores = gene_scores_raw.to_numpy()  # shape: (n_cells, n_classes)
    col_names = list(gene_scores_raw.columns) if hasattr(gene_scores_raw, "columns") else None

    align_score = 0.0
    dominant_labels = []
    cluster_sizes = []

    for class_name in unique_class:
        idx = np.where(leiden_class == class_name)[0]  # index of cells in this cluster

        if len(idx) == 0:
            continue

        cluster_scores = gene_scores[idx]  # shape: (n_cells_in_cluster, n_classes)

        sum_scores = cluster_scores.sum(axis=0)  # shape: (n_classes,)

        # max over class
        max_sum = sum_scores.max()
        max_idx = int(sum_scores.argmax())
        dominant_labels.append(max_idx)
        cluster_sizes.append(len(idx))
        dom_name = col_names[max_idx] if col_names is not None else str(max_idx)

        align_score += max_sum

    if len(dominant_labels) > 0:
        dom_series = pd.Series(dominant_labels).value_counts().sort_index()
        if col_names is not None:
            dom_series.index = [col_names[i] for i in dom_series.index]

    align_score = align_score / len(leiden_class)
    align_score = round(align_score, 4)

    return align_score


def get_entropy_align_score(adata, res, label="Level2"):
    """
    Calculate alignment score using 1 - entropy(sum_scores) per cluster.

    For each Leiden cluster, compute sum_scores across soft labels and use
    1 - entropy(sum_scores) as the cluster purity. The final score is
    weighted by cluster size and normalized by total cells.
    """
    leiden_class = adata.obs[f'leiden_{res}'].to_numpy()
    unique_class = np.unique(leiden_class)
    gene_scores = adata.obsm[label].to_numpy()  # shape: (n_cells, n_classes)

    align_score = 0.0
    total_cells = len(leiden_class)
    if total_cells == 0:
        return 0.0

    for class_name in unique_class:
        idx = np.where(leiden_class == class_name)[0]
        if len(idx) == 0:
            continue
        cluster_scores = gene_scores[idx]
        sum_scores = cluster_scores.sum(axis=0)
        if np.all(sum_scores == 0):
            continue
        ent = scipy.stats.entropy(sum_scores)
        cluster_score = 1.0 - ent
        align_score += cluster_score * len(idx)

    align_score = align_score / total_cells
    align_score = round(float(align_score), 4)
    return align_score


def clustering_metrics(adata, pred_label_key, true_label_key):
    pred_labels = adata.obs[pred_label_key].values
    true_labels = adata.obs[true_label_key].values

    true_labels = pd.Categorical(true_labels).codes
    pred_labels = pd.Categorical(pred_labels).codes

    # print(len(np.unique(pred_labels)), len(np.unique(true_labels)))
    ari = adjusted_rand_score(true_labels, pred_labels)
    nmi = normalized_mutual_info_score(true_labels, pred_labels)
    # print(f"ARI: {ari:.4f}, NMI: {nmi:.4f}", len(np.unique(pred_labels)), len(np.unique(true_labels)))
    return ari, nmi


def get_align_score(adata, res, label = "Level2"):
    leiden_class = adata.obs[f'leiden_{res}']
    unique_class = np.unique(leiden_class)

    align_score = 0
    for class_name in unique_class:
        class_label = adata.obs[label][leiden_class == class_name]
        max_label_num = class_label.value_counts().max()
        align_score += max_label_num
        
    align_score = align_score / len(leiden_class)
    align_score = align_score.round(4)

    return align_score
