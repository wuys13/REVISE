import numpy as np
import pandas as pd
import squidpy as sq
import scanpy as sc
import scipy

from tqdm import tqdm


def get_spatial_score(adata, res):

    nn_graph_space = adata.obsp["spatial_connectivities"]

    labels = adata.obs[f"leiden_{res}"].to_numpy()


    same_label_matrix = (labels[:, None] == labels[None, :]).astype(int)
    same_label_matrix = scipy.sparse.csr_matrix(same_label_matrix)
    spatial_count = (nn_graph_space.multiply(same_label_matrix)).sum(axis=1).A1
    # spatial_count = (nn_graph_space.dot(labels) == labels).sum(axis=1)
    
    adata.obs[f"spatial_score_{res}"] = spatial_count

    return adata


def get_align_score(adata, res, label = "Level2"):
    leiden_class = adata.obs[f'leiden_{res}']
    unique_class = np.unique(leiden_class)
    # 计算每个class数目最多label的数量
    align_score = 0
    for class_name in unique_class:
        class_label = adata.obs[label][leiden_class == class_name]
        max_label_num = class_label.value_counts().max()
        align_score += max_label_num
        
    align_score = align_score / len(leiden_class)
    align_score = align_score.round(4)

    return align_score

import numpy as np

def get_weighted_align_score(adata, res, label="Level2"):
    # leiden cluster
    leiden_class = adata.obs[f'leiden_{res}'].to_numpy()
    unique_class = np.unique(leiden_class)

    # soft label
    gene_scores = adata.obsm[label].to_numpy()  # shape: (n_cells, n_classes)

    align_score = 0.0

    for class_name in unique_class:
        idx = np.where(leiden_class == class_name)[0]  # index of cells in this cluster

        if len(idx) == 0:
            continue

        cluster_scores = gene_scores[idx]  # shape: (n_cells_in_cluster, n_classes)

        sum_scores = cluster_scores.sum(axis=0)  # shape: (n_classes,)

        # max over class
        max_sum = sum_scores.max()

        align_score += max_sum

    align_score = align_score / len(leiden_class)
    align_score = round(align_score, 4)

    return align_score




def get_best_cluster(adata, neighbors_method = "pca", alpha = 0.2, resolutions = [0.1, 0.3, 0.5, 0.7, 0.9], label = "Level2"):
    adata = adata.copy()
    adata_raw = adata.copy()

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

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=100)
    adata = adata[:, adata.var['highly_variable']]
    sc.pp.pca(adata, n_comps=30)
    print(
        "Preprocess: n_obs=%s n_vars=%s n_hvg=%s n_pcs=%s"
        % (adata.n_obs, adata.n_vars, int(adata.var['highly_variable'].sum()), 30)
    )

    sc.pp.neighbors(adata, n_pcs=30)
    nn_graph_genes = adata.obsp["connectivities"]
    # spatial proximity graph
    sq.gr.spatial_neighbors(adata)
    nn_graph_space = adata.obsp["spatial_connectivities"]
    joint_graph = (1 - alpha) * nn_graph_genes + alpha * nn_graph_space

    if neighbors_method == "pca":
        adjacency_graph = nn_graph_genes
    elif neighbors_method == "spatial":
        adjacency_graph = nn_graph_space
    elif neighbors_method == "joint":
        adjacency_graph = joint_graph
    else:
        raise ValueError("neighbors_method must be pca or spatial")

    print(
        "Graph config: method=%s alpha=%s resolutions=%s"
        % (neighbors_method, alpha, resolutions)
    )
    _log_graph_stats("nn_graph_genes", nn_graph_genes)
    _log_graph_stats("nn_graph_space", nn_graph_space)
    _log_graph_stats("adjacency_graph", adjacency_graph)
    _log_graph_diff("diff_genes_space", nn_graph_genes, nn_graph_space)
    _log_graph_diff("diff_adj_genes", adjacency_graph, nn_graph_genes)
    _log_graph_diff("diff_adj_space", adjacency_graph, nn_graph_space)
    
    merge_df = pd.DataFrame()
    # count = 0
    # previous_align_score = 0
    for res in tqdm(resolutions, desc = "leiden"):
        sc.tl.leiden(adata, adjacency=adjacency_graph, resolution=res, key_added=f"leiden_{res}" )
    
        adata = get_spatial_score(adata, res = res)
        # align_score = get_align_score(adata, res = res, label = label)
        align_score = get_weighted_align_score(adata, res = res, label = label)
        mean_score = np.mean(adata.obs[f"spatial_score_{res}"])
        cluster_num = len(np.unique(adata.obs[f"leiden_{res}"]))
        # weighted_score = mean_score * np.log(cluster_num)
        df = pd.DataFrame({
            "resolution": res,
            "cluster_num": cluster_num,
            "mean_score": mean_score,
            "align_score": align_score,
        }, index = [0])
        
    
        print(f"Resolution {res}: {cluster_num} clusters mean spatial score: {mean_score:.4f} {align_score}...")

        sc.pl.scatter(adata, x = "x", y = "y", color = [f"leiden_{res}", f"spatial_score_{res}"] )
        
        # increase_rate = (align_score - previous_align_score) / align_score
        # if increase_rate <= 0.01:
        #     count  = count + 1

        # previous_align_score = align_score
        # if count >= 2:
        #     print(f"Breaking at resolution {res}")
        #     break
        
        merge_df = pd.concat([merge_df, df], axis = 0)

    
    best_res = merge_df[merge_df["align_score"] == merge_df["align_score"].max()]["resolution"].values[-1]
    print(f"Best resolution: {best_res}")
    merge_df.reset_index(drop = True, inplace=True)
    adata_raw.obs = adata.obs.copy()
    
    return adata_raw, merge_df, best_res

import matplotlib.pyplot as plt
def plot_sc_SVC(adata, color, title = None, file_name = None):

    plt.figure(figsize=(10, 8*len(color)))
    sc.pl.scatter(adata, x="x", y="y", 
            color = color,
            title=title, show = False,
            )
    plt.savefig(file_name, dpi = 300)
    plt.close()
    


def split_adata_evenly(adata_sp, split_size = 20000, seed = 42):
    np.random.seed(seed)
    total_cells = adata_sp.n_obs
    num_samples = total_cells // split_size
    if num_samples == 0:
        num_samples = 1

    shuffled_indices = np.random.permutation(adata_sp.obs.index)

    sample_indices = np.array_split(shuffled_indices, num_samples)

    adata_samples = [adata_sp[sample_index].copy() for sample_index in sample_indices]
    print(f"Split {total_cells} cells into {num_samples} samples")

    return adata_samples


from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

def get_cm(adata, row, col):
    grouped = adata.obs.groupby([row, col]).size()
    cm_df = grouped.unstack(fill_value=0)
    
    return cm_df
def plot_cm(cm_df, save_dir=None):
        
    plt.figure(figsize=(10, 8))
    
    sns.heatmap(cm_df, 
                    annot=True, 
                    fmt='',
                    cmap='Reds',
                    )

    
    plt.title('Confusion Matrix')
    plt.xlabel('Expert anno')
    plt.ylabel('sc_SVC cluster')
    
    if save_dir is not None:
        plt.savefig(f"{save_dir}/cm.pdf", bbox_inches='tight')
        plt.close()
    else:
        plt.show()
