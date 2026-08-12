import scanpy as sc
import squidpy as sq


def get_adjacency_graph(
        adata,
        data_type,
        neighbors_method="pca",
        alpha=0.2,
        gene_neighbor_num=30,
        spatial_neighbor_num=30):
    """
    Construct adjacency graph based on gene expression or spatial connectivity.
    
    This function builds a k-nearest neighbor graph using either:
    - Gene expression similarity (PCA-based)
    - Spatial proximity
    - A combination of both (joint)
    
    Args:
        adata: Input AnnData object. The adjacency matrix will be computed
            and stored in adata.obsp["connectivities"] or adata.obsp["spatial_connectivities"]
        data_type: Data type to determine hyperparameters ("sp" or "sc")
            - "sp": Uses 200 top genes, 30 PCs
            - "sc": Uses 2000 top genes, 50 PCs, includes scaling
        neighbors_method: Method for neighbor computation
            - "pca": Gene expression-based neighbors using PCA
            - "spatial": Spatial proximity-based neighbors
            - "joint": Weighted combination of both (1-alpha)*gene + alpha*spatial
        alpha: Weight for spatial connectivity in joint method (0-1)
        gene_neighbor_num: Number of neighbors for gene expression graph
        spatial_neighbor_num: Number of neighbors for spatial connectivity graph

    Returns:
        scipy.sparse.csr_matrix: Adjacency matrix of shape (n_obs, n_obs)
            representing connectivity between spots/cells
    """
    if data_type == "sp":
        n_top_genes = 200
        n_comps = 30
    elif data_type == "sc":
        n_top_genes = 2000
        n_comps = 50
    elif data_type == "sc_app":
        n_top_genes = 100
        n_comps = 30
    else:
        raise NotImplementedError('data_type must be "sp" or "sc"')

    adata = adata.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    if adata.var.shape[0] > n_top_genes:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
        adata = adata[:, adata.var['highly_variable']]
    print(f"{data_type} Finial gene number: {adata.shape[1]}")

    if data_type == "sc":
        sc.pp.scale(adata, max_value=10)

    sc.pp.pca(adata, n_comps=n_comps)
    sc.pp.neighbors(adata, n_pcs=n_comps, use_rep='X_pca', n_neighbors=gene_neighbor_num)

    nn_graph_genes = adata.obsp["connectivities"]

    if neighbors_method == "pca":
        adjacency_graph = nn_graph_genes
    elif neighbors_method == "spatial":
        sq.gr.spatial_neighbors(adata, n_neighbors=spatial_neighbor_num)
        nn_graph_space = adata.obsp["spatial_connectivities"]
        adjacency_graph = nn_graph_space
    elif neighbors_method == "joint":
        sq.gr.spatial_neighbors(adata)
        nn_graph_space = adata.obsp["spatial_connectivities"]
        joint_graph = (1 - alpha) * nn_graph_genes + alpha * nn_graph_space
        adjacency_graph = joint_graph
        print(f"Setting joint graph with alpha = {alpha}")
    else:
        raise ValueError("neighbors_method must be pca/spatial/joint")

    return adjacency_graph
