import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from tqdm import tqdm

from revise.backend.ops.topology import get_adjacency_graph


def get_sc_obs(spot_names, all_cells_in_spot, spatial_coords):
    """
    Create DataFrame mapping spots to their constituent cells, with spatial coordinates.

    Args:
        spot_names: Array or list of spot names
        all_cells_in_spot: Dictionary mapping spot names to lists of cell IDs
        spatial_coords: numpy array or DataFrame of shape (n_spots, 2), giving x and y for each spot,
                        with order corresponding to spot_names

    Returns:
        pd.DataFrame: DataFrame with columns 'spot_name', 'cell_id', 'x', 'y'
            containing all spot-cell mappings and their spatial coordinates
    """
    sc_obs = pd.DataFrame()
    # handle spatial_coords as np array or pd.DataFrame
    if isinstance(spatial_coords, pd.DataFrame):
        coords_df = spatial_coords
    else:
        coords_df = pd.DataFrame(spatial_coords, columns=['x', 'y'], index=spot_names)

    for idx, i in enumerate(spot_names):
        sc_ids = all_cells_in_spot[i]
        # Get x, y for the spot
        if i in coords_df.index:
            x, y = coords_df.loc[i, 'x'], coords_df.loc[i, 'y']
        else:  # fallback to int index if not indexable by name
            x, y = coords_df.iloc[idx]['x'], coords_df.iloc[idx]['y']
        spot_sc_obs = pd.DataFrame({
            'spot_name': [i] * len(sc_ids),
            'cell_id': sc_ids,
            'x': [x] * len(sc_ids),
            'y': [y] * len(sc_ids),
        })
        sc_obs = pd.concat([sc_obs, spot_sc_obs], axis=0)
    sc_obs.reset_index(drop=True, inplace=True)
    return sc_obs

def get_true_cell_type(SVC_obs, adata_sc):
    """
    Extract true cell type labels and spatial coordinates for SVC cells.
    
    Args:
        SVC_obs: DataFrame with 'cell_id' column
        adata_sc: Single-cell AnnData object containing true labels and coordinates
        
    Returns:
        pd.DataFrame: Updated SVC_obs with added columns:
            - 'true_cell_type': True cell type labels
            - 'x', 'y': Spatial coordinates
    """
    true_cell_type_df = pd.DataFrame(adata_sc.obs)
    true_cell_type_df['cell_id'] = true_cell_type_df['cell_id'].astype(str)
    true_cell_type_df.set_index('cell_id', inplace=True)
    if 'clusters' not in true_cell_type_df.columns:
        true_cell_type_df['clusters'] = true_cell_type_df['Level1']

    SVC_obs['cell_id'] = SVC_obs['cell_id'].astype(str)
    aligned = true_cell_type_df.reindex(SVC_obs['cell_id'])
    SVC_obs['true_cell_type'] = aligned['clusters'].fillna('Unknown').to_numpy()
    known = aligned[['x', 'y']].notna().all(axis=1).to_numpy()
    if np.any(known):
        SVC_obs.loc[known, ["x", "y"]] = aligned.loc[known, ["x", "y"]].to_numpy()

    return SVC_obs

def construct_sc_ref(adata_sc, key_type, type_list=None):
    """
    Construct single-cell reference profiles by averaging expressions per cell type.
    
    Args:
        adata_sc: Single-cell AnnData object
        key_type: Column name in adata_sc.obs containing cell type labels
        type_list: List of cell types to include. If None, uses all unique
            values in adata_sc.obs[key_type]
            
    Returns:
        pd.DataFrame: Reference expression matrix of shape (n_types, n_genes)
            with cell types as index and genes as columns. Each row contains
            the mean expression profile for that cell type.
    """
    n_gene, n_type = adata_sc.shape[1], len(type_list)
    sc_ref = np.zeros((n_type, n_gene))
    X = adata_sc.X if isinstance(adata_sc.X, np.ndarray) else adata_sc.X.toarray()
    for i, cell_type in tqdm(enumerate(type_list)):
        sc_ref[i] = np.mean(X[adata_sc.obs[key_type] == cell_type], axis=0)
    sc_ref = pd.DataFrame(sc_ref, index=type_list, columns=adata_sc.var_names)

    return sc_ref


def get_subcluster(adata_sc, compare_df, sc_resolutions=(1, 2, 3, 4, 5), celltype_col ='Level1'):
    """
    Generate subclustered single-cell data for in-panel vs all-panel comparison.
    
    This function performs Leiden clustering at multiple resolutions for each
    cell type, creating two versions:
    - all_panel: Uses all genes for clustering
    - in_panel: Uses only selected genes (from compare_df) for clustering
    
    Args:
        adata_sc: Single-cell AnnData object
        compare_df: DataFrame with gene names as index (used to select in_panel genes)
        sc_resolutions: Tuple of resolution values for Leiden clustering
        celltype_col: Column name in adata_sc.obs containing cell type labels
            
    Returns:
        tuple: (adata_sc_all_panel, adata_sc_in_panel)
            - adata_sc_all_panel: AnnData with subclusters based on all genes
            - adata_sc_in_panel: AnnData with subclusters based on selected genes
    """
    select_genes = compare_df.index.copy()

    cts = list(adata_sc.obs[celltype_col].unique())
    adata_sc_all_panel = []
    adata_sc_in_panel = []
    for select_ct in tqdm(cts, desc="Cell types"):
        ct_adata_sc = adata_sc[adata_sc.obs[celltype_col] == select_ct]
        ct_adata_sc_all_panel = ct_adata_sc.copy()
        ct_adata_sc_in_panel = ct_adata_sc.copy()
        ct_adata_sc_in_panel = ct_adata_sc_in_panel[:,select_genes]

        adjacency_graph_sc = get_adjacency_graph(ct_adata_sc_all_panel, data_type = "sc", neighbors_method = "pca", alpha = 0)
        for res in tqdm(sc_resolutions, desc = "all panel leiden"):
            sc.tl.leiden(ct_adata_sc_all_panel, adjacency=adjacency_graph_sc, resolution=res, key_added=f"leiden_{res}" )
            ct_adata_sc_all_panel.obs[f'leiden_{res}'] = [f"{select_ct}_{i}" for i in ct_adata_sc_all_panel.obs[f'leiden_{res}']]
        adata_sc_all_panel.append(ct_adata_sc_all_panel)

        adjacency_graph_sp = get_adjacency_graph(ct_adata_sc_in_panel, data_type = "sp", neighbors_method = "pca", alpha = 0)
        for res in tqdm(sc_resolutions, desc = "in panel leiden"):
            sc.tl.leiden(ct_adata_sc_in_panel, adjacency=adjacency_graph_sp, resolution=res, key_added=f"leiden_{res}")
            ct_adata_sc_in_panel.obs[f'leiden_{res}'] = [f"{select_ct}_{i}" for i in ct_adata_sc_in_panel.obs[f'leiden_{res}']]
        adata_sc_in_panel.append(ct_adata_sc_in_panel)

    adata_sc_all_panel = sc.concat(adata_sc_all_panel)
    adata_sc_in_panel = sc.concat(adata_sc_in_panel)

    return adata_sc_all_panel, adata_sc_in_panel


def merge_subcluster(ct_adata_sc, subcluster, mode="mean"):
    """Merge cells inside a subcluster into pseudo-bulk profiles.

    Parameters
    ----------
    ct_adata_sc : AnnData
        Original single-cell data (usually sparse matrix).
    subcluster : str
        Column in ``adata.obs`` used for grouping cells.
    mode : str, default: ``"mean"``
        Aggregation strategy: ``"mean"``, ``"max"``, or ``"median"``.

    Returns
    -------
    AnnData
        New AnnData object containing the merged pseudo-bulk data.
    """
    adata = ct_adata_sc.copy()

    subclusters = adata.obs[subcluster].unique()
    n_genes = adata.n_vars
    n_groups = len(subclusters)

    merged_matrix = np.zeros((n_groups, n_genes))
    group_labels = []

    # Check if input is sparse matrix
    is_sparse_input = issparse(adata.X)

    for i, cluster in enumerate(subclusters):
        # Get cell indices for current subcluster
        cluster_mask = adata.obs[subcluster] == cluster
        if cluster_mask.sum() == 0:
            continue

        cluster_cells = adata[cluster_mask, :]
        group_labels.append(cluster)

        # Extract expression matrix
        X_cluster = cluster_cells.X

        # Calculate expression values based on mode (optimized for sparse matrix processing)
        if mode == "mean":
            if is_sparse_input:
                # For sparse matrix, use mean calculation (avoid converting to dense matrix)
                cluster_expr = np.array(X_cluster.mean(axis=0)).flatten()
            else:
                cluster_expr = X_cluster.mean(axis=0)

        elif mode == "max":
            if is_sparse_input:
                # For sparse matrix, calculate maximum of each column
                cluster_expr = X_cluster.max(axis=0).toarray().flatten()
            else:
                cluster_expr = X_cluster.max(axis=0)

        elif mode == "median":
            if is_sparse_input:
                cluster_expr = np.median(X_cluster.toarray(), axis=0)
            else:
                cluster_expr = np.median(X_cluster, axis=0)

        else:
            raise ValueError(f"Unsupported mode: {mode}. Use 'mean', 'max', or 'median'.")

        # Ensure cluster_expr is a 1D array
        if hasattr(cluster_expr, 'A1'):  # Handle sparse matrix output
            merged_matrix[i, :] = cluster_expr.A1
        else:
            merged_matrix[i, :] = cluster_expr

    # Create new AnnData object
    merged_adata = sc.AnnData(X=merged_matrix)
    merged_adata.var_names = adata.var_names.copy()
    merged_adata.obs[subcluster] = group_labels
    merged_adata.obs_names = [f"{subcluster}_{label}" for label in group_labels]

    # Copy var information from original data
    merged_adata.var = adata.var.copy()

    return merged_adata


def get_prune_adata(adata):
    """
    Prune zero-variance genes by converting them to integer type.
    
    This function identifies genes with zero variance and converts their
    expression values to integers, which can help with downstream processing
    and memory efficiency.
    
    Args:
        adata: AnnData object to prune
        
    Returns:
        AnnData: Copy of adata with zero-variance genes converted to integers
    """
    adata = adata.copy()
    adata_int = adata.copy()
    adata_int.X = adata_int.X.astype(int)
    zero_var_mask = np.asarray(adata_int.X.var(axis=0) == 0).ravel()
    adata.X[:, zero_var_mask] = adata_int.X[:, zero_var_mask]
    return adata
