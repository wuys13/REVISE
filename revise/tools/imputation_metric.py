import squidpy as sq
import pandas as pd
import numpy as np
from tqdm import tqdm
from scipy.sparse import diags

def compute_conditional_MoranI(adata_sp, cell_type_col="Level1", save_dir=None):
    '''
    Diffuse expression leads to counterintuitively elevated spatial autocorrelation, requiring cell-type-specific analysis.
    '''

    adata_sp = adata_sp.copy()

    sc.pp.normalize_total(adata_sp, target_sum=1e4)
    sc.pp.log1p(adata_sp)

    all_moranI_df = pd.DataFrame()
    
    # Compute spatial neighbors
    sq.gr.spatial_neighbors(adata_sp)
    
    # Get the connectivity matrix
    W = adata_sp.obsp['spatial_connectivities'].tocoo()
    
    # Get cell types
    ctype = adata_sp.obs[cell_type_col].values
    
    # Create masks
    same_mask = (ctype[W.row] == ctype[W.col])
    diff_mask = (ctype[W.row] != ctype[W.col])
    
    # Create copies of the connectivity matrix for same and different cell types
    W_same = W.copy()
    W_diff = W.copy()
    
    # Apply masks to the connectivity matrices
    W_same.data[~same_mask] = 0
    W_diff.data[~diff_mask] = 0
    
    # Convert back to CSR format and add to adata
    adata_sp.obsp['spatial_connectivities_same'] = W_same.tocsr()
    adata_sp.obsp['spatial_connectivities_diff'] = W_diff.tocsr()
    
    # Compute Moran's I for same cell type
    sq.gr.spatial_autocorr(
        adata_sp,
        connectivity_key="spatial_connectivities_same",
        mode="moran",
        genes=adata_sp.var_names,
        n_perms=1,
        n_jobs=1,
    )
    df_same = adata_sp.uns["moranI"]["I"].loc[adata_sp.var_names]
    df_same[df_same >= 1] = np.nan
    df_same[df_same <= -1] = np.nan
    df_same = df_same.rename("I_same")
    
    # Compute Moran's I for different cell types
    sq.gr.spatial_autocorr(
        adata_sp,
        connectivity_key="spatial_connectivities_diff",
        mode="moran",
        genes=adata_sp.var_names,
        n_perms=1,
        n_jobs=1,
    )
    df_diff = adata_sp.uns["moranI"]["I"].loc[adata_sp.var_names]
    df_diff[df_diff >= 1] = np.nan
    df_diff[df_diff <= -1] = np.nan
    df_diff = df_diff.rename("I_diff")
    
    # Combine results
    all_moranI_df = pd.concat([df_same, df_diff], axis=1)
    
    return all_moranI_df
