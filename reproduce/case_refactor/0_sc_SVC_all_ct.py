from tqdm import tqdm
import sys
sys.path.append("src")

# input data -------------------
import os
output_dir = "results"

patient_id = "P2CRC"
data_type = "Xenium"
output_dir = f"{output_dir}/{patient_id}_{data_type}"
os.makedirs(output_dir, exist_ok=True)

sc_ref_file="adata_sc_all_reanno.h5ad"




alpha = 0.2

raw_data_path = "./raw_data"

raw_file_name = f"{raw_data_path}/{patient_id}_{data_type}.h5ad"
sc_file_name = f"{raw_data_path}/adata_sc_all_reanno.h5ad"




# read and process data -------------------
import scanpy as sc
import pandas as pd

adata_sp = sc.read(raw_file_name)
# adata_sp.obs['total_counts'] = adata_sp.X.sum(axis = 1)
adata_sp = adata_sp[adata_sp.obs['transcript_counts']>=60,:]
sc.pp.filter_genes(adata_sp, min_cells=100)

adata_sp.obs = adata_sp.obs[['cell_id','x','y']]
# adata_sp.obs_names = adata_sp.obs['cell_id']

adata_sc = sc.read(sc_file_name)
adata_sc = adata_sc[adata_sc.obs['Patient'] == patient_id, :]
adata_sc.obs = adata_sc.obs[['Level1','Level2']]
sc.pp.filter_genes(adata_sc, min_cells=100)
adata_sc.obs['Level1'].replace({"Mono/Macro": "Mono_Macro"}, inplace=True)
adata_sc_raw = adata_sc.copy()

overlap_genes = adata_sp.var_names.intersection(adata_sc.var_names)
adata_sp = adata_sp[:, overlap_genes]

adata_sp



# revise cluster mode -------------------
from revise.sc_anno import annotate
def get_sc_SVC_adata(ct_adata_sp, ct_adata_sc):

    ct_adata_sp = annotate(ct_adata_sp, ct_adata_sc,
                                ct_name = "Level2")


    from revise.sc_SVC_cluster_mode import get_best_cluster, split_adata_evenly
    resolutions = [0.6, 0.7, 0.8]
    # resolutions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    neighbors_method = "joint"
    adata_samples = split_adata_evenly(ct_adata_sp)

    sc_SVC_adata, merge_df, best_res = get_best_cluster(adata_samples[0], neighbors_method = neighbors_method, alpha = alpha, resolutions = resolutions)

    sc_SVC_adata


    sc_SVC_adata.obs['SVC_cluster'] = sc_SVC_adata.obs[f'leiden_{best_res}']
    sp_cluster_num = merge_df.loc[merge_df['resolution'] == best_res, 'cluster_num'].values[0]
    sp_cluster_num

    ct_adata_sc = annotate(ct_adata_sc, sc_SVC_adata,
                                ct_name = 'SVC_cluster',
                                )
    ct_adata_sc


    # from revise.sc_SVC_cluster_mode import get_cm, plot_cm
    # cm_df = get_cm(ct_adata_sc, 'SVC_cluster', 'Level2')
    # plot_cm(cm_df, save_dir = output_dir)
    # cm_df

    return sc_SVC_adata, ct_adata_sc



adata_sp = annotate(adata_sp, adata_sc,
                            ct_name = "Level1"
                            )
cts = adata_sc.obs['Level1'].unique()

for select_ct in tqdm(cts, desc = "Cell types"):
    print(select_ct)
    ct_adata_sc = adata_sc[adata_sc.obs['Level1'] == select_ct]
    ct_adata_sp = adata_sp[adata_sp.obs['Level1'] == select_ct]

    ct_sc_SVC_spatial, ct_sc_SVC_expr = get_sc_SVC_adata(ct_adata_sp, ct_adata_sc)
    ct_sc_SVC_spatial.write(f"{output_dir}/{select_ct}_spatial.h5ad")
    ct_sc_SVC_expr.write(f"{output_dir}/{select_ct}_expr.h5ad")