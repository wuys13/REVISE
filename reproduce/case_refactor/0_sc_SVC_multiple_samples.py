# 实验：Harmony 之后能不能做OT
# Or 单独做完合并的时候，怎么通过Harmony 合并？
select_ct = "T"

import os
output_dir = "results/sc_SVC_case"

patient_id = "P2CRC"
data_type = "Xenium"
output_dir = f"{output_dir}/{patient_id}_{data_type}"
sc_ref_file="adata_sc_all_reanno.h5ad"


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


import scanpy as sc
import pandas as pd

adata_sp = sc.read(raw_file_name)
# adata_sp.obs['total_counts'] = adata_sp.X.sum(axis = 1)
adata_sp = adata_sp[adata_sp.obs['transcript_counts']>=60,:]
sc.pp.filter_genes(adata_sp, min_cells=100)

# adata_sp.obs = adata_sp.obs[['cell_id','x','y']]
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


adata_sp = annotate(adata_sp, adata_sc,
                            ct_name = "Level1", 
                            by_patient = True, patient_key = "Patient" # 这里要实现一下，内部 adata_sp, adata_sc 分别要按照 patient_key 分组，之后再合并
                            )


# select cell type for further analysis
ct_adata_sp = adata_sp[adata_sp.obs['Level1'] == select_ct]
ct_adata_sc = adata_sc[adata_sc.obs['Level1'] == select_ct]


ct_sc_SVC_spatial, ct_sc_SVC_expr = get_sc_SVC_adata(
    ct_adata_sp, ct_adata_sc,  
    by_patient = True, batch_key = "Patient") # 这里要实现一下，内部 adata_sp, adata_sc 分别要按照 patient_key 分组，之后再合并


import scanpy as sc
import harmonypy as hm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def run_harmony(adata, batch_key, n_pcs=20, covariates = None):
    """
    对AnnData对象运行Harmony批次矫正
    :param adata: AnnData对象
    :param batch_key: 批次信息的obs列名
    :param n_pcs: 用于Harmony的主成分数
    :param covariates: 协变量
    :return: 矫正后的AnnData对象（.obsm['X_pca_harmony']），Harmony结果对象
    """
    # 预处理
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    if adata.shape[0] > 2000:
        n_top_genes = 2000
    else:
        n_top_genes = 200
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    adata = adata[:, adata.var['highly_variable']]
    # sc.pp.scale(adata)
    sc.tl.pca(adata, n_comps=n_pcs)
    # Harmony
    ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, batch_key,
    covariates = covariates)
    adata.obsm['X_pca_harmony'] = ho.Z_corr.T
    return adata, ho


ct_sc_SVC_expr = run_harmony(ct_sc_SVC_expr, batch_key = "Patient", covariates = ["SVC_cluster"]) # 将之前每个数据集的SVC_cluster 作为协变量