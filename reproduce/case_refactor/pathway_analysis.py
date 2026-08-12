import scanpy as sc
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np

def read_gmt(file_path):
    """
    读取 GMT 文件并返回一个字典，键为通路名称，值为基因列表。
    """
    pathway_dict = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            pathway_name = parts[0]
            genes = parts[2:]  # 跳过前两列（通路名称和描述）
            pathway_dict[pathway_name] = genes
    return pathway_dict


def plot_sp(adata_sp, color, size = 40, save_path = None):
    
    plt.figure(figsize=(12, 10))
    sc.pl.scatter(adata_sp, x="x", y="y", color=color, size = size)
    if save_path:
        plt.savefig(save_path)

def get_sampled_adata(adata_sp, n_samples = 30000, seed = 42):
    """
    从 adata_sp 中随机抽取 n_samples 个样本，并返回一个新的 AnnData 对象。
    """
    if n_samples is not None:
        np.random.seed(seed)
        sampled_adata = adata_sp[np.random.choice(adata_sp.n_obs, n_samples, replace=False), :].copy()
    return sampled_adata


def get_select_adata(file_name, plot_data_type, patient_id, data_type, n_samples, seed):
    adata_sp = sc.read_h5ad(file_name)
    
    if plot_data_type == "original":
        print(f"Processing original data for {patient_id}_{data_type}")
    else:
        adata_sp.X = adata_sp.layers["ot_smooth"]
        print(f"Processing sp_SVC data for {patient_id}_{data_type}")
    
    print(adata_sp)

    # print(adata_sp.obs['Level1'].value_counts())
    
    select_adata = get_sampled_adata(adata_sp, n_samples = n_samples, seed = seed)
    # sc.pp.normalize_total(select_adata, target_sum=1e4)
    # sc.pp.log1p(select_adata)

    return select_adata

import omicverse as ov
import os
from tqdm import tqdm
def plot_pathway(select_adata, pathway_dict, score_method, save_dir, plot_data_type, save_file_format = "png"):
    if score_method == "AUC":
        save_path = f"{save_dir}/{plot_data_type}/pathway_{score_method}"
        os.makedirs(save_path, exist_ok=True)

        for pathway_name, genes in tqdm(pathway_dict.items(), desc = "AUC Plot"):
            try:
                valid_genes = [gene for gene in genes if gene in select_adata.var_names]
                print(len(valid_genes), len(genes))
                
                # sc.pp.normalize_total(select_adata, target_sum=1e4)  # 标准化
                # sc.pp.log1p(select_adata)  # 对数转换

                ov.single.geneset_aucell(
                    select_adata,
                    geneset_name=pathway_name,
                    geneset=valid_genes
                )
                color = f"{pathway_name}_aucell"
                file_path = f"{save_path}/{pathway_name}.{save_file_format}"
                plot_sp(select_adata, color, save_path = file_path,
                        #  size = 60
                        )
            except Exception as e:
                print(f"Error processing {pathway_name}: {e}")

    elif score_method == "score_genes":
        sc.pp.scale(select_adata)
        save_path = f"{save_dir}/{plot_data_type}/pathway_{score_method}"
        os.makedirs(save_path, exist_ok=True)

        for pathway_name, genes in tqdm(pathway_dict.items(), desc = "Score Plot"):
            valid_genes = [gene for gene in genes if gene in select_adata.var_names]
            print(len(valid_genes), len(genes))
            sc.tl.score_genes(
                select_adata,
                gene_list=valid_genes,
                score_name=pathway_name,
                use_raw=False,  
            )
            color = f"{pathway_name}"
            file_path = f"{save_path}/{pathway_name}.{save_file_format}"
            plot_sp(select_adata, color, save_path = file_path)
    else:
        print("Unknown score method")

def get_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw_data_path", type=str, default="/home/wys/Sim2Real-ST/REVISE_data_process/raw_data")
    parser.add_argument("--svc_data_path", type=str, default="../REVISE/results/HD")
    parser.add_argument("--patient_id", type=str, default="P1CRC")
    parser.add_argument("--data_type", type=str, default="HD") # HD_sp_SVC
    parser.add_argument("--score_method", type=str, default="AUC")
    parser.add_argument("--n_samples", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--save_file_format", type=str, default="png")
    parser.add_argument("--gmt_selection", type=str, default="select_5")

    parser.add_argument("--save_dir", type=str, default="./output/sp_SVC_case")
    return parser.parse_args()

def main():
    args = get_args()
    raw_data_path = args.raw_data_path
    svc_data_path = args.svc_data_path
    patient_id = args.patient_id
    data_type = args.data_type
    score_method = args.score_method
    n_samples = args.n_samples
    seed = args.seed
    save_file_format = args.save_file_format
    save_dir = args.save_dir

    if args.gmt_selection == "select_5":
        gmt_file_path = './pathway/select_5.gmt'
    elif args.gmt_selection == "all":
        gmt_file_path = './pathway/h.all.v2025.1.Hs.symbols.gmt'
    pathway_dict = read_gmt(gmt_file_path)
    print(len(pathway_dict))


    save_dir = f"{save_dir}/{patient_id}_{data_type}"
    raw_file_name = f"{raw_data_path}/{patient_id}_{data_type}.h5ad"
    sp_SVC_file_name = f"{svc_data_path}/{patient_id}/{patient_id}_{data_type}_pot_REVISE.h5ad"
    
    for plot_data_type in tqdm(["original", "sp_SVC"], desc = "plot_data_type"):
        if plot_data_type == "sp_SVC":
            file_name = sp_SVC_file_name
        elif plot_data_type == "original":
            file_name = raw_file_name
        else:
            print("Unknown plot data type")
            continue
        
        select_adata = get_select_adata(file_name, plot_data_type, patient_id, data_type, n_samples, seed)
        plot_pathway(select_adata, pathway_dict, score_method, save_dir, plot_data_type, save_file_format)

    print("Finish")

if __name__ == "__main__":
    main()