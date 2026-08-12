import os
import argparse
import scanpy as sc

import revise.application as application
from revise.tools.log import Logger

from revise.conf.application_sp_conf import ApplicationSpConf

def get_args():
    parser = argparse.ArgumentParser("REVISE application sp-SVC")
    parser.add_argument("--raw_data_path", type=str, help="data root path")
    parser.add_argument("--sample_name", type=str, help="sample name of datasets")
    parser.add_argument("--st_file", type=str, help="st file name")
    parser.add_argument("--sc_ref_file", type=str, help="single cell reference file name")
    parser.add_argument("--patient_key", default="Patient", type=str, help="patient key in sc_ref_file")
    parser.add_argument("--sample_size", default=None, type=int, help="subsample for test")

    parser.add_argument("--save_path", default="results/sp_SVC_case", type=str, help="save path")
    parser.add_argument("--log_path", default="logs", type=str, help="save path")

    args_cls = parser.parse_args()
    return args_cls


def preprocess(adata, data_type):
    adata = adata.copy()
    if data_type == "st":
        sc.pp.filter_cells(adata, min_counts=20)
        sc.pp.filter_genes(adata, min_cells=30)
    elif data_type == "sc":
        sc.pp.filter_cells(adata, min_genes=20)
        sc.pp.filter_genes(adata, min_cells=50)

        replace_columns = {k: k.replace("/", "_") for k in adata.obs["Level1"].unique().tolist() if "/" in k}
        adata.obs["Level1"].replace(replace_columns, inplace=True)
        replace_columns = {k: k.replace("/", "_") for k in adata.obs["Level2"].unique().tolist() if "/" in k}
        adata.obs["Level2"].replace(replace_columns, inplace=True)
    return adata


def main(args_cls):
    config = ApplicationSpConf(
        sample_name=args_cls.sample_name,
        annotate_mode="pot",
        raw_data_path=args_cls.raw_data_path,
        result_root_path=args_cls.save_path,
        cell_type_col="Level1",
        confidence_col="Confidence",
        unknown_key="Unknown",
        st_file=args_cls.st_file,
        sc_ref_file=args_cls.sc_ref_file
    )
    os.makedirs(config.result_dir, exist_ok=True)
    logger = Logger(name="application_sp", log_file=os.path.join(config.result_dir, "application_sp.log")).get_logger()
    logger.info(f"begin run REVISE application sp")
    logger.info(f"read st_file: {config.st_file_path}, sc_ref_file: {config.sc_ref_file_path}")
    adata_st = sc.read_h5ad(config.st_file_path)
    
    # prevent OOM or accelerate the running process
    if args_cls.sample_size is not None:
        sc.pp.subsample(adata_st, n_obs=args_cls.sample_size)
    adata_st = preprocess(adata_st, "st")

    adata_sc_ref = sc.read_h5ad(config.sc_ref_file_path)
    if args_cls.patient_key is not None:
        adata_sc_ref = adata_sc_ref[adata_sc_ref.obs[args_cls.patient_key] == args_cls.sample_name, :]
    adata_sc_ref = preprocess(adata_sc_ref, "sc")
    
    logger.info(f"adata_st shape: {adata_st.shape}, adata_sc_ref shape: {adata_sc_ref.shape}")
    
    sp_svc = application.SpSVC(adata_st, adata_sc_ref, config, logger)
    sp_svc.global_anchoring()
    sp_svc.local_refinement()
    for key, svc_adata in sp_svc.svc.items():
        # svc_adata.write_h5ad(os.path.join(config.result_dir, f"{key}.h5ad"))
        svc_adata.write_h5ad(os.path.join(config.result_dir, "sp_SVC.h5ad"))


if __name__ == "__main__":
    args = get_args()
    main(args)
    print(f"Finish patient_id: {args.sample_name} .....")
