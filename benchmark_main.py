import argparse
import os

import scanpy as sc

import revise.benchmark as benchmark
from revise.conf.benchmark_impute_conf import BenchmarkImputeConf
from revise.conf.benchmark_seg_conf import BenchmarkSegConf
from revise.conf.benchmark_sr_conf import BenchmarkSrConf
from revise.tools.log import Logger

from tqdm import tqdm


# base settings based on Sim2Real-ST benchmark: can be modified according to your needs
seg_methods = ["seg_1", "seg_2", "seg_3", "seg_4"]
bin2cell_methods = ["bin2cell"]
batch_nums = [1, 2, 3, 4]
spot_sizes = [20, 50, 100, 200]


def get_args():
    parser = argparse.ArgumentParser("REVISE benchmark")
    parser.add_argument("--cf", type=str, help="confounding factors, segmentation/bin2cell/batch_effect/spot_size/gene_panel/gene_dropout")
    parser.add_argument("--raw_data_path", type=str, help="data root path")
    parser.add_argument("--sample_name", type=str, help="sample name of datasets")
    parser.add_argument("--task", type=str, help="task of benchmark")

    parser.add_argument("--st_file", type=str, help="st file name")
    parser.add_argument("--gt_svc_file", type=str, help="ground truth svc file name: raw Xenium")
    parser.add_argument("--sc_ref_file", type=str, help="single cell reference file name")
    parser.add_argument("--save_path", default="results", type=str, help="save path")
    parser.add_argument("--log_path", default="logs", type=str, help="save path")
    args_cls = parser.parse_args()
    return args_cls


def run_benchmark(cf, config, logger):

    logger.info(f"read st_file: {config.st_file_path}, gt_svc_file: {config.gt_svc_file_path}, sc_ref_file: {config.sc_ref_file_path}")

    adata_st = sc.read_h5ad(config.st_file_path)
    adata_real_st = sc.read_h5ad(config.gt_svc_file_path)
    adata_sc_ref = sc.read_h5ad(config.sc_ref_file_path)

    if cf == "segmentation" or cf == "bin2cell": # sp-SVC
        svc = benchmark.SpSVC(adata_st, adata_sc_ref, config, adata_real_st, logger)
    elif cf == "batch_effect" or cf == "spot_size": # sc-SVC
        svc = benchmark.ScSVCSr(adata_st, adata_sc_ref, config, adata_real_st, logger)
    elif cf == "gene_panel" or cf == "gene_dropout": # sc-SVC
        svc = benchmark.ScSVCImpute(adata_st, adata_sc_ref, config, adata_real_st, logger)
    else:
        raise ValueError(f"{cf} is not supported")

    logger.info(f"Starting {cf} benchmark ------------------------------------------")
    benchmark.main(svc)
    logger.info(f"Finished {cf} benchmark ------------------------------------------")


def get_logger(args_cls, log_prefix):
    log_dir = os.path.join(args_cls.log_path, args_cls.task, args_cls.sample_name)
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(str(log_dir), f'{log_prefix}.log')
    logger = Logger(name=args_cls.task, log_file=log_file).get_logger()
    return logger


def main(args_cls):
    if args_cls.cf == "segmentation":
        if args_cls.st_file is None or args_cls.gt_svc_file is None or args_cls.sc_ref_file is None:
            raise ValueError(f"cf '{args_cls.cf}' must specify st_file and gt_svc_file or sc_ref_file")
        for seg_method in tqdm(seg_methods, desc="Segmentation Methods"):
            config = BenchmarkSegConf(
                sample_name=args_cls.sample_name,
                annotate_mode="pot",
                raw_data_path=str(os.path.join(args_cls.raw_data_path, args_cls.task)),
                result_root_path=str(os.path.join(args_cls.save_path, args_cls.task)),
                cell_type_col="Level1",
                confidence_col="Confidence",
                unknown_key="Unknown",
                st_file=args_cls.st_file,
                gt_svc_file=args_cls.gt_svc_file,
                sc_ref_file=args_cls.sc_ref_file,
                seg_method=seg_method,
            )
            logger = get_logger(args_cls, seg_method)
            run_benchmark(args_cls.cf, config, logger)

    elif args_cls.cf == "bin2cell":
        if args_cls.st_file is None or args_cls.gt_svc_file is None or args_cls.sc_ref_file is None:
            raise ValueError(f"cf '{args_cls.cf}' must specify st_file and gt_svc_file or sc_ref_file")
        for seg_method in tqdm(bin2cell_methods, desc="Bin2Cell Methods"):
            config = BenchmarkSegConf(
                sample_name=args_cls.sample_name,
                annotate_mode="pot",
                raw_data_path=str(os.path.join(args_cls.raw_data_path, args_cls.task)),
                result_root_path=str(os.path.join(args_cls.save_path, args_cls.task)),
                cell_type_col="Level1",
                confidence_col="Confidence",
                unknown_key="Unknown",
                st_file=args_cls.st_file,
                gt_svc_file=args_cls.gt_svc_file,
                sc_ref_file=args_cls.sc_ref_file,
                seg_method=seg_method,
            )
            logger = get_logger(args_cls, seg_method)
            run_benchmark(args_cls.cf, config, logger)
            
    elif args_cls.cf == "batch_effect":
        spot_size = 50 # can be set to any parameter in global spot_size_list
        for batch_num in tqdm(batch_nums, desc="Batch Effect Settings"):
            logger = get_logger(args_cls, f"{spot_size}_{batch_num}")
            if batch_num == 0: # simulation
                logger.info("simulation")
                gt_svc_file = "simulated_xenium.h5ad"
                sc_ref_file = "simulated_xenium.h5ad"
            elif batch_num == 1: # real data, same ref
                logger.info(f"batch num {batch_num}, real data, same ref")
                gt_svc_file = "selected_xenium.h5ad"
                sc_ref_file = "selected_xenium.h5ad"
            elif batch_num == 2: # real data, pair ref
                logger.info(f"batch num {batch_num}, real data, pair patient with part cell types")
                gt_svc_file = "selected_xenium.h5ad"
                sc_ref_file = "real_sc_ref_part.h5ad"
            elif batch_num == 3: # real data, other patient ref
                logger.info(f"batch num {batch_num}, real data, pair ref with all sc")
                gt_svc_file = "selected_xenium.h5ad"
                sc_ref_file = "real_sc_ref_all.h5ad"
            elif batch_num == 4: # real data, other patient ref
                logger.info(f"batch num {batch_num}, real data, other patient ref")
                gt_svc_file = "selected_xenium.h5ad"
                sc_ref_file = "real_sc_ref_others.h5ad"
            else:
                raise NotImplementedError(f"batch_num {batch_num} not implemented")
            config = BenchmarkSrConf(
                sample_name=args_cls.sample_name,
                annotate_mode="pot",
                raw_data_path=str(os.path.join(args_cls.raw_data_path, args_cls.task)),
                result_root_path=str(os.path.join(args_cls.save_path, args_cls.task)),
                cell_type_col="Level1",
                confidence_col="Confidence",
                unknown_key="Unknown",
                st_file="xenium_spot.h5ad",
                gt_svc_file=gt_svc_file,
                sc_ref_file=sc_ref_file,
                spot_size=spot_size
            )
            run_benchmark(args_cls.cf, config, logger)

    elif args_cls.cf == "spot_size":
        batch_num = 3 # the most realistic setting
        gt_svc_file = "selected_xenium.h5ad"
        sc_ref_file = "real_sc_ref_all.h5ad"

        for spot_size in tqdm(spot_sizes, desc="Spot Sizes"):
            config = BenchmarkSrConf(
                sample_name=args_cls.sample_name,
                annotate_mode="pot",
                raw_data_path=str(os.path.join(args_cls.raw_data_path, args_cls.task)),
                result_root_path=str(os.path.join(args_cls.save_path, args_cls.task)),
                cell_type_col="Level1",
                confidence_col="Confidence",
                unknown_key="Unknown",
                st_file="xenium_spot.h5ad",
                gt_svc_file=gt_svc_file,
                sc_ref_file=sc_ref_file,
                spot_size=spot_size
            )
            # log_dir = os.path.join(args_cls.log_path, args_cls.task, args_cls.sample_name)
            # os.makedirs(log_dir, exist_ok=True)
            # log_file = os.path.join(log_dir, f'{spot_size}_{batch_num}.log')
            logger = get_logger(args_cls, f"{spot_size}_{batch_num}")
            run_benchmark(args_cls.cf, config, logger)

    elif args_cls.cf == "gene_panel":

        st_file = "selected_xenium.h5ad"
        gt_svc_file = "selected_xenium.h5ad"
        sc_ref_file = "real_sc_ref.h5ad"

        config = BenchmarkImputeConf(
            sample_name=args_cls.sample_name,
            annotate_mode="pot",
            raw_data_path=str(os.path.join(args_cls.raw_data_path, args_cls.task)),
            result_root_path=str(os.path.join(args_cls.save_path, args_cls.task)),
            cell_type_col="Level1",
            confidence_col="Confidence",
            unknown_key="Unknown",
            st_file=st_file,
            gt_svc_file=gt_svc_file,
            sc_ref_file=sc_ref_file,
        )
        # log_dir = os.path.join(args_cls.log_path, args_cls.task, args_cls.sample_name)
        # os.makedirs(log_dir, exist_ok=True)
        # log_file = os.path.join(log_dir, f'{args_cls.cf}.log')
        logger = get_logger(args_cls, f"{args_cls.cf}")
        run_benchmark(args_cls.cf, config, logger)
    
    elif args_cls.cf == "gene_dropout":

        st_file = "selected_xenium.h5ad"
        gt_svc_file = "selected_xenium.h5ad"
        sc_ref_file = "real_sc_ref.h5ad"

        config = BenchmarkImputeConf(
            sample_name=args_cls.sample_name,
            annotate_mode="pot",
            raw_data_path=str(os.path.join(args_cls.raw_data_path, args_cls.task)),
            result_root_path=str(os.path.join(args_cls.save_path, args_cls.task)),
            cell_type_col="Level1",
            confidence_col="Confidence",
            unknown_key="Unknown",
            st_file=st_file,
            gt_svc_file=gt_svc_file,
            sc_ref_file=sc_ref_file,
        )
        # log_dir = os.path.join(args_cls.log_path, args_cls.task, args_cls.sample_name)
        # os.makedirs(log_dir, exist_ok=True)
        # log_file = os.path.join(log_dir, f'{args_cls.cf}.log')
        logger = get_logger(args_cls, f"{args_cls.cf}")
        run_benchmark(args_cls.cf, config, logger)

    else:
        raise NotImplementedError(f"config confounding factor {args_cls.cf} not implemented ...")


if __name__ == "__main__":
    args = get_args()
    print(args)
    main(args)
    print("Done: ", args.cf)
