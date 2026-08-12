import os
import argparse
import numpy as np
import scanpy as sc

from revise.tools.log import Logger
from revise.application import ScSVC
from revise.conf.application_sc_conf import ApplicationScConf


def _parse_resolutions(value):
    chunks = [chunk.strip() for chunk in value.replace(";", ",").split(",") if chunk.strip()]
    if not chunks:
        raise argparse.ArgumentTypeError("resolutions must contain at least one float value")
    try:
        return [float(chunk) for chunk in chunks]
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"invalid resolution list: {value}") from exc


def _build_output_dir(save_path, sample_name, data_type, select_ct):
    return os.path.join(save_path, f"{sample_name}_{data_type}", select_ct)


def main(args):
    output_dir = _build_output_dir(args.save_path, args.sample_name, args.data_type, args.select_ct)
    os.makedirs(output_dir, exist_ok=True)

    config = ApplicationScConf(
        sample_name=args.sample_name,
        raw_data_path=args.raw_data_path,
        result_root_path=output_dir,
        cell_type_col=args.cell_type_col,
        confidence_col="Confidence",
        unknown_key="Unknown",
        st_file=f"{args.data_type}.h5ad",
        sc_ref_file=args.sc_ref_file,
        annotate_mode=args.annotate_mode,
        annotate_pot_reg=args.annotate_pot_reg,
        annotate_pot_reg_m=args.annotate_pot_reg_m,
        annotate_pot_reg_type=args.annotate_pot_reg_type,
        prep_st_min_counts=args.st_min_transcripts,
        prep_st_min_cells=args.st_min_cells,
        prep_sc_min_cells=args.sc_min_cells,
    )

    logger = Logger(name="run.log", log_file=os.path.join(output_dir, "application_sc.log")).get_logger()
    logger.info("Loading spatial data from %s", config.st_file_path)
    adata_sp = sc.read(config.st_file_path)
    if "transcript_counts" not in adata_sp.obs:
        x_mat = adata_sp.X
        transcript_counts = np.asarray(x_mat.sum(axis=1)).ravel()
        adata_sp.obs["transcript_counts"] = transcript_counts
    adata_sp = adata_sp[adata_sp.obs["transcript_counts"] >= args.st_min_transcripts, :]
    if adata_sp.n_obs == 0:
        raise ValueError("No spatial cells remain after transcript filtering")
    sc.pp.filter_genes(adata_sp, min_cells=args.st_min_cells)
    if adata_sp.n_vars == 0:
        raise ValueError("No spatial genes remain after min_cells filtering")

    logger.info(f"Loading sc reference data from {config.sc_ref_file_path}")
    adata_sc = sc.read(config.sc_ref_file_path)
    if "Patient" not in adata_sc.obs:
        raise KeyError("Patient not found in sc reference AnnData obs")
    adata_sc = adata_sc[adata_sc.obs["Patient"] == args.sample_name, :]
    if adata_sc.n_obs == 0:
        raise ValueError(f"No sc reference cells found for sample_name={args.sample_name}")
    required_cols = ["Level1", "Level2"]
    missing = [col for col in required_cols if col not in adata_sc.obs]
    if missing:
        raise KeyError(f"Missing required sc reference obs columns: {missing}")
    adata_sc.obs = adata_sc.obs.loc[:, required_cols].copy()
    sc.pp.filter_genes(adata_sc, min_cells=args.sc_min_cells)
    if adata_sc.n_vars == 0:
        raise ValueError("No sc reference genes remain after min_cells filtering")
    adata_sc.obs["Level1"] = adata_sc.obs["Level1"].replace({"Mono/Macro": "Mono_Macro"})

    logger.info("Aligning genes between spatial and sc reference data")
    overlap_genes = adata_sp.var_names.intersection(adata_sc.var_names)
    if overlap_genes.empty:
        raise ValueError("No overlapping genes between spatial and sc reference data")
    adata_sp = adata_sp[:, overlap_genes]

    sc_svc = ScSVC(adata_sp, adata_sc, config, logger)
    sc_svc.global_anchoring()

    resolutions = (
        [float(args.select_resolution)]
        if args.select_resolution is not None
        else list(args.resolutions)
    )

    logger.info("Begin local cluster: %s", args.select_ct)
    sc_svc_spatial, sc_svc_expr = sc_svc.local_refinement(
        args.select_ct,
        args.sub_cell_type_col,
        resolutions,
        select_res=args.select_resolution,
    )

    if args.save_adata_flag:
        sc_svc_spatial.write(os.path.join(output_dir, "sc_SVC_spatial.h5ad"))
        sc_svc_expr.write(os.path.join(output_dir, "sc_SVC_expr.h5ad"))
        logger.info(f"Done. Output adata saved to {output_dir}")
    else:
        print("No data saved. You can refer to the output directory for the reconstructed AnnData.")


def get_args():
    parser = argparse.ArgumentParser(description="Run sc-SVC reconstruction pipeline.")
    parser.add_argument("--sample_name", help="Patient/sample ID used to subset sc reference data.")
    parser.add_argument("--data_type")
    parser.add_argument("--raw_data_path", help="Directory containing `{sample_name}_{data_type}.h5ad` and sc reference file.")
    parser.add_argument("--save_path", help="Root directory for output.")
    parser.add_argument("--sc_ref_file", help="sc reference h5ad filename.")
    parser.add_argument("--select_ct", help="Cell type (Level1) for local refinement.")
    parser.add_argument("--cell_type_col", default="Level1", help="Cell type column name for global annotation.")
    parser.add_argument("--sub_cell_type_col", default="Level2", help="Sub-cell type column name for local refinement.")
    parser.add_argument("--st_min_transcripts", type=int, default=60, help="Minimum transcript counts per spatial cell.")
    parser.add_argument("--st_min_cells", type=int, default=100, help="Minimum spatial cells per gene.")
    parser.add_argument("--sc_min_cells", type=int, default=100, help="Minimum sc reference cells per gene.")
    parser.add_argument(
        "--annotate_mode",
        default="tacco",
        choices=["tacco", "pot"],
        help="Annotation method for global/local anchoring.",
    )
    parser.add_argument(
        "--annotate_pot_reg",
        type=float,
        default=0.1,
        help="POT regularization parameter for annotation.",
    )
    parser.add_argument(
        "--annotate_pot_reg_m",
        type=float,
        default=0.0,
        help="POT marginal regularization for annotation.",
    )
    parser.add_argument(
        "--annotate_pot_reg_type",
        type=str,
        default="entropy",
        help="POT regularization type for annotation.",
    )
    parser.add_argument(
        "--resolutions",
        type=_parse_resolutions,
        default="0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0",
        help="Leiden resolutions for local clustering, comma-separated (e.g., 0.6,0.7,0.8).",
    )
    parser.add_argument("--select_resolution", default=None, type=float, help="Specific leiden resolutions for local clustering")
    parser.add_argument("--save_adata_flag", action="store_true", help="Save reconstructed AnnData to disk.")

    args_cls = parser.parse_args()
    return args_cls


if __name__ == "__main__":
    args = get_args()
    main(args)
    print(f"Finish sample_name: {args.sample_name} .....")
