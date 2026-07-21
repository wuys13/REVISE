#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import scanpy as sc

from revise.analysis.biological_metrics import compute_conditional_moran_i
from revise.analysis.biological_metrics import compute_identity_metrics
from revise.analysis.biological_metrics import compute_local_label_entropy
from revise.analysis.biological_metrics import compute_tmp_mer
from revise.analysis.biological_metrics import make_cell_type_mean_baseline
from revise.analysis.biological_metrics import summarize_conditional_moran_i


def _parse_genes(raw: str | None) -> list[str] | None:
    if raw is None:
        return None
    genes = [token.strip() for token in raw.split(",") if token.strip()]
    return genes or None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compute biology-facing post-reconstruction metrics from an AnnData file: "
            "MISC/MIDC, local label entropy, ASW, and optional TMP/MER."
        )
    )
    parser.add_argument("--input-h5ad", required=True, type=Path, help="Input reconstructed or raw H5AD.")
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory for metric CSV outputs.")
    parser.add_argument("--label-col", default="Level1", help="Cell-type label column in adata.obs.")
    parser.add_argument(
        "--connectivity-key",
        default="spatial_connectivities",
        help="Spatial connectivity matrix key in adata.obsp; built from obsm['spatial'] if absent.",
    )
    parser.add_argument(
        "--genes",
        default=None,
        help="Optional comma-separated gene list for MISC/MIDC. Defaults to all genes.",
    )
    parser.add_argument(
        "--marker-yaml",
        type=Path,
        default=None,
        help="Optional marker YAML for TMP/MER. Supports selected_marker_genes or plain mapping.",
    )
    parser.add_argument(
        "--embedding-key",
        default="X_pca",
        help="Embedding key for ASW. Falls back to adata.X if missing.",
    )
    parser.add_argument(
        "--true-label-col",
        default=None,
        help="Optional reference label column for ARI/NMI identity metrics.",
    )
    parser.add_argument(
        "--pred-label-col",
        default=None,
        help="Optional predicted label column for ARI/NMI identity metrics.",
    )
    parser.add_argument(
        "--no-normalize",
        action="store_true",
        help="Do not normalize/log1p before MISC/MIDC.",
    )
    parser.add_argument(
        "--include-self-entropy",
        action="store_true",
        help="Include the focal cell in each local label entropy neighborhood.",
    )
    parser.add_argument(
        "--cell-type-mean-baseline-h5ad",
        type=Path,
        default=None,
        help="Optional output H5AD for a cell-type-mean negative-control baseline.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.input_h5ad)
    genes = _parse_genes(args.genes)
    normalize = not args.no_normalize

    conditional = compute_conditional_moran_i(
        adata,
        cell_type_col=args.label_col,
        genes=genes,
        connectivity_key=args.connectivity_key,
        normalize=normalize,
        log1p=normalize,
    )
    conditional.to_csv(args.output_dir / "conditional_moran.csv", index=False, float_format="%.8f")
    summarize_conditional_moran_i(conditional).to_csv(
        args.output_dir / "conditional_moran_summary.csv",
        index=False,
        float_format="%.8f",
    )

    entropy, entropy_summary = compute_local_label_entropy(
        adata,
        label_col=args.label_col,
        connectivity_key=args.connectivity_key,
        include_self=args.include_self_entropy,
        normalize=True,
    )
    entropy.to_csv(args.output_dir / "local_label_entropy.csv", index=False, float_format="%.8f")
    entropy_summary.to_csv(
        args.output_dir / "local_label_entropy_summary.csv",
        index=False,
        float_format="%.8f",
    )

    identity = compute_identity_metrics(
        adata,
        label_col=args.label_col,
        embedding_key=args.embedding_key,
        true_label_col=args.true_label_col,
        pred_label_col=args.pred_label_col,
    )
    identity.to_csv(args.output_dir / "identity_metrics.csv", index=False, float_format="%.8f")

    if args.marker_yaml is not None:
        tmp_mer, tmp_mer_macro = compute_tmp_mer(
            adata,
            label_col=args.label_col,
            marker_yaml=args.marker_yaml,
        )
        tmp_mer.to_csv(args.output_dir / "tmp_mer.csv", index=False, float_format="%.8f")
        tmp_mer_macro.to_csv(args.output_dir / "tmp_mer_macro.csv", index=False, float_format="%.8f")

    if args.cell_type_mean_baseline_h5ad is not None:
        baseline = make_cell_type_mean_baseline(adata, label_col=args.label_col)
        args.cell_type_mean_baseline_h5ad.parent.mkdir(parents=True, exist_ok=True)
        baseline.write_h5ad(args.cell_type_mean_baseline_h5ad)


if __name__ == "__main__":
    main()
