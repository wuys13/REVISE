from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import scanpy as sc
from anndata import AnnData


def read_gmt(path: str | Path) -> dict[str, list[str]]:
    """Read gene sets from a GMT file while preserving row and gene order."""
    gene_sets: dict[str, list[str]] = {}
    with Path(path).open("r", encoding="utf-8") as handle:
        for row_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3 or not fields[0] or not any(fields[2:]):
                raise ValueError(f"Malformed GMT row {row_number}")
            name = fields[0]
            if name in gene_sets:
                raise ValueError(f"Duplicate gene set {name!r}")
            gene_sets[name] = fields[2:]
    return gene_sets


def score_genes(
    adata: AnnData,
    genes: Sequence[str],
    *,
    score_name: str,
) -> tuple[AnnData, str]:
    """Run Scanpy gene-set scoring on a copy using genes present in AnnData."""
    overlapping = [gene for gene in genes if gene in adata.var_names]
    if not overlapping:
        raise ValueError("Gene set has no overlap with adata.var_names")

    work = adata.copy()
    sc.tl.score_genes(
        work,
        gene_list=overlapping,
        score_name=score_name,
        use_raw=False,
    )
    return work, score_name
