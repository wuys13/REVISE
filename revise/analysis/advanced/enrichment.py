from __future__ import annotations

import pandas as pd


ENRICHMENT_RESULT_COLUMNS = (
    "Gene_set",
    "Term",
    "Overlap",
    "P-value",
    "Adjusted P-value",
    "Old P-value",
    "Old Adjusted P-value",
    "Odds Ratio",
    "Combined Score",
    "Genes",
)


class EnrichmentExecutionError(RuntimeError):
    """Raised when the enrichment provider fails during execution."""


def _empty_result() -> pd.DataFrame:
    return pd.DataFrame(columns=ENRICHMENT_RESULT_COLUMNS)


def _require_gseapy():
    try:
        import gseapy
    except ImportError as exc:
        raise ImportError(
            "gseapy is required for pathway enrichment features; "
            'install it with `python -m pip install "revise-svc[pathway]"`.'
        ) from exc
    return gseapy


def get_enrichment(deg_genes, geneset_file, cutoff=0.05):
    """Run Enrichr while preserving dependency and execution failure states."""
    if not deg_genes:
        return _empty_result()

    gp = _require_gseapy()
    try:
        result = gp.enrichr(
            gene_list=deg_genes,
            gene_sets=geneset_file,
            organism="human",
            cutoff=cutoff,
        ).results
    except Exception as exc:
        raise EnrichmentExecutionError("Enrichr execution failed") from exc

    if result.empty:
        return _empty_result()
    return result


def get_enrichment_local(
    deg_genes, geneset_file, cutoff=0.05, background=None
):
    """Run GSEApy enrichment against a local gene-set resource."""
    if not deg_genes:
        return _empty_result()

    gp = _require_gseapy()
    try:
        result = gp.enrich(
            gene_list=deg_genes,
            gene_sets=geneset_file,
            background=background,
            cutoff=cutoff,
            no_plot=True,
        ).results
    except Exception as exc:
        raise EnrichmentExecutionError("GSEApy local enrichment execution failed") from exc

    if result.empty:
        return _empty_result()
    return result


__all__ = [
    "ENRICHMENT_RESULT_COLUMNS",
    "EnrichmentExecutionError",
    "get_enrichment",
    "get_enrichment_local",
]
