from __future__ import annotations

import scanpy as sc
from anndata import AnnData


def filter_reference(
    adata: AnnData,
    filter_column: str | None = None,
    filter_value: str | None = None,
) -> AnnData:
    if filter_column is None and filter_value is None:
        return adata.copy()
    filtered = adata[adata.obs[filter_column] == filter_value, :].copy()
    if filtered.n_obs == 0:
        raise ValueError(
            f"Reference filter {filter_column!r} == {filter_value!r} matched no rows"
        )
    return filtered


def preprocess_spatial(
    adata: AnnData,
    min_transcript_counts: int | None = 60,
    min_cell_counts: int = 100,
    min_counts: int | None = None,
) -> AnnData:
    result = adata.copy()
    if min_counts is not None:
        sc.pp.filter_cells(result, min_counts=min_counts)
    if min_transcript_counts is not None:
        result = result[
            result.obs["transcript_counts"] >= min_transcript_counts,
            :,
        ].copy()
    sc.pp.filter_genes(result, min_cells=min_cell_counts)
    return result


def preprocess_reference(
    adata: AnnData,
    min_transcript_counts: int | None = None,
    min_cell_counts: int = 100,
    min_genes: int | None = None,
) -> AnnData:
    result = adata.copy()
    if min_genes is not None:
        sc.pp.filter_cells(result, min_genes=min_genes)
    if min_transcript_counts is not None:
        result = result[
            result.obs["transcript_counts"] >= min_transcript_counts,
            :,
        ].copy()
    sc.pp.filter_genes(result, min_cells=min_cell_counts)
    return result


__all__ = ["filter_reference", "preprocess_reference", "preprocess_spatial"]
