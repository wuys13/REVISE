from __future__ import annotations

from typing import Iterable

import pandas as pd
import scanpy as sc
from anndata import AnnData

from revise.utils.labels import normalize_cell_type_label

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


def normalize_reference_labels(
    adata: AnnData,
    columns: Iterable[str | None],
    *,
    trim: bool = False,
) -> AnnData:
    """Normalize route-visible reference labels without changing gene axes."""
    result = adata.copy()
    for column in columns:
        if column is None or column not in result.obs:
            continue
        labels = result.obs[column].astype(str)
        if trim:
            normalized = labels.map(normalize_cell_type_label)
        else:
            if not labels.str.contains("/", regex=False).any():
                continue
            normalized = labels.str.replace("/", "_", regex=False)
        pairs = pd.DataFrame({"original": labels, "normalized": normalized}).drop_duplicates()
        collisions = pairs.groupby("normalized", sort=False)["original"].nunique()
        if (collisions > 1).any():
            names = collisions[collisions > 1].index.tolist()
            raise ValueError(
                f"Reference labels in {column!r} collide after slash normalization: "
                f"{names[:5]}"
            )
        result.obs[column] = normalized
    return result


def prepare_sc_svc_pair(
    spatial: AnnData,
    reference: AnnData,
    *,
    broad_column: str,
    subtype_column: str | None = None,
) -> tuple[AnnData, AnnData]:
    """Prepare the standard sc-SVC pair after generic preprocessing."""
    required_columns = list(dict.fromkeys([broad_column, subtype_column]))
    required_columns = [column for column in required_columns if column is not None]
    missing = [column for column in required_columns if column not in reference.obs]
    if missing:
        raise KeyError(f"Missing required columns in sc reference: {missing}")
    result_reference = reference.copy()
    result_reference.obs = result_reference.obs.loc[:, required_columns].copy()
    result_reference = normalize_reference_labels(
        result_reference,
        required_columns,
        trim=True,
    )
    overlap_genes = spatial.var_names.intersection(result_reference.var_names)
    if overlap_genes.empty:
        raise ValueError("No overlapping genes between spatial and sc reference data")
    return spatial[:, overlap_genes].copy(), result_reference


__all__ = [
    "filter_reference",
    "normalize_reference_labels",
    "prepare_sc_svc_pair",
    "preprocess_reference",
    "preprocess_spatial",
]
