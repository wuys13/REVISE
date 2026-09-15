"""Descriptive Moran comparison on an explicitly supplied shared spatial graph."""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse


def _side_moran(adata: AnnData, weights: sparse.csr_matrix, support_reason: str) -> dict:
    results = {}
    weight_sum = float(weights.sum())
    if adata.X is None:
        raise ValueError("Expression matrix X is missing")
    # Bound dense working memory by a gene block rather than the full gene axis.
    for start in range(0, adata.n_vars, 64):
        block = adata.X[:, start:start + 64]
        values = np.asarray(block.toarray() if sparse.issparse(block) else block, dtype=float)
        if not np.isfinite(values).all():
            raise ValueError("Expression values must be finite")
        centered = values - values.mean(axis=0)
        denominator = np.einsum("ij,ij->j", centered, centered)
        numerator = np.einsum("ij,ij->j", centered, weights @ centered)
        for column, gene in enumerate(adata.var_names[start:start + 64]):
            if support_reason:
                results[gene] = (np.nan, "insufficient_support", support_reason)
            elif denominator[column] == 0:
                results[gene] = (np.nan, "not_computable", "constant_expression")
            else:
                value = adata.n_obs * numerator[column] / (weight_sum * denominator[column])
                if not np.isfinite(value):
                    raise ValueError("Moran computation produced a nonfinite value")
                results[gene] = (float(value), "computed", "")
    return results


def compare_moran(
    raw: AnnData,
    reconstruction: AnnData,
    weights: sparse.spmatrix,
    *,
    min_units: int,
    min_edges: int,
) -> pd.DataFrame:
    """Compare all provided genes without changing expression or final weights.

    The caller supplies prepared expression and the final shared graph in exact
    observation order. Edges count nonzero directed weights; isolates remain in
    the declared cohort. No normalization, graph construction or inference is
    performed here.
    """
    if type(min_units) is not int or min_units < 1 or type(min_edges) is not int or min_edges < 1:
        raise ValueError("min_units and min_edges must be positive integers")
    if (not raw.obs_names.is_unique or not reconstruction.obs_names.is_unique
            or not raw.obs_names.equals(reconstruction.obs_names) or raw.n_obs == 0):
        raise ValueError("Raw and reconstruction observation axes must be unique, nonempty and identical")
    for adata in (raw, reconstruction):
        if not adata.var_names.is_unique or adata.n_vars == 0:
            raise ValueError("Gene axes must be unique and nonempty")
    graph = sparse.csr_matrix(weights, dtype=float, copy=True)
    graph.sum_duplicates()
    graph.eliminate_zeros()
    if graph.shape != (raw.n_obs, raw.n_obs):
        raise ValueError("Graph shape must match the observation axis")
    if not np.isfinite(graph.data).all() or (graph.data < 0).any() or graph.diagonal().any():
        raise ValueError("Graph must have finite nonnegative weights and no self edges")
    reason = "no_edges" if graph.nnz == 0 else (
        "too_few_units" if raw.n_obs < min_units else
        "too_few_edges" if graph.nnz < min_edges else ""
    )
    sides = [_side_moran(adata, graph, reason) for adata in (raw, reconstruction)]
    rows = []
    for gene in dict.fromkeys([*raw.var_names, *reconstruction.var_names]):
        row = {"gene_id": gene, "n_units": raw.n_obs, "n_edges": graph.nnz}
        for name, values in zip(("raw", "reconstruction"), sides):
            value, status, why = values.get(gene, (np.nan, "unmeasured", "gene_not_provided"))
            row.update({f"{name}_moran": value, f"{name}_status": status, f"{name}_reason": why})
        computed = row["raw_status"] == row["reconstruction_status"] == "computed"
        row["delta"] = row["reconstruction_moran"] - row["raw_moran"] if computed else np.nan
        row["comparison_status"] = "computed" if computed else "unavailable"
        rows.append(row)
    return pd.DataFrame(rows)
