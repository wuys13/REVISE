"""Observation-paired partition comparison helpers for reconstruction analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Hashable, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.optimize import linear_sum_assignment
from scipy.sparse import issparse
from sklearn.metrics import (
    adjusted_mutual_info_score,
    adjusted_rand_score,
    f1_score,
    mutual_info_score,
    normalized_mutual_info_score,
)


@dataclass(frozen=True)
class PartitionComparison:
    """Canonical outputs from one directed partition-comparison edge."""

    summary: pd.DataFrame
    mapping: pd.DataFrame
    contingency: pd.DataFrame
    assignments: pd.DataFrame


def select_level1_resolution(sweep: pd.DataFrame) -> pd.Series:
    """Select the highest-Level1-ARI resolution, breaking ties toward lower values."""
    required = {"resolution", "ARI"}
    missing = required - set(sweep.columns)
    if missing:
        raise KeyError(f"Resolution sweep is missing columns: {sorted(missing)}")

    candidates = sweep.loc[:, list(sweep.columns)].copy()
    candidates["resolution"] = pd.to_numeric(candidates["resolution"], errors="coerce")
    candidates["ARI"] = pd.to_numeric(candidates["ARI"], errors="coerce")
    candidates = candidates.loc[
        np.isfinite(candidates["resolution"])
        & (candidates["resolution"] > 0)
        & np.isfinite(candidates["ARI"])
    ]
    if candidates.empty:
        raise ValueError("Resolution sweep has no finite ARI candidates")
    return candidates.sort_values(["ARI", "resolution"], ascending=[False, True]).iloc[0]


def align_observation_pairs(
    raw: AnnData,
    reconstructed: AnnData,
) -> tuple[AnnData, AnnData, dict[str, int | bool]]:
    """Copy, strictly pair and reorder Raw/Reconstructed observations and genes."""
    if not raw.obs_names.is_unique or not reconstructed.obs_names.is_unique:
        raise ValueError("Raw and reconstructed observation IDs must be unique")
    if set(raw.obs_names) != set(reconstructed.obs_names):
        raise ValueError("Raw and reconstructed inputs must contain the same observation IDs")

    raw_copy = raw.copy()
    reordered = not raw_copy.obs_names.equals(reconstructed.obs_names)
    recon_copy = reconstructed[raw_copy.obs_names, :].copy()
    shared_genes = raw_copy.var_names.intersection(recon_copy.var_names, sort=False)
    if len(shared_genes) == 0:
        raise ValueError("Raw and reconstructed inputs must share at least one gene")

    raw_copy = raw_copy[:, shared_genes].copy()
    recon_copy = recon_copy[:, shared_genes].copy()
    return raw_copy, recon_copy, {
        "n_units": int(raw_copy.n_obs),
        "n_shared_genes": int(len(shared_genes)),
        "reordered_reconstruction": reordered,
    }


def _feature_variance(matrix) -> np.ndarray:
    if issparse(matrix):
        mean = np.asarray(matrix.mean(axis=0)).ravel()
        mean_square = np.asarray(matrix.multiply(matrix).mean(axis=0)).ravel()
    else:
        values = np.asarray(matrix, dtype=float)
        mean = values.mean(axis=0)
        mean_square = np.square(values).mean(axis=0)
    return np.maximum(mean_square - np.square(mean), 0.0)


def select_shared_feature_names(
    raw: AnnData,
    reconstructed: AnnData,
    *,
    n_top_genes: int = 2000,
) -> list[str]:
    """Choose one deterministic shared feature set from paired expression matrices."""
    if raw.var_names.tolist() != reconstructed.var_names.tolist():
        raise ValueError("Raw and reconstructed inputs must have identical ordered genes")
    if n_top_genes < 1:
        raise ValueError("n_top_genes must be at least one")
    if raw.n_vars <= n_top_genes:
        return raw.var_names.astype(str).tolist()
    scores = (_feature_variance(raw.X) + _feature_variance(reconstructed.X)) / 2
    order = np.argsort(-scores, kind="stable")[:n_top_genes]
    selected = set(raw.var_names[order].astype(str))
    return [str(gene) for gene in raw.var_names if str(gene) in selected]


def prepare_leiden_graph(
    adata: AnnData,
    *,
    feature_names: Sequence[str],
    n_neighbors: int = 10,
    n_pcs: int = 40,
) -> AnnData:
    """Build an expression graph without dropping paired observations or features."""
    import scanpy as sc

    if not feature_names:
        raise ValueError("feature_names must not be empty")
    if not set(feature_names) <= set(adata.var_names):
        raise ValueError("All feature_names must occur in adata.var_names")
    work = adata[:, list(feature_names)].copy()
    if work.n_obs < 3 or work.n_vars < 2:
        raise ValueError("Leiden partitioning requires at least three units and two features")
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    usable_pcs = min(n_pcs, work.n_obs - 1, work.n_vars - 1)
    if usable_pcs < 1:
        raise ValueError("Leiden partitioning has no usable principal components")
    sc.tl.pca(work, n_comps=usable_pcs, svd_solver="arpack")
    sc.pp.neighbors(work, n_neighbors=min(n_neighbors, work.n_obs - 1), n_pcs=usable_pcs)
    return work


def leiden_labels(
    graph_adata: AnnData,
    *,
    resolution: float,
    random_state: int = 42,
) -> pd.Series:
    """Return deterministic Leiden labels from a prepared expression graph."""
    import scanpy as sc

    if not np.isfinite(resolution) or resolution <= 0:
        raise ValueError("resolution must be positive and finite")
    work = graph_adata.copy()
    key = "reconstruction_impact_leiden"
    sc.tl.leiden(work, resolution=float(resolution), key_added=key, random_state=random_state)
    return work.obs[key].astype(str).rename(key)


def _as_label_series(values: pd.Series | Sequence[Hashable], name: str) -> pd.Series:
    series = values.copy() if isinstance(values, pd.Series) else pd.Series(values)
    if not series.index.is_unique:
        raise ValueError(f"{name} labels must have unique unit IDs")
    if series.isna().any():
        raise ValueError(f"{name} labels must not contain missing values")
    return series.astype(str)


def _variation_of_information(left: pd.Series, right: pd.Series) -> float:
    left_probs = left.value_counts(normalize=True)
    right_probs = right.value_counts(normalize=True)
    entropy_left = float(-(left_probs * np.log(left_probs)).sum())
    entropy_right = float(-(right_probs * np.log(right_probs)).sum())
    return float(entropy_left + entropy_right - 2 * mutual_info_score(left, right))


def compare_partitions(
    raw_labels: pd.Series | Sequence[Hashable],
    reconstructed_labels: pd.Series | Sequence[Hashable],
    *,
    comparison_edge: str = "raw_to_recon_expression",
) -> PartitionComparison:
    """Compare paired partitions using maximum-overlap Hungarian label matching."""
    raw = _as_label_series(raw_labels, "Raw")
    recon = _as_label_series(reconstructed_labels, "Reconstructed")
    if set(raw.index) != set(recon.index):
        raise ValueError("Raw and reconstructed labels must contain the same unit IDs")
    recon = recon.reindex(raw.index)

    raw_categories = sorted(raw.unique().tolist())
    recon_categories = sorted(recon.unique().tolist())
    contingency = pd.crosstab(raw, recon).reindex(
        index=raw_categories, columns=recon_categories, fill_value=0
    )
    row_indices, col_indices = linear_sum_assignment(-contingency.to_numpy())
    matched_pairs = {
        raw_categories[row_index]: recon_categories[col_index]
        for row_index, col_index in zip(row_indices, col_indices)
    }

    mapping_rows: list[dict[str, object]] = []
    matched_recon = set(matched_pairs.values())
    for raw_cluster in raw_categories:
        recon_cluster = matched_pairs.get(raw_cluster)
        raw_size = int(contingency.loc[raw_cluster].sum())
        overlap = int(contingency.loc[raw_cluster, recon_cluster]) if recon_cluster else 0
        recon_size = int(contingency[recon_cluster].sum()) if recon_cluster else 0
        precision = overlap / recon_size if recon_size else np.nan
        recall = overlap / raw_size if raw_size else np.nan
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        mapping_rows.append(
            {
                "raw_cluster": raw_cluster,
                "recon_cluster": recon_cluster,
                "matched": recon_cluster is not None,
                "overlap_n": overlap,
                "raw_cluster_n": raw_size,
                "recon_cluster_n": recon_size,
                "precision": precision,
                "recall": recall,
                "f1": f1,
            }
        )
    for recon_cluster in recon_categories:
        if recon_cluster not in matched_recon:
            mapping_rows.append(
                {
                    "raw_cluster": pd.NA,
                    "recon_cluster": recon_cluster,
                    "matched": False,
                    "overlap_n": 0,
                    "raw_cluster_n": 0,
                    "recon_cluster_n": int(contingency[recon_cluster].sum()),
                    "precision": 0.0,
                    "recall": np.nan,
                    "f1": 0.0,
                }
            )
    mapping = pd.DataFrame(mapping_rows)

    mapped_raw = raw.map(matched_pairs).fillna("__unmatched_raw_cluster__")
    unit_changed = mapped_raw != recon
    assignments = pd.DataFrame(
        {
            "unit_id": raw.index.astype(str),
            "raw_cluster": raw.to_numpy(),
            "matched_raw_cluster": mapped_raw.to_numpy(),
            "recon_cluster": recon.to_numpy(),
            "unit_changed": unit_changed.to_numpy(),
            "comparison_edge": comparison_edge,
        },
        index=raw.index,
    )

    labels_for_f1 = sorted(set(mapped_raw) | set(recon))
    matched_accuracy = float(1 - unit_changed.mean())
    matched_macro_f1 = float(
        f1_score(recon, mapped_raw, labels=labels_for_f1, average="macro", zero_division=0)
    )
    summary = pd.DataFrame(
        [
            {
                "comparison_edge": comparison_edge,
                "n_units": int(raw.shape[0]),
                "n_raw_clusters": len(raw_categories),
                "n_recon_clusters": len(recon_categories),
                "matched_accuracy": matched_accuracy,
                "st_unit_change_fraction": float(1 - matched_accuracy),
                "matched_macro_f1": matched_macro_f1,
                "balanced_cluster_change": float(1 - matched_macro_f1),
                "ARI": float(adjusted_rand_score(raw, recon)),
                "AMI": float(adjusted_mutual_info_score(raw, recon)),
                "NMI": float(normalized_mutual_info_score(raw, recon)),
                "VI": _variation_of_information(raw, recon),
                "status": "ok",
            }
        ]
    )
    return PartitionComparison(summary, mapping, contingency, assignments)
