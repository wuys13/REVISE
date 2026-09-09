"""Observation-paired partition comparison helpers for reconstruction analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Hashable, Sequence
import warnings

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.optimize import linear_sum_assignment
from scipy.sparse import SparseEfficiencyWarning, issparse
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


def select_resolution_for_target_cluster_count(
    sweep: pd.DataFrame,
    *,
    target_cluster_count: int,
) -> pd.Series:
    """Choose the lowest resolution with the closest observed cluster count."""
    if target_cluster_count < 1:
        raise ValueError("target_cluster_count must be at least one")
    required = {"resolution", "n_clusters"}
    missing = required - set(sweep.columns)
    if missing:
        raise KeyError(f"Cluster-count sweep is missing columns: {sorted(missing)}")
    candidates = sweep.copy()
    candidates["resolution"] = pd.to_numeric(candidates["resolution"], errors="coerce")
    candidates["n_clusters"] = pd.to_numeric(candidates["n_clusters"], errors="coerce")
    candidates = candidates.loc[
        np.isfinite(candidates["resolution"])
        & (candidates["resolution"] > 0)
        & np.isfinite(candidates["n_clusters"])
        & (candidates["n_clusters"] >= 1)
    ].copy()
    if candidates.empty:
        raise ValueError("Cluster-count sweep has no finite candidates")
    candidates["n_clusters"] = candidates["n_clusters"].astype(int)
    candidates["cluster_count_gap"] = (candidates["n_clusters"] - target_cluster_count).abs()
    selected = candidates.sort_values(
        ["cluster_count_gap", "resolution"], ascending=[True, True]
    ).iloc[0].copy()
    selected["target_cluster_count"] = int(target_cluster_count)
    selected["status"] = (
        "ok" if int(selected["cluster_count_gap"]) <= 1 else "unmatched_cluster_complexity"
    )
    return selected


def summarize_change_by_level1(
    assignments: pd.DataFrame,
    level1_labels: pd.Series,
    *,
    min_report_n: int = 30,
    z_value: float = 1.96,
) -> pd.DataFrame:
    """Summarize an already-global Hungarian change call by Raw Level1 labels."""
    if "unit_changed" not in assignments:
        raise KeyError("assignments must contain unit_changed")
    if min_report_n < 1:
        raise ValueError("min_report_n must be at least one")
    if not assignments.index.is_unique or not level1_labels.index.is_unique:
        raise ValueError("assignments and Level1 labels must have unique unit IDs")
    if set(assignments.index) != set(level1_labels.index):
        raise ValueError("assignments and Level1 labels must contain the same unit IDs")

    work = pd.DataFrame(
        {
            "level1": level1_labels.reindex(assignments.index).astype(str),
            "unit_changed": assignments["unit_changed"].astype(bool),
        }
    )
    grouped = work.groupby("level1", sort=True)["unit_changed"].agg(
        total_units="size", changed_units="sum"
    )
    grouped.loc["Overall"] = [int(work.shape[0]), int(work["unit_changed"].sum())]
    grouped = grouped.reset_index()
    grouped["changed_units"] = grouped["changed_units"].astype(int)
    grouped["change_fraction"] = grouped["changed_units"] / grouped["total_units"]

    n = grouped["total_units"].astype(float)
    p = grouped["change_fraction"].astype(float)
    denominator = 1 + z_value**2 / n
    center = (p + z_value**2 / (2 * n)) / denominator
    half_width = z_value * np.sqrt((p * (1 - p) + z_value**2 / (4 * n)) / n) / denominator
    grouped["wilson_ci_lower"] = center - half_width
    grouped["wilson_ci_upper"] = center + half_width
    grouped["low_sample_size"] = grouped["total_units"] < min_report_n
    return grouped


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


def _nonzero_counts(matrix, *, axis: int) -> np.ndarray:
    if issparse(matrix):
        return np.asarray((matrix != 0).sum(axis=axis)).ravel()
    return np.count_nonzero(np.asarray(matrix), axis=axis)


def filter_paired_sp_svc_inputs(
    raw: AnnData,
    reconstructed: AnnData,
    *,
    min_genes: int = 50,
    min_cells: int = 3,
) -> tuple[AnnData, AnnData, dict[str, int | bool]]:
    """Define QC on Raw counts and apply the retained IDs and genes to both carriers."""
    if min_genes < 1 or min_cells < 1:
        raise ValueError("min_genes and min_cells must be at least one")
    raw_work, recon_work, audit = align_observation_pairs(raw, reconstructed)
    input_units = int(raw_work.n_obs)
    input_genes = int(raw_work.n_vars)
    retained_units = _nonzero_counts(raw_work.X, axis=1) >= min_genes
    raw_work = raw_work[retained_units].copy()
    recon_work = recon_work[retained_units].copy()
    after_cell_qc_units = int(raw_work.n_obs)
    if raw_work.n_obs < 3:
        raise ValueError("Raw QC retained fewer than three paired observations")
    retained_genes = (
        (_nonzero_counts(raw_work.X, axis=0) >= min_cells)
        & ~raw_work.var_names.str.startswith("MT-")
    )
    raw_work = raw_work[:, retained_genes].copy()
    recon_work = recon_work[:, retained_genes].copy()
    if raw_work.n_vars < 2:
        raise ValueError("Raw QC retained fewer than two shared genes")
    retained_after_gene_filter = _nonzero_counts(raw_work.X, axis=1) > 0
    raw_work = raw_work[retained_after_gene_filter].copy()
    recon_work = recon_work[retained_after_gene_filter].copy()
    if raw_work.n_obs < 3:
        raise ValueError("Raw QC retained fewer than three paired observations")
    audit.update(
        {
            "input_units": input_units,
            "excluded_raw_qc_units": input_units - int(raw_work.n_obs),
            "excluded_post_gene_filter_units": after_cell_qc_units - int(raw_work.n_obs),
            "n_units": int(raw_work.n_obs),
            "input_shared_genes": input_genes,
            "excluded_raw_qc_genes": input_genes - int(raw_work.n_vars),
            "n_shared_genes": int(raw_work.n_vars),
        }
    )
    return raw_work, recon_work, audit


def select_raw_hvg_feature_names(
    raw: AnnData,
    *,
    n_top_genes: int = 2000,
) -> list[str]:
    """Select one canonical-compatible Raw HVG set for both sp-SVC carriers."""
    if n_top_genes < 1:
        raise ValueError("n_top_genes must be at least one")
    if raw.n_vars <= n_top_genes:
        return raw.var_names.astype(str).tolist()
    import scanpy as sc

    work = raw.copy()
    if issparse(work.X):
        work.X = work.X.tocsc()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Some cells have zero counts")
        sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=SparseEfficiencyWarning)
        sc.pp.highly_variable_genes(
            work,
            n_top_genes=n_top_genes,
            flavor="seurat_v3",
            check_values=False,
        )
    selected = set(work.var_names[work.var["highly_variable"]].astype(str))
    return [str(gene) for gene in raw.var_names if str(gene) in selected]


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
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Some cells have zero counts")
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
    if "connectivities" not in graph_adata.obsp:
        raise KeyError("graph_adata must contain a precomputed connectivities graph")
    # Leiden uses the neighbor graph only.  Copying ``graph_adata`` here would
    # duplicate the full expression matrix once per candidate resolution.
    work = AnnData(obs=graph_adata.obs.copy())
    key = "reconstruction_impact_leiden"
    sc.tl.leiden(
        work,
        adjacency=graph_adata.obsp["connectivities"],
        resolution=float(resolution),
        key_added=key,
        random_state=random_state,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
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


def _label_sort_key(value: str) -> tuple[int, float | str, str]:
    try:
        numeric = float(value)
    except ValueError:
        return (1, value, value)
    return (0, numeric, value)


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

    raw_categories = sorted(raw.unique().tolist(), key=_label_sort_key)
    recon_categories = sorted(recon.unique().tolist(), key=_label_sort_key)
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
