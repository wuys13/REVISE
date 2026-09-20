from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from anndata import AnnData
from scipy import sparse
from sklearn.metrics import (
    adjusted_mutual_info_score,
    adjusted_rand_score,
    normalized_mutual_info_score,
    silhouette_score,
)

ResultTable = pd.DataFrame
EvaluationResults = dict[str, ResultTable | dict[str, object]]

_RESULT_FILES = {
    "identity_metrics": "identity_metrics.csv",
    "tmp_mer": "tmp_mer_by_cell_type.csv",
    "tmp_mer_summary": "tmp_mer_summary.csv",
    "conditional_moran": "misc_midc_by_gene.csv",
    "conditional_moran_summary": "misc_midc_summary.csv",
    "global_moran": "moran_all_by_gene.csv",
    "cell_type_moran": "moran_by_cell_type_and_gene.csv",
}


def _matrix_sum_mean(matrix) -> tuple[float, float]:
    if sparse.issparse(matrix):
        return float(matrix.sum()), float(matrix.mean())
    values = np.asarray(matrix)
    return float(values.sum()), float(values.mean())


def _deduplicate(items: Sequence[str]) -> list[str]:
    result: list[str] = []
    seen: set[str] = set()
    for item in items:
        token = str(item)
        folded = token.casefold()
        if folded not in seen:
            result.append(token)
            seen.add(folded)
    return result


def _normalize_label(value: object) -> str:
    return "".join(char for char in str(value).lower() if char.isalnum())


def _resolve_label(
    target: str,
    labels: np.ndarray,
    marker_aliases: Mapping[str, Sequence[str]] | None,
) -> tuple[np.ndarray, str]:
    direct = labels == target
    if np.any(direct):
        return direct, target

    target_normalized = _normalize_label(target)
    unique_labels = pd.unique(labels)
    normalized_matches = [
        str(label)
        for label in unique_labels
        if _normalize_label(label) == target_normalized
    ]
    if len(normalized_matches) == 1:
        matched = normalized_matches[0]
        return labels == matched, matched

    aliases = marker_aliases.get(target, ()) if marker_aliases else ()
    for alias in aliases:
        alias_matches = labels == str(alias)
        if np.any(alias_matches):
            return alias_matches, str(alias)
    for alias in aliases:
        alias_normalized = _normalize_label(alias)
        normalized_matches = [
            str(label)
            for label in unique_labels
            if _normalize_label(label) == alias_normalized
        ]
        if len(normalized_matches) == 1:
            matched = normalized_matches[0]
            return labels == matched, matched
    return direct, target


def _select_genes(adata: AnnData, genes: Sequence[str] | None) -> list[str]:
    if genes is None:
        return [str(gene) for gene in adata.var_names]
    selected = _deduplicate([str(gene) for gene in genes])
    missing = [gene for gene in selected if gene not in adata.var_names]
    if missing:
        raise KeyError(f"Genes not found in adata.var_names: {missing}")
    return selected


def _prepare_spatial_expression(adata: AnnData, genes: Sequence[str]) -> AnnData:
    work = adata.copy()
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    return work[:, list(genes)].copy()


def _connectivity_prefix(connectivity_key: str) -> str:
    suffix = "_connectivities"
    if not connectivity_key.endswith(suffix):
        raise KeyError(
            f"Cannot build missing connectivity key {connectivity_key!r}; "
            f"generated keys must end with {suffix!r}"
        )
    return connectivity_key[: -len(suffix)]


def _get_or_build_connectivity(
    adata: AnnData,
    *,
    spatial_key: str,
    connectivity_key: str,
) -> sparse.csr_matrix:
    if connectivity_key not in adata.obsp:
        if spatial_key not in adata.obsm:
            raise KeyError(f"Missing spatial coordinates: adata.obsm[{spatial_key!r}]")
        sq.gr.spatial_neighbors(
            adata,
            spatial_key=spatial_key,
            key_added=_connectivity_prefix(connectivity_key),
        )
    return adata.obsp[connectivity_key].tocsr()


def _squidpy_moran_i(
    adata: AnnData,
    *,
    connectivity_key: str,
    genes: Sequence[str],
) -> np.ndarray:
    result = sq.gr.spatial_autocorr(
        adata,
        connectivity_key=connectivity_key,
        genes=list(genes),
        mode="moran",
        transformation=True,
        n_perms=None,
        corr_method=None,
        use_raw=False,
        copy=True,
    )
    return result.loc[list(genes), "I"].to_numpy(dtype=float)


def compute_identity_metrics(
    adata: AnnData,
    *,
    cell_type_col: str,
    pred_label_col: str,
    embedding_key: str,
) -> pd.DataFrame:
    """Compute ARI, AMI, NMI and ASW from existing labels and an embedding."""
    missing = [
        column for column in (cell_type_col, pred_label_col) if column not in adata.obs
    ]
    if missing:
        raise KeyError(f"Columns not found in adata.obs: {missing}")
    if embedding_key not in adata.obsm:
        raise KeyError(f"Embedding not found: adata.obsm[{embedding_key!r}]")

    reference = adata.obs[cell_type_col].astype(str).to_numpy()
    predicted = adata.obs[pred_label_col].astype(str).to_numpy()
    embedding = np.asarray(adata.obsm[embedding_key], dtype=float)
    unique_reference = np.unique(reference)
    if len(unique_reference) < 2 or len(unique_reference) >= adata.n_obs:
        asw = np.nan
    else:
        asw = float(silhouette_score(embedding, reference, metric="euclidean"))

    values = {
        "ARI": adjusted_rand_score(reference, predicted),
        "AMI": adjusted_mutual_info_score(reference, predicted),
        "NMI": normalized_mutual_info_score(reference, predicted),
        "ASW": asw,
    }
    return pd.DataFrame(
        {
            "metric": list(values),
            "value": [float(value) for value in values.values()],
            "reference_label_col": cell_type_col,
            "pred_label_col": pred_label_col,
            "embedding_key": embedding_key,
        }
    )


def compute_tmp_mer(
    adata: AnnData,
    *,
    cell_type_col: str,
    marker_map: Mapping[str, Sequence[str]],
    marker_aliases: Mapping[str, Sequence[str]] | None = None,
    eps: float = 1e-8,
) -> pd.DataFrame:
    """Compute cell-type-level TMP and MER from ``adata.X`` as supplied."""
    if cell_type_col not in adata.obs:
        raise KeyError(f"Column not found in adata.obs: {cell_type_col!r}")
    if not marker_map:
        raise ValueError("marker_map must not be empty")
    if not math.isfinite(eps) or eps <= 0:
        raise ValueError("eps must be a positive finite number")

    normalized_map = {
        str(target): _deduplicate([str(gene) for gene in markers])
        for target, markers in marker_map.items()
    }
    labels = adata.obs[cell_type_col].astype(str).to_numpy()
    gene_lookup = {str(gene).casefold(): str(gene) for gene in adata.var_names}
    matched_labels: dict[str, tuple[np.ndarray, str]] = {}
    label_owners: dict[str, str] = {}
    for target in normalized_map:
        mask, label_used = _resolve_label(target, labels, marker_aliases)
        if np.any(mask):
            if label_used in label_owners:
                raise ValueError(
                    f"Marker targets {label_owners[label_used]!r} and {target!r} "
                    f"both resolve to {label_used!r}"
                )
            label_owners[label_used] = target
        matched_labels[target] = mask, label_used

    rows: list[dict[str, object]] = []
    for target, markers in normalized_map.items():
        on_defined = markers
        off_defined = _deduplicate(
            [
                gene
                for other_target, other_markers in normalized_map.items()
                if other_target != target
                for gene in other_markers
            ]
        )
        on_folded = {gene.casefold() for gene in on_defined}
        off_defined = [gene for gene in off_defined if gene.casefold() not in on_folded]
        on_used = [
            gene_lookup[gene.casefold()]
            for gene in on_defined
            if gene.casefold() in gene_lookup
        ]
        off_used = [
            gene_lookup[gene.casefold()]
            for gene in off_defined
            if gene.casefold() in gene_lookup
        ]
        mask, label_used = matched_labels[target]
        n_cells = int(mask.sum())
        status = "ok"
        if n_cells == 0:
            status = "no_cells_for_label"
        elif not on_used or not off_used:
            status = "marker_missing"

        row: dict[str, object] = {
            "target_cell_type": target,
            "label_used": label_used,
            "metric_status": status,
            "n_cells": n_cells,
            "n_on_markers_defined": len(on_defined),
            "n_off_markers_defined": len(off_defined),
            "n_on_markers_used": len(on_used),
            "n_off_markers_used": len(off_used),
            "marker_coverage_on": len(on_used) / len(on_defined)
            if on_defined
            else np.nan,
            "marker_coverage_off": len(off_used) / len(off_defined)
            if off_defined
            else np.nan,
            "on_markers_used": ",".join(on_used),
            "off_markers_used": ",".join(off_used),
        }
        if status == "ok":
            sum_on, mean_on = _matrix_sum_mean(adata[mask, on_used].X)
            sum_off, mean_off = _matrix_sum_mean(adata[mask, off_used].X)
            row.update(
                {
                    "TMP": float(sum_on / (sum_on + sum_off + eps)),
                    "MER": float(mean_on / (mean_off + eps)),
                    "log2MER": float(math.log2((mean_on + eps) / (mean_off + eps))),
                    "off_target_contamination": mean_off,
                    "mean_on": mean_on,
                    "mean_off": mean_off,
                    "sum_on": sum_on,
                    "sum_off": sum_off,
                }
            )
        else:
            row.update(
                {
                    key: np.nan
                    for key in (
                        "TMP",
                        "MER",
                        "log2MER",
                        "off_target_contamination",
                        "mean_on",
                        "mean_off",
                        "sum_on",
                        "sum_off",
                    )
                }
            )
        rows.append(row)
    return pd.DataFrame(rows)


def compute_conditional_moran_i(
    adata: AnnData,
    *,
    cell_type_col: str,
    genes: Sequence[str] | None = None,
    spatial_key: str = "spatial",
    connectivity_key: str = "spatial_connectivities",
) -> pd.DataFrame:
    """Compute per-gene MISC and MIDC on same- and different-type edges."""
    if cell_type_col not in adata.obs:
        raise KeyError(f"Column not found in adata.obs: {cell_type_col!r}")
    selected_genes = _select_genes(adata, genes)
    work = _prepare_spatial_expression(adata, selected_genes)
    connectivity = _get_or_build_connectivity(
        work,
        spatial_key=spatial_key,
        connectivity_key=connectivity_key,
    )
    labels = work.obs[cell_type_col].astype(str).to_numpy()
    edges = connectivity.tocoo()
    same_mask = labels[edges.row] == labels[edges.col]
    same = sparse.coo_matrix(
        (edges.data[same_mask], (edges.row[same_mask], edges.col[same_mask])),
        shape=edges.shape,
    ).tocsr()
    different = sparse.coo_matrix(
        (edges.data[~same_mask], (edges.row[~same_mask], edges.col[~same_mask])),
        shape=edges.shape,
    ).tocsr()
    work.obsp[connectivity_key] = same
    misc = _squidpy_moran_i(
        work,
        connectivity_key=connectivity_key,
        genes=selected_genes,
    )
    work.obsp[connectivity_key] = different
    midc = _squidpy_moran_i(
        work,
        connectivity_key=connectivity_key,
        genes=selected_genes,
    )
    return pd.DataFrame(
        {
            "Gene": selected_genes,
            "MISC": misc,
            "MIDC": midc,
            "n_same_edges": int(same.nnz),
            "n_diff_edges": int(different.nnz),
        }
    )


def compute_global_moran_i(
    adata: AnnData,
    *,
    genes: Sequence[str] | None = None,
    spatial_key: str = "spatial",
    connectivity_key: str = "spatial_connectivities",
) -> pd.DataFrame:
    """Compute per-gene Moran's I across all spatial units."""
    selected_genes = _select_genes(adata, genes)
    work = _prepare_spatial_expression(adata, selected_genes)
    connectivity = _get_or_build_connectivity(
        work,
        spatial_key=spatial_key,
        connectivity_key=connectivity_key,
    )
    return pd.DataFrame(
        {
            "Gene": selected_genes,
            "group": "All",
            "MoranI": _squidpy_moran_i(
                work,
                connectivity_key=connectivity_key,
                genes=selected_genes,
            ),
            "n_units": int(work.n_obs),
            "n_edges": int(connectivity.nnz),
        }
    )


def compute_cell_type_moran_i(
    adata: AnnData,
    *,
    cell_type_col: str,
    genes: Sequence[str] | None = None,
    spatial_key: str = "spatial",
    connectivity_key: str = "spatial_connectivities",
    min_cell_type_size: int = 10,
) -> pd.DataFrame:
    """Compute per-gene Moran's I after independently graphing each cell type."""
    if cell_type_col not in adata.obs:
        raise KeyError(f"Column not found in adata.obs: {cell_type_col!r}")
    if min_cell_type_size < 1:
        raise ValueError("min_cell_type_size must be at least 1")
    selected_genes = _select_genes(adata, genes)
    labels = adata.obs[cell_type_col].astype(str)
    rows: list[pd.DataFrame] = []
    for cell_type in sorted(pd.unique(labels)):
        selected = labels == cell_type
        n_units = int(selected.sum())
        if n_units < min_cell_type_size:
            continue
        work = _prepare_spatial_expression(adata[selected].copy(), selected_genes)
        if connectivity_key in work.obsp:
            del work.obsp[connectivity_key]
        connectivity = _get_or_build_connectivity(
            work,
            spatial_key=spatial_key,
            connectivity_key=connectivity_key,
        )
        rows.append(
            pd.DataFrame(
                {
                    "Gene": selected_genes,
                    "group": str(cell_type),
                    "MoranI": _squidpy_moran_i(
                        work,
                        connectivity_key=connectivity_key,
                        genes=selected_genes,
                    ),
                    "n_units": n_units,
                    "n_edges": int(connectivity.nnz),
                }
            )
        )
    columns = ["Gene", "group", "MoranI", "n_units", "n_edges"]
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=columns)


def evaluate_adata(
    adata: AnnData,
    *,
    cell_type_col: str,
    pred_label_col: str,
    embedding_key: str,
    marker_map: Mapping[str, Sequence[str]],
    marker_aliases: Mapping[str, Sequence[str]] | None = None,
    genes: Sequence[str] | None = None,
    spatial_key: str = "spatial",
    connectivity_key: str = "spatial_connectivities",
    min_cell_type_size: int = 10,
) -> EvaluationResults:
    """Compute all standalone Figure 3 biological-recovery metrics."""
    selected_genes = _select_genes(adata, genes)
    identity = compute_identity_metrics(
        adata,
        cell_type_col=cell_type_col,
        pred_label_col=pred_label_col,
        embedding_key=embedding_key,
    )
    tmp_mer = compute_tmp_mer(
        adata,
        cell_type_col=cell_type_col,
        marker_map=marker_map,
        marker_aliases=marker_aliases,
    )
    conditional = compute_conditional_moran_i(
        adata,
        cell_type_col=cell_type_col,
        genes=selected_genes,
        spatial_key=spatial_key,
        connectivity_key=connectivity_key,
    )
    global_moran = compute_global_moran_i(
        adata,
        genes=selected_genes,
        spatial_key=spatial_key,
        connectivity_key=connectivity_key,
    )
    cell_type_moran = compute_cell_type_moran_i(
        adata,
        cell_type_col=cell_type_col,
        genes=selected_genes,
        spatial_key=spatial_key,
        connectivity_key=connectivity_key,
        min_cell_type_size=min_cell_type_size,
    )

    successful_targets = tmp_mer[tmp_mer["metric_status"] == "ok"]
    tmp_mer_summary = pd.DataFrame(
        {
            "n_targets_ok": [int(successful_targets.shape[0])],
            "TMP_macro": [float(successful_targets["TMP"].mean())],
            "MER_macro": [float(successful_targets["MER"].mean())],
            "log2MER_macro": [float(successful_targets["log2MER"].mean())],
            "off_target_contamination_macro": [
                float(successful_targets["off_target_contamination"].mean())
            ],
        }
    )
    conditional_summary = pd.DataFrame(
        {
            "n_genes": [int(conditional.shape[0])],
            "MISC_mean": [float(conditional["MISC"].mean())],
            "MIDC_mean": [float(conditional["MIDC"].mean())],
            "MISC_median": [float(conditional["MISC"].median())],
            "MIDC_median": [float(conditional["MIDC"].median())],
        }
    )
    cell_type_counts = {
        str(label): int(count)
        for label, count in adata.obs[cell_type_col].astype(str).value_counts().items()
    }
    skipped = sorted(
        label for label, count in cell_type_counts.items() if count < min_cell_type_size
    )
    run_metadata: dict[str, object] = {
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "cell_type_col": cell_type_col,
        "pred_label_col": pred_label_col,
        "embedding_key": embedding_key,
        "spatial_key": spatial_key,
        "connectivity_key": connectivity_key,
        "min_cell_type_size": int(min_cell_type_size),
        "cell_type_counts": cell_type_counts,
        "skipped_cell_types": skipped,
        "marker_map": {
            str(target): [str(gene) for gene in markers]
            for target, markers in marker_map.items()
        },
        "marker_aliases": {
            str(target): [str(alias) for alias in aliases]
            for target, aliases in (marker_aliases or {}).items()
        },
        "genes": selected_genes,
    }
    return {
        "identity_metrics": identity,
        "tmp_mer": tmp_mer,
        "tmp_mer_summary": tmp_mer_summary,
        "conditional_moran": conditional,
        "conditional_moran_summary": conditional_summary,
        "global_moran": global_moran,
        "cell_type_moran": cell_type_moran,
        "run_metadata": run_metadata,
    }


def save_evaluation_results(
    results: Mapping[str, ResultTable | Mapping[str, object]],
    output_dir: str | Path,
    *,
    metadata: Mapping[str, object] | None = None,
) -> None:
    """Save detailed metric tables and concise human-readable run metadata."""
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    for key, filename in _RESULT_FILES.items():
        table = results[key]
        if not isinstance(table, pd.DataFrame):
            raise TypeError(f"results[{key!r}] must be a pandas DataFrame")
        table.to_csv(destination / filename, index=False, float_format="%.8f")

    run_metadata = dict(results.get("run_metadata", {}))
    if metadata:
        run_metadata.update(metadata)
    with (destination / "run_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(run_metadata, handle, ensure_ascii=False, indent=2, default=str)
        handle.write("\n")
