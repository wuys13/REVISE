from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy import sparse
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import normalized_mutual_info_score
from sklearn.metrics import silhouette_samples


def _dedup_keep_order(items: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for item in items:
        if item not in seen:
            out.append(item)
            seen.add(item)
    return out


def _to_numpy_matrix(x) -> np.ndarray:
    if sparse.issparse(x):
        arr = x.toarray()
    elif isinstance(x, pd.DataFrame):
        arr = x.values
    else:
        arr = np.asarray(x)
    arr = arr.astype(float, copy=False)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    return arr


def _matrix_sum_mean(x) -> tuple[float, float]:
    if x is None:
        return 0.0, 0.0
    if sparse.issparse(x):
        if x.shape[0] == 0 or x.shape[1] == 0:
            return 0.0, 0.0
        return float(x.sum()), float(x.mean())
    arr = np.asarray(x)
    if arr.size == 0:
        return 0.0, 0.0
    return float(arr.sum()), float(arr.mean())


def _normalize_label_token(value: str) -> str:
    return "".join(ch for ch in str(value).lower() if ch.isalnum())


def _match_target_label(
    target: str,
    labels: np.ndarray,
    aliases: Mapping[str, Sequence[str]] | None = None,
) -> tuple[np.ndarray, str]:
    direct = labels == target
    if np.any(direct):
        return direct, target

    target_norm = _normalize_label_token(target)
    unique_labels = pd.unique(labels)
    norm_hits = [lab for lab in unique_labels if _normalize_label_token(str(lab)) == target_norm]
    if len(norm_hits) == 1:
        chosen = str(norm_hits[0])
        return labels == chosen, chosen

    alias_table: dict[str, Sequence[str]] = {
        "monocytes": ["Mono/Macro", "Mono_Macro", "Mono_macro", "Monocyte", "Monocytes"],
        "macrophage": ["Mono/Macro", "Mono_Macro", "Mono_macro", "Macrophage", "Macrophages"],
    }
    if aliases:
        alias_table.update({_normalize_label_token(k): v for k, v in aliases.items()})

    for alias in alias_table.get(target_norm, []):
        mask = labels == alias
        if np.any(mask):
            return mask, str(alias)
    return direct, target


def load_marker_map(marker_yaml: str | Path) -> dict[str, list[str]]:
    """Load a cell-type marker map from a YAML file.

    The file may either be a plain ``{cell_type: [genes...]}`` mapping or use
    the historical ``selected_marker_genes`` wrapper key.
    """
    try:
        import yaml
    except ImportError as exc:  # pragma: no cover
        raise ImportError("PyYAML is required to read marker YAML files") from exc

    path = Path(marker_yaml)
    with path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle) or {}
    raw = cfg.get("selected_marker_genes", cfg)
    if not isinstance(raw, Mapping):
        raise ValueError(f"Invalid marker YAML format: {path}")

    marker_map: dict[str, list[str]] = {}
    for key, value in raw.items():
        if isinstance(value, (list, tuple)):
            marker_map[str(key)] = [str(gene) for gene in value]
    if not marker_map:
        raise ValueError(f"No marker genes found in marker YAML: {path}")
    return marker_map


def _prepare_connectivity(
    adata: AnnData,
    connectivity_key: str,
    *,
    copy: bool,
) -> tuple[AnnData, sparse.csr_matrix]:
    work = adata.copy() if copy else adata
    if connectivity_key not in work.obsp:
        import squidpy as sq

        sq.gr.spatial_neighbors(work)
    if connectivity_key not in work.obsp:
        raise KeyError(f"Missing spatial connectivity matrix: adata.obsp[{connectivity_key!r}]")
    return work, work.obsp[connectivity_key].tocsr()


def _moran_i_for_matrix(x, connectivity: sparse.spmatrix) -> np.ndarray:
    w = connectivity.tocsr().astype(float)
    n_obs = w.shape[0]
    s0 = float(w.sum())
    x_arr = _to_numpy_matrix(x)
    if n_obs == 0 or s0 == 0 or x_arr.shape[0] != n_obs:
        return np.full(x_arr.shape[1], np.nan)

    centered = x_arr - np.nanmean(x_arr, axis=0, keepdims=True)
    denom = np.nansum(centered * centered, axis=0)
    weighted = w @ centered
    numer = np.nansum(centered * weighted, axis=0)
    values = (n_obs / s0) * numer / denom
    values[~np.isfinite(values)] = np.nan
    return values


def compute_conditional_moran_i(
    adata: AnnData,
    cell_type_col: str = "Level1",
    genes: Sequence[str] | None = None,
    connectivity_key: str = "spatial_connectivities",
    normalize: bool = True,
    log1p: bool = True,
) -> pd.DataFrame:
    """Compute MISC and MIDC from same-type and cross-type spatial edges.

    MISC is Moran's I on the spatial graph restricted to same-cell-type neighbor
    pairs. MIDC is Moran's I on the graph restricted to different-cell-type
    neighbor pairs.
    """
    if cell_type_col not in adata.obs:
        raise KeyError(f"Column {cell_type_col!r} not found in adata.obs")

    work, connectivity = _prepare_connectivity(adata, connectivity_key, copy=True)
    if normalize:
        sc.pp.normalize_total(work, target_sum=1e4)
    if log1p:
        sc.pp.log1p(work)

    if genes is None:
        selected_genes = list(map(str, work.var_names))
    else:
        selected_genes = [gene for gene in genes if gene in work.var_names]
        if not selected_genes:
            raise ValueError("None of the requested genes are present in adata.var_names")
        work = work[:, selected_genes].copy()

    labels = work.obs[cell_type_col].astype(str).to_numpy()
    coo = connectivity.tocoo()
    same_mask = labels[coo.row] == labels[coo.col]
    diff_mask = ~same_mask
    same = sparse.coo_matrix(
        (coo.data[same_mask], (coo.row[same_mask], coo.col[same_mask])),
        shape=coo.shape,
    ).tocsr()
    diff = sparse.coo_matrix(
        (coo.data[diff_mask], (coo.row[diff_mask], coo.col[diff_mask])),
        shape=coo.shape,
    ).tocsr()

    misc = _moran_i_for_matrix(work.X, same)
    midc = _moran_i_for_matrix(work.X, diff)
    return pd.DataFrame(
        {
            "Gene": selected_genes,
            "MISC": misc,
            "MIDC": midc,
            "n_same_edges": int(same.nnz),
            "n_diff_edges": int(diff.nnz),
        }
    )


def summarize_conditional_moran_i(moran_df: pd.DataFrame) -> pd.DataFrame:
    """Summarize per-gene MISC/MIDC values into a one-row table."""
    return pd.DataFrame(
        {
            "n_genes": [int(moran_df.shape[0])],
            "MISC_mean": [float(pd.to_numeric(moran_df["MISC"], errors="coerce").mean())],
            "MIDC_mean": [float(pd.to_numeric(moran_df["MIDC"], errors="coerce").mean())],
            "MISC_median": [float(pd.to_numeric(moran_df["MISC"], errors="coerce").median())],
            "MIDC_median": [float(pd.to_numeric(moran_df["MIDC"], errors="coerce").median())],
        }
    )


def compute_tmp_mer(
    adata: AnnData,
    label_col: str,
    marker_map: Mapping[str, Sequence[str]] | None = None,
    marker_yaml: str | Path | None = None,
    eps: float = 1e-8,
    aliases: Mapping[str, Sequence[str]] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute Target Marker Purity and Marker Enrichment Ratio.

    TMP compares on-target marker signal with all off-target marker signal. MER
    compares mean on-target expression with mean off-target expression.
    """
    if label_col not in adata.obs:
        raise KeyError(f"Column {label_col!r} not found in adata.obs")
    if marker_map is None:
        if marker_yaml is None:
            raise ValueError("Either marker_map or marker_yaml must be provided")
        marker_map = load_marker_map(marker_yaml)

    labels = adata.obs[label_col].astype(str).to_numpy()
    lower_to_actual: dict[str, str] = {}
    for gene in map(str, adata.var_names):
        lower_to_actual.setdefault(gene.lower(), gene)

    rows: list[dict[str, object]] = []
    for target, markers in marker_map.items():
        on_defined = _dedup_keep_order(str(gene) for gene in markers)
        off_defined = _dedup_keep_order(
            str(gene)
            for other_target, other_markers in marker_map.items()
            if other_target != target
            for gene in other_markers
        )
        on_set = {gene.lower() for gene in on_defined}
        off_defined = [gene for gene in off_defined if gene.lower() not in on_set]
        on_used = [lower_to_actual[gene.lower()] for gene in on_defined if gene.lower() in lower_to_actual]
        off_used = [lower_to_actual[gene.lower()] for gene in off_defined if gene.lower() in lower_to_actual]

        mask, matched_label = _match_target_label(str(target), labels, aliases=aliases)
        n_cells = int(mask.sum())
        row: dict[str, object] = {
            "target_cell_type": str(target),
            "label_used": matched_label,
            "metric_status": "ok",
            "n_cells": n_cells,
            "n_on_markers_defined": len(on_defined),
            "n_off_markers_defined": len(off_defined),
            "n_on_markers_used": len(on_used),
            "n_off_markers_used": len(off_used),
            "marker_coverage_on": len(on_used) / len(on_defined) if on_defined else np.nan,
            "marker_coverage_off": len(off_used) / len(off_defined) if off_defined else np.nan,
            "on_markers_used": ",".join(on_used),
            "off_markers_used": ",".join(off_used),
        }

        if n_cells == 0:
            row["metric_status"] = "no_cells_for_label"
        elif len(on_used) == 0 or len(off_used) == 0:
            row["metric_status"] = "marker_missing"

        if row["metric_status"] != "ok":
            row.update(
                {
                    "TMP": np.nan,
                    "MER": np.nan,
                    "log2MER": np.nan,
                    "off_target_contamination": np.nan,
                    "mean_on": np.nan,
                    "mean_off": np.nan,
                    "sum_on": np.nan,
                    "sum_off": np.nan,
                }
            )
            rows.append(row)
            continue

        x_on = adata[mask, on_used].X
        x_off = adata[mask, off_used].X
        sum_on, mean_on = _matrix_sum_mean(x_on)
        sum_off, mean_off = _matrix_sum_mean(x_off)
        row.update(
            {
                "TMP": float(sum_on / (sum_on + sum_off + eps)),
                "MER": float(mean_on / (mean_off + eps)),
                "log2MER": float(math.log2((mean_on + eps) / (mean_off + eps))),
                "off_target_contamination": float(mean_off),
                "mean_on": float(mean_on),
                "mean_off": float(mean_off),
                "sum_on": float(sum_on),
                "sum_off": float(sum_off),
            }
        )
        rows.append(row)

    per_target = pd.DataFrame(rows)
    return per_target, summarize_tmp_mer(per_target)


def summarize_tmp_mer(tmp_mer_df: pd.DataFrame) -> pd.DataFrame:
    """Macro-average successful TMP/MER rows."""
    ok = tmp_mer_df[tmp_mer_df["metric_status"] == "ok"].copy()
    return pd.DataFrame(
        {
            "n_targets_ok": [int(ok.shape[0])],
            "TMP_macro": [float(ok["TMP"].mean()) if ok.shape[0] else np.nan],
            "MER_macro": [float(ok["MER"].mean()) if ok.shape[0] else np.nan],
            "log2MER_macro": [float(ok["log2MER"].mean()) if ok.shape[0] else np.nan],
            "off_target_contamination_macro": [
                float(ok["off_target_contamination"].mean()) if ok.shape[0] else np.nan
            ],
        }
    )


def shannon_entropy_from_labels(labels: Sequence[object], normalize: bool = False) -> float:
    """Compute Shannon entropy for a label vector."""
    arr = np.asarray(labels)
    arr = arr[pd.notna(arr)]
    if arr.size == 0:
        return float("nan")
    _, counts = np.unique(arr.astype(str), return_counts=True)
    probs = counts / counts.sum()
    probs = probs[probs > 0]
    entropy = float(-(probs * np.log(probs)).sum())
    if normalize and counts.size > 1:
        entropy /= float(np.log(counts.size))
    return entropy


def compute_local_label_entropy(
    adata: AnnData,
    label_col: str,
    connectivity_key: str = "spatial_connectivities",
    include_self: bool = False,
    normalize: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute local neighborhood label entropy for each observation."""
    if label_col not in adata.obs:
        raise KeyError(f"Column {label_col!r} not found in adata.obs")
    work, connectivity = _prepare_connectivity(adata, connectivity_key, copy=True)
    labels = work.obs[label_col].astype(str).to_numpy()
    w = connectivity.tocsr()
    rows: list[dict[str, object]] = []
    for i, obs_name in enumerate(map(str, work.obs_names)):
        start, end = w.indptr[i], w.indptr[i + 1]
        neighbor_idx = list(map(int, w.indices[start:end]))
        if include_self:
            neighbor_idx.append(i)
        neighbor_labels = labels[neighbor_idx] if neighbor_idx else np.asarray([], dtype=str)
        rows.append(
            {
                "obs_name": obs_name,
                "label": labels[i],
                "n_neighbors": len(neighbor_idx),
                "local_label_entropy": shannon_entropy_from_labels(neighbor_labels, normalize=normalize),
            }
        )
    per_cell = pd.DataFrame(rows)
    summary = pd.DataFrame(
        {
            "n_obs": [int(per_cell.shape[0])],
            "local_label_entropy_mean": [
                float(pd.to_numeric(per_cell["local_label_entropy"], errors="coerce").mean())
            ],
            "local_label_entropy_median": [
                float(pd.to_numeric(per_cell["local_label_entropy"], errors="coerce").median())
            ],
        }
    )
    return per_cell, summary


def compute_asw(
    adata: AnnData,
    label_col: str,
    embedding_key: str = "X_pca",
) -> float:
    """Compute average silhouette width for labels in an embedding."""
    if label_col not in adata.obs:
        raise KeyError(f"Column {label_col!r} not found in adata.obs")
    labels = adata.obs[label_col].astype(str).to_numpy()
    if embedding_key in adata.obsm:
        embedding = np.asarray(adata.obsm[embedding_key], dtype=float)
    else:
        embedding = _to_numpy_matrix(adata.X)
    n_obs = embedding.shape[0]
    n_labels = len(np.unique(labels))
    if n_obs < 3 or n_labels < 2 or n_labels >= n_obs:
        return float("nan")
    return float(np.mean(silhouette_samples(embedding, labels, metric="euclidean")))


def compute_identity_metrics(
    adata: AnnData,
    label_col: str,
    embedding_key: str = "X_pca",
    true_label_col: str | None = None,
    pred_label_col: str | None = None,
) -> pd.DataFrame:
    """Compute identity-structure metrics for reconstructed observations.

    ASW only needs one label column. ARI/NMI are optional because they require
    paired predicted and reference labels in ``adata.obs``.
    """
    rows: list[dict[str, object]] = [
        {
            "metric": "ASW",
            "value": compute_asw(adata, label_col=label_col, embedding_key=embedding_key),
            "label_col": label_col,
        }
    ]
    if true_label_col is None and pred_label_col is None:
        return pd.DataFrame(rows)
    if true_label_col is None or pred_label_col is None:
        raise ValueError("true_label_col and pred_label_col must be provided together")
    missing = [col for col in (true_label_col, pred_label_col) if col not in adata.obs]
    if missing:
        raise KeyError(f"Column(s) not found in adata.obs: {missing}")

    true_labels = pd.Categorical(adata.obs[true_label_col].astype(str)).codes
    pred_labels = pd.Categorical(adata.obs[pred_label_col].astype(str)).codes
    rows.extend(
        [
            {
                "metric": "ARI",
                "value": float(adjusted_rand_score(true_labels, pred_labels)),
                "label_col": pred_label_col,
                "reference_label_col": true_label_col,
            },
            {
                "metric": "NMI",
                "value": float(normalized_mutual_info_score(true_labels, pred_labels)),
                "label_col": pred_label_col,
                "reference_label_col": true_label_col,
            },
        ]
    )
    return pd.DataFrame(rows)


def make_cell_type_mean_baseline(
    adata: AnnData,
    label_col: str,
    reference_adata: AnnData | None = None,
    reference_label_col: str | None = None,
) -> AnnData:
    """Create a negative-control baseline using cell-type mean expression."""
    if label_col not in adata.obs:
        raise KeyError(f"Column {label_col!r} not found in adata.obs")
    ref = adata if reference_adata is None else reference_adata
    ref_label_col = label_col if reference_label_col is None else reference_label_col
    if ref_label_col not in ref.obs:
        raise KeyError(f"Column {ref_label_col!r} not found in reference_adata.obs")

    genes = list(adata.var_names.intersection(ref.var_names))
    if not genes:
        raise ValueError("No shared genes between target AnnData and reference AnnData")

    target = adata[:, genes].copy()
    ref_sub = ref[:, genes].copy()
    ref_x = _to_numpy_matrix(ref_sub.X)
    ref_labels = ref_sub.obs[ref_label_col].astype(str).to_numpy()
    target_labels = target.obs[label_col].astype(str).to_numpy()
    global_mean = ref_x.mean(axis=0)
    means = {
        label: ref_x[ref_labels == label].mean(axis=0)
        for label in pd.unique(ref_labels)
        if np.any(ref_labels == label)
    }
    baseline_x = np.vstack([means.get(label, global_mean) for label in target_labels])
    target.X = baseline_x.astype(float, copy=False)
    target.uns["cell_type_mean_baseline"] = {
        "label_col": label_col,
        "reference_label_col": ref_label_col,
        "n_reference_obs": int(ref_sub.n_obs),
        "n_genes": int(len(genes)),
    }
    return target
