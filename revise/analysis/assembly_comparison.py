"""Fair, read-only comparison of independently generated assembly H5AD files.

The routines in this module deliberately do not reconstruct data.  They align
already generated method outputs on real observation identifiers and gene names,
then run the same expression-only preprocessing and Leiden sweep independently
for each method.  An original spatial object is used only for labels and plot
coordinates after clustering.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData, read_h5ad
from scipy import sparse
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

DEFAULT_METHODS = ("mean", "random", "within_cluster", "outside_cluster")
DEFAULT_RESOLUTIONS = (0.6, 0.7, 0.8)
INPUT_ASSUMPTIONS = {
    "generated_expression": "finite, nonnegative, unlogged linear expression declared by caller",
    "alignment": "exact obs_names and var_names; no positional pairing or zero fill",
    "working_expression": "fresh copy normalized to 1e4 per observation, then log1p",
    "baseline_role": "labels and spatial coordinates only; never expression clustering",
    "selection": "no automatic best resolution or winning method",
}


@dataclass
class LoadedAssemblyInputs:
    """Native inputs plus the hashes used to prove they were not overwritten."""

    methods: dict[str, AnnData]
    baseline: AnnData
    paths: dict[str, Path]
    hashes: dict[str, str]


@dataclass
class TypeComparison:
    """Intermediate and final objects for one broad cell type."""

    broad_type: str
    status: str
    issues: list[str]
    shared_ids: pd.Index
    shared_genes: pd.Index
    baseline_labels: pd.Series
    aligned: dict[str, AnnData] = field(default_factory=dict)
    clustered: dict[str, AnnData] = field(default_factory=dict)
    metrics: pd.DataFrame = field(default_factory=pd.DataFrame)
    contingencies: dict[tuple[str, float], pd.DataFrame] = field(default_factory=dict)


@dataclass
class AssemblyComparison:
    """Comparison result without any automatic ranking or winner selection."""

    coverage: pd.DataFrame
    metrics: pd.DataFrame
    by_type: dict[str, TypeComparison]
    assumptions: dict[str, str]


def file_sha256(path: str | Path, *, chunk_size: int = 1024 * 1024) -> str:
    """Return a streaming SHA-256 digest for a file."""

    digest = sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_assembly_inputs(
    method_paths: Mapping[str, str | Path],
    baseline_path: str | Path,
    *,
    expected_methods: Sequence[str] | None = DEFAULT_METHODS,
) -> LoadedAssemblyInputs:
    """Read native H5AD inputs and record immutable-source hashes.

    The function only reads the supplied files.  It does not write normalized or
    clustered data back to the source paths.
    """

    if not method_paths:
        raise ValueError("method_paths must contain at least one generated H5AD")
    if expected_methods is not None:
        missing = sorted(set(expected_methods).difference(method_paths))
        extra = sorted(set(method_paths).difference(expected_methods))
        if missing or extra:
            raise ValueError(
                f"method paths differ from expected methods; missing={missing}, extra={extra}"
            )

    paths = {
        name: Path(path).expanduser().resolve() for name, path in method_paths.items()
    }
    paths["baseline"] = Path(baseline_path).expanduser().resolve()
    absent = [f"{name}: {path}" for name, path in paths.items() if not path.is_file()]
    if absent:
        raise FileNotFoundError("Missing H5AD input(s): " + "; ".join(absent))

    hashes = {name: file_sha256(path) for name, path in paths.items()}
    methods = {name: read_h5ad(paths[name]) for name in method_paths}
    baseline = read_h5ad(paths["baseline"])
    return LoadedAssemblyInputs(methods, baseline, paths, hashes)


def assert_input_hashes_unchanged(inputs: LoadedAssemblyInputs) -> None:
    """Raise if any native H5AD changed after analysis."""

    changed = {}
    for name, path in inputs.paths.items():
        current = file_sha256(path)
        if current != inputs.hashes[name]:
            changed[name] = {"before": inputs.hashes[name], "after": current}
    if changed:
        raise RuntimeError(f"Native H5AD input hash changed: {changed}")


def _matrix_data(x: Any) -> np.ndarray:
    if sparse.issparse(x):
        return np.asarray(x.data)
    return np.asarray(x)


def validate_generated_expression(adata: AnnData, *, method: str) -> None:
    """Validate the caller-declared unlogged linear-expression boundary."""

    if not adata.obs_names.is_unique:
        raise ValueError(f"{method}: obs_names must be unique real identifiers")
    if not adata.var_names.is_unique:
        raise ValueError(f"{method}: var_names must be unique gene identifiers")
    values = _matrix_data(adata.X)
    if not np.issubdtype(values.dtype, np.number):
        raise TypeError(f"{method}: X must be numeric")
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{method}: X must contain only finite values")
    if np.any(values < 0):
        raise ValueError(f"{method}: X must be nonnegative unlogged linear expression")


def resolve_broad_column(
    adata: AnnData,
    requested: str | None,
    *,
    fallback: str = "revise_Level1",
    object_name: str = "object",
) -> tuple[str, bool]:
    """Resolve a broad-label column while explicitly reporting fallback use."""

    if requested is not None and requested in adata.obs:
        return requested, False
    if fallback in adata.obs:
        return fallback, requested != fallback
    raise KeyError(
        f"{object_name}: broad label column {requested!r} is unavailable and "
        f"explicit fallback {fallback!r} is also unavailable"
    )


def _alias_mask(values: pd.Series, aliases: Sequence[str]) -> pd.Series:
    cleaned_aliases = {str(value).strip().casefold() for value in aliases}
    cleaned_values = values.astype("string").str.strip().str.casefold()
    return cleaned_values.isin(cleaned_aliases).fillna(False)


def _empty_metrics() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "broad_type",
            "method",
            "resolution",
            "status",
            "matched_label_ids",
            "missing_baseline_labels",
            "baseline_classes",
            "clusters",
            "ARI",
            "NMI",
        ]
    )


def _resolution_key(resolution: float) -> str:
    return f"leiden_{format(float(resolution), '.12g')}"


def _validate_resolutions(resolutions: Sequence[float]) -> tuple[float, ...]:
    values = tuple(float(value) for value in resolutions)
    if not values or any(not np.isfinite(value) or value <= 0 for value in values):
        raise ValueError("resolutions must contain positive finite values")
    if len(set(values)) != len(values):
        raise ValueError("resolutions must not contain duplicates")
    return values


def _coverage_row(
    *,
    broad_type: str,
    object_name: str,
    original_ids: pd.Index,
    shared_ids: pd.Index,
    broad_column: str,
    fallback_used: bool,
    original_genes: pd.Index | None,
    shared_genes: pd.Index | None,
) -> dict[str, Any]:
    excluded_ids = tuple(sorted(set(original_ids).difference(shared_ids)))
    if original_genes is None or shared_genes is None:
        excluded_genes: Any = pd.NA
        original_gene_count: Any = pd.NA
        shared_gene_count: Any = pd.NA
    else:
        excluded_genes = tuple(sorted(set(original_genes).difference(shared_genes)))
        original_gene_count = len(original_genes)
        shared_gene_count = len(shared_genes)
    return {
        "broad_type": broad_type,
        "object": object_name,
        "broad_column": broad_column,
        "fallback_used": fallback_used,
        "missing_type": len(original_ids) == 0,
        "original_type_ids": len(original_ids),
        "shared_ids": len(shared_ids),
        "excluded_ids": excluded_ids,
        "original_genes": original_gene_count,
        "shared_genes": shared_gene_count,
        "excluded_genes": excluded_genes,
    }


def _prepare_type_scope(
    methods: Mapping[str, AnnData],
    baseline: AnnData,
    *,
    broad_type: str,
    aliases: Sequence[str],
    broad_column: str | None,
    broad_fallback: str,
    baseline_subtype_column: str,
) -> tuple[TypeComparison, list[dict[str, Any]]]:
    selected: dict[str, AnnData] = {}
    selected_ids: dict[str, pd.Index] = {}
    resolved_columns: dict[str, tuple[str, bool]] = {}

    for method, adata in methods.items():
        column, used_fallback = resolve_broad_column(
            adata,
            broad_column,
            fallback=broad_fallback,
            object_name=method,
        )
        resolved_columns[method] = (column, used_fallback)
        subset = adata[_alias_mask(adata.obs[column], aliases).to_numpy(), :]
        selected[method] = subset
        selected_ids[method] = subset.obs_names

    baseline_column, baseline_fallback_used = resolve_broad_column(
        baseline,
        broad_column,
        fallback=broad_fallback,
        object_name="baseline",
    )
    if baseline_subtype_column not in baseline.obs:
        raise KeyError(
            f"baseline subtype column {baseline_subtype_column!r} not found in baseline.obs"
        )
    baseline_subset = baseline[
        _alias_mask(baseline.obs[baseline_column], aliases).to_numpy(), :
    ]
    selected_ids["baseline"] = baseline_subset.obs_names

    id_sets = [set(index.astype(str)) for index in selected_ids.values()]
    shared_ids = (
        pd.Index(sorted(set.intersection(*id_sets))) if id_sets else pd.Index([])
    )
    gene_sets = [set(adata.var_names.astype(str)) for adata in methods.values()]
    shared_genes = (
        pd.Index(sorted(set.intersection(*gene_sets))) if gene_sets else pd.Index([])
    )

    coverage = []
    for method, adata in methods.items():
        column, used_fallback = resolved_columns[method]
        coverage.append(
            _coverage_row(
                broad_type=broad_type,
                object_name=method,
                original_ids=selected_ids[method],
                shared_ids=shared_ids,
                broad_column=column,
                fallback_used=used_fallback,
                original_genes=adata.var_names,
                shared_genes=shared_genes,
            )
        )
    coverage.append(
        _coverage_row(
            broad_type=broad_type,
            object_name="baseline",
            original_ids=selected_ids["baseline"],
            shared_ids=shared_ids,
            broad_column=baseline_column,
            fallback_used=baseline_fallback_used,
            original_genes=None,
            shared_genes=None,
        )
    )

    issues = []
    missing_types = [name for name, ids in selected_ids.items() if len(ids) == 0]
    if missing_types:
        issues.append("broad type absent from: " + ", ".join(missing_types))
    if len(shared_ids) < 3:
        issues.append(
            f"only {len(shared_ids)} shared real IDs; at least 3 are required"
        )
    if len(shared_genes) < 2:
        issues.append(f"only {len(shared_genes)} shared genes; at least 2 are required")

    if issues:
        baseline_labels = pd.Series(
            pd.array([], dtype="string"),
            index=pd.Index([], dtype="object"),
            name=baseline_subtype_column,
        )
        return (
            TypeComparison(
                broad_type,
                "unavailable",
                issues,
                shared_ids,
                shared_genes,
                baseline_labels,
            ),
            coverage,
        )

    aligned = {
        method: selected[method][shared_ids, shared_genes].copy() for method in methods
    }
    baseline_labels = baseline.obs.loc[shared_ids, baseline_subtype_column].copy()
    baseline_labels.name = baseline_subtype_column
    return (
        TypeComparison(
            broad_type,
            "ready",
            [],
            shared_ids,
            shared_genes,
            baseline_labels,
            aligned=aligned,
        ),
        coverage,
    )


def cluster_aligned_methods(
    aligned: Mapping[str, AnnData],
    *,
    resolutions: Sequence[float] = DEFAULT_RESOLUTIONS,
    seed: int = 42,
    target_sum: float = 1e4,
    max_pcs: int = 40,
    max_neighbors: int = 10,
) -> dict[str, AnnData]:
    """Run an equal, independent expression pipeline for every method."""

    resolution_values = _validate_resolutions(resolutions)
    if max_pcs < 1 or max_neighbors < 2:
        raise ValueError("max_pcs must be >=1 and max_neighbors must be >=2")
    clustered: dict[str, AnnData] = {}
    for method, native_scope in aligned.items():
        validate_generated_expression(native_scope, method=method)
        if native_scope.n_obs < 3 or native_scope.n_vars < 2:
            raise ValueError(
                f"{method}: clustering requires at least 3 IDs and 2 genes"
            )

        work = native_scope.copy()
        sc.pp.normalize_total(work, target_sum=target_sum)
        sc.pp.log1p(work)
        n_pcs = min(max_pcs, work.n_obs - 1, work.n_vars - 1)
        sc.pp.pca(work, n_comps=n_pcs, svd_solver="arpack", random_state=seed)
        n_neighbors = min(max_neighbors, work.n_obs - 1)
        sc.pp.neighbors(work, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=seed)
        for resolution in resolution_values:
            sc.tl.leiden(
                work,
                resolution=resolution,
                key_added=_resolution_key(resolution),
                random_state=seed,
                flavor="igraph",
                n_iterations=2,
                directed=False,
            )
        work.uns["assembly_comparison"] = {
            "assumptions": dict(INPUT_ASSUMPTIONS),
            "seed": int(seed),
            "resolutions": list(resolution_values),
            "n_pcs": int(n_pcs),
            "n_neighbors": int(n_neighbors),
            "shared_feature_order": "lexicographic exact-name intersection",
        }
        clustered[method] = work
    return clustered


def evaluate_against_baseline(
    clustered: Mapping[str, AnnData],
    baseline_labels: pd.Series,
    *,
    broad_type: str,
    resolutions: Sequence[float] = DEFAULT_RESOLUTIONS,
) -> tuple[pd.DataFrame, dict[tuple[str, float], pd.DataFrame]]:
    """Evaluate cluster labels on exact matched IDs with explicit missing handling."""

    rows: list[dict[str, Any]] = []
    contingencies: dict[tuple[str, float], pd.DataFrame] = {}
    for method, work in clustered.items():
        for resolution in _validate_resolutions(resolutions):
            cluster_key = _resolution_key(resolution)
            joined = pd.concat(
                [
                    baseline_labels.rename("baseline"),
                    work.obs[cluster_key].rename("cluster"),
                ],
                axis=1,
                join="inner",
            )
            valid = joined["baseline"].notna() & joined["cluster"].notna()
            scored = joined.loc[valid]
            contingency = pd.crosstab(
                scored["baseline"],
                scored["cluster"],
                rownames=[baseline_labels.name or "baseline"],
                colnames=[cluster_key],
                dropna=False,
            )
            contingencies[(method, resolution)] = contingency
            row: dict[str, Any] = {
                "broad_type": broad_type,
                "method": method,
                "resolution": resolution,
                "status": "ok" if len(scored) >= 2 else "unavailable",
                "matched_label_ids": len(scored),
                "missing_baseline_labels": len(joined) - len(scored),
                "baseline_classes": scored["baseline"].nunique(dropna=True),
                "clusters": scored["cluster"].nunique(dropna=True),
                "ARI": np.nan,
                "NMI": np.nan,
            }
            if len(scored) >= 2:
                row["ARI"] = adjusted_rand_score(scored["baseline"], scored["cluster"])
                row["NMI"] = normalized_mutual_info_score(
                    scored["baseline"], scored["cluster"]
                )
            rows.append(row)
    return pd.DataFrame(rows, columns=_empty_metrics().columns), contingencies


def compare_assembly_methods(
    methods: Mapping[str, AnnData],
    baseline: AnnData,
    *,
    type_aliases: Mapping[str, Sequence[str]],
    broad_column: str | None = "revise_Level1",
    broad_fallback: str = "revise_Level1",
    baseline_subtype_column: str = "SVC_cluster",
    resolutions: Sequence[float] = DEFAULT_RESOLUTIONS,
    seed: int = 42,
) -> AssemblyComparison:
    """Compare independent assemblies without ranking methods or resolutions."""

    if not methods:
        raise ValueError("methods must not be empty")
    if not type_aliases:
        raise ValueError("type_aliases must not be empty")
    resolution_values = _validate_resolutions(resolutions)
    for method, adata in methods.items():
        validate_generated_expression(adata, method=method)
    if not baseline.obs_names.is_unique:
        raise ValueError("baseline: obs_names must be unique real identifiers")

    coverage_rows: list[dict[str, Any]] = []
    metric_frames: list[pd.DataFrame] = []
    by_type: dict[str, TypeComparison] = {}
    for broad_type, aliases in type_aliases.items():
        if not aliases:
            raise ValueError(f"{broad_type}: aliases must not be empty")
        result, type_coverage = _prepare_type_scope(
            methods,
            baseline,
            broad_type=broad_type,
            aliases=aliases,
            broad_column=broad_column,
            broad_fallback=broad_fallback,
            baseline_subtype_column=baseline_subtype_column,
        )
        coverage_rows.extend(type_coverage)
        if result.status == "ready":
            result.clustered = cluster_aligned_methods(
                result.aligned,
                resolutions=resolution_values,
                seed=seed,
            )
            result.metrics, result.contingencies = evaluate_against_baseline(
                result.clustered,
                result.baseline_labels,
                broad_type=broad_type,
                resolutions=resolution_values,
            )
            result.status = "ok"
            metric_frames.append(result.metrics)
        by_type[broad_type] = result

    coverage = pd.DataFrame(coverage_rows)
    metrics = (
        pd.concat(metric_frames, ignore_index=True)
        if metric_frames
        else _empty_metrics()
    )
    return AssemblyComparison(coverage, metrics, by_type, dict(INPUT_ASSUMPTIONS))


def spatial_plot_frame(
    result: TypeComparison,
    baseline: AnnData,
    *,
    method: str,
    resolution: float,
    spatial_key: str = "spatial",
) -> pd.DataFrame:
    """Build baseline/new long-form plot data on exact shared IDs."""

    if result.status != "ok":
        raise ValueError(f"{result.broad_type}: comparison is {result.status}")
    if method not in result.clustered:
        raise KeyError(f"unknown method {method!r}")
    if spatial_key not in baseline.obsm:
        raise KeyError(f"baseline.obsm lacks spatial key {spatial_key!r}")
    coordinates = np.asarray(baseline[result.shared_ids].obsm[spatial_key])
    if coordinates.ndim != 2 or coordinates.shape[1] < 2:
        raise ValueError(
            f"baseline.obsm[{spatial_key!r}] must have at least two columns"
        )
    key = _resolution_key(resolution)
    if key not in result.clustered[method].obs:
        raise KeyError(f"resolution {resolution} was not clustered for {method}")
    base = pd.DataFrame(
        {
            "id": result.shared_ids,
            "x": coordinates[:, 0],
            "y": coordinates[:, 1],
            "panel": "baseline",
            "label": result.baseline_labels.reindex(result.shared_ids)
            .astype("string")
            .to_numpy(),
        }
    )
    new = base[["id", "x", "y"]].copy()
    new["panel"] = f"{method} | r={format(float(resolution), '.12g')}"
    new["label"] = (
        result.clustered[method]
        .obs.loc[result.shared_ids, key]
        .astype("string")
        .to_numpy()
    )
    return pd.concat([base, new], ignore_index=True)


def plot_spatial_comparison(
    result: TypeComparison,
    baseline: AnnData,
    *,
    method: str,
    resolution: float,
    spatial_key: str = "spatial",
    point_size: float = 8.0,
) -> tuple[plt.Figure, np.ndarray]:
    """Plot original subtype labels beside independently inferred clusters."""

    frame = spatial_plot_frame(
        result,
        baseline,
        method=method,
        resolution=resolution,
        spatial_key=spatial_key,
    )
    panels = list(frame["panel"].drop_duplicates())
    figure, axes = plt.subplots(
        1, len(panels), figsize=(5 * len(panels), 4), squeeze=False
    )
    for axis, panel in zip(axes[0], panels):
        subset = frame.loc[frame["panel"] == panel]
        labels = subset["label"].fillna("<missing>").astype(str)
        categories = sorted(labels.unique())
        palette = plt.get_cmap("tab20", max(len(categories), 1))
        for index, label in enumerate(categories):
            points = subset.loc[labels == label]
            axis.scatter(
                points["x"],
                points["y"],
                s=point_size,
                color="#bdbdbd" if label == "<missing>" else palette(index),
                label=label,
                linewidths=0,
            )
        axis.set_title(panel)
        axis.set_aspect("equal")
        axis.invert_yaxis()
        axis.set_xlabel("spatial x")
        axis.set_ylabel("spatial y")
        axis.legend(markerscale=2, fontsize=7, frameon=False)
    figure.suptitle(
        f"{result.broad_type}: baseline label versus new expression clustering"
    )
    figure.tight_layout()
    return figure, axes
