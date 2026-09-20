"""Fair, read-only comparison of independently generated assembly H5AD files.

This module consumes already generated H5AD files. It never reconstructs an
input and it never writes back to a native input. A comparison scope is made
from exact observation IDs and exact gene names, then each method receives its
own normalization, PCA, neighbor graph, and Leiden sweep.

The comparison is deliberately explicit about the three current broad types.
The only label normalization performed here is the historical slash to
underscore spelling normalization on the configured broad-label column.
Missing labels, IDs, genes, and coordinates remain visible in the result
rather than being silently replaced.
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
DEFAULT_CELL_TYPES = ("T", "Mono_Macro", "Fibroblast")
DEFAULT_BROAD_COLUMN = "Level1"
DEFAULT_BASELINE_SUBTYPE_COLUMN = "SVC_cluster"
DEFAULT_SPATIAL_KEY = "spatial"
DEFAULT_RESOLUTIONS = (0.6, 0.7, 0.8)

INPUT_ASSUMPTIONS = {
    "generated_expression": "finite, nonnegative, unlogged linear expression declared by caller",
    "alignment": "exact real obs_names and exact var_names; no positional pairing or zero fill",
    "broad_labels": "explicit Level1 cell types; normalize slash to underscore and preserve missing values",
    "working_expression": "fresh copy normalized to 1e4 per observation, then log1p",
    "baseline_role": "SVC_cluster labels and spatial coordinates only; never expression clustering",
    "spatial_coordinates": (
        "each method obsm[spatial] is aligned by real ID and must be finite, shape-equal, "
        "and exactly equal to the baseline; no rescaling; source and unit are inherited "
        "from formal sample.yaml and baseline provenance"
    ),
    "selection": "all requested resolutions and methods are retained; no automatic winner",
}


@dataclass
class LoadedAssemblyInputs:
    """Native inputs plus hashes used to prove they were not overwritten."""

    methods: dict[str, AnnData]
    baseline: AnnData
    paths: dict[str, Path]
    hashes: dict[str, str]


@dataclass
class TypeComparison:
    """Intermediate and final objects for one explicit broad cell type."""

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
    coordinate_checks: dict[str, dict[str, Any]] = field(default_factory=dict)


@dataclass
class AssemblyComparison:
    """Comparison result without automatic ranking or winner selection."""

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

    The function only reads the supplied files. It does not write normalized
    or clustered data back to the source paths.
    """

    if not method_paths:
        raise ValueError("method_paths must contain at least one generated H5AD")
    if expected_methods is not None:
        expected = set(expected_methods)
        supplied = set(method_paths)
        missing = sorted(expected.difference(supplied))
        extra = sorted(supplied.difference(expected))
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


def input_summary_frame(inputs: LoadedAssemblyInputs) -> pd.DataFrame:
    """Return compact native shape/path/hash evidence for a notebook report."""

    rows = []
    for name in [*inputs.methods, "baseline"]:
        adata = inputs.baseline if name == "baseline" else inputs.methods[name]
        rows.append(
            {
                "object": name,
                "path": str(inputs.paths[name]),
                "n_obs": int(adata.n_obs),
                "n_vars": int(adata.n_vars),
                "sha256": inputs.hashes[name],
            }
        )
    return pd.DataFrame(rows)


def _matrix_data(x: Any) -> np.ndarray:
    if sparse.issparse(x):
        return np.asarray(x.data)
    return np.asarray(x)


def _validate_ids(adata: AnnData, *, object_name: str) -> None:
    ids = pd.Index(adata.obs_names)
    if not ids.is_unique:
        raise ValueError(f"{object_name}: obs_names must be unique real identifiers")
    blank = np.asarray(ids.isna()) | (np.asarray(ids.astype(str).str.strip()) == "")
    if bool(np.any(blank)):
        raise ValueError(f"{object_name}: obs_names must not contain empty identifiers")


def validate_generated_expression(adata: AnnData, *, method: str) -> None:
    """Validate the caller-declared unlogged linear-expression boundary."""

    _validate_ids(adata, object_name=method)
    if not adata.var_names.is_unique:
        raise ValueError(f"{method}: var_names must be unique gene identifiers")
    values = _matrix_data(adata.X)
    if not np.issubdtype(values.dtype, np.number):
        raise TypeError(f"{method}: X must be numeric")
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{method}: X must contain only finite values")
    if np.any(values < 0):
        raise ValueError(f"{method}: X must be nonnegative unlogged linear expression")


def _normalize_label(value: Any) -> Any:
    """Normalize only slash spelling while leaving NA and IDs untouched."""

    if value is None:
        return value
    try:
        if bool(pd.isna(value)):
            return value
    except (TypeError, ValueError):
        pass
    return value.replace("/", "_") if isinstance(value, str) else value


def normalize_broad_labels(values: pd.Series) -> pd.Series:
    """Return a copy with slash replaced by underscore and NA values preserved."""

    result = values.astype(object).copy()
    return result.map(_normalize_label)


def _require_column(adata: AnnData, column: str, *, object_name: str) -> None:
    if column not in adata.obs:
        raise KeyError(f"{object_name}: required broad label column {column!r} not found")


def _select_cell_type(adata: AnnData, *, broad_column: str, broad_type: str) -> np.ndarray:
    values = normalize_broad_labels(adata.obs[broad_column])
    target = str(broad_type).replace("/", "_")
    return values.eq(target).fillna(False).to_numpy(dtype=bool)


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


def _as_coordinate_array(value: Any, *, object_name: str, spatial_key: str) -> np.ndarray:
    if sparse.issparse(value):
        value = value.toarray()
    coordinates = np.asarray(value)
    if coordinates.ndim != 2 or coordinates.shape[1] < 2:
        raise ValueError(
            f"{object_name}.obsm[{spatial_key!r}] must be a 2D array with at least two columns"
        )
    if not np.issubdtype(coordinates.dtype, np.number):
        raise TypeError(f"{object_name}.obsm[{spatial_key!r}] must be numeric")
    return coordinates


def _spatial_array(
    adata: AnnData, *, spatial_key: str, object_name: str
) -> np.ndarray:
    """Validate spatial structure without requiring unused rows to be finite."""

    if spatial_key not in adata.obsm:
        raise KeyError(f"{object_name}.obsm lacks spatial key {spatial_key!r}")
    coordinates = _as_coordinate_array(
        adata.obsm[spatial_key], object_name=object_name, spatial_key=spatial_key
    )
    if coordinates.shape[0] != adata.n_obs:
        raise ValueError(
            f"{object_name}.obsm[{spatial_key!r}] has {coordinates.shape[0]} rows; "
            f"expected {adata.n_obs}"
        )
    return coordinates


def _coordinates_for_ids(
    adata: AnnData,
    ids: pd.Index,
    *,
    spatial_key: str,
    object_name: str,
) -> np.ndarray:
    coordinates = _spatial_array(
        adata, spatial_key=spatial_key, object_name=object_name
    )
    positions = adata.obs_names.get_indexer(ids)
    if np.any(positions < 0):
        missing = ids[np.flatnonzero(positions < 0)].tolist()
        raise KeyError(f"{object_name}: missing spatial rows for IDs {missing[:5]}")
    selected = coordinates[positions]
    if not np.all(np.isfinite(selected)):
        raise ValueError(
            f"{object_name}.obsm[{spatial_key!r}] must be finite for compared IDs"
        )
    return selected


def _coordinate_check(
    method: AnnData,
    baseline: AnnData,
    ids: pd.Index,
    *,
    method_name: str,
    spatial_key: str,
) -> dict[str, Any]:
    """Compare one method's own coordinates with baseline coordinates by ID."""

    check: dict[str, Any] = {
        "status": "ok",
        "method": method_name,
        "spatial_key": spatial_key,
        "n_ids": int(len(ids)),
        "method_shape": None,
        "baseline_shape": None,
        "method_finite": False,
        "baseline_finite": False,
        "exact_equal": False,
        "difference_count": None,
        "max_abs_difference": None,
        "issue": None,
    }
    try:
        method_coordinates = _coordinates_for_ids(
            method,
            ids,
            spatial_key=spatial_key,
            object_name=method_name,
        )
        baseline_coordinates = _coordinates_for_ids(
            baseline,
            ids,
            spatial_key=spatial_key,
            object_name="baseline",
        )
    except (KeyError, TypeError, ValueError) as error:
        check.update(status="unavailable", issue=str(error))
        return check

    check["method_shape"] = tuple(method_coordinates.shape)
    check["baseline_shape"] = tuple(baseline_coordinates.shape)
    check["method_finite"] = bool(np.all(np.isfinite(method_coordinates)))
    check["baseline_finite"] = bool(np.all(np.isfinite(baseline_coordinates)))
    if method_coordinates.shape != baseline_coordinates.shape:
        check.update(
            status="unavailable",
            issue=(
                f"spatial shape differs from baseline: method={method_coordinates.shape}, "
                f"baseline={baseline_coordinates.shape}"
            ),
        )
        return check
    equal = np.array_equal(method_coordinates, baseline_coordinates)
    check["exact_equal"] = bool(equal)
    if not equal:
        difference = np.abs(method_coordinates - baseline_coordinates)
        check["difference_count"] = int(np.count_nonzero(difference))
        check["max_abs_difference"] = float(np.max(difference))
        check.update(
            status="unavailable",
            issue=(
                f"spatial coordinates differ from baseline for {check['difference_count']} "
                f"values (max_abs_difference={check['max_abs_difference']})"
            ),
        )
    return check


def _coverage_row(
    *,
    broad_type: str,
    object_name: str,
    original_ids: pd.Index,
    shared_ids: pd.Index,
    broad_column: str,
    original_genes: pd.Index | None,
    shared_genes: pd.Index | None,
    coordinate_check: Mapping[str, Any] | None = None,
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
    check = dict(coordinate_check or {})
    return {
        "broad_type": broad_type,
        "object": object_name,
        "broad_column": broad_column,
        "missing_type": len(original_ids) == 0,
        "original_type_ids": len(original_ids),
        "shared_ids": len(shared_ids),
        "excluded_ids": excluded_ids,
        "original_genes": original_gene_count,
        "shared_genes": shared_gene_count,
        "excluded_genes": excluded_genes,
        "spatial_key": check.get("spatial_key", pd.NA),
        "spatial_status": check.get("status", "not_checked"),
        "spatial_shape": check.get("method_shape", pd.NA),
        "baseline_spatial_shape": check.get("baseline_shape", pd.NA),
        "spatial_finite": check.get("method_finite", pd.NA),
        "baseline_spatial_finite": check.get("baseline_finite", pd.NA),
        "spatial_exact_equal": check.get("exact_equal", pd.NA),
        "spatial_difference_count": check.get("difference_count", pd.NA),
        "spatial_max_abs_difference": check.get("max_abs_difference", pd.NA),
        "spatial_issue": check.get("issue", pd.NA),
    }


def _prepare_type_scope(
    methods: Mapping[str, AnnData],
    baseline: AnnData,
    *,
    broad_type: str,
    broad_column: str,
    baseline_subtype_column: str,
    spatial_key: str,
) -> tuple[TypeComparison, list[dict[str, Any]]]:
    selected: dict[str, AnnData] = {}
    selected_ids: dict[str, pd.Index] = {}
    for method, adata in methods.items():
        subset = adata[
            _select_cell_type(adata, broad_column=broad_column, broad_type=broad_type), :
        ]
        selected[method] = subset
        selected_ids[method] = pd.Index(subset.obs_names)

    baseline_subset = baseline[
        _select_cell_type(baseline, broad_column=broad_column, broad_type=broad_type), :
    ]
    selected_ids["baseline"] = pd.Index(baseline_subset.obs_names)

    id_sets = [set(index) for index in selected_ids.values()]
    shared_ids = pd.Index(sorted(set.intersection(*id_sets))) if id_sets else pd.Index([])
    gene_sets = [set(adata.var_names) for adata in methods.values()]
    shared_genes = pd.Index(sorted(set.intersection(*gene_sets))) if gene_sets else pd.Index([])

    coordinate_checks: dict[str, dict[str, Any]] = {}
    if len(shared_ids):
        for method, adata in selected.items():
            coordinate_checks[method] = _coordinate_check(
                adata,
                baseline,
                shared_ids,
                method_name=method,
                spatial_key=spatial_key,
            )

    coverage = []
    for method, adata in methods.items():
        coverage.append(
            _coverage_row(
                broad_type=broad_type,
                object_name=method,
                original_ids=selected_ids[method],
                shared_ids=shared_ids,
                broad_column=broad_column,
                original_genes=adata.var_names,
                shared_genes=shared_genes,
                coordinate_check=coordinate_checks.get(method),
            )
        )

    baseline_check: dict[str, Any] = {}
    if len(shared_ids):
        try:
            baseline_coordinates = _coordinates_for_ids(
                baseline,
                shared_ids,
                spatial_key=spatial_key,
                object_name="baseline",
            )
            baseline_check = {
                "status": "ok",
                "spatial_key": spatial_key,
                "baseline_shape": tuple(baseline_coordinates.shape),
                "method_shape": tuple(baseline_coordinates.shape),
                "baseline_finite": bool(np.all(np.isfinite(baseline_coordinates))),
                "method_finite": bool(np.all(np.isfinite(baseline_coordinates))),
                "exact_equal": True,
            }
        except (KeyError, TypeError, ValueError) as error:
            baseline_check = {
                "status": "unavailable",
                "spatial_key": spatial_key,
                "issue": str(error),
            }
    coverage.append(
        _coverage_row(
            broad_type=broad_type,
            object_name="baseline",
            original_ids=selected_ids["baseline"],
            shared_ids=shared_ids,
            broad_column=broad_column,
            original_genes=None,
            shared_genes=None,
            coordinate_check=baseline_check,
        )
    )

    issues: list[str] = []
    missing_types = [name for name, ids in selected_ids.items() if len(ids) == 0]
    if missing_types:
        issues.append("broad type absent from: " + ", ".join(missing_types))
    if len(shared_ids) < 3:
        issues.append(f"only {len(shared_ids)} shared real IDs; at least 3 are required")
    if len(shared_genes) < 2:
        issues.append(f"only {len(shared_genes)} shared genes; at least 2 are required")
    for method, check in coordinate_checks.items():
        if check.get("status") != "ok":
            issues.append(f"{method}: {check.get('issue', 'spatial coordinate check failed')}")
    if baseline_check and baseline_check.get("status") != "ok":
        issues.append(
            f"baseline: {baseline_check.get('issue', 'spatial coordinate check failed')}"
        )

    if issues:
        baseline_labels = pd.Series(
            pd.array([], dtype="string"),
            index=pd.Index([], dtype="object"),
            name=baseline_subtype_column,
        )
        return (
            TypeComparison(
                broad_type=broad_type,
                status="unavailable",
                issues=issues,
                shared_ids=shared_ids,
                shared_genes=shared_genes,
                baseline_labels=baseline_labels,
                coordinate_checks=coordinate_checks,
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
            broad_type=broad_type,
            status="ready",
            issues=[],
            shared_ids=shared_ids,
            shared_genes=shared_genes,
            baseline_labels=baseline_labels,
            aligned=aligned,
            coordinate_checks=coordinate_checks,
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
            raise ValueError(f"{method}: clustering requires at least 3 IDs and 2 genes")

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
    """Evaluate labels on exact matched IDs with explicit missing handling."""

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


def _validate_cell_types(cell_types: Sequence[str]) -> tuple[str, ...]:
    if isinstance(cell_types, (str, bytes)):
        raise TypeError("cell_types must be a sequence of explicit type names")
    values = tuple(str(value).replace("/", "_") for value in cell_types)
    if not values:
        raise ValueError("cell_types must not be empty")
    if len(set(values)) != len(values):
        raise ValueError("cell_types must not contain duplicate canonical names")
    return values


def compare_assembly_methods(
    methods: Mapping[str, AnnData],
    baseline: AnnData,
    *,
    cell_types: Sequence[str] = DEFAULT_CELL_TYPES,
    broad_column: str = DEFAULT_BROAD_COLUMN,
    baseline_subtype_column: str = DEFAULT_BASELINE_SUBTYPE_COLUMN,
    spatial_key: str = DEFAULT_SPATIAL_KEY,
    resolutions: Sequence[float] = DEFAULT_RESOLUTIONS,
    seed: int = 42,
) -> AssemblyComparison:
    """Compare four assemblies for explicit broad cell types.

    cell_types is intentionally an exact-value list. Historical broad labels
    are normalized only by replacing slash with underscore; the original
    AnnData objects, observation IDs, and SVC_cluster values are untouched.
    """

    if not methods:
        raise ValueError("methods must not be empty")
    resolution_values = _validate_resolutions(resolutions)
    type_values = _validate_cell_types(cell_types)
    for method, adata in methods.items():
        validate_generated_expression(adata, method=method)
        _require_column(adata, broad_column, object_name=method)
    _validate_ids(baseline, object_name="baseline")
    _require_column(baseline, broad_column, object_name="baseline")
    if baseline_subtype_column not in baseline.obs:
        raise KeyError(
            f"baseline subtype column {baseline_subtype_column!r} not found in baseline.obs"
        )
    if spatial_key not in baseline.obsm:
        raise KeyError(f"baseline.obsm lacks spatial key {spatial_key!r}")
    _spatial_array(baseline, spatial_key=spatial_key, object_name="baseline")

    coverage_rows: list[dict[str, Any]] = []
    metric_frames: list[pd.DataFrame] = []
    by_type: dict[str, TypeComparison] = {}
    for broad_type in type_values:
        result, type_coverage = _prepare_type_scope(
            methods,
            baseline,
            broad_type=broad_type,
            broad_column=broad_column,
            baseline_subtype_column=baseline_subtype_column,
            spatial_key=spatial_key,
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


def _assert_plot_coordinates(
    result: TypeComparison,
    baseline: AnnData,
    *,
    method: str,
    spatial_key: str,
) -> np.ndarray:
    if result.status != "ok":
        raise ValueError(f"{result.broad_type}: comparison is {result.status}")
    if method not in result.clustered or method not in result.aligned:
        raise KeyError(f"unknown method {method!r}")
    check = _coordinate_check(
        result.aligned[method],
        baseline,
        result.shared_ids,
        method_name=method,
        spatial_key=spatial_key,
    )
    if check.get("status") != "ok":
        raise ValueError(f"{result.broad_type} | {method}: {check.get('issue')}")
    return _coordinates_for_ids(
        result.aligned[method],
        result.shared_ids,
        spatial_key=spatial_key,
        object_name=method,
    )


def spatial_plot_frame(
    result: TypeComparison,
    baseline: AnnData,
    *,
    method: str,
    resolution: float,
    spatial_key: str = DEFAULT_SPATIAL_KEY,
) -> pd.DataFrame:
    """Build baseline/new long-form plot data after the coordinate gate."""

    coordinates = _assert_plot_coordinates(
        result, baseline, method=method, spatial_key=spatial_key
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
        result.clustered[method].obs.loc[result.shared_ids, key]
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
    spatial_key: str = DEFAULT_SPATIAL_KEY,
    point_size: float = 8.0,
) -> tuple[plt.Figure, np.ndarray]:
    """Plot original subtype labels beside independent expression clusters."""

    frame = spatial_plot_frame(
        result,
        baseline,
        method=method,
        resolution=resolution,
        spatial_key=spatial_key,
    )
    panels = list(frame["panel"].drop_duplicates())
    figure, axes = plt.subplots(1, len(panels), figsize=(5 * len(panels), 4), squeeze=False)
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
    figure.suptitle(f"{result.broad_type}: baseline label versus new expression clustering")
    figure.tight_layout()
    return figure, axes
