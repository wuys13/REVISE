"""Deterministic long-table outputs for reconstruction-impact calculations.

This module deliberately only serializes already-computed analysis objects.  It
does not choose thresholds, relabel clusters, or recompute scientific metrics.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import gzip
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd

from revise.analysis.reconstruction_impact import (
    AnatomyRegionAnalysis,
    PartitionAnalysis,
    RawLevel2Mapping,
    SpatialImpactAnalysis,
)


_COMPARISON_COLUMNS = [
    "comparison_id",
    "scope_comparison_id",
    "sample_id",
    "task_cell_type",
    "scope",
    "comparison_edge",
    "raw_view",
    "reconstruction_view",
    "label_source",
    "observation_basis",
    "gene_rule",
    "normalization",
]


@dataclass(frozen=True)
class ImpactScope:
    """Already-computed objects belonging to one published comparison scope."""

    comparison_id: str
    sample_id: str
    task_cell_type: str | None
    scope: str
    comparison: Mapping[str, Any] = field(default_factory=dict)
    partition: PartitionAnalysis | None = None
    spatial: SpatialImpactAnalysis | None = None
    anatomy: AnatomyRegionAnalysis | None = None
    raw_level2: RawLevel2Mapping | None = None


def _scope_from_value(value: ImpactScope | Mapping[str, Any]) -> ImpactScope:
    if isinstance(value, ImpactScope):
        return value
    if not isinstance(value, Mapping):
        raise TypeError("Each impact scope must be an ImpactScope or mapping")
    return ImpactScope(**value)


def _json_safe(value: Any) -> Any:
    """Convert numpy/pandas values to JSON values without emitting NaN."""
    if value is None or value is pd.NA or value is pd.NaT:
        return None
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, pd.DataFrame):
        return _json_safe(value.to_dict(orient="records"))
    if isinstance(value, pd.Series):
        return _json_safe(value.to_dict())
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, (list, tuple, set)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if math.isfinite(number) else None
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (pd.Timestamp,)):
        return value.isoformat()
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    return value


def _write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(_json_safe(value), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _write_csv(path: Path, frame: pd.DataFrame, *, compressed: bool = False) -> None:
    """Write a stable CSV with explicit columns and reproducible gzip bytes."""
    path.parent.mkdir(parents=True, exist_ok=True)
    work = frame.copy()
    work.columns = [str(column) for column in work.columns]
    text = work.to_csv(
        index=False,
        na_rep="",
        float_format="%.17g",
        lineterminator="\n",
    )
    if not compressed:
        path.write_text(text, encoding="utf-8", newline="")
        return
    with path.open("wb") as handle:
        with gzip.GzipFile(fileobj=handle, mode="wb", mtime=0) as stream:
            stream.write(text.encode("utf-8"))


def _as_frame(value: Any) -> pd.DataFrame:
    if value is None:
        return pd.DataFrame()
    if isinstance(value, pd.DataFrame):
        return value.copy()
    return pd.DataFrame(value)


def _with_id(frame: Any, column: str) -> pd.DataFrame:
    work = _as_frame(frame)
    if column not in work.columns:
        work.insert(0, column, work.index.map(str))
    else:
        work[column] = work[column].map(lambda value: None if pd.isna(value) else str(value))
    return work.reset_index(drop=True)


def _context(scope: ImpactScope, comparison_id: str | None = None) -> dict[str, Any]:
    return {
        "comparison_id": str(comparison_id or scope.comparison_id),
        "scope_comparison_id": str(scope.comparison_id),
        "sample_id": str(scope.sample_id),
        "task_cell_type": scope.task_cell_type,
        "scope": str(scope.scope),
    }


def _add_context(
    frame: Any,
    scope: ImpactScope,
    *,
    edge: str | None = None,
    comparison_id: str | None = None,
    extra: Mapping[str, Any] | None = None,
) -> pd.DataFrame:
    work = _as_frame(frame)
    values = _context(scope, comparison_id=comparison_id)
    if edge is not None:
        values["comparison_edge"] = str(edge)
    if extra:
        values.update(extra)
    for column, value in values.items():
        work[column] = value
    ordered = list(values) + [column for column in work.columns if column not in values]
    return work.loc[:, ordered].reset_index(drop=True)


def _concat(frames: Iterable[pd.DataFrame], columns: list[str] | None = None) -> pd.DataFrame:
    material = [frame for frame in frames if frame is not None]
    if not material:
        return pd.DataFrame(columns=columns or [])
    known_columns: list[str] = []
    for frame in material:
        for column in frame.columns:
            if column not in known_columns:
                known_columns.append(column)
    if columns:
        known_columns = columns + [column for column in known_columns if column not in columns]
    material = [frame for frame in material if not frame.empty]
    if not material:
        return pd.DataFrame(columns=known_columns)
    # Avoid pandas' all-NA dtype warning while restoring those columns below.
    material = [
        frame.loc[:, [
            column for column in frame.columns
            if not frame[column].isna().all()
        ]]
        for frame in material
    ]
    result = pd.concat(material, ignore_index=True, sort=False)
    return result.reindex(columns=known_columns)


def _partition_edges(analysis: PartitionAnalysis) -> list[tuple[str, Any, str]]:
    """Return standard and complexity edges while preserving their source group."""
    entries: list[tuple[str, Any, str]] = []
    for group, comparisons in (
        ("comparison", analysis.comparisons),
        ("complexity", analysis.complexity_comparisons),
    ):
        for edge, comparison in comparisons.items():
            entries.append((str(edge), comparison, group))
    return entries


def _partition_comparison_id(scope: ImpactScope, edge: str) -> str:
    return f"{scope.comparison_id}::partition::{edge}"


def _contingency_long(
    contingency: Any,
    scope: ImpactScope,
    edge: str,
    comparison_id: str,
) -> pd.DataFrame:
    matrix = _as_frame(contingency)
    rows: list[dict[str, Any]] = []
    for raw_cluster, values in matrix.iterrows():
        numeric = pd.to_numeric(values, errors="coerce")
        row_total = float(numeric.sum(min_count=1)) if numeric.notna().any() else np.nan
        for recon_cluster, value in numeric.items():
            count = None if pd.isna(value) else float(value)
            if count is not None and count.is_integer():
                count = int(count)
            rows.append(
                {
                    **_context(scope, comparison_id=comparison_id),
                    "comparison_edge": edge,
                    "raw_cluster": str(raw_cluster),
                    "recon_cluster": str(recon_cluster),
                    "normalization": "absolute",
                    "value": count,
                    "count": count,
                    "denominator": None,
                }
            )
            normalized = None
            if count is not None and math.isfinite(row_total) and row_total != 0:
                normalized = count / row_total
            rows.append(
                {
                    **_context(scope, comparison_id=comparison_id),
                    "comparison_edge": edge,
                    "raw_cluster": str(raw_cluster),
                    "recon_cluster": str(recon_cluster),
                    "normalization": "raw_cluster_fraction",
                    "value": normalized,
                    "count": count,
                    "denominator": None if not math.isfinite(row_total) else row_total,
                }
            )
    return pd.DataFrame(rows)


def _write_partition(
    destination: Path,
    scopes: list[ImpactScope],
    register: Any,
) -> None:
    summaries: list[pd.DataFrame] = []
    assignments: list[pd.DataFrame] = []
    mappings: list[pd.DataFrame] = []
    contingencies: list[pd.DataFrame] = []
    resolution_sweeps: list[pd.DataFrame] = []
    complexity_sweeps: list[pd.DataFrame] = []
    level1_summaries: list[pd.DataFrame] = []
    audit_rows: list[dict[str, Any]] = []

    for scope in scopes:
        analysis = scope.partition
        if analysis is None:
            continue
        edges = _partition_edges(analysis)
        for edge, comparison, group in edges:
            edge_comparison_id = _partition_comparison_id(scope, edge)
            headline_eligible = (
                group == "comparison"
                and analysis.matched_cluster_status != "unmatched_cluster_complexity"
            )
            summary = _as_frame(comparison.summary)
            summaries.append(
                _add_context(
                    summary,
                    scope,
                    edge=edge,
                    comparison_id=edge_comparison_id,
                    extra={
                        "comparison_group": group,
                        "matched_cluster_status": analysis.matched_cluster_status,
                        "headline_eligible": headline_eligible,
                        "resolution": analysis.resolution if group == "comparison" else analysis.complexity_resolution,
                        "resolution_source": analysis.resolution_source,
                    },
                )
            )
            mappings.append(
                _add_context(
                    _as_frame(comparison.mapping),
                    scope,
                    edge=edge,
                    comparison_id=edge_comparison_id,
                    extra={"comparison_group": group},
                )
            )
            assignments.append(
                _add_context(
                    _with_id(comparison.assignments, "unit_id"),
                    scope,
                    edge=edge,
                    comparison_id=edge_comparison_id,
                    extra={"comparison_group": group},
                )
            )
            contingencies.append(
                _contingency_long(comparison.contingency, scope, edge, edge_comparison_id)
            )

        if not edges:
            continue
        main_edge = edges[0][0]
        sweep = _as_frame(analysis.sweep)
        if not sweep.empty:
            sweep["selected"] = np.isclose(
                pd.to_numeric(sweep.get("resolution"), errors="coerce"),
                float(analysis.resolution),
                equal_nan=False,
            )
        resolution_sweeps.append(
            _add_context(
                sweep,
                scope,
                edge=main_edge,
                comparison_id=_partition_comparison_id(scope, main_edge),
                extra={"selection_status": analysis.matched_cluster_status},
            )
        )
        complexity = _as_frame(analysis.complexity_sweep)
        if not complexity.empty:
            complexity["selected"] = np.isclose(
                pd.to_numeric(complexity.get("resolution"), errors="coerce"),
                float(analysis.complexity_resolution),
                equal_nan=False,
            )
        complexity_sweeps.append(
            _add_context(
                complexity,
                scope,
                edge=main_edge,
                comparison_id=_partition_comparison_id(scope, main_edge),
                extra={"selection_status": analysis.matched_cluster_status},
            )
        )
        level1_summaries.append(
            _add_context(
                analysis.change_by_level1,
                scope,
                edge=main_edge,
                comparison_id=_partition_comparison_id(scope, main_edge),
            )
        )
        audit_rows.append(
            {
                **_context(scope),
                "resolution": analysis.resolution,
                "resolution_source": analysis.resolution_source,
                "complexity_resolution": analysis.complexity_resolution,
                "matched_cluster_status": analysis.matched_cluster_status,
                "feature_names": list(analysis.feature_names),
                "audit": analysis.audit,
                "representation_audit": analysis.representation_audit,
                "edges": [{"edge": edge, "group": group} for edge, _, group in edges],
            }
        )

    if not audit_rows:
        return
    base = destination / "partition"
    files = [
        ("partition_summary", "summary.csv", _concat(summaries)),
        ("partition_assignments", "assignments.csv.gz", _concat(assignments), True),
        ("partition_mapping", "mapping.csv", _concat(mappings)),
        ("partition_contingency", "contingency.csv", _concat(contingencies)),
        ("partition_resolution_sweep", "resolution_sweep.csv", _concat(resolution_sweeps)),
        ("partition_complexity_sweep", "complexity_sweep.csv", _concat(complexity_sweeps)),
        ("partition_level1_summary", "level1_summary.csv", _concat(level1_summaries)),
    ]
    for item in files:
        role, name, frame, *compressed = item
        path = base / name
        _write_csv(path, frame, compressed=bool(compressed))
        register(role, path, f"Reconstruction-impact partition {name} long table")
    audit_path = base / "audit.json"
    _write_json(audit_path, {"scopes": audit_rows})
    register("partition_audit", audit_path, "Partition comparison and resolution audit")


def _window_metric_long(
    frame: Any,
    scope: ImpactScope,
    scale: Any,
    n_draws: Any,
) -> pd.DataFrame:
    source = _as_frame(frame)
    if source.empty:
        return pd.DataFrame()
    if {"metric", "baseline", "delta"} <= set(source.columns):
        return _add_context(source, scope)
    rows: list[dict[str, Any]] = []
    for metric in ("k_obs", "entropy", "neff", "evenness"):
        recon_col = f"{metric}_recon"
        if recon_col not in source.columns:
            continue
        for baseline, raw_col, delta_col, raw_sd_col, delta_sd_col in (
            (
                "raw_leiden",
                f"{metric}_raw",
                f"delta_{metric}_vs_raw_leiden",
                f"{metric}_sd_raw",
                f"delta_{metric}_vs_raw_leiden_sd",
            ),
            (
                "raw_level2",
                f"{metric}_level2",
                f"delta_{metric}_vs_raw_level2",
                f"{metric}_sd_level2",
                f"delta_{metric}_vs_raw_level2_sd",
            ),
        ):
            if raw_col not in source.columns and delta_col not in source.columns:
                continue
            for _, row in source.iterrows():
                valid = row.get("valid_window")
                if valid is None or pd.isna(valid):
                    support_status = None
                elif bool(valid):
                    support_status = "computed"
                else:
                    support_status = "insufficient_support"
                raw_value = row.get(raw_col)
                delta_value = row.get(delta_col)
                baseline_present = not all(
                    value is None or pd.isna(value)
                    for value in (raw_value, delta_value)
                )
                status = support_status
                if support_status == "computed" and not baseline_present:
                    status = "unmeasured"
                rows.append(
                    {
                        **_context(scope),
                        "window_id": row.get("window_id"),
                        "window_x": row.get("window_x"),
                        "window_y": row.get("window_y"),
                        "window_x_index": row.get("window_x_index"),
                        "window_y_index": row.get("window_y_index"),
                        "level1_region": row.get("level1_region"),
                        "n_units": row.get("n_units"),
                        "valid_window": valid,
                        "status": status,
                        "baseline_status": status,
                        "scale": row.get("scale", scale),
                        "baseline": baseline,
                        "metric": metric,
                        "raw": raw_value,
                        "reconstruction": row.get(recon_col),
                        "delta": delta_value,
                        "raw_sd": row.get(raw_sd_col),
                        "reconstruction_sd": row.get(f"{metric}_sd_recon"),
                        "delta_sd": row.get(delta_sd_col),
                        "n_draws": row.get("n_draws", n_draws),
                        "n_changed_units": row.get("n_changed_units"),
                        "unit_change_fraction": row.get("unit_change_fraction"),
                        "in_state_region": row.get("in_state_region"),
                        "in_gain_region": row.get("in_gain_region"),
                    }
                )
    return pd.DataFrame(rows)


def _diversity_long(frame: Any, scope: ImpactScope) -> pd.DataFrame:
    source = _as_frame(frame)
    if source.empty:
        return pd.DataFrame()
    rows: list[dict[str, Any]] = []
    for metric in ("k_obs", "entropy", "neff", "evenness"):
        for baseline, raw_suffix, delta_suffix in (
            ("raw_leiden", "raw", "delta"),
            ("raw_level2", "level2", "delta_vs_raw_level2"),
        ):
            prefix = {
                "median": "median",
                "q1": "q1",
                "q3": "q3",
            }
            value_columns = {}
            for statistic, statistic_prefix in prefix.items():
                value_columns[f"raw_{statistic}"] = f"{statistic_prefix}_{metric}_{raw_suffix}"
                value_columns[f"reconstruction_{statistic}"] = f"{statistic_prefix}_{metric}_recon"
                value_columns[f"delta_{statistic}"] = (
                    f"{statistic_prefix}_delta_{metric}"
                    if delta_suffix == "delta"
                    else f"{statistic_prefix}_delta_{metric}_vs_raw_level2"
                )
            if not any(column in source.columns for column in value_columns.values()):
                continue
            for _, row in source.iterrows():
                payload = row.to_dict()
                payload.update(
                    {
                        "baseline": baseline,
                        "metric": metric,
                        "raw": row.get(value_columns["raw_median"]),
                        "reconstruction": row.get(value_columns["reconstruction_median"]),
                        "delta": row.get(value_columns["delta_median"]),
                    }
                )
                for output_column, source_column in value_columns.items():
                    payload[output_column] = row.get(source_column)
                rows.append(payload)
    return _add_context(pd.DataFrame(rows), scope)


def _scale_sensitivity_long(frame: Any, scope: ImpactScope) -> pd.DataFrame:
    source = _as_frame(frame)
    if source.empty:
        return pd.DataFrame()
    if {"baseline", "metric"} <= set(source.columns):
        result = _add_context(source, scope)
        result["sensitivity_type"] = "scale"
        return result
    rows: list[dict[str, Any]] = []
    for metric in ("k_obs", "entropy", "neff", "evenness"):
        for baseline, suffix in (
            ("raw_leiden", "vs_raw_leiden"),
            ("raw_level2", "vs_raw_level2"),
        ):
            column = f"median_delta_{metric}_{suffix}"
            if column not in source.columns:
                continue
            for _, row in source.iterrows():
                payload = row.to_dict()
                payload.update(
                    {
                        "sensitivity_type": "scale",
                        "baseline": baseline,
                        "metric": metric,
                        "delta": row.get(column),
                    }
                )
                rows.append(payload)
    return _add_context(pd.DataFrame(rows), scope)


def _bootstrap_frame(frame: Any, threshold: Mapping[str, Any], region_type: str, scope: ImpactScope) -> pd.DataFrame:
    source = _as_frame(frame)
    if source.empty:
        source = pd.DataFrame(
            [{"bootstrap_draw": None, "bootstrap_threshold": None, "status": threshold.get("status")}]
        )
    elif "bootstrap_draw" not in source.columns:
        source.insert(0, "bootstrap_draw", np.arange(len(source), dtype=int))
    source["region_type"] = region_type
    source["threshold_status"] = threshold.get("status")
    return _add_context(source, scope)


def _write_spatial(destination: Path, scopes: list[ImpactScope], register: Any) -> None:
    metrics: list[pd.DataFrame] = []
    assignments: list[pd.DataFrame] = []
    diversity: list[pd.DataFrame] = []
    cluster_change: list[pd.DataFrame] = []
    extent: list[pd.DataFrame] = []
    sensitivity: list[pd.DataFrame] = []
    bootstrap: list[pd.DataFrame] = []
    audit_rows: list[dict[str, Any]] = []

    for scope in scopes:
        analysis = scope.spatial
        if analysis is None:
            continue
        scale = analysis.scale_audit.get("main_window_side_um")
        n_draws = analysis.scale_audit.get("rarefaction_draws")
        metrics.append(_window_metric_long(analysis.window_metrics, scope, scale, n_draws))
        assignments.append(_add_context(_with_id(analysis.unit_assignments, "unit_id"), scope, extra={"scale": scale}))
        diversity.append(_diversity_long(analysis.diversity_by_anatomy, scope))
        cluster_change.append(_add_context(analysis.cluster_change_by_anatomy, scope))
        state_extent = _as_frame(analysis.region_extent_by_anatomy)
        gain_extent = _as_frame(analysis.gain_region_extent_by_anatomy)
        state_extent["region_type"] = "state"
        gain_extent["region_type"] = "gain"
        for frame in (state_extent, gain_extent):
            if "baseline" not in frame.columns:
                frame["baseline"] = "raw_leiden"
            extent.append(_add_context(frame, scope, extra={"scale": scale}))

        support = _add_context(analysis.support_sensitivity, scope, extra={"sensitivity_type": "support"})
        scale_frame = _scale_sensitivity_long(analysis.scale_sensitivity, scope)
        threshold = _add_context(analysis.threshold_sensitivity, scope, extra={"sensitivity_type": "threshold"})
        sensitivity.extend([support, scale_frame, threshold])
        bootstrap.append(_bootstrap_frame(analysis.state_threshold_bootstrap, analysis.state_threshold, "state", scope))
        bootstrap.append(_bootstrap_frame(analysis.gain_threshold_bootstrap, analysis.gain_threshold, "gain", scope))
        audit_rows.append(
            {
                **_context(scope),
                "scale_audit": analysis.scale_audit,
                "support_selection": analysis.support_selection,
                "support_sensitivity": analysis.support_sensitivity,
                "scale_sensitivity": analysis.scale_sensitivity,
                "threshold_sensitivity": analysis.threshold_sensitivity,
                "state_threshold": analysis.state_threshold,
                "gain_threshold": analysis.gain_threshold,
            }
        )

    if not audit_rows:
        return
    base = destination / "spatial"
    files = [
        ("spatial_window_metrics", "window_metrics.csv", _concat(metrics)),
        ("spatial_unit_window_assignments", "unit_window_assignments.csv.gz", _concat(assignments), True),
        ("spatial_diversity_by_anatomy", "diversity_by_anatomy.csv", _concat(diversity)),
        ("spatial_cluster_change_by_anatomy", "cluster_change_by_anatomy.csv", _concat(cluster_change)),
        ("spatial_region_extent_by_anatomy", "region_extent_by_anatomy.csv", _concat(extent)),
        ("spatial_sensitivity", "sensitivity.csv", _concat(sensitivity)),
        ("spatial_threshold_bootstrap", "threshold_bootstrap.csv", _concat(bootstrap)),
    ]
    for item in files:
        role, name, frame, *compressed = item
        path = base / name
        _write_csv(path, frame, compressed=bool(compressed))
        register(role, path, f"Reconstruction-impact spatial {name} long table")
    audit_path = base / "audit.json"
    audit_payload: dict[str, Any] = {"scopes": audit_rows}
    if len(audit_rows) == 1:
        audit_payload.update(
            {
                key: value
                for key, value in audit_rows[0].items()
                if key not in {"comparison_id", "sample_id", "task_cell_type", "scope"}
            }
        )
    _write_json(audit_path, audit_payload)
    register("spatial_audit", audit_path, "Spatial scale, support, threshold and status audit")


def _write_anatomy(destination: Path, scopes: list[ImpactScope], register: Any) -> None:
    assignments: list[pd.DataFrame] = []
    maps: list[pd.DataFrame] = []
    summaries: list[pd.DataFrame] = []
    audits: list[dict[str, Any]] = []
    for scope in scopes:
        analysis = scope.anatomy
        if analysis is None:
            continue
        assignments.append(_add_context(_with_id(analysis.anatomy_unit_assignments, "unit_id"), scope))
        maps.append(_add_context(analysis.anatomy_windows, scope))
        summaries.append(_add_context(analysis.anatomy_context_summary, scope))
        audits.append(
            {
                **_context(scope),
                "scale_audit": analysis.scale_audit,
                "support_selection": analysis.support_selection,
                "support_sensitivity": analysis.support_sensitivity,
            }
        )
    if not audits:
        return
    base = destination / "anatomy"
    files = [
        ("anatomy_assignments", "assignments.csv.gz", _concat(assignments), True),
        ("anatomy_window_map", "window_map.csv", _concat(maps)),
        ("anatomy_summary", "summary.csv", _concat(summaries)),
    ]
    for item in files:
        role, name, frame, *compressed = item
        path = base / name
        _write_csv(path, frame, compressed=bool(compressed))
        register(role, path, f"Reconstruction-impact anatomy {name} long table")
    audit_path = base / "audit.json"
    _write_json(audit_path, {"scopes": audits})
    register("anatomy_audit", audit_path, "Anatomy scale and support audit")


def _posterior_long(posterior: Any, scope: ImpactScope) -> pd.DataFrame:
    source = _with_id(posterior, "unit_id")
    if source.empty:
        return pd.DataFrame(columns=[*list(_context(scope)), "unit_id", "raw_level2", "posterior"])
    result = source.melt(id_vars=["unit_id"], var_name="raw_level2", value_name="posterior")
    return _add_context(result, scope)


def _write_raw_level2(destination: Path, scopes: list[ImpactScope], register: Any) -> None:
    assignments: list[pd.DataFrame] = []
    posterior: list[pd.DataFrame] = []
    audits: list[dict[str, Any]] = []
    for scope in scopes:
        mapping = scope.raw_level2
        if mapping is None:
            continue
        assignments.append(_add_context(_with_id(mapping.assignments, "unit_id"), scope))
        posterior.append(_posterior_long(mapping.posterior, scope))
        audits.append({**_context(scope), "audit": mapping.audit})
    if not audits:
        return
    base = destination / "raw_level2"
    files = [
        ("raw_level2_assignments", "assignments.csv.gz", _concat(assignments), True),
        ("raw_level2_posterior", "posterior.csv.gz", _concat(posterior), True),
    ]
    for item in files:
        role, name, frame, *compressed = item
        path = base / name
        _write_csv(path, frame, compressed=bool(compressed))
        register(role, path, f"Raw Level2 {name} long table")
    audit_path = base / "audit.json"
    _write_json(audit_path, {"scopes": audits})
    register("raw_level2_audit", audit_path, "Raw Level2 mapping and posterior audit")


def write_impact_outputs(
    output_dir: str | Path,
    scopes: Iterable[ImpactScope | Mapping[str, Any]],
    *,
    audit: Mapping[str, Any] | None = None,
    calculation: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Write all available impact aspects and return a batch-compatible manifest."""
    resolved_scopes = [_scope_from_value(scope) for scope in scopes]
    identifiers = [str(scope.comparison_id) for scope in resolved_scopes]
    if any(not identifier.strip() for identifier in identifiers):
        raise ValueError("Impact scopes require a nonempty comparison_id")
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("Impact scopes require unique comparison_id values")

    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    artifacts: dict[str, dict[str, str]] = {}

    def register(role: str, path: Path, description: str) -> None:
        relative = path.relative_to(destination).as_posix()
        if role in artifacts:
            raise ValueError(f"Duplicate impact artifact role: {role}")
        artifacts[role] = {"path": relative, "description": description}

    comparison_rows: list[dict[str, Any]] = []
    comparison_extra: set[str] = set()
    for scope in resolved_scopes:
        metadata = dict(scope.comparison)
        metadata.pop("comparison_id", None)
        base_row = {
            **_context(scope),
            **metadata,
            "comparison_id": str(scope.comparison_id),
            "comparison_edge": "spatial_context",
        }
        comparison_rows.append(base_row)
        if scope.partition is not None:
            for edge, _, _ in _partition_edges(scope.partition):
                edge_id = _partition_comparison_id(scope, edge)
                comparison_rows.append(
                    {
                        **_context(scope, comparison_id=edge_id),
                        **metadata,
                        "comparison_id": edge_id,
                        "comparison_edge": edge,
                    }
                )
        comparison_extra.update(metadata)
    comparison_columns = _COMPARISON_COLUMNS + [
        column for column in sorted(comparison_extra) if column not in _COMPARISON_COLUMNS
    ]
    comparisons = pd.DataFrame(comparison_rows).reindex(columns=comparison_columns)
    comparisons_path = destination / "comparisons.csv"
    _write_csv(comparisons_path, comparisons)
    register("comparisons", comparisons_path, "Published reconstruction-impact comparison definitions")

    _write_partition(destination, resolved_scopes, register)
    _write_spatial(destination, resolved_scopes, register)
    _write_anatomy(destination, resolved_scopes, register)
    _write_raw_level2(destination, resolved_scopes, register)

    root_audit = dict(audit or {})
    root_audit.setdefault("schema_version", 1)
    root_audit["comparison_ids"] = identifiers
    root_audit["n_scopes"] = len(resolved_scopes)
    audit_path = destination / "audit.json"
    _write_json(audit_path, root_audit)
    register("audit", audit_path, "Reconstruction-impact root provenance and scope audit")

    resolved_calculation = {
        "input_view": "reconstruction_impact",
        "parameters": {"n_scopes": len(resolved_scopes)},
        "comparison_basis": "paired_reconstruction_impact",
    }
    if calculation is not None:
        resolved_calculation.update(dict(calculation))
    if not isinstance(resolved_calculation.get("input_view"), str) or not resolved_calculation["input_view"].strip():
        raise ValueError("Impact calculation.input_view must be a nonempty string")
    if not isinstance(resolved_calculation.get("parameters"), Mapping):
        raise ValueError("Impact calculation.parameters must be a mapping")
    if not isinstance(resolved_calculation.get("comparison_basis"), str) or not resolved_calculation["comparison_basis"].strip():
        raise ValueError("Impact calculation.comparison_basis must be a nonempty string")
    resolved_calculation["parameters"] = dict(resolved_calculation["parameters"])
    return {"artifacts": artifacts, "calculation": resolved_calculation}
