"""Configuration and artifact helpers for reconstruction-impact analyses."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import gc
from pathlib import Path
from typing import Any, Mapping

import pandas as pd
import yaml
from anndata import AnnData
from sklearn.metrics import adjusted_rand_score

from revise.analysis.basic.partition_change import (
    PartitionComparison,
    align_observation_pairs,
    compare_partitions,
    leiden_labels,
    prepare_leiden_graph,
    select_level1_resolution,
    select_shared_feature_names,
)
from revise.analysis.basic.spatial_region import (
    assign_anatomy_candidates,
    assign_square_windows,
    compute_window_diversity,
    convert_coordinates_to_microns,
    select_min_parent_support,
    summarize_region,
    summarize_region_by_anatomy,
)


DEFAULTS: dict[str, Any] = {
    "partition_change": {
        "level1_resolution_candidates": [0.3, 0.5, 0.8],
        "within_level1_resolution": 0.5,
        "random_state": 42,
        "comparisons": [],
    },
    "spatial_region": {
        "coordinate_unit": "pixel",
        "cell_equivalent_um": 8.0,
        "main_window_multiplier": 5,
        "window_multipliers": [1, 2, 3, 5, 7, 10],
        "neff_threshold": 2.0,
        "neff_thresholds": [1.5, 2.0, 2.5, 3.0, 3.5],
        "min_units_per_window": 1,
        "anatomy_region": {
            "tumor_label": "Tumor",
            "normal_source_label": "Intestinal Epithelial",
        },
    },
}


@dataclass(frozen=True)
class PartitionAnalysis:
    """Resolved partition nodes and directed comparison edges for one scope."""

    resolution: float
    resolution_source: str
    sweep: pd.DataFrame
    audit: dict[str, int | bool]
    comparisons: dict[str, PartitionComparison]
    representation_audit: dict[str, int | bool]
    feature_names: list[str]


@dataclass(frozen=True)
class SpatialImpactAnalysis:
    """Spatial tables resolved for one route at the declared physical scale."""

    scale_audit: dict[str, float | int | str | None]
    support_selection: dict[str, int | str | None]
    support_sensitivity: pd.DataFrame
    unit_assignments: pd.DataFrame
    anatomy_unit_assignments: pd.DataFrame
    window_metrics: pd.DataFrame
    anatomy_windows: pd.DataFrame
    anatomy_summary: pd.DataFrame
    scale_sensitivity: pd.DataFrame
    threshold_sensitivity: pd.DataFrame


def run_partition_analysis(
    raw: AnnData,
    reconstructed: AnnData,
    *,
    level1_col: str,
    final_cluster_key: str | None = None,
    route_kind: str = "sp_svc",
    resolution_mode: str | None = None,
    resolution_candidates: list[float] | None = None,
    within_level1_resolution: float = 0.5,
    random_state: int = 42,
    n_top_genes: int = 2000,
    feature_names: list[str] | None = None,
) -> PartitionAnalysis:
    """Build the one scientifically applicable comparison edge for a route.

    A multi-Level1 scope calibrates one resolution on Raw's Level1 ARI.  A
    single-Level1 scope uses the fixed internal resolution.  Both expression
    nodes share the same observation IDs and feature names before their own
    PCA/neighbor graphs are constructed.
    """
    if route_kind not in {"sp_svc", "sc_svc"}:
        raise ValueError("route_kind must be sp_svc or sc_svc")
    if level1_col not in raw.obs:
        raise KeyError(f"Raw input is missing {level1_col!r}")
    raw_work, recon_work, audit = align_observation_pairs(raw, reconstructed)
    if level1_col not in recon_work.obs:
        raise KeyError(f"Reconstructed input is missing {level1_col!r}")
    if feature_names is None:
        features = select_shared_feature_names(raw_work, recon_work, n_top_genes=n_top_genes)
    else:
        features = [str(feature) for feature in feature_names]
        if not features or not set(features) <= set(raw_work.var_names):
            raise ValueError("feature_names must be a non-empty shared ordered gene list")
    raw_graph = prepare_leiden_graph(raw_work, feature_names=features)
    level1 = raw_work.obs[level1_col].astype(str)
    candidates = resolution_candidates or [0.3, 0.5, 0.8]

    use_level1_ari = (
        resolution_mode == "level1_ari"
        or (resolution_mode is None and level1.nunique() > 1)
    )
    if resolution_mode not in {None, "level1_ari", "fixed_within_level1"}:
        raise ValueError("resolution_mode must be level1_ari or fixed_within_level1")
    if use_level1_ari:
        sweep_rows: list[dict[str, float]] = []
        labels_by_resolution: dict[float, pd.Series] = {}
        for resolution in candidates:
            labels = leiden_labels(raw_graph, resolution=resolution, random_state=random_state)
            labels_by_resolution[float(resolution)] = labels
            sweep_rows.append(
                {
                    "resolution": float(resolution),
                    "ARI": float(adjusted_rand_score(level1, labels)),
                }
            )
        sweep = pd.DataFrame(sweep_rows)
        selected = select_level1_resolution(sweep)
        resolution = float(selected["resolution"])
        raw_labels = labels_by_resolution[resolution]
        resolution_source = "raw_level1_ari"
    else:
        resolution = float(within_level1_resolution)
        raw_labels = leiden_labels(raw_graph, resolution=resolution, random_state=random_state)
        sweep = pd.DataFrame(columns=["resolution", "ARI"])
        resolution_source = "fixed_within_level1"

    representation_audit: dict[str, int | bool] = {}
    if route_kind == "sp_svc":
        recon_graph = prepare_leiden_graph(recon_work, feature_names=features)
        recon_expression_labels = leiden_labels(
            recon_graph, resolution=resolution, random_state=random_state
        )
        comparisons = {
            "raw_to_recon_expression": compare_partitions(
                raw_labels,
                recon_expression_labels,
                comparison_edge="raw_to_recon_expression",
            )
        }
    else:
        if final_cluster_key is None:
            raise ValueError("sc_svc requires final_cluster_key")
        if final_cluster_key not in recon_work.obs:
            raise KeyError(f"Reconstructed input is missing {final_cluster_key!r}")
        final_labels_source = recon_work.obs[final_cluster_key]
        if final_labels_source.isna().any():
            raise ValueError(f"{final_cluster_key!r} must not contain missing labels")
        final_labels = final_labels_source.astype(str)
        differences = raw_work[:, features].X != recon_work[:, features].X
        difference_nnz = int(differences.nnz) if hasattr(differences, "nnz") else int(differences.sum())
        representation_audit = {
            "spatial_expression_identical": difference_nnz == 0,
            "expression_difference_nnz": difference_nnz,
        }
        comparisons = {
            "raw_to_final_svc": compare_partitions(
                raw_labels,
                final_labels,
                comparison_edge="raw_to_final_svc",
            )
        }
    return PartitionAnalysis(
        resolution, resolution_source, sweep, audit, comparisons, representation_audit, features
    )


def compute_spatial_impact(
    *,
    full_coordinates: pd.DataFrame,
    full_level1_labels: pd.Series,
    paired_coordinates: pd.DataFrame,
    raw_labels: pd.Series,
    reconstructed_labels: pd.Series,
    unit_changed: pd.Series,
    microns_per_coordinate: float,
    cell_equivalent_um: float = 8.0,
    main_window_multiplier: int = 5,
    window_multipliers: list[int] | None = None,
    neff_threshold: float = 2.0,
    neff_thresholds: list[float] | None = None,
    tumor_label: str = "Tumor",
    normal_source_label: str = "Intestinal Epithelial",
) -> SpatialImpactAnalysis:
    """Compute fixed-anatomy context and diversity changes on physical windows."""
    if cell_equivalent_um <= 0 or main_window_multiplier <= 0:
        raise ValueError("cell_equivalent_um and main_window_multiplier must be positive")
    paired_ids = paired_coordinates.index
    for name, labels in (("Raw", raw_labels), ("Reconstructed", reconstructed_labels), ("unit_changed", unit_changed)):
        if not labels.index.is_unique or set(labels.index) != set(paired_ids):
            raise ValueError(f"{name} labels must contain exactly the paired coordinate IDs")
    if not full_level1_labels.index.is_unique or set(full_level1_labels.index) != set(full_coordinates.index):
        raise ValueError("Full Level1 labels must contain exactly the full coordinate IDs")

    full_um = convert_coordinates_to_microns(full_coordinates, microns_per_coordinate=microns_per_coordinate)
    paired_um = convert_coordinates_to_microns(paired_coordinates, microns_per_coordinate=microns_per_coordinate)
    origin = (float(full_um["x"].min()), float(full_um["y"].min()))
    main_side = float(cell_equivalent_um * main_window_multiplier)
    anatomy_assignments = assign_square_windows(full_um, window_side_length=main_side, origin=origin)
    anatomy_windows = assign_anatomy_candidates(
        anatomy_assignments,
        full_level1_labels,
        tumor_label=tumor_label,
        normal_source_label=normal_source_label,
    )
    paired_windows = assign_square_windows(paired_um, window_side_length=main_side, origin=origin)
    support_selection, support_sensitivity = select_min_parent_support(paired_windows)
    minimum = support_selection["min_parent_units"]
    if minimum is None:
        minimum = int(paired_windows.shape[0]) + 1
    metrics = compute_window_diversity(
        paired_windows,
        raw_labels,
        reconstructed_labels,
        min_units_per_window=int(minimum),
    )
    changed = unit_changed.astype(bool).reindex(paired_windows.index)
    window_change = (
        pd.DataFrame({"window_id": paired_windows["window_id"], "unit_changed": changed})
        .groupby("window_id", sort=True)["unit_changed"]
        .agg(n_changed_units="sum", unit_change_fraction="mean")
        .reset_index()
    )
    metrics = metrics.merge(window_change, on="window_id", how="left", validate="one_to_one")
    anatomy_context = anatomy_windows.drop(columns=["window_x", "window_y"])
    metrics = metrics.merge(anatomy_context, on="window_id", how="left", validate="one_to_one")
    if metrics["level1_region"].isna().any():
        raise ValueError("Every paired window must map to full-cohort anatomy")
    summary, metrics = summarize_region(metrics, threshold=neff_threshold, window_side_length=main_side)
    anatomy_summary = summarize_region_by_anatomy(metrics, anatomy_windows, window_side_length=main_side)
    unit_assignments = paired_windows.loc[:, ["x", "y", "window_id"]].copy()
    unit_assignments["raw_cluster"] = raw_labels.reindex(unit_assignments.index).astype(str)
    unit_assignments["reconstructed_cluster"] = reconstructed_labels.reindex(unit_assignments.index).astype(str)
    unit_assignments["unit_changed"] = changed
    unit_assignments = unit_assignments.merge(
        anatomy_windows.loc[:, ["window_id", "level1_region"]], on="window_id", how="left", validate="many_to_one"
    )
    unit_assignments.index = paired_windows.index

    threshold_rows = []
    for threshold in neff_thresholds or [1.5, 2.0, 2.5, 3.0, 3.5]:
        threshold_summary, _ = summarize_region(metrics, threshold=float(threshold), window_side_length=main_side)
        threshold_rows.append(threshold_summary)
    threshold_sensitivity = pd.concat(threshold_rows, ignore_index=True)

    scale_rows = []
    for multiplier in window_multipliers or [1, 2, 3, 5, 7, 10]:
        side = float(cell_equivalent_um * multiplier)
        assignments = assign_square_windows(paired_um, window_side_length=side, origin=origin)
        scale_metrics = compute_window_diversity(
            assignments,
            raw_labels,
            reconstructed_labels,
            min_units_per_window=int(minimum),
        )
        scale_summary, _ = summarize_region(scale_metrics, threshold=neff_threshold, window_side_length=side)
        scale_summary["window_multiplier"] = int(multiplier)
        scale_rows.append(scale_summary)
        del assignments, scale_metrics, scale_summary
        gc.collect()
    scale_sensitivity = pd.concat(scale_rows, ignore_index=True)
    scale_audit = {
        "coordinate_unit": "coordinate",
        "microns_per_coordinate": float(microns_per_coordinate),
        "cell_equivalent_um": float(cell_equivalent_um),
        "main_window_multiplier": int(main_window_multiplier),
        "main_window_side_um": main_side,
        "main_window_area_um2": main_side**2,
        "origin_x_um": origin[0],
        "origin_y_um": origin[1],
    }
    return SpatialImpactAnalysis(
        scale_audit=scale_audit,
        support_selection=support_selection,
        support_sensitivity=support_sensitivity,
        unit_assignments=unit_assignments,
        anatomy_unit_assignments=anatomy_assignments,
        window_metrics=metrics,
        anatomy_windows=anatomy_windows,
        anatomy_summary=anatomy_summary,
        scale_sensitivity=scale_sensitivity,
        threshold_sensitivity=threshold_sensitivity,
    )


def _merge_defaults(config: dict[str, Any], defaults: Mapping[str, Any]) -> dict[str, Any]:
    merged = deepcopy(config)
    for key, default in defaults.items():
        if key not in merged:
            merged[key] = deepcopy(default)
        elif isinstance(default, Mapping) and isinstance(merged[key], Mapping):
            merged[key] = _merge_defaults(dict(merged[key]), default)
    return merged


def file_sha256(path: str | Path, *, block_size: int = 1024 * 1024) -> str:
    """Return a streaming SHA-256 digest for an analysis input artifact."""
    import hashlib

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def load_reconstruction_impact_config(path: str | Path) -> dict[str, Any]:
    """Load post-reconstruction analysis config and reject invalid carrier roles."""
    config_path = Path(path)
    with config_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    if not isinstance(config, dict):
        raise ValueError("Reconstruction-impact config must be a mapping")
    if config.get("schema_version") != 1:
        raise ValueError("Reconstruction-impact config requires schema_version: 1")
    sample = config.get("sample")
    if not isinstance(sample, Mapping) or not sample.get("id"):
        raise ValueError("Reconstruction-impact config requires sample.id")
    config = _merge_defaults(config, DEFAULTS)
    comparisons = config["partition_change"].get("comparisons")
    if not isinstance(comparisons, list):
        raise ValueError("partition_change.comparisons must be a list")
    for comparison in comparisons:
        if not isinstance(comparison, Mapping):
            raise ValueError("Each partition comparison must be a mapping")
        for key, value in comparison.items():
            if "reconstructed" in str(key) and isinstance(value, str) and value.endswith("/expr.h5ad"):
                raise ValueError(
                    "The expression-side expr.h5ad is not observation-paired and cannot be a spatial comparison input"
                )
    output = config.setdefault("output", {})
    output.setdefault("dir", f"output/reconstruction_impact/{sample['id']}")
    return config


def write_analysis_artifacts(
    output_dir: str | Path,
    *,
    config: Mapping[str, Any],
    manifest: Mapping[str, Any],
    input_audit: pd.DataFrame,
) -> Path:
    """Persist the common resolved-config, manifest and input-audit contract."""
    import json

    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    with (destination / "resolved_config.yaml").open("w", encoding="utf-8") as handle:
        yaml.safe_dump(dict(config), handle, sort_keys=False)
    (destination / "manifest.json").write_text(
        json.dumps(dict(manifest), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    input_audit.to_csv(destination / "input_audit.csv", index=False)
    return destination


def write_partition_artifacts(output_dir: str | Path, analysis: PartitionAnalysis) -> Path:
    """Write one route's comparison-edge tables without notebook-local schema logic."""
    destination = Path(output_dir) / "st_unit"
    destination.mkdir(parents=True, exist_ok=True)
    summaries = []
    for edge, comparison in analysis.comparisons.items():
        summary = comparison.summary.copy()
        summary["resolution"] = analysis.resolution
        summary["resolution_source"] = analysis.resolution_source
        summaries.append(summary)
        comparison.mapping.to_csv(destination / f"{edge}_mapping.csv", index=False)
        comparison.contingency.to_csv(destination / f"{edge}_contingency_absolute.csv")
        normalized = comparison.contingency.div(comparison.contingency.sum(axis=1), axis=0).fillna(0.0)
        normalized.to_csv(destination / f"{edge}_contingency_normalized.csv")
        comparison.assignments.to_csv(destination / f"{edge}_unit_assignments.csv.gz", index=False, compression="gzip")
    pd.concat(summaries, ignore_index=True).to_csv(destination / "partition_summary.csv", index=False)
    analysis.sweep.to_csv(destination / "resolution_sweep.csv", index=False)
    pd.DataFrame([analysis.audit | analysis.representation_audit]).to_csv(
        destination / "partition_audit.csv", index=False
    )
    return destination


def write_spatial_artifacts(output_dir: str | Path, analysis: SpatialImpactAnalysis) -> Path:
    """Write route-independent spatial tables resolved by :func:`compute_spatial_impact`."""
    import json

    destination = Path(output_dir) / "spatial"
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "scale_audit.json").write_text(
        json.dumps(analysis.scale_audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (destination / "support_selection.json").write_text(
        json.dumps(analysis.support_selection, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    analysis.support_sensitivity.to_csv(destination / "support_sensitivity.csv", index=False)
    analysis.unit_assignments.to_csv(destination / "unit_window_assignments.csv.gz", compression="gzip")
    analysis.anatomy_unit_assignments.to_csv(destination / "anatomy_unit_window_assignments.csv.gz", compression="gzip")
    analysis.anatomy_windows.to_csv(destination / "anatomy_candidate_map.csv", index=False)
    analysis.window_metrics.to_csv(destination / "window_metrics.csv", index=False)
    analysis.anatomy_summary.to_csv(destination / "anatomy_region_summary.csv", index=False)
    analysis.scale_sensitivity.to_csv(destination / "scale_sensitivity.csv", index=False)
    analysis.threshold_sensitivity.to_csv(destination / "threshold_sensitivity.csv", index=False)
    return destination
