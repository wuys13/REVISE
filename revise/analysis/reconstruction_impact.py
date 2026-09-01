"""Configuration and artifact helpers for reconstruction-impact analyses."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import gc
from pathlib import Path
from typing import Any, Mapping

import numpy as np
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
    select_resolution_for_target_cluster_count,
    select_shared_feature_names,
    summarize_change_by_level1,
)
from revise.analysis.basic.spatial_region import (
    assign_anatomy_candidates,
    assign_square_windows,
    compute_window_diversity,
    compute_rarefied_window_diversity,
    convert_coordinates_to_microns,
    select_region_threshold,
    select_window_scale,
    summarize_anatomy_context,
    summarize_cluster_change_by_anatomy,
    summarize_diversity_by_anatomy,
    summarize_region_extent_by_anatomy,
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
        "candidate_window_sides_um": [16.0, 24.0, 32.0, 40.0, 56.0, 80.0],
        "min_parent_units": 4,
        "rarefaction_draws": 200,
        "threshold_bootstraps": 500,
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
    complexity_comparisons: dict[str, PartitionComparison]
    complexity_sweep: pd.DataFrame
    complexity_resolution: float
    change_by_level1: pd.DataFrame
    matched_cluster_status: str


@dataclass(frozen=True)
class SpatialImpactAnalysis:
    """Spatial tables resolved for one route at the declared physical scale."""

    scale_audit: dict[str, float | int | str | None]
    support_selection: dict[str, float | int | str | None]
    support_sensitivity: pd.DataFrame
    unit_assignments: pd.DataFrame
    anatomy_unit_assignments: pd.DataFrame
    window_metrics: pd.DataFrame
    anatomy_windows: pd.DataFrame
    anatomy_context_summary: pd.DataFrame
    cluster_change_by_anatomy: pd.DataFrame
    diversity_by_anatomy: pd.DataFrame
    region_extent_by_anatomy: pd.DataFrame
    scale_sensitivity: pd.DataFrame
    threshold_sensitivity: pd.DataFrame
    gain_region_extent_by_anatomy: pd.DataFrame
    state_threshold: dict[str, float | int | str | None]
    gain_threshold: dict[str, float | int | str | None]
    state_threshold_bootstrap: pd.DataFrame
    gain_threshold_bootstrap: pd.DataFrame


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
    """Separate same-resolution complexity diagnostics from matched-K change."""
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
    else:
        resolution = float(within_level1_resolution)
        raw_labels = leiden_labels(raw_graph, resolution=resolution, random_state=random_state)
        sweep = pd.DataFrame(columns=["resolution", "ARI"])

    representation_audit: dict[str, int | bool] = {}
    complexity_comparisons: dict[str, PartitionComparison]
    complexity_sweep = pd.DataFrame()

    def matched_sweep(graph: AnnData, target: int) -> tuple[pd.Series, pd.DataFrame, pd.Series]:
        coarse = np.round(np.arange(0.1, 1.5 + 0.001, 0.1), 2).tolist()
        labels_by_resolution: dict[float, pd.Series] = {}
        rows: list[dict[str, float | int | str]] = []
        for value in coarse:
            labels = leiden_labels(graph, resolution=value, random_state=random_state)
            labels_by_resolution[value] = labels
            rows.append({"resolution": value, "n_clusters": int(labels.nunique()), "stage": "coarse"})
        coarse_selected = select_resolution_for_target_cluster_count(
            pd.DataFrame(rows), target_cluster_count=target
        )
        lower = max(0.02, float(coarse_selected["resolution"]) - 0.1)
        upper = float(coarse_selected["resolution"]) + 0.1
        fine = np.round(np.arange(lower, upper + 0.001, 0.02), 2).tolist()
        for value in fine:
            if value in labels_by_resolution:
                continue
            labels = leiden_labels(graph, resolution=value, random_state=random_state)
            labels_by_resolution[value] = labels
            rows.append({"resolution": value, "n_clusters": int(labels.nunique()), "stage": "fine"})
        sweep_frame = pd.DataFrame(rows).sort_values("resolution", kind="stable").reset_index(drop=True)
        selected = select_resolution_for_target_cluster_count(
            sweep_frame, target_cluster_count=target
        )
        sweep_frame["target_cluster_count"] = int(target)
        sweep_frame["cluster_count_gap"] = (sweep_frame["n_clusters"] - int(target)).abs()
        return labels_by_resolution[float(selected["resolution"])], sweep_frame, selected

    if route_kind == "sp_svc":
        recon_graph = prepare_leiden_graph(recon_work, feature_names=features)
        recon_complexity_labels = leiden_labels(
            recon_graph, resolution=resolution, random_state=random_state
        )
        complexity_comparisons = {
            "raw_to_recon_same_resolution": compare_partitions(
                raw_labels,
                recon_complexity_labels,
                comparison_edge="raw_to_recon_same_resolution",
            )
        }
        recon_expression_labels, complexity_sweep, matched = matched_sweep(
            recon_graph, int(raw_labels.nunique())
        )
        comparisons = {
            "raw_to_recon_expression": compare_partitions(
                raw_labels,
                recon_expression_labels,
                comparison_edge="raw_to_recon_expression",
            )
        }
        matched_resolution = float(matched["resolution"])
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
        complexity_comparisons = {
            "raw_to_final_svc_complexity_diagnostic": compare_partitions(
                raw_labels,
                final_labels,
                comparison_edge="raw_to_final_svc_complexity_diagnostic",
            )
        }
        raw_labels, complexity_sweep, matched = matched_sweep(raw_graph, int(final_labels.nunique()))
        comparisons = {
            "raw_to_final_svc": compare_partitions(
                raw_labels,
                final_labels,
                comparison_edge="raw_to_final_svc",
            )
        }
        matched_resolution = float(matched["resolution"])
    matched_comparison = next(iter(comparisons.values()))
    return PartitionAnalysis(
        matched_resolution,
        "matched_cluster_count",
        complexity_sweep,
        audit,
        comparisons,
        representation_audit,
        features,
        complexity_comparisons,
        sweep,
        resolution,
        summarize_change_by_level1(matched_comparison.assignments, level1),
        str(matched["status"]),
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
    candidate_window_sides_um: list[float] | None = None,
    min_parent_units: int = 4,
    rarefaction_draws: int = 200,
    threshold_bootstraps: int = 500,
    tumor_label: str = "Tumor",
    normal_source_label: str = "Intestinal Epithelial",
) -> SpatialImpactAnalysis:
    """Compute parent-specific rarefied diversity and data-driven State/Gain Regions."""
    paired_ids = paired_coordinates.index
    for name, labels in (("Raw", raw_labels), ("Reconstructed", reconstructed_labels), ("unit_changed", unit_changed)):
        if not labels.index.is_unique or set(labels.index) != set(paired_ids):
            raise ValueError(f"{name} labels must contain exactly the paired coordinate IDs")
    if not full_level1_labels.index.is_unique or set(full_level1_labels.index) != set(full_coordinates.index):
        raise ValueError("Full Level1 labels must contain exactly the full coordinate IDs")

    full_um = convert_coordinates_to_microns(full_coordinates, microns_per_coordinate=microns_per_coordinate)
    paired_um = convert_coordinates_to_microns(paired_coordinates, microns_per_coordinate=microns_per_coordinate)
    origin = (float(full_um["x"].min()), float(full_um["y"].min()))
    candidates = candidate_window_sides_um or [16.0, 24.0, 32.0, 40.0, 56.0, 80.0]
    support_selection, support_sensitivity = select_window_scale(
        paired_um,
        candidate_window_sides=candidates,
        min_parent_units=min_parent_units,
        origin=origin,
    )
    main_side = support_selection["window_side_length"]
    if main_side is None:
        main_side = float(min(candidates))
    main_side = float(main_side)
    anatomy_assignments = assign_square_windows(full_um, window_side_length=main_side, origin=origin)
    anatomy_windows = assign_anatomy_candidates(
        anatomy_assignments,
        full_level1_labels,
        tumor_label=tumor_label,
        normal_source_label=normal_source_label,
    )
    paired_windows = assign_square_windows(paired_um, window_side_length=main_side, origin=origin)
    metrics = compute_rarefied_window_diversity(
        paired_windows,
        raw_labels,
        reconstructed_labels,
        min_parent_units=min_parent_units,
        n_draws=rarefaction_draws,
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
    valid_metrics = metrics.loc[metrics["valid_window"]]
    state_threshold, state_bootstrap = select_region_threshold(
        valid_metrics["neff_recon"].to_numpy(),
        n_bootstrap=threshold_bootstraps,
    )
    gain_threshold, gain_bootstrap = select_region_threshold(
        valid_metrics.loc[valid_metrics["delta_neff"] > 0, "delta_neff"].to_numpy(),
        n_bootstrap=threshold_bootstraps,
    )
    metrics["in_state_region"] = pd.Series(pd.NA, index=metrics.index, dtype="boolean")
    metrics["in_gain_region"] = pd.Series(pd.NA, index=metrics.index, dtype="boolean")
    if state_threshold["status"] == "ok":
        metrics.loc[:, "in_state_region"] = (
            metrics["valid_window"] & (metrics["neff_recon"] >= float(state_threshold["threshold"]))
        )
    if gain_threshold["status"] == "ok":
        metrics.loc[:, "in_gain_region"] = (
            metrics["valid_window"] & (metrics["delta_neff"] >= float(gain_threshold["threshold"]))
        )
    unit_assignments = paired_windows.loc[:, ["x", "y", "window_id"]].copy()
    unit_assignments["raw_cluster"] = raw_labels.reindex(unit_assignments.index).astype(str)
    unit_assignments["reconstructed_cluster"] = reconstructed_labels.reindex(unit_assignments.index).astype(str)
    unit_assignments["unit_changed"] = changed
    unit_assignments = unit_assignments.merge(
        anatomy_windows.loc[:, ["window_id", "level1_region"]], on="window_id", how="left", validate="many_to_one"
    )
    unit_assignments.index = paired_windows.index
    anatomy_context_summary = summarize_anatomy_context(
        anatomy_assignments,
        anatomy_windows,
        window_side_length=main_side,
    )
    cluster_change_by_anatomy = summarize_cluster_change_by_anatomy(
        unit_assignments,
        metrics,
    )
    diversity_by_anatomy = summarize_diversity_by_anatomy(metrics)
    metrics["in_region"] = metrics["in_state_region"]
    region_extent_by_anatomy = summarize_region_extent_by_anatomy(
        metrics,
        window_side_length=main_side,
    )
    metrics["in_region"] = metrics["in_gain_region"]
    gain_region_extent_by_anatomy = summarize_region_extent_by_anatomy(
        metrics,
        window_side_length=main_side,
    )
    metrics = metrics.drop(columns=["in_region"])
    threshold_sensitivity = pd.concat(
        [
            pd.DataFrame([state_threshold | {"region_type": "state"}]),
            pd.DataFrame([gain_threshold | {"region_type": "gain"}]),
        ],
        ignore_index=True,
    )
    scale_rows = []
    for side in candidates:
        assignments = assign_square_windows(paired_um, window_side_length=side, origin=origin)
        scale_metrics = compute_window_diversity(
            assignments,
            raw_labels,
            reconstructed_labels,
            min_units_per_window=min_parent_units,
        )
        valid = scale_metrics.loc[scale_metrics["valid_window"]]
        scale_rows.append(
            {
                "window_side_length": float(side),
                "n_valid_windows": int(valid.shape[0]),
                "median_neff_recon": float(valid["neff_recon"].median()) if not valid.empty else np.nan,
                "median_delta_neff": float(valid["delta_neff"].median()) if not valid.empty else np.nan,
            }
        )
        del assignments, scale_metrics
        gc.collect()
    scale_sensitivity = pd.DataFrame(scale_rows)
    scale_audit = {
        "coordinate_unit": "coordinate",
        "microns_per_coordinate": float(microns_per_coordinate),
        "main_window_side_um": main_side,
        "main_window_area_um2": main_side**2,
        "min_parent_units": int(min_parent_units),
        "rarefaction_draws": int(rarefaction_draws),
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
        anatomy_context_summary=anatomy_context_summary,
        cluster_change_by_anatomy=cluster_change_by_anatomy,
        diversity_by_anatomy=diversity_by_anatomy,
        region_extent_by_anatomy=region_extent_by_anatomy,
        scale_sensitivity=scale_sensitivity,
        threshold_sensitivity=threshold_sensitivity,
        gain_region_extent_by_anatomy=gain_region_extent_by_anatomy,
        state_threshold=state_threshold,
        gain_threshold=gain_threshold,
        state_threshold_bootstrap=state_bootstrap,
        gain_threshold_bootstrap=gain_bootstrap,
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
    analysis.complexity_sweep.to_csv(destination / "complexity_raw_resolution_sweep.csv", index=False)
    analysis.sweep.to_csv(destination / "matched_k_resolution_sweep.csv", index=False)
    analysis.change_by_level1.to_csv(destination / "change_by_level1.csv", index=False)
    for edge, comparison in analysis.complexity_comparisons.items():
        comparison.summary.to_csv(destination / f"{edge}_summary.csv", index=False)
        comparison.mapping.to_csv(destination / f"{edge}_mapping.csv", index=False)
        comparison.contingency.to_csv(destination / f"{edge}_contingency_absolute.csv")
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
    analysis.anatomy_context_summary.to_csv(destination / "anatomy_context_summary.csv", index=False)
    analysis.cluster_change_by_anatomy.to_csv(destination / "cluster_change_by_anatomy.csv", index=False)
    analysis.diversity_by_anatomy.to_csv(destination / "diversity_by_anatomy.csv", index=False)
    analysis.region_extent_by_anatomy.to_csv(destination / "region_extent_by_anatomy.csv", index=False)
    analysis.gain_region_extent_by_anatomy.to_csv(
        destination / "gain_region_extent_by_anatomy.csv", index=False
    )
    analysis.scale_sensitivity.to_csv(destination / "scale_sensitivity.csv", index=False)
    analysis.threshold_sensitivity.to_csv(destination / "threshold_sensitivity.csv", index=False)
    analysis.state_threshold_bootstrap.to_csv(
        destination / "state_threshold_bootstrap.csv", index=False
    )
    analysis.gain_threshold_bootstrap.to_csv(
        destination / "gain_threshold_bootstrap.csv", index=False
    )
    return destination
