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

from revise.backend.kernels.ot import OTKernel

from revise.analysis.basic.partition_change import (
    PartitionComparison,
    align_observation_pairs,
    compare_partitions,
    filter_paired_sp_svc_inputs,
    leiden_labels,
    prepare_leiden_graph,
    select_level1_resolution,
    select_raw_hvg_feature_names,
    select_resolution_for_target_cluster_count,
    select_shared_feature_names,
    summarize_change_by_level1,
)
from revise.analysis.basic.spatial_region import (
    assign_anatomy_candidates,
    assign_square_windows,
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

    @property
    def matched_k_sweep(self) -> pd.DataFrame:
        """Return the sweep used to choose the matched-K comparison resolution."""
        return self.sweep

    @property
    def raw_complexity_sweep(self) -> pd.DataFrame:
        """Return the Raw-side resolution/Level1-ARI complexity sweep."""
        return self.complexity_sweep


@dataclass(frozen=True)
class AnatomyRegionAnalysis:
    """One route-level Level1 anatomy Region definition."""

    scale_audit: dict[str, float | int | str | None]
    support_selection: dict[str, float | int | str | None]
    support_sensitivity: pd.DataFrame
    anatomy_unit_assignments: pd.DataFrame
    anatomy_windows: pd.DataFrame
    anatomy_context_summary: pd.DataFrame


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


@dataclass(frozen=True)
class RawLevel2Mapping:
    """Reference-derived Level2 labels assigned from one Raw parent expression cohort."""

    labels: pd.Series
    posterior: pd.DataFrame
    assignments: pd.DataFrame
    audit: dict[str, Any]


def select_raw_level1_parent_cohort(
    raw_level1_labels: pd.Series,
    carrier_ids: pd.Index,
    *,
    parent_value: str,
) -> tuple[pd.Index, dict[str, int | str]]:
    """Restrict a carrier to its original Raw Level1 parent with explicit audit."""
    if not raw_level1_labels.index.is_unique or not carrier_ids.is_unique:
        raise ValueError("Raw Level1 labels and carrier IDs must be unique")
    missing = carrier_ids.difference(raw_level1_labels.index)
    if not missing.empty:
        raise ValueError(f"Carrier contains {len(missing)} IDs absent from Raw")
    normalized_parent = str(parent_value).replace("/", "_")
    carrier_level1 = (
        raw_level1_labels.reindex(carrier_ids)
        .astype(str)
        .str.replace("/", "_", regex=False)
    )
    selected = carrier_ids[carrier_level1.eq(normalized_parent).to_numpy()]
    if selected.empty:
        raise ValueError(f"Carrier has no original Raw Level1 units for {normalized_parent!r}")
    audit: dict[str, int | str] = {
        "parent": normalized_parent,
        "carrier_units": int(len(carrier_ids)),
        "raw_level1_parent_units": int(len(selected)),
        "excluded_raw_level1_mismatch": int(len(carrier_ids) - len(selected)),
    }
    return selected, audit


def map_raw_level2_labels(
    raw_parent: AnnData,
    reference: AnnData,
    *,
    parent_value: str,
    method: str,
    level1_col: str = "Level1",
    level2_col: str = "Level2",
    reference_filter_column: str | None = None,
    reference_filter_value: str | None = None,
    pot_reg: float = 0.1,
    pot_reg_m: float = 0.0,
    pot_reg_type: str = "entropy",
    tacco_multi_center: int | None = None,
    tacco_lamb: float | None = None,
) -> RawLevel2Mapping:
    """Map reference Level2 labels onto Raw expression after Level1 parent selection."""
    if method not in {"pot", "tacco"}:
        raise ValueError("method must be pot or tacco")
    if level1_col not in raw_parent.obs or level1_col not in reference.obs:
        raise KeyError(f"Raw and reference inputs must contain {level1_col!r}")
    if level2_col not in reference.obs:
        raise KeyError(f"Reference input must contain {level2_col!r}")

    normalized_parent = str(parent_value).replace("/", "_")
    raw = raw_parent.copy()
    raw_level1 = raw.obs[level1_col].astype(str).str.replace("/", "_", regex=False)
    if not raw_level1.eq(normalized_parent).all():
        raise ValueError("Raw Level2 mapping input must contain exactly one requested Level1 parent")
    raw.obs[level1_col] = raw_level1

    selected_reference = reference.copy()
    if reference_filter_column is not None:
        if reference_filter_column not in selected_reference.obs:
            raise KeyError(f"Reference input is missing {reference_filter_column!r}")
        selected_reference = selected_reference[
            selected_reference.obs[reference_filter_column].astype(str).eq(str(reference_filter_value))
        ].copy()
    selected_reference.obs[level1_col] = (
        selected_reference.obs[level1_col].astype(str).str.replace("/", "_", regex=False)
    )
    selected_reference = selected_reference[
        selected_reference.obs[level1_col].eq(normalized_parent)
    ].copy()
    if selected_reference.n_obs == 0:
        raise ValueError(f"Reference has no cells for Level1 parent {normalized_parent!r}")
    if selected_reference.obs[level2_col].isna().any():
        raise ValueError("Reference Level2 labels must not be missing")
    if raw.var_names.intersection(selected_reference.var_names).empty:
        raise ValueError("Raw and reference inputs must share genes for Level2 mapping")

    confidence_col = f"{level2_col}_confidence"
    annotated = OTKernel.annotate(
        raw,
        selected_reference,
        method=method,
        annotation_key=level2_col,
        confidence_key=confidence_col,
        pot_reg=pot_reg,
        pot_reg_m=pot_reg_m,
        pot_reg_type=pot_reg_type,
        pot_verbose=False,
        multi_center=tacco_multi_center,
        lamb=tacco_lamb,
    )
    if level2_col not in annotated.obs or annotated.obs[level2_col].isna().any():
        raise ValueError("Raw Level2 mapping returned missing assignments")
    if confidence_col not in annotated.obs or level2_col not in annotated.obsm:
        raise ValueError("Raw Level2 mapping must return confidence and posterior outputs")
    labels = annotated.obs[level2_col].astype(str).rename("raw_level2")
    posterior_source = annotated.obsm[level2_col]
    posterior = (
        posterior_source.copy()
        if isinstance(posterior_source, pd.DataFrame)
        else pd.DataFrame(posterior_source, index=annotated.obs_names)
    )
    assignments = pd.DataFrame(
        {
            "unit_id": annotated.obs_names.astype(str),
            "raw_level2": labels.to_numpy(),
            "confidence": annotated.obs[confidence_col].to_numpy(dtype=float),
            "mapping_method": method,
            "source": "raw_expression_reference_level2",
        },
        index=annotated.obs_names,
    )
    audit = {
        "source": "raw_expression_reference_level2",
        "method": method,
        "parent": normalized_parent,
        "n_raw_units": int(raw.n_obs),
        "n_reference_cells": int(selected_reference.n_obs),
        "n_reference_level2": int(selected_reference.obs[level2_col].nunique()),
        "n_mapped_level2": int(labels.nunique()),
        "n_missing": int(labels.isna().sum()),
        "reference_filter_column": reference_filter_column,
        "reference_filter_value": reference_filter_value,
    }
    return RawLevel2Mapping(labels, posterior, assignments, audit)


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
    raw_qc_min_genes: int = 50,
    raw_qc_min_cells: int = 3,
) -> PartitionAnalysis:
    """Separate same-resolution complexity diagnostics from matched-K change."""
    if route_kind not in {"sp_svc", "sc_svc"}:
        raise ValueError("route_kind must be sp_svc or sc_svc")
    if level1_col not in raw.obs:
        raise KeyError(f"Raw input is missing {level1_col!r}")
    if route_kind == "sp_svc":
        raw_work, recon_work, audit = filter_paired_sp_svc_inputs(
            raw,
            reconstructed,
            min_genes=raw_qc_min_genes,
            min_cells=raw_qc_min_cells,
        )
    else:
        raw_work, recon_work, audit = align_observation_pairs(raw, reconstructed)
    if level1_col not in recon_work.obs:
        raise KeyError(f"Reconstructed input is missing {level1_col!r}")
    if feature_names is None:
        if route_kind == "sp_svc":
            features = select_raw_hvg_feature_names(raw_work, n_top_genes=n_top_genes)
            audit["feature_selection"] = "raw_canonical_seurat_v3_shared"
        else:
            features = select_shared_feature_names(raw_work, recon_work, n_top_genes=n_top_genes)
            audit["feature_selection"] = "paired_average_variance"
    else:
        features = [str(feature) for feature in feature_names]
        if not features or not set(features) <= set(raw_work.var_names):
            raise ValueError("feature_names must be a non-empty shared ordered gene list")
        audit["feature_selection"] = "preselected_shared_features"
    audit["n_features"] = int(len(features))
    audit["leiden_backend"] = "igraph"
    audit["leiden_n_iterations"] = 2
    audit["leiden_random_state"] = int(random_state)
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
    matched_k_sweep = pd.DataFrame()

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
        recon_expression_labels, matched_k_sweep, matched = matched_sweep(
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
        raw_labels, matched_k_sweep, matched = matched_sweep(raw_graph, int(final_labels.nunique()))
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
        resolution=matched_resolution,
        resolution_source="matched_cluster_count",
        sweep=matched_k_sweep,
        audit=audit,
        comparisons=comparisons,
        representation_audit=representation_audit,
        feature_names=features,
        complexity_comparisons=complexity_comparisons,
        complexity_sweep=sweep,
        complexity_resolution=resolution,
        change_by_level1=summarize_change_by_level1(matched_comparison.assignments, level1),
        matched_cluster_status=str(matched["status"]),
    )


def compute_anatomy_regions(
    *,
    full_coordinates: pd.DataFrame,
    full_level1_labels: pd.Series,
    microns_per_coordinate: float,
    candidate_window_sides_um: list[float] | None = None,
    min_parent_units: int = 4,
    cell_equivalent_um: float = 8.0,
    tumor_label: str = "Tumor",
    normal_source_label: str = "Intestinal Epithelial",
) -> AnatomyRegionAnalysis:
    """Define one full-tissue Level1 anatomy Region grid for a route."""
    if not full_level1_labels.index.is_unique or set(full_level1_labels.index) != set(full_coordinates.index):
        raise ValueError("Full Level1 labels must contain exactly the full coordinate IDs")
    if cell_equivalent_um <= 0:
        raise ValueError("cell_equivalent_um must be positive")
    full_um = convert_coordinates_to_microns(full_coordinates, microns_per_coordinate=microns_per_coordinate)
    origin = (float(full_um["x"].min()), float(full_um["y"].min()))
    candidates = candidate_window_sides_um or [16.0, 24.0, 32.0, 40.0, 56.0, 80.0]
    support_selection, support_sensitivity = select_window_scale(
        full_um,
        candidate_window_sides=candidates,
        min_parent_units=min_parent_units,
        origin=origin,
    )
    main_side = support_selection["window_side_length"]
    if main_side is None:
        main_side = float(min(candidates))
    main_side = float(main_side)
    assignments = assign_square_windows(full_um, window_side_length=main_side, origin=origin)
    windows = assign_anatomy_candidates(
        assignments,
        full_level1_labels,
        tumor_label=tumor_label,
        normal_source_label=normal_source_label,
    )
    summary = summarize_anatomy_context(assignments, windows, window_side_length=main_side)
    return AnatomyRegionAnalysis(
        scale_audit={
            "coordinate_unit": "coordinate",
            "microns_per_coordinate": float(microns_per_coordinate),
            "cell_equivalent_um": float(cell_equivalent_um),
            "main_window_side_um": main_side,
            "main_window_cells_per_side": float(main_side / cell_equivalent_um),
            "main_window_area_um2": main_side**2,
            "min_parent_units": int(min_parent_units),
            "origin_x_um": origin[0],
            "origin_y_um": origin[1],
        },
        support_selection=support_selection,
        support_sensitivity=support_sensitivity,
        anatomy_unit_assignments=assignments,
        anatomy_windows=windows,
        anatomy_context_summary=summary,
    )


def _map_to_anatomy_regions(
    coordinates: pd.DataFrame,
    anatomy: AnatomyRegionAnalysis,
) -> pd.DataFrame:
    """Map points to the already-defined route-level anatomy grid."""
    assignments = assign_square_windows(
        coordinates,
        window_side_length=float(anatomy.scale_audit["main_window_side_um"]),
        origin=(
            float(anatomy.scale_audit["origin_x_um"]),
            float(anatomy.scale_audit["origin_y_um"]),
        ),
    )
    regions = anatomy.anatomy_windows.set_index("window_id")["level1_region"]
    result = assignments.loc[:, ["window_id"]].rename(columns={"window_id": "anatomy_window_id"})
    result["level1_region"] = result["anatomy_window_id"].map(regions)
    if result["level1_region"].isna().any():
        raise ValueError("Every parent coordinate must map to route-level anatomy")
    return result


def _summarize_parent_window_anatomy(unit_assignments: pd.DataFrame) -> pd.DataFrame:
    """Assign each parent window from its observed macro-anatomy units.

    Parent and anatomy grids intentionally use independent, occupancy-selected
    scales. A parent-window center can therefore lie in an empty anatomy cell;
    its observed units are the only valid context assignment.
    """
    counts = (
        unit_assignments.groupby(
            ["window_id", "anatomy_window_id", "level1_region"], sort=True
        )
        .size()
        .rename("n_units")
        .reset_index()
    )
    return (
        counts.sort_values(
            ["window_id", "n_units", "anatomy_window_id"],
            ascending=[True, False, True],
            kind="stable",
        )
        .drop_duplicates("window_id")
        .loc[:, ["window_id", "anatomy_window_id", "level1_region"]]
    )


def compute_spatial_impact(
    *,
    full_coordinates: pd.DataFrame,
    full_level1_labels: pd.Series,
    paired_coordinates: pd.DataFrame,
    raw_labels: pd.Series,
    raw_level2_labels: pd.Series | None = None,
    reconstructed_labels: pd.Series,
    unit_changed: pd.Series,
    microns_per_coordinate: float,
    candidate_window_sides_um: list[float] | None = None,
    min_parent_units: int = 4,
    rarefaction_draws: int = 200,
    threshold_bootstraps: int = 500,
    random_state: int = 42,
    cell_equivalent_um: float = 8.0,
    anatomy_analysis: AnatomyRegionAnalysis | None = None,
    tumor_label: str = "Tumor",
    normal_source_label: str = "Intestinal Epithelial",
) -> SpatialImpactAnalysis:
    """Compute parent-specific diversity against fixed route-level anatomy Regions."""
    paired_ids = paired_coordinates.index
    for name, labels in (("Raw", raw_labels), ("Reconstructed", reconstructed_labels), ("unit_changed", unit_changed)):
        if not labels.index.is_unique or set(labels.index) != set(paired_ids):
            raise ValueError(f"{name} labels must contain exactly the paired coordinate IDs")
    if raw_level2_labels is not None and (
        not raw_level2_labels.index.is_unique
        or set(raw_level2_labels.index) != set(paired_ids)
    ):
        raise ValueError("Raw Level2 labels must contain exactly the paired coordinate IDs")
    if not full_level1_labels.index.is_unique or set(full_level1_labels.index) != set(full_coordinates.index):
        raise ValueError("Full Level1 labels must contain exactly the full coordinate IDs")

    paired_um = convert_coordinates_to_microns(paired_coordinates, microns_per_coordinate=microns_per_coordinate)
    anatomy = anatomy_analysis or compute_anatomy_regions(
        full_coordinates=full_coordinates,
        full_level1_labels=full_level1_labels,
        microns_per_coordinate=microns_per_coordinate,
        candidate_window_sides_um=candidate_window_sides_um,
        min_parent_units=min_parent_units,
        cell_equivalent_um=cell_equivalent_um,
        tumor_label=tumor_label,
        normal_source_label=normal_source_label,
    )
    origin = (
        float(anatomy.scale_audit["origin_x_um"]),
        float(anatomy.scale_audit["origin_y_um"]),
    )
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
    paired_windows = assign_square_windows(paired_um, window_side_length=main_side, origin=origin)
    metrics = compute_rarefied_window_diversity(
        paired_windows,
        raw_labels,
        reconstructed_labels,
        raw_level2_labels=raw_level2_labels,
        min_parent_units=min_parent_units,
        n_draws=rarefaction_draws,
        random_state=random_state,
    )
    changed = unit_changed.astype(bool).reindex(paired_windows.index)
    window_change = (
        pd.DataFrame({"window_id": paired_windows["window_id"], "unit_changed": changed})
        .groupby("window_id", sort=True)["unit_changed"]
        .agg(n_changed_units="sum", unit_change_fraction="mean")
        .reset_index()
    )
    metrics = metrics.merge(window_change, on="window_id", how="left", validate="one_to_one")
    unit_assignments = paired_windows.loc[
        :, ["x", "y", "window_id", "window_x_index", "window_y_index"]
    ].copy()
    unit_assignments["raw_cluster"] = raw_labels.reindex(unit_assignments.index).astype(str)
    unit_assignments["reconstructed_cluster"] = reconstructed_labels.reindex(unit_assignments.index).astype(str)
    unit_assignments["unit_changed"] = changed
    unit_context = _map_to_anatomy_regions(unit_assignments.loc[:, ["x", "y"]], anatomy)
    unit_assignments["anatomy_window_id"] = unit_context["anatomy_window_id"].to_numpy()
    unit_assignments["level1_region"] = unit_context["level1_region"].to_numpy()
    unit_assignments.index = paired_windows.index
    metric_context = _summarize_parent_window_anatomy(unit_assignments)
    metrics = metrics.merge(metric_context, on="window_id", how="left", validate="one_to_one")
    if metrics["level1_region"].isna().any():
        raise ValueError("Every paired window must contain mapped anatomy units")
    valid_metrics = metrics.loc[metrics["valid_window"]]
    state_threshold, state_bootstrap = select_region_threshold(
        valid_metrics["neff_recon"].to_numpy(),
        n_bootstrap=threshold_bootstraps,
        random_state=random_state,
    )
    gain_threshold, gain_bootstrap = select_region_threshold(
        valid_metrics.loc[valid_metrics["delta_neff"] > 0, "delta_neff"].to_numpy(),
        n_bootstrap=threshold_bootstraps,
        random_state=random_state,
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
    threshold_sensitivity = pd.DataFrame(
        [
            state_threshold | {"region_type": "state"},
            gain_threshold | {"region_type": "gain"},
        ]
    )
    scale_rows = []
    for side in candidates:
        assignments = assign_square_windows(paired_um, window_side_length=side, origin=origin)
        scale_metrics = compute_rarefied_window_diversity(
            assignments,
            raw_labels,
            reconstructed_labels,
            raw_level2_labels=raw_level2_labels,
            min_parent_units=min_parent_units,
            n_draws=rarefaction_draws,
            random_state=random_state,
        )
        valid = scale_metrics.loc[scale_metrics["valid_window"]]
        scale_row: dict[str, float | int] = {
            "window_side_length": float(side),
            "n_valid_windows": int(valid.shape[0]),
            "rarefaction_draws": int(rarefaction_draws),
            "random_state": int(random_state),
        }
        for metric in ("k_obs", "entropy", "neff", "evenness"):
            for column in (
                f"{metric}_raw",
                f"{metric}_recon",
                f"delta_{metric}",
                f"{metric}_level2",
                f"delta_{metric}_vs_raw_leiden",
                f"delta_{metric}_vs_raw_level2",
            ):
                if column not in valid:
                    continue
                values = pd.to_numeric(valid[column], errors="coerce")
                scale_row[f"median_{column}"] = (
                    float(values.median()) if values.notna().any() else np.nan
                )
        scale_rows.append(scale_row)
        del assignments, scale_metrics
        gc.collect()
    scale_sensitivity = pd.DataFrame(scale_rows)
    scale_audit = {
        "coordinate_unit": "coordinate",
        "microns_per_coordinate": float(microns_per_coordinate),
        "main_window_side_um": main_side,
        "main_window_cells_per_side": float(main_side / cell_equivalent_um),
        "main_window_area_um2": main_side**2,
        "min_parent_units": int(min_parent_units),
        "rarefaction_draws": int(rarefaction_draws),
        "rarefaction_random_state": int(random_state),
        "threshold_random_state": int(random_state),
        "origin_x_um": origin[0],
        "origin_y_um": origin[1],
    }
    return SpatialImpactAnalysis(
        scale_audit=scale_audit,
        support_selection=support_selection,
        support_sensitivity=support_sensitivity,
        unit_assignments=unit_assignments,
        anatomy_unit_assignments=anatomy.anatomy_unit_assignments,
        window_metrics=metrics,
        anatomy_windows=anatomy.anatomy_windows,
        anatomy_context_summary=anatomy.anatomy_context_summary,
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
        summary["matched_cluster_status"] = analysis.matched_cluster_status
        summary["headline_eligible"] = analysis.matched_cluster_status != "unmatched_cluster_complexity"
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


def write_anatomy_artifacts(output_dir: str | Path, analysis: AnatomyRegionAnalysis) -> Path:
    """Write one route-level anatomy Region definition and scale decision."""
    import json

    destination = Path(output_dir) / "anatomy"
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "scale_audit.json").write_text(
        json.dumps(analysis.scale_audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (destination / "support_selection.json").write_text(
        json.dumps(analysis.support_selection, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    analysis.support_sensitivity.to_csv(destination / "scale_decision.csv", index=False)
    analysis.anatomy_unit_assignments.to_csv(
        destination / "window_assignments.csv.gz", compression="gzip"
    )
    analysis.anatomy_windows.to_csv(destination / "anatomy_candidate_map.csv", index=False)
    analysis.anatomy_context_summary.to_csv(destination / "region_summary.csv", index=False)
    return destination


def write_raw_level2_artifacts(
    output_dir: str | Path, mapping: RawLevel2Mapping
) -> Path:
    """Persist reference-derived Raw Level2 labels and their provenance."""
    import json

    destination = Path(output_dir) / "raw_level2"
    destination.mkdir(parents=True, exist_ok=True)
    mapping.assignments.to_csv(
        destination / "assignments.csv.gz", index=False, compression="gzip"
    )
    mapping.posterior.to_csv(destination / "posterior.csv.gz", compression="gzip")
    (destination / "audit.json").write_text(
        json.dumps(mapping.audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
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
    (destination / "state_threshold.json").write_text(
        json.dumps(analysis.state_threshold, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (destination / "gain_threshold.json").write_text(
        json.dumps(analysis.gain_threshold, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    analysis.state_threshold_bootstrap.to_csv(
        destination / "state_threshold_bootstrap.csv", index=False
    )
    analysis.gain_threshold_bootstrap.to_csv(
        destination / "gain_threshold_bootstrap.csv", index=False
    )
    return destination
