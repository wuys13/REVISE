"""Batch adapter for the route-scoped reconstruction-impact calculations."""

from __future__ import annotations

from collections.abc import Mapping
import json
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

from revise.analysis.reconstruction_impact import (
    DEFAULTS,
    RawLevel2Mapping,
    compute_anatomy_regions,
    compute_spatial_impact,
    file_sha256,
    map_raw_level2_labels,
    run_partition_analysis,
)


CODE_DEPENDENCIES = (
    "revise.analysis.reconstruction_impact",
    "revise.analysis.basic.partition_change",
    "revise.analysis.basic.spatial_region",
    "revise.backend.kernels.ot",
    "revise.backend.ops.assignment",
    "revise.backend.ops.distance",
    "revise.backend.ops.tacco_runtime",
    "revise.analysis.impact_outputs",
)


def _section(parameters: Mapping[str, Any], name: str) -> dict[str, Any]:
    value = parameters.get(name, {})
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise ValueError(f"parameters.{name} must be a mapping")
    return dict(value)


def _positive_limit(value: Any, name: str) -> int | None:
    if value is None:
        return None
    if type(value) is not int:
        raise ValueError(f"{name} must be a positive integer or null")
    if value < 1:
        raise ValueError(f"{name} must be a positive integer or null")
    return value


def _strict_integer(value: Any, name: str, *, default: int | None = None, minimum: int | None = None) -> int:
    if value is None:
        if default is None:
            raise ValueError(f"{name} must be an explicit integer")
        value = default
    if type(value) is not int:
        raise ValueError(f"{name} must be an explicit integer")
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return value


def _normalise_parent_value(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"parameters.parent_labels.{name} must be a nonempty exact raw label")
    return value


def _parent_specs(reconstruction: Mapping[str, Any], parameters: Mapping[str, Any]) -> list[tuple[str, str]]:
    values = parameters.get("parent_labels")
    if not isinstance(values, Mapping) or not values:
        raise ValueError("parameters.parent_labels must map task names to exact Raw Level1 labels")
    specs = []
    for name, value in values.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError("parameters.parent_labels keys must be nonempty task names")
        if reconstruction.get("modality") == "hST" and name == "global":
            raise ValueError("parent task name 'global' is reserved for the hST global scope")
        specs.append((name, _normalise_parent_value(value, name)))
    raw_labels = [value for _, value in specs]
    if len(set(raw_labels)) != len(raw_labels):
        raise ValueError("parameters.parent_labels must not map multiple tasks to one raw label")
    if reconstruction.get("modality") == "iST":
        task = reconstruction.get("cell_type")
        if task not in dict(specs):
            raise ValueError(f"No exact parent_labels entry for iST task {task!r}")
        return [(task, dict(specs)[task])]
    return specs


def _exact_ids(adata: ad.AnnData, level1_col: str, raw_label: str) -> pd.Index:
    if level1_col not in adata.obs:
        raise KeyError(f"Raw input is missing {level1_col!r}")
    values = adata.obs[level1_col].astype(str)
    selected = adata.obs_names[values.eq(raw_label).to_numpy()]
    if selected.empty:
        raise ValueError(f"Raw input has no exact {level1_col} label {raw_label!r}")
    return selected


def _sample_ids(ids: pd.Index, limit: int | None, seed: int) -> pd.Index:
    ids = pd.Index(ids)
    if limit is None or limit >= len(ids):
        return ids.copy()
    if limit < 1:
        raise ValueError("sampling limits must be positive")
    generator = np.random.default_rng(seed)
    positions = np.sort(generator.choice(len(ids), size=limit, replace=False))
    return ids[positions]


def _coordinates(view) -> pd.DataFrame:
    values = view.coordinates
    if values is None:
        raise ValueError(f"{view.source} has no spatial coordinates")
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[1] < 2:
        raise ValueError(f"{view.source} spatial coordinates must be two-dimensional")
    return pd.DataFrame(values[:, :2], index=view.observation_ids, columns=["x", "y"])


def _paired_carriers(raw: ad.AnnData, spatial: ad.AnnData, ids: pd.Index, level1_col: str):
    raw_scope = raw[ids, :].copy()
    reconstructed = spatial[ids, :].copy()
    labels = raw_scope.obs[level1_col].copy()
    if level1_col not in reconstructed.obs or not reconstructed.obs[level1_col].astype(str).equals(
        labels.astype(str)
    ):
        reconstructed.obs[level1_col] = labels.to_numpy()
    return raw_scope, reconstructed


def _partition_options(
    partition: Mapping[str, Any],
    *,
    modality: str,
    level1_col: str,
    route_kind: str,
    random_state: int,
) -> dict[str, Any]:
    resolution_mode = partition.get("mode")
    resolution_candidates = partition.get(
        "level1_resolution_candidates",
        DEFAULTS["partition_change"]["level1_resolution_candidates"],
    )
    options = {
        "level1_col": level1_col,
        "route_kind": route_kind,
        "resolution_mode": resolution_mode,
        "resolution_candidates": [float(value) for value in resolution_candidates],
        "within_level1_resolution": float(
            partition.get("within_level1_resolution", DEFAULTS["partition_change"]["within_level1_resolution"])
        ),
        "random_state": _strict_integer(random_state, "partition_change.random_state"),
        "n_top_genes": _strict_integer(partition.get("n_top_genes"), "partition_change.n_top_genes", default=2000, minimum=1),
        "raw_qc_min_genes": _strict_integer(
            partition.get("raw_qc_min_genes"), "partition_change.raw_qc_min_genes", default=50, minimum=1
        ),
        "raw_qc_min_cells": _strict_integer(
            partition.get("raw_qc_min_cells"), "partition_change.raw_qc_min_cells", default=3, minimum=1
        ),
    }
    feature_names = partition.get("feature_names")
    if feature_names is not None:
        options["feature_names"] = [str(value) for value in feature_names]
    if route_kind == "sc_svc":
        options["final_cluster_key"] = partition.get("reconstructed_cluster_key", "SVC_cluster")
    if modality not in {"hST", "iST"}:
        raise ValueError(f"reconstruction-impact does not support {modality!r}")
    return options


def _mapping_options(
    mapping: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    method = mapping.get("method")
    if method not in {"pot", "tacco"}:
        raise ValueError("parameters.raw_level2_mapping.method must be pot or tacco")
    mapping_level1 = mapping.get("level1_column")
    level2_col = mapping.get("level2_column")
    if not isinstance(mapping_level1, str) or not mapping_level1.strip():
        raise ValueError("raw_level2_mapping.level1_column must be nonempty")
    if not isinstance(level2_col, str) or not level2_col.strip():
        raise ValueError("raw_level2_mapping.level2_column must be nonempty")
    kwargs: dict[str, Any] = {
        "method": method,
        "level1_col": mapping_level1,
        "level2_col": level2_col,
    }
    filter_config = mapping.get("reference_filter")
    if filter_config is not None:
        if not isinstance(filter_config, Mapping):
            raise ValueError("raw_level2_mapping.reference_filter must be a mapping")
        kwargs["reference_filter_column"] = filter_config.get("column")
        kwargs["reference_filter_value"] = filter_config.get("value")
    method_config = mapping.get(method, {})
    if method_config is not None:
        if not isinstance(method_config, Mapping):
            raise ValueError(f"raw_level2_mapping.{method} must be a mapping")
        if method == "pot":
            for key, target in (("reg", "pot_reg"), ("reg_m", "pot_reg_m"), ("reg_type", "pot_reg_type")):
                if key in method_config:
                    kwargs[target] = method_config[key]
        else:
            for key, target in (("multi_center", "tacco_multi_center"), ("lamb", "tacco_lamb")):
                if key in method_config:
                    kwargs[target] = method_config[key]
    return kwargs, {"method": method, "level1_col": mapping_level1, "level2_col": level2_col}


def _map_level2_for_scope(
    raw_scope: ad.AnnData,
    reference: ad.AnnData,
    parent_values: list[str],
    mapping: Mapping[str, Any],
    *,
    level1_col: str,
) -> RawLevel2Mapping:
    kwargs, audit_base = _mapping_options(mapping)
    mappings = []
    for raw_label in parent_values:
        parent_ids = _exact_ids(raw_scope, level1_col, raw_label)
        raw_parent = raw_scope[parent_ids, :].copy()
        mappings.append(
            map_raw_level2_labels(raw_parent, reference, parent_value=raw_label, **kwargs)
        )
    labels = pd.concat([item.labels for item in mappings]).reindex(raw_scope.obs_names)
    if labels.isna().any():
        raise ValueError("Raw Level2 mapping did not cover every scoped Raw observation")
    assignments = pd.concat([item.assignments for item in mappings], axis=0).reindex(raw_scope.obs_names)
    posterior = pd.concat([item.posterior for item in mappings], axis=0, sort=False).reindex(raw_scope.obs_names)
    return RawLevel2Mapping(
        labels=labels,
        posterior=posterior,
        assignments=assignments,
        audit={
            **audit_base,
            "parents": list(parent_values),
            "n_raw_units": int(raw_scope.n_obs),
            "n_reference_cells": int(reference.n_obs),
            "submappings": [item.audit for item in mappings],
        },
    )


def _comparison_edge(partition, route_kind: str):
    expected = "raw_to_final_svc" if route_kind == "sc_svc" else "raw_to_recon_expression"
    if expected not in partition.comparisons:
        raise KeyError(f"Partition analysis did not return expected comparison edge {expected!r}")
    return expected, partition.comparisons[expected]


def _partition_gene_audit(raw_scope: ad.AnnData, partition) -> list[dict[str, str]]:
    """Join Raw-QC eligibility with the subsequent Raw-HVG choice."""
    raw_qc = partition.audit.get("raw_qc_gene_reasons")
    if not isinstance(raw_qc, list):
        return []
    raw_qc_by_gene = {
        item.get("gene_id"): item.get("raw_qc_status")
        for item in raw_qc
        if isinstance(item, Mapping)
        and isinstance(item.get("gene_id"), str)
        and isinstance(item.get("raw_qc_status"), str)
    }
    features = set(pd.Index(partition.feature_names).astype(str))
    audit = []
    for gene_id in raw_scope.var_names.astype(str):
        raw_qc_status = raw_qc_by_gene.get(gene_id)
        if raw_qc_status is None:
            continue
        audit.append(
            {
                "gene_id": gene_id,
                "raw_qc_status": raw_qc_status,
                "partition_feature_status": (
                    "selected_raw_hvg"
                    if raw_qc_status == "retained" and gene_id in features
                    else "not_selected_raw_hvg"
                    if raw_qc_status == "retained"
                    else "not_eligible_raw_qc"
                ),
            }
        )
    return audit


def _scope_analysis(
    *,
    raw_scope: ad.AnnData,
    reconstructed_scope: ad.AnnData,
    full_coordinates: pd.DataFrame,
    full_level1_labels: pd.Series,
    raw_level2: RawLevel2Mapping | None,
    anatomy,
    partition_options: Mapping[str, Any],
    spatial_options: Mapping[str, Any],
    route_kind: str,
    scope_audit: dict[str, Any],
):
    partition = run_partition_analysis(raw_scope, reconstructed_scope, **dict(partition_options))
    edge, comparison = _comparison_edge(partition, route_kind)
    assignments = comparison.assignments
    paired_ids = pd.Index(assignments.index)
    raw_labels = assignments["raw_cluster"].astype(str).reindex(paired_ids)
    reconstructed_labels = assignments["recon_cluster"].astype(str)
    # The explicit reindex below keeps Raw IDs as the sole spatial observation axis.
    paired_coordinates = full_coordinates.reindex(paired_ids).copy()
    if paired_coordinates.isna().any().any():
        raise ValueError("Partition assignments contain IDs absent from the frozen Raw coordinates")
    spatial = None
    if raw_level2 is not None:
        scoped_level2 = raw_level2.labels.reindex(paired_ids)
        spatial = compute_spatial_impact(
            full_coordinates=full_coordinates,
            full_level1_labels=full_level1_labels,
            paired_coordinates=paired_coordinates,
            raw_labels=raw_labels,
            raw_level2_labels=scoped_level2,
            reconstructed_labels=reconstructed_labels,
            unit_changed=assignments["unit_changed"].astype(bool),
            anatomy_analysis=anatomy,
            **dict(spatial_options),
        )
    retained_ids = pd.Index(assignments.index)
    excluded_ids = raw_scope.obs_names[~raw_scope.obs_names.isin(retained_ids)]
    feature_names = pd.Index(partition.feature_names).astype(str)
    excluded_features = raw_scope.var_names.astype(str)[~raw_scope.var_names.astype(str).isin(feature_names)]
    scope_audit = dict(scope_audit)
    scope_audit.update(
        {
            "partition_edge": edge,
            "partition_audit": dict(partition.audit),
            "matched_cluster_status": partition.matched_cluster_status,
            "retained_observation_ids": retained_ids.astype(str).tolist(),
            "excluded_observation_ids": excluded_ids.astype(str).tolist(),
            "feature_ids": feature_names.tolist(),
            "excluded_feature_ids": excluded_features.tolist(),
            "n_retained_observations": int(len(retained_ids)),
            "n_excluded_observations": int(len(excluded_ids)),
            "gene_audit": _partition_gene_audit(raw_scope, partition),
        }
    )
    return {
        "partition": partition,
        "spatial": spatial,
        "raw_level2": raw_level2,
        "audit": scope_audit,
    }


def _write_outputs(output_dir: Path, bundle: Mapping[str, Any]) -> dict[str, Any]:
    """Call the paired writer using the adapter's stable result-bundle seam."""
    from revise.analysis import impact_outputs

    writer = getattr(impact_outputs, "write_impact_outputs", None)
    scope_type = getattr(impact_outputs, "ImpactScope", None)
    if writer is None or scope_type is None:
        raise AttributeError("impact_outputs.write_impact_outputs and ImpactScope are required")
    identity = bundle["common_inputs"]["record_identity"]
    sample_id = str(identity.get("sample_id") or "unknown_sample")
    task_cell_type = identity.get("cell_type")
    scopes = []
    for scope_name, value in bundle["scopes"].items():
        scope_audit = value["audit"]
        comparison_id = json.dumps(
            {"sample_id": sample_id, "task_cell_type": task_cell_type, "scope": scope_name},
            sort_keys=True, ensure_ascii=False, separators=(",", ":"),
        )
        scopes.append(
            scope_type(
                comparison_id=comparison_id,
                sample_id=sample_id,
                task_cell_type=task_cell_type,
                scope=scope_name,
                comparison={
                    "comparison_edge": scope_audit.get("partition_edge", "multiple"),
                    "raw_view": "raw_leiden",
                    "reconstruction_view": "spatial_carrier",
                    "label_source": bundle["level1_col"],
                    "observation_basis": "covered_raw_observation_ids",
                    "gene_rule": "raw_defined_features",
                    "normalization": "helper_defined",
                },
                partition=value["partition"],
                spatial=value["spatial"],
                anatomy=bundle["anatomy"],
                raw_level2=value["raw_level2"],
            )
        )
    audit = {
        "schema_version": 1,
        "common_inputs": bundle["common_inputs"],
        "parameters": bundle["parameters"],
        "resources": bundle["resources"],
        "scope_audits": {name: value["audit"] for name, value in bundle["scopes"].items()},
        "anatomy_scale_audit": bundle["anatomy"].scale_audit,
    }
    calculation = {
        "input_view": "raw+aligned_raw+spatial_carrier",
        "parameters": bundle["parameters"],
        "comparison_basis": "Raw versus reconstructed spatial carrier with Raw Level2 reference baseline",
        "route_kind": bundle["route_kind"],
        "scopes": list(bundle["scopes"]),
        "resource": bundle["resources"]["raw_level2_reference"],
    }
    result = writer(output_dir, scopes=scopes, audit=audit, calculation=calculation)
    if not isinstance(result, Mapping):
        raise ValueError("impact output writer must return an artifact mapping")
    return dict(result)


def run(context):
    """Run route-scoped impact calculations and hand their bundle to the writer."""
    parameters = context.parameters
    reconstruction = context.reconstruction
    modality = reconstruction.get("modality")
    if modality not in {"hST", "iST"}:
        raise ValueError("reconstruction-impact supports hST and iST only")
    if modality == "iST" and reconstruction.get("ist_mapping") != "paired":
        raise ValueError("reconstruction-impact requires paired iST spatial output")

    level1_col = parameters.get("level1_column")
    if not isinstance(level1_col, str) or not level1_col.strip():
        raise ValueError("parameters.level1_column must be a nonempty Raw Level1 column")
    parent_specs = _parent_specs(reconstruction, parameters)
    partition = _section(parameters, "partition_change")
    spatial = _section(parameters, "spatial_region")
    mapping = _section(parameters, "raw_level2_mapping")
    route_kind = parameters.get("route_kind")
    if route_kind is None:
        route_kind = "sc_svc" if modality == "iST" else "sp_svc"
    if route_kind not in {"sp_svc", "sc_svc"}:
        raise ValueError("route_kind must be sp_svc or sc_svc")
    partition_seed = _strict_integer(
        partition.get("random_state"), "partition_change.random_state", default=42
    )
    spatial_seed = _strict_integer(
        spatial.get("random_state"), "spatial_region.random_state", default=partition_seed
    )
    sample_limit = _positive_limit(
        partition.get("sample_n_units"), "partition_change.sample_n_units"
    )
    parent_limit = _positive_limit(
        partition.get("parent_sample_n_units"), "partition_change.parent_sample_n_units"
    )

    reference_path = context.resources.get("raw_level2_reference")
    if reference_path is None:
        raise ValueError("resources.raw_level2_reference is required for Raw Level2 mapping")
    reference_path = Path(reference_path)
    if not reference_path.is_file():
        raise FileNotFoundError(f"Raw Level2 reference does not exist: {reference_path}")

    inputs = context.inputs
    raw_view = inputs.raw()
    aligned_view = inputs.aligned_raw()
    spatial_view = inputs.spatial()
    full_raw = raw_view.adata
    if level1_col not in full_raw.obs:
        raise KeyError(f"Raw input is missing {level1_col!r}")
    full_labels = full_raw.obs[level1_col].astype(str)
    full_coordinates = _coordinates(raw_view)
    covered_raw = aligned_view.adata
    if level1_col not in covered_raw.obs:
        raise KeyError(f"Covered Raw input is missing {level1_col!r}")
    spatial_ids = pd.Index(spatial_view.adata.obs_names)
    if modality == "hST":
        if not spatial_ids.equals(covered_raw.obs_names):
            raise ValueError("hST spatial carrier must cover the aligned Raw observation axis exactly")
        full_reconstructed = spatial_view.adata
    else:
        _, task_raw_label = parent_specs[0]
        task_ids = _exact_ids(covered_raw, level1_col, task_raw_label)
        if not spatial_ids.isin(task_ids).all():
            raise ValueError(
                f"iST spatial carrier contains observations outside exact parent {task_raw_label!r}"
            )
        full_reconstructed = spatial_view.adata
    spatial_microns = spatial.get("microns_per_coordinate")
    if spatial_microns is None:
        raise ValueError("spatial_region.microns_per_coordinate is required")
    anatomy_region = spatial.get("anatomy_region", {})
    if not isinstance(anatomy_region, Mapping):
        raise ValueError("spatial_region.anatomy_region must be a mapping")
    spatial_defaults = DEFAULTS["spatial_region"]
    spatial_options = {
        "microns_per_coordinate": float(spatial_microns),
        "candidate_window_sides_um": [
            float(value)
            for value in spatial.get(
                "candidate_window_sides_um", spatial_defaults["candidate_window_sides_um"]
            )
        ],
        "min_parent_units": _strict_integer(
            spatial.get("min_parent_units"),
            "spatial_region.min_parent_units",
            default=spatial_defaults["min_parent_units"],
            minimum=1,
        ),
        "rarefaction_draws": _strict_integer(
            spatial.get("rarefaction_draws"),
            "spatial_region.rarefaction_draws",
            default=spatial_defaults["rarefaction_draws"],
            minimum=1,
        ),
        "threshold_bootstraps": _strict_integer(
            spatial.get("threshold_bootstraps"),
            "spatial_region.threshold_bootstraps",
            default=spatial_defaults["threshold_bootstraps"],
            minimum=1,
        ),
        "random_state": spatial_seed,
        "cell_equivalent_um": float(
            spatial.get("cell_equivalent_um", spatial_defaults["cell_equivalent_um"])
        ),
        "tumor_label": anatomy_region.get(
            "tumor_label", spatial_defaults["anatomy_region"]["tumor_label"]
        ),
        "normal_source_label": anatomy_region.get(
            "normal_source_label", spatial_defaults["anatomy_region"]["normal_source_label"]
        ),
    }
    anatomy = compute_anatomy_regions(
        full_coordinates=full_coordinates,
        full_level1_labels=full_labels,
        microns_per_coordinate=spatial_options["microns_per_coordinate"],
        candidate_window_sides_um=spatial_options["candidate_window_sides_um"],
        min_parent_units=spatial_options["min_parent_units"],
        cell_equivalent_um=spatial_options["cell_equivalent_um"],
        tumor_label=spatial_options["tumor_label"],
        normal_source_label=spatial_options["normal_source_label"],
    )
    reference = ad.read_h5ad(reference_path)
    partition_options = _partition_options(
        partition,
        modality=modality,
        level1_col=level1_col,
        route_kind=route_kind,
        random_state=partition_seed,
    )

    scopes: list[tuple[str, pd.Index, str, int]] = []
    if modality == "hST":
        global_ids = _sample_ids(covered_raw.obs_names, sample_limit, partition_seed)
        scopes.append(("global", global_ids, "global", partition_seed))
        parent_seed_base = partition_seed + 1
        for index, (task_name, raw_label) in enumerate(parent_specs):
            parent_ids = _exact_ids(covered_raw, level1_col, raw_label)
            scopes.append(
                (
                    task_name,
                    _sample_ids(parent_ids, parent_limit, parent_seed_base + index),
                    raw_label,
                    parent_seed_base + index,
                )
            )
    else:
        task_name, raw_label = parent_specs[0]
        parent_ids = _exact_ids(covered_raw, level1_col, raw_label)
        scopes.append((task_name, _sample_ids(parent_ids, parent_limit, partition_seed), raw_label, partition_seed))

    scope_results: dict[str, dict[str, Any]] = {}
    for scope_name, scope_ids, raw_label, sampling_seed in scopes:
        raw_scope, reconstructed_scope = _paired_carriers(
            covered_raw, full_reconstructed, scope_ids, level1_col
        )
        scope_partition_options = dict(partition_options)
        if modality == "hST" and scope_name != "global":
            scope_partition_options["resolution_mode"] = "fixed_within_level1"
        global_partition_only = modality == "hST" and scope_name == "global"
        raw_level2 = None
        if not global_partition_only:
            raw_level2 = _map_level2_for_scope(
                raw_scope,
                reference,
                [raw_label],
                mapping,
                level1_col=level1_col,
            )
        scope_results[scope_name] = _scope_analysis(
            raw_scope=raw_scope,
            reconstructed_scope=reconstructed_scope,
            full_coordinates=full_coordinates,
            full_level1_labels=full_labels,
            raw_level2=raw_level2,
            anatomy=anatomy,
            partition_options=scope_partition_options,
            spatial_options=spatial_options,
            route_kind=route_kind,
            scope_audit={
                "scope": scope_name,
                "raw_level1_label": raw_label,
                "sampling_seed": int(sampling_seed),
                "sample_limit": sample_limit if scope_name == "global" else parent_limit,
                "sampled_observation_ids": scope_ids.astype(str).tolist(),
                "full_covered_observation_ids": covered_raw.obs_names.astype(str).tolist(),
                "full_raw_anatomy_observation_ids": full_raw.obs_names.astype(str).tolist(),
                "spatial_status": (
                    "not_applicable_parent_baseline" if global_partition_only else "computed"
                ),
                "raw_level2_status": (
                    "not_applicable_parent_baseline" if global_partition_only else "computed"
                ),
            },
        )

    bundle = {
        "modality": modality,
        "route_kind": route_kind,
        "level1_col": level1_col,
        "common_inputs": {
            "raw_source": raw_view.source,
            "aligned_raw_source": aligned_view.source,
            "spatial_source": spatial_view.source,
            "view_provenance": {
                "raw": dict(raw_view.provenance),
                "aligned_raw": dict(aligned_view.provenance),
                "spatial": dict(spatial_view.provenance),
            },
            "raw_source_observation_ids": raw_view.observation_ids.astype(str).tolist(),
            "covered_observation_ids": aligned_view.observation_ids.astype(str).tolist(),
            "raw_observation_ids": covered_raw.obs_names.astype(str).tolist(),
            "raw_gene_ids": covered_raw.var_names.astype(str).tolist(),
            "coordinate_audit": {
                "unit": reconstruction.get("coordinates", {}).get("unit"),
                "microns_per_coordinate": spatial_options["microns_per_coordinate"],
                "observation_ids": full_coordinates.index.astype(str).tolist(),
                "origin_x": float(full_coordinates["x"].min()),
                "origin_y": float(full_coordinates["y"].min()),
            },
            "record_identity": {
                key: reconstruction.get(key)
                for key in ("sample_id", "cell_type", "modality", "ist_mapping", "fingerprint")
            },
        },
        "anatomy": anatomy,
        "scopes": scope_results,
        "parameters": {
            "parent_labels": dict(parent_specs),
            "sample_n_units": sample_limit,
            "parent_sample_n_units": parent_limit,
            "partition_change": dict(partition),
            "spatial_region": dict(spatial),
            "raw_level2_mapping": dict(mapping),
        },
        "resources": {
            "raw_level2_reference": str(reference_path),
            "raw_level2_reference_identity": {
                "path": str(reference_path),
                "sha256": file_sha256(reference_path),
            },
        },
    }
    written = _write_outputs(context.output_dir, bundle)
    artifacts = written.get("artifacts", written)
    if not isinstance(artifacts, Mapping) or not artifacts:
        raise ValueError("impact output writer must return a nonempty artifact mapping")
    calculation = written.get("calculation")
    if calculation is None:
        calculation = {
            "input_view": "raw+aligned_raw+spatial_carrier",
            "parameters": bundle["parameters"],
            "comparison_basis": "Raw versus reconstructed spatial carrier with Raw Level2 reference baseline",
            "route_kind": route_kind,
            "scopes": list(scope_results),
            "resource": str(reference_path),
        }
    elif not isinstance(calculation, Mapping):
        raise ValueError("impact output writer calculation must be a mapping")
    else:
        calculation = dict(calculation)
    return {"artifacts": dict(artifacts), "calculation": calculation}
