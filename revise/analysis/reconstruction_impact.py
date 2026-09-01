"""Configuration and artifact helpers for reconstruction-impact analyses."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
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


DEFAULTS: dict[str, Any] = {
    "partition_change": {
        "level1_resolution_candidates": [0.3, 0.5, 0.8],
        "within_level1_resolution": 0.5,
        "random_state": 42,
        "comparisons": [],
    },
    "spatial_region": {
        "coordinate_unit": "um",
        "base_unit_um": 8.0,
        "window_side_length_um": 40.0,
        "window_side_lengths_um": [8.0, 16.0, 24.0, 40.0, 56.0, 80.0],
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


def run_partition_analysis(
    raw: AnnData,
    reconstructed: AnnData,
    *,
    level1_col: str,
    final_cluster_key: str | None = None,
    resolution_mode: str | None = None,
    resolution_candidates: list[float] | None = None,
    within_level1_resolution: float = 0.5,
    random_state: int = 42,
    n_top_genes: int = 2000,
) -> PartitionAnalysis:
    """Build comparable expression partitions and the applicable three-node edges.

    A multi-Level1 scope calibrates one resolution on Raw's Level1 ARI.  A
    single-Level1 scope uses the fixed internal resolution.  Both expression
    nodes share the same observation IDs and feature names before their own
    PCA/neighbor graphs are constructed.
    """
    if level1_col not in raw.obs:
        raise KeyError(f"Raw input is missing {level1_col!r}")
    raw_work, recon_work, audit = align_observation_pairs(raw, reconstructed)
    if level1_col not in recon_work.obs:
        raise KeyError(f"Reconstructed input is missing {level1_col!r}")
    features = select_shared_feature_names(raw_work, recon_work, n_top_genes=n_top_genes)
    raw_graph = prepare_leiden_graph(raw_work, feature_names=features)
    recon_graph = prepare_leiden_graph(recon_work, feature_names=features)
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
    if final_cluster_key is not None:
        if final_cluster_key not in recon_work.obs:
            raise KeyError(f"Reconstructed input is missing {final_cluster_key!r}")
        final_labels_source = recon_work.obs[final_cluster_key]
        if final_labels_source.isna().any():
            raise ValueError(f"{final_cluster_key!r} must not contain missing labels")
        final_labels = final_labels_source.astype(str)
        comparisons["recon_expression_to_final_svc"] = compare_partitions(
            recon_expression_labels,
            final_labels,
            comparison_edge="recon_expression_to_final_svc",
        )
        comparisons["raw_to_final_svc"] = compare_partitions(
            raw_labels,
            final_labels,
            comparison_edge="raw_to_final_svc",
        )
    return PartitionAnalysis(resolution, resolution_source, sweep, audit, comparisons)


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
