"""Compute one reconstruction-impact route in bounded stages for notebook execution.

Each stage writes the same internal workflow artifacts as the source notebook.
The temporary checkpoint is used only to let a constrained Jupyter execution
render the real results without repeating the expensive workflow in one process.
"""

# ruff: noqa: E402 -- the repository root must be added before local imports.

from __future__ import annotations

import argparse
import pickle
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import anndata as ad
import numpy as np
import pandas as pd

from revise.analysis.reconstruction_impact import (
    compute_anatomy_regions,
    compute_spatial_impact,
    file_sha256,
    load_reconstruction_impact_config,
    map_raw_level2_labels,
    run_partition_analysis,
    select_raw_level1_parent_cohort,
    write_analysis_artifacts,
    write_anatomy_artifacts,
    write_partition_artifacts,
    write_raw_level2_artifacts,
    write_spatial_artifacts,
)


ROUTES = {
    "visiumhd": "configs/analysis/reconstruction_impact_visiumhd_p1crc.yaml",
    "xenium": "configs/analysis/reconstruction_impact_xenium_p2crc_fibroblast.yaml",
}
PARENT_SOURCE = {"Fibroblast": "Fibroblast", "Mono_Macro": "Mono/Macro", "T": "T"}


def coordinates(adata_obj: ad.AnnData) -> pd.DataFrame:
    return pd.DataFrame(adata_obj.obsm["spatial"], index=adata_obj.obs_names, columns=["x", "y"])


def checkpoint_path(config: dict, override: str | None) -> Path:
    output = Path(override or config["output"]["dir"])
    return output / ".runtime" / "notebook_state.pkl"


def load_state(path: Path) -> dict:
    with path.open("rb") as handle:
        return pickle.load(handle)


def save_state(path: Path, state: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        pickle.dump(state, handle, protocol=pickle.HIGHEST_PROTOCOL)


def raw_level2_mapping(raw_parent: ad.AnnData, config: dict, parent: str):
    mapping_config = config["raw_level2_mapping"]
    level1 = mapping_config["level1_column"]
    reference_context = ad.read_h5ad(ROOT / mapping_config["reference_h5ad"], backed="r")
    reference_mask = (
        reference_context.obs[level1]
        .astype(str)
        .str.replace("/", "_", regex=False)
        .eq(parent)
    )
    reference_filter = mapping_config.get("reference_filter")
    if reference_filter:
        reference_mask &= reference_context.obs[reference_filter["column"]].astype(str).eq(
            str(reference_filter["value"])
        )
    reference_parent = reference_context[reference_mask].to_memory()
    reference_context.file.close()
    pot = mapping_config.get("pot", {})
    tacco = mapping_config.get("tacco", {})
    return map_raw_level2_labels(
        raw_parent,
        reference_parent,
        parent_value=parent,
        method=mapping_config["method"],
        level1_col=level1,
        level2_col=mapping_config["level2_column"],
        reference_filter_column=(reference_filter or {}).get("column"),
        reference_filter_value=(reference_filter or {}).get("value"),
        pot_reg=pot.get("reg", 0.1),
        pot_reg_m=pot.get("reg_m", 0.0),
        pot_reg_type=pot.get("reg_type", "entropy"),
        tacco_multi_center=tacco.get("multi_center"),
        tacco_lamb=tacco.get("lamb"),
    )


def prepare_state(route: str, config: dict, output_root: Path) -> dict:
    raw_path = ROOT / config["context"]["h5ad"]
    level1 = config["context"]["level1_column"]
    raw = ad.read_h5ad(raw_path, backed="r")
    full_coordinates = coordinates(raw)
    full_level1 = raw.obs[level1].astype(str).copy()
    region = config["spatial_region"]
    tissue_um = full_coordinates * region["microns_per_coordinate"]
    xlim = tuple(tissue_um["x"].agg(["min", "max"]))
    ylim = tuple(tissue_um["y"].agg(["min", "max"]))
    anatomy = compute_anatomy_regions(
        full_coordinates=full_coordinates,
        full_level1_labels=full_level1,
        microns_per_coordinate=region["microns_per_coordinate"],
        candidate_window_sides_um=region["candidate_window_sides_um"],
        min_parent_units=region["min_parent_units"],
        cell_equivalent_um=8.0,
        **region["anatomy_region"],
    )
    write_anatomy_artifacts(output_root, anatomy)
    raw.file.close()
    return {
        "CONFIG": config,
        "OUTPUT_DIR": output_root,
        "CELL_EQUIVALENT_UM": 8.0,
        "LEVEL1": level1,
        "full_coordinates": full_coordinates,
        "full_level1": full_level1,
        "TISSUE_COORDINATES_UM": tissue_um,
        "TISSUE_XLIM": xlim,
        "TISSUE_YLIM": ylim,
        "ORIGIN_UM": (xlim[0], ylim[0]),
        "ANATOMY": anatomy,
        "RUNS": {},
        "route": route,
    }


def visium_scope_ids(config: dict, state: dict, scope: str, reconstructed_ids: pd.Index) -> pd.Index:
    limit = config["partition_change"]["sample_n_units"]
    use_full = False
    rng = np.random.default_rng(config["partition_change"]["random_state"])

    def sample(ids: pd.Index) -> pd.Index:
        return ids if use_full or len(ids) <= limit else pd.Index(rng.choice(ids.to_numpy(), size=limit, replace=False))

    global_ids = sample(reconstructed_ids)
    if scope == "global":
        return global_ids
    for parent, source_label in PARENT_SOURCE.items():
        candidate = reconstructed_ids[state["full_level1"].reindex(reconstructed_ids).eq(source_label)]
        parent_ids = sample(candidate)
        if parent == scope:
            return parent_ids
    raise KeyError(scope)


def compute_visium_scope(config: dict, state: dict, scope: str) -> None:
    comparison_config = config["partition_change"]["comparisons"][0]
    raw_path = ROOT / comparison_config["raw_h5ad"]
    recon_path = ROOT / comparison_config["reconstructed_spatial_h5ad"]
    raw_context = ad.read_h5ad(raw_path, backed="r")
    recon_context = ad.read_h5ad(recon_path, backed="r")
    ids = visium_scope_ids(config, state, scope, recon_context.obs_names.copy())
    raw_scope = raw_context[ids].to_memory()
    recon_scope = recon_context[ids].to_memory()
    resolution_mode = "level1_ari" if scope == "global" else "fixed_within_level1"
    partition = run_partition_analysis(
        raw_scope,
        recon_scope,
        level1_col=state["LEVEL1"],
        route_kind="sp_svc",
        resolution_mode=resolution_mode,
        resolution_candidates=config["partition_change"]["level1_resolution_candidates"],
        within_level1_resolution=config["partition_change"]["within_level1_resolution"],
        random_state=config["partition_change"]["random_state"],
        n_top_genes=config["partition_change"]["n_top_genes"],
        raw_qc_min_genes=config["partition_change"]["raw_qc_min_genes"],
        raw_qc_min_cells=config["partition_change"]["raw_qc_min_cells"],
    )
    if scope == "global":
        state["GLOBAL_PARTITION"] = partition
        write_partition_artifacts(state["OUTPUT_DIR"] / "global", partition)
    else:
        comparison = next(iter(partition.comparisons.values()))
        retained_ids = comparison.assignments.index
        raw_scope = raw_scope[retained_ids].copy()
        recon_scope = recon_scope[retained_ids].copy()
        assignments = comparison.assignments.set_index("unit_id")
        level2_mapping = raw_level2_mapping(raw_scope, config, scope)
        impact = compute_spatial_impact(
            full_coordinates=state["full_coordinates"],
            full_level1_labels=state["full_level1"],
            paired_coordinates=coordinates(raw_scope),
            raw_labels=assignments["raw_cluster"],
            raw_level2_labels=level2_mapping.labels,
            reconstructed_labels=assignments["recon_cluster"],
            unit_changed=assignments["unit_changed"],
            microns_per_coordinate=config["spatial_region"]["microns_per_coordinate"],
            candidate_window_sides_um=config["spatial_region"]["candidate_window_sides_um"],
            min_parent_units=config["spatial_region"]["min_parent_units"],
            rarefaction_draws=config["spatial_region"]["rarefaction_draws"],
            threshold_bootstraps=config["spatial_region"]["threshold_bootstraps"],
            cell_equivalent_um=8.0,
            anatomy_analysis=state["ANATOMY"],
            **config["spatial_region"]["anatomy_region"],
        )
        state["RUNS"][scope] = {
            "partition": partition,
            "comparison": comparison,
            "impact": impact,
            "raw_level2": level2_mapping,
            "available": int(state["full_level1"].reindex(recon_context.obs_names).eq(PARENT_SOURCE[scope]).sum()),
            "sampled": len(ids),
            "used": raw_scope.n_obs,
        }
        write_partition_artifacts(state["OUTPUT_DIR"] / scope, partition)
        write_raw_level2_artifacts(state["OUTPUT_DIR"] / scope, level2_mapping)
        write_spatial_artifacts(state["OUTPUT_DIR"] / scope, impact)
    raw_context.file.close()
    recon_context.file.close()


def compute_xenium_parent(config: dict, state: dict, parent: str) -> None:
    comparison_config = next(item for item in config["partition_change"]["comparisons"] if item["name"] == parent)
    raw_context = ad.read_h5ad(ROOT / comparison_config["raw_h5ad"], backed="r")
    recon_parent = ad.read_h5ad(ROOT / comparison_config["reconstructed_spatial_h5ad"])
    parent_ids, cohort_audit = select_raw_level1_parent_cohort(
        raw_context.obs[state["LEVEL1"]], recon_parent.obs_names, parent_value=parent
    )
    recon_parent = recon_parent[parent_ids].copy()
    raw_parent = raw_context[parent_ids].to_memory()
    partition = run_partition_analysis(
        raw_parent,
        recon_parent,
        level1_col=state["LEVEL1"],
        final_cluster_key=comparison_config["reconstructed_cluster_key"],
        route_kind="sc_svc",
        resolution_mode="fixed_within_level1",
        within_level1_resolution=config["partition_change"]["within_level1_resolution"],
        random_state=config["partition_change"]["random_state"],
        n_top_genes=config["partition_change"]["n_top_genes"],
    )
    comparison = next(iter(partition.comparisons.values()))
    assignments = comparison.assignments.set_index("unit_id")
    level2_mapping = raw_level2_mapping(raw_parent, config, parent)
    impact = compute_spatial_impact(
        full_coordinates=state["full_coordinates"],
        full_level1_labels=state["full_level1"],
        paired_coordinates=coordinates(raw_parent),
        raw_labels=assignments["raw_cluster"],
        raw_level2_labels=level2_mapping.labels,
        reconstructed_labels=assignments["recon_cluster"],
        unit_changed=assignments["unit_changed"],
        microns_per_coordinate=config["spatial_region"]["microns_per_coordinate"],
        candidate_window_sides_um=config["spatial_region"]["candidate_window_sides_um"],
        min_parent_units=config["spatial_region"]["min_parent_units"],
        rarefaction_draws=config["spatial_region"]["rarefaction_draws"],
        threshold_bootstraps=config["spatial_region"]["threshold_bootstraps"],
        cell_equivalent_um=8.0,
        anatomy_analysis=state["ANATOMY"],
        **config["spatial_region"]["anatomy_region"],
    )
    state["RUNS"][parent] = {
        "partition": partition,
        "comparison": comparison,
        "impact": impact,
        "raw_level2": level2_mapping,
        "cohort_audit": cohort_audit,
        "available": cohort_audit["carrier_units"],
        "used": raw_parent.n_obs,
        "expression_h5ad": comparison_config["expression_h5ad"],
    }
    write_partition_artifacts(state["OUTPUT_DIR"] / parent, partition)
    write_raw_level2_artifacts(state["OUTPUT_DIR"] / parent, level2_mapping)
    write_spatial_artifacts(state["OUTPUT_DIR"] / parent, impact)
    raw_context.file.close()


def finalize(route: str, config: dict, state: dict) -> None:
    output = state["OUTPUT_DIR"]
    if route == "visiumhd":
        comparison = config["partition_change"]["comparisons"][0]
        state["manifest"] = {
            "analysis_contract_version": 4,
            "route": "sp_svc",
            "sampling": {
                "global_input": state["GLOBAL_PARTITION"].audit["input_units"],
                "global_raw_qc_retained": state["GLOBAL_PARTITION"].audit["n_units"],
                "full_switch": False,
            },
            "partition_preprocessing": {
                "raw_qc_min_genes": config["partition_change"]["raw_qc_min_genes"],
                "raw_qc_min_cells": config["partition_change"]["raw_qc_min_cells"],
                "feature_selection": config["partition_change"]["feature_selection"],
                "leiden_backend": "igraph",
                "leiden_random_state": config["partition_change"]["random_state"],
            },
            "anatomy_scale": state["ANATOMY"].scale_audit,
            "inputs": {
                "raw": {"path": comparison["raw_h5ad"], "sha256": file_sha256(ROOT / comparison["raw_h5ad"])},
                "reconstructed": {"path": comparison["reconstructed_spatial_h5ad"], "sha256": file_sha256(ROOT / comparison["reconstructed_spatial_h5ad"])},
                "level2_reference": {
                    "path": config["raw_level2_mapping"]["reference_h5ad"],
                    "sha256": file_sha256(ROOT / config["raw_level2_mapping"]["reference_h5ad"]),
                },
            },
            "raw_level2_mapping": {
                name: run["raw_level2"].audit for name, run in state["RUNS"].items()
            },
        }
        state["input_audit"] = pd.DataFrame([
            {"role": "Raw full Level1 context", "units": len(state["full_level1"]), "paired": True},
            {"role": "Reconstructed spatial carrier", "units": len(state["full_coordinates"]), "paired": True},
        ])
    else:
        raw_path = ROOT / config["context"]["h5ad"]
        state["manifest"] = {
            "analysis_contract_version": 4,
            "route": "sc_svc",
            "anatomy_scale": state["ANATOMY"].scale_audit,
            "inputs": {"raw": {"path": config["context"]["h5ad"], "sha256": file_sha256(raw_path)}},
            "level2_reference": {
                "path": config["raw_level2_mapping"]["reference_h5ad"],
                "sha256": file_sha256(ROOT / config["raw_level2_mapping"]["reference_h5ad"]),
            },
            "parents": {
                name: {
                    "expression_h5ad": run["expression_h5ad"],
                    "spatial_expression_identical": run["partition"].representation_audit.get("spatial_expression_identical"),
                    "raw_level2_mapping": run["raw_level2"].audit,
                    "cohort_audit": run["cohort_audit"],
                }
                for name, run in state["RUNS"].items()
            },
        }
        state["input_audit"] = pd.DataFrame([
            {"role": "Raw full Level1 context", "units": len(state["full_level1"]), "paired": True},
            *[
                {
                    "role": f"{name} spatial carrier after Raw Level1 audit",
                    "units": run["used"],
                    "excluded_units": run["cohort_audit"]["excluded_raw_level1_mismatch"],
                    "paired": True,
                }
                for name, run in state["RUNS"].items()
            ],
            {
                "role": "Raw Level2 single-cell reference",
                "units": sum(run["raw_level2"].audit["n_reference_cells"] for run in state["RUNS"].values()),
                "paired": False,
            },
        ])
    write_analysis_artifacts(output, config=config, manifest=state["manifest"], input_audit=state["input_audit"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("route", choices=ROUTES)
    parser.add_argument("stage", choices=["anatomy", "global", "Fibroblast", "Mono_Macro", "T", "finalize"])
    parser.add_argument("--output-root")
    args = parser.parse_args()
    config = load_reconstruction_impact_config(ROOT / ROUTES[args.route])
    output = Path(args.output_root or config["output"]["dir"])
    if not output.is_absolute():
        output = ROOT / output
    output = output.resolve()
    checkpoint = checkpoint_path(config, str(output))
    if args.stage == "anatomy":
        state = prepare_state(args.route, config, output)
    else:
        state = load_state(checkpoint)
        state["OUTPUT_DIR"] = output
        if args.route == "visiumhd" and args.stage in {"global", *PARENT_SOURCE}:
            compute_visium_scope(config, state, args.stage)
        elif args.route == "xenium" and args.stage in PARENT_SOURCE:
            compute_xenium_parent(config, state, args.stage)
        elif args.stage == "finalize":
            finalize(args.route, config, state)
        else:
            raise ValueError(f"Stage {args.stage!r} is not valid for {args.route}")
    save_state(checkpoint, state)
    print(checkpoint)


if __name__ == "__main__":
    main()
