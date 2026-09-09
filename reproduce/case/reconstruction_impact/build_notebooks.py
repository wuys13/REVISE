"""Build the two executable reconstruction-impact route notebooks.

The analysis algorithms live in :mod:`revise.analysis`; this small renderer only
keeps the two notebooks structurally identical and gives each result its own
evidence cell.
"""

from __future__ import annotations

from pathlib import Path

import nbformat as nbf


ROOT = Path(__file__).resolve().parents[3]
NOTEBOOK_DIR = Path(__file__).resolve().parent


HELPERS = r'''import os
import sys
import warnings
from pathlib import Path

os.environ.setdefault("KMP_WARNINGS", "0")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("KMP_USE_SHM", "0")
os.environ.setdefault("KMP_AFFINITY", "disabled")
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message="If reg_type = entropy, then the matrix c is overwritten")
warnings.filterwarnings("ignore", message="Changing the sparsity structure of a csr_matrix is expensive.*")

PROJECT_ROOT = Path(os.environ.get("REVISE_REPOSITORY_ROOT", Path.cwd())).resolve()
if not (PROJECT_ROOT / "configs").is_dir():
    PROJECT_ROOT = Path.cwd().resolve()
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap, TwoSlopeNorm
import numpy as np
import pandas as pd
import seaborn as sns
from IPython.display import Markdown, display

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

ROOT = PROJECT_ROOT

ANATOMY_ORDER = ["Tumor", "Normal", "Interface", "Other"]
ANATOMY_COLORS = {"Tumor": "#d73027", "Normal": "#4575b4", "Interface": "#984ea3", "Other": "#d9d9d9"}

def coordinates(adata_obj):
    return pd.DataFrame(adata_obj.obsm["spatial"], index=adata_obj.obs_names, columns=["x", "y"])

def compact_iqr(values):
    values = pd.Series(values).dropna()
    if values.empty:
        return "NA"
    return f"{values.median():.3g} [{values.quantile(.25):.3g}, {values.quantile(.75):.3g}]"

def save_figure(fig, name):
    figure_dir = OUTPUT_DIR / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_dir / f"{name}.png", dpi=180, bbox_inches="tight")

def set_spatial_axes(ax):
    ax.set(xlim=TISSUE_XLIM, ylim=TISSUE_YLIM, xlabel="x (µm)", ylabel="y (µm)")
    ax.set_aspect("equal", adjustable="box")

def window_field(ax, frame, value, *, side_um, cmap="viridis", norm=None, title="", categorical=False):
    valid = frame.dropna(subset=[value])
    if valid.empty:
        set_spatial_axes(ax); ax.set_title(title + " (no valid windows)"); return None
    x_indices = np.arange(valid["window_x_index"].min(), valid["window_x_index"].max() + 1)
    y_indices = np.arange(valid["window_y_index"].min(), valid["window_y_index"].max() + 1)
    grid = np.full((len(y_indices), len(x_indices)), np.nan)
    for row in valid.itertuples():
        grid[int(row.window_y_index - y_indices[0]), int(row.window_x_index - x_indices[0])] = getattr(row, value)
    x_edges = ORIGIN_UM[0] + np.arange(x_indices[0], x_indices[-1] + 2) * side_um
    y_edges = ORIGIN_UM[1] + np.arange(y_indices[0], y_indices[-1] + 2) * side_um
    image = ax.pcolormesh(x_edges, y_edges, grid, shading="flat", cmap=cmap, norm=norm)
    set_spatial_axes(ax); ax.set_title(title)
    return image

def anatomy_map(ax, focus=None):
    regions = ANATOMY.anatomy_windows.copy()
    codes = regions["level1_region"].map({name: i for i, name in enumerate(ANATOMY_ORDER)})
    if focus is not None:
        codes = np.where(regions["level1_region"].eq(focus), ANATOMY_ORDER.index(focus), ANATOMY_ORDER.index("Other"))
    frame = regions.assign(anatomy_code=codes)
    cmap = ListedColormap([ANATOMY_COLORS[name] for name in ANATOMY_ORDER])
    image = window_field(ax, frame, "anatomy_code", side_um=float(ANATOMY.scale_audit["main_window_side_um"]), cmap=cmap, norm=BoundaryNorm(np.arange(-.5, 4.5), 4), title=focus or "Tumor / Normal / Interface / Other")
    if image is not None and focus is None:
        colorbar = plt.colorbar(image, ax=ax, ticks=range(4), fraction=.046, pad=.04)
        colorbar.ax.set_yticklabels(ANATOMY_ORDER)

def plot_support_curve(impact, title):
    table = impact.support_sensitivity
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.4))
    cells = table["window_side_length"] / CELL_EQUIVALENT_UM
    axes[0].plot(cells, table["retained_parent_unit_fraction"], marker="o")
    axes[0].axvline(impact.scale_audit["main_window_cells_per_side"], color="black", ls="--")
    axes[0].set(xlabel="cells per side", ylabel="retained-unit fraction", title="retention")
    axes[1].plot(cells, table["valid_window_fraction"], marker="o")
    axes[1].axvline(impact.scale_audit["main_window_cells_per_side"], color="black", ls="--")
    axes[1].set(xlabel="cells per side", ylabel="valid-window fraction", title="valid windows")
    fig.suptitle(title); fig.tight_layout(); return fig

def plot_contingency(comparisons, title):
    fig, axes = plt.subplots(1, len(comparisons), figsize=(4.5 * len(comparisons), 3.8), squeeze=False)
    for ax, (scope, comparison) in zip(axes[0], comparisons.items()):
        table = comparison.contingency.div(comparison.contingency.sum(axis=1), axis=0).fillna(0)
        sns.heatmap(table, ax=ax, cmap="mako", vmin=0, vmax=1, cbar=True)
        ax.set(title=scope, xlabel="Reconstructed / Final", ylabel="Raw")
    fig.suptitle(title); fig.tight_layout(); return fig

def plot_changed_units(runs, title):
    fig, axes = plt.subplots(1, len(runs), figsize=(4.5 * len(runs), 4), squeeze=False)
    for ax, (scope, run) in zip(axes[0], runs.items()):
        frame = run["impact"].unit_assignments
        colors = np.where(frame["unit_changed"], "#d73027", "#d9d9d9")
        ax.scatter(frame["x"], frame["y"], c=colors, s=1, linewidths=0, rasterized=True)
        set_spatial_axes(ax); ax.set_title(scope)
    fig.suptitle(title); fig.tight_layout(); return fig

def plot_cluster_pair(run, title):
    frame = run["impact"].unit_assignments
    mapping = run["comparison"].mapping.set_index("recon_cluster")["raw_cluster"].to_dict()
    raw_values = sorted(frame["raw_cluster"].unique())
    palette = {value: color for value, color in zip(raw_values, sns.color_palette("tab20", len(raw_values)))}
    recon_colors = [palette.get(mapping.get(value), "#111111") for value in frame["reconstructed_cluster"]]
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    axes[0].scatter(frame["x"], frame["y"], c=frame["raw_cluster"].map(palette), s=1, linewidths=0, rasterized=True)
    axes[1].scatter(frame["x"], frame["y"], c=recon_colors, s=1, linewidths=0, rasterized=True)
    for ax, label in zip(axes, ["Raw Leiden", "Reconstructed / Final"]):
        set_spatial_axes(ax); ax.set_title(label)
    fig.suptitle(title + " — matched clusters share colors"); fig.tight_layout(); return fig

def plot_reconstructed_clusters(run, title):
    frame = run["impact"].unit_assignments
    values = sorted(frame["reconstructed_cluster"].unique())
    palette = {value: color for value, color in zip(values, sns.color_palette("tab20", len(values)))}
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    ax.scatter(frame["x"], frame["y"], c=frame["reconstructed_cluster"].map(palette), s=1, linewidths=0, rasterized=True)
    set_spatial_axes(ax); ax.set_title(title); fig.tight_layout(); return fig

def plot_metric_comparison(run, metric, title):
    metrics = run["impact"].window_metrics.loc[lambda x: x["valid_window"]].copy()
    side = float(run["impact"].scale_audit["main_window_side_um"])
    state_columns = [f"{metric}_raw", f"{metric}_recon", f"{metric}_level2"]
    state_max = max(float(metrics[state_columns].max().max()), 1.0)
    state_min = 0.0 if metric == "evenness" else 1.0
    state_max = 1.0 if metric == "evenness" else state_max
    delta_columns = [f"delta_{metric}_vs_raw_leiden", f"delta_{metric}_vs_raw_level2"]
    delta_max = max(float(metrics[delta_columns].abs().max().max()), 0.01)
    state_norm = plt.Normalize(state_min, state_max)
    delta_norm = TwoSlopeNorm(vcenter=0, vmin=-delta_max, vmax=delta_max)
    panels = [
        (f"{metric}_raw", "Raw Leiden", "viridis", state_norm),
        (f"{metric}_recon", "Reconstructed", "viridis", state_norm),
        (f"delta_{metric}_vs_raw_leiden", "Recon − Raw Leiden", "coolwarm", delta_norm),
        (f"{metric}_level2", "Raw Level2", "viridis", state_norm),
        (f"{metric}_recon", "Reconstructed", "viridis", state_norm),
        (f"delta_{metric}_vs_raw_level2", "Recon − Raw Level2", "coolwarm", delta_norm),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.2), sharex=True, sharey=True)
    for ax, (column, label, cmap, norm) in zip(axes.flat, panels):
        image = window_field(ax, metrics, column, side_um=side, cmap=cmap, norm=norm, title=label)
        if image is not None:
            plt.colorbar(image, ax=ax, fraction=.046, pad=.04)
    fig.suptitle(title); fig.tight_layout(); return fig

def plot_high_diversity(run, title):
    metrics = run["impact"].window_metrics.copy()
    side = float(run["impact"].scale_audit["main_window_side_um"])
    valid = metrics.loc[metrics["valid_window"]].copy()
    maximum = max(float(valid["neff_recon"].max()), 1.0)
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), sharex=True, sharey=True)
    image = window_field(axes[0], valid, "neff_recon", side_um=side, cmap="viridis", norm=plt.Normalize(1, maximum), title="Reconstructed Neff")
    if image is not None:
        plt.colorbar(image, ax=axes[0], fraction=.046, pad=.04)
    threshold_ok = run["impact"].state_threshold.get("status") == "ok"
    mask = valid.assign(region_mask=valid["in_state_region"].fillna(False).astype(float))
    if threshold_ok:
        image = window_field(axes[1], mask, "region_mask", side_um=side, cmap=ListedColormap(["#efefef", "#fdae61"]), norm=BoundaryNorm([-0.5, 0.5, 1.5], 2), title="High-diversity Region")
        if image is not None:
            colorbar = plt.colorbar(image, ax=axes[1], ticks=[0, 1], fraction=.046, pad=.04)
            colorbar.ax.set_yticklabels(["outside", "inside"])
    else:
        window_field(axes[1], mask, "region_mask", side_um=side, cmap=ListedColormap(["#efefef"]), norm=BoundaryNorm([-0.5, 0.5], 1), title="No stable Region threshold")
    fig.suptitle(title); fig.tight_layout(); return fig

def raw_level2_mapping(raw_parent, parent):
    mapping_config = CONFIG["raw_level2_mapping"]
    level1 = mapping_config["level1_column"]
    reference_context = ad.read_h5ad(ROOT / mapping_config["reference_h5ad"], backed="r")
    reference_mask = reference_context.obs[level1].astype(str).str.replace("/", "_", regex=False).eq(parent)
    reference_filter = mapping_config.get("reference_filter")
    if reference_filter:
        reference_mask &= reference_context.obs[reference_filter["column"]].astype(str).eq(str(reference_filter["value"]))
    reference_parent = reference_context[reference_mask].to_memory()
    reference_context.file.close()
    pot = mapping_config.get("pot", {}); tacco = mapping_config.get("tacco", {})
    return map_raw_level2_labels(raw_parent, reference_parent, parent_value=parent, method=mapping_config["method"], level1_col=level1, level2_col=mapping_config["level2_column"], reference_filter_column=(reference_filter or {}).get("column"), reference_filter_value=(reference_filter or {}).get("value"), pot_reg=pot.get("reg", .1), pot_reg_m=pot.get("reg_m", 0.), pot_reg_type=pot.get("reg_type", "entropy"), tacco_multi_center=tacco.get("multi_center"), tacco_lamb=tacco.get("lamb"))

def extent_table(impact, column):
    source = impact.region_extent_by_anatomy if column == "in_state_region" else impact.gain_region_extent_by_anatomy
    return source.loc[source["level1_region"].isin(["Overall", "Tumor", "Normal", "Interface"])]
'''


def setup_cell(kind: str) -> str:
    config = "reconstruction_impact_visiumhd_p1crc.yaml" if kind == "visium" else "reconstruction_impact_xenium_p2crc_fibroblast.yaml"
    if kind == "visium":
        return with_cached_state(f'''CONFIG = load_reconstruction_impact_config(ROOT / "configs" / "analysis" / "{config}")
OUTPUT_DIR = Path(os.environ.get("REVISE_ANALYSIS_OUTPUT_ROOT", CONFIG["output"]["dir"]))
USE_FULL_VISIUMHD_COHORT = False
VISIUMHD_SAMPLE_N_UNITS = 30_000
CELL_EQUIVALENT_UM = 8.0
LEVEL1 = CONFIG["context"]["level1_column"]
comparison_config = CONFIG["partition_change"]["comparisons"][0]
raw_context = ad.read_h5ad(ROOT / comparison_config["raw_h5ad"], backed="r")
recon_context = ad.read_h5ad(ROOT / comparison_config["reconstructed_spatial_h5ad"], backed="r")
full_coordinates = coordinates(raw_context)
full_level1 = raw_context.obs[LEVEL1].astype(str).copy()
TISSUE_COORDINATES_UM = full_coordinates * CONFIG["spatial_region"]["microns_per_coordinate"]
TISSUE_XLIM = tuple(TISSUE_COORDINATES_UM["x"].agg(["min", "max"]))
TISSUE_YLIM = tuple(TISSUE_COORDINATES_UM["y"].agg(["min", "max"]))
ORIGIN_UM = (TISSUE_XLIM[0], TISSUE_YLIM[0])
ANATOMY = compute_anatomy_regions(full_coordinates=full_coordinates, full_level1_labels=full_level1, microns_per_coordinate=CONFIG["spatial_region"]["microns_per_coordinate"], candidate_window_sides_um=CONFIG["spatial_region"]["candidate_window_sides_um"], min_parent_units=CONFIG["spatial_region"]["min_parent_units"], cell_equivalent_um=CELL_EQUIVALENT_UM, **CONFIG["spatial_region"]["anatomy_region"])
write_anatomy_artifacts(OUTPUT_DIR, ANATOMY)
reconstructed_ids = recon_context.obs_names.copy()
if not set(reconstructed_ids) <= set(raw_context.obs_names): raise ValueError("Every reconstructed observation must occur in Raw.")
rng = np.random.default_rng(CONFIG["partition_change"]["random_state"])
def deterministic_same_id_sample(ids, limit):
    ids = pd.Index(ids)
    return ids if USE_FULL_VISIUMHD_COHORT or len(ids) <= limit else pd.Index(rng.choice(ids.to_numpy(), size=limit, replace=False))
global_ids = deterministic_same_id_sample(reconstructed_ids, VISIUMHD_SAMPLE_N_UNITS)
raw_global = raw_context[global_ids].to_memory(); recon_global = recon_context[global_ids].to_memory()
GLOBAL_PARTITION = run_partition_analysis(raw_global, recon_global, level1_col=LEVEL1, route_kind="sp_svc", resolution_mode="level1_ari", resolution_candidates=CONFIG["partition_change"]["level1_resolution_candidates"], within_level1_resolution=CONFIG["partition_change"]["within_level1_resolution"], random_state=CONFIG["partition_change"]["random_state"], n_top_genes=CONFIG["partition_change"]["n_top_genes"], raw_qc_min_genes=CONFIG["partition_change"]["raw_qc_min_genes"], raw_qc_min_cells=CONFIG["partition_change"]["raw_qc_min_cells"])
PARENT_SOURCE = {{"Fibroblast": "Fibroblast", "Mono_Macro": "Mono/Macro", "T": "T"}}
RUNS = {{}}
for parent, source_label in PARENT_SOURCE.items():
    parent_ids = deterministic_same_id_sample(reconstructed_ids[full_level1.reindex(reconstructed_ids).eq(source_label)], VISIUMHD_SAMPLE_N_UNITS)
    raw_parent = raw_context[parent_ids].to_memory(); recon_parent = recon_context[parent_ids].to_memory()
    partition = run_partition_analysis(raw_parent, recon_parent, level1_col=LEVEL1, route_kind="sp_svc", resolution_mode="fixed_within_level1", within_level1_resolution=CONFIG["partition_change"]["within_level1_resolution"], random_state=CONFIG["partition_change"]["random_state"], n_top_genes=CONFIG["partition_change"]["n_top_genes"], raw_qc_min_genes=CONFIG["partition_change"]["raw_qc_min_genes"], raw_qc_min_cells=CONFIG["partition_change"]["raw_qc_min_cells"])
    comparison = next(iter(partition.comparisons.values())); retained_ids = comparison.assignments.index; raw_parent = raw_parent[retained_ids].copy(); recon_parent = recon_parent[retained_ids].copy(); assignments = comparison.assignments.set_index("unit_id")
    level2_mapping = raw_level2_mapping(raw_parent, parent)
    impact = compute_spatial_impact(full_coordinates=full_coordinates, full_level1_labels=full_level1, paired_coordinates=coordinates(raw_parent), raw_labels=assignments["raw_cluster"], raw_level2_labels=level2_mapping.labels, reconstructed_labels=assignments["recon_cluster"], unit_changed=assignments["unit_changed"], microns_per_coordinate=CONFIG["spatial_region"]["microns_per_coordinate"], candidate_window_sides_um=CONFIG["spatial_region"]["candidate_window_sides_um"], min_parent_units=CONFIG["spatial_region"]["min_parent_units"], rarefaction_draws=CONFIG["spatial_region"]["rarefaction_draws"], threshold_bootstraps=CONFIG["spatial_region"]["threshold_bootstraps"], cell_equivalent_um=CELL_EQUIVALENT_UM, anatomy_analysis=ANATOMY, **CONFIG["spatial_region"]["anatomy_region"])
    RUNS[parent] = {{"partition": partition, "comparison": comparison, "impact": impact, "raw_level2": level2_mapping, "available": int(full_level1.reindex(reconstructed_ids).eq(source_label).sum()), "sampled": len(parent_ids), "used": raw_parent.n_obs}}
for name, run in RUNS.items(): write_partition_artifacts(OUTPUT_DIR / name, run["partition"]); write_raw_level2_artifacts(OUTPUT_DIR / name, run["raw_level2"]); write_spatial_artifacts(OUTPUT_DIR / name, run["impact"])
manifest = {{"analysis_contract_version": 4, "route": "sp_svc", "sampling": {{"global_input": len(global_ids), "global_raw_qc_retained": GLOBAL_PARTITION.audit["n_units"], "full_switch": USE_FULL_VISIUMHD_COHORT}}, "partition_preprocessing": {{"raw_qc_min_genes": CONFIG["partition_change"]["raw_qc_min_genes"], "raw_qc_min_cells": CONFIG["partition_change"]["raw_qc_min_cells"], "feature_selection": CONFIG["partition_change"]["feature_selection"], "leiden_backend": "igraph", "leiden_random_state": CONFIG["partition_change"]["random_state"]}}, "anatomy_scale": ANATOMY.scale_audit, "raw_level2_mapping": {{name: run["raw_level2"].audit for name, run in RUNS.items()}}, "inputs": {{"raw": {{"path": comparison_config["raw_h5ad"], "sha256": file_sha256(ROOT / comparison_config["raw_h5ad"])}}, "reconstructed": {{"path": comparison_config["reconstructed_spatial_h5ad"], "sha256": file_sha256(ROOT / comparison_config["reconstructed_spatial_h5ad"])}}, "level2_reference": {{"path": CONFIG["raw_level2_mapping"]["reference_h5ad"], "sha256": file_sha256(ROOT / CONFIG["raw_level2_mapping"]["reference_h5ad"])}}}}}}
input_audit = pd.DataFrame([{{"role": "Raw full Level1 context", "units": raw_context.n_obs, "paired": True}}, {{"role": "Reconstructed spatial carrier", "units": recon_context.n_obs, "paired": True}}])
write_analysis_artifacts(OUTPUT_DIR, config=CONFIG, manifest=manifest, input_audit=input_audit)
raw_context.file.close(); recon_context.file.close()''')
    return with_cached_state(f'''CONFIG = load_reconstruction_impact_config(ROOT / "configs" / "analysis" / "{config}")
OUTPUT_DIR = Path(os.environ.get("REVISE_ANALYSIS_OUTPUT_ROOT", CONFIG["output"]["dir"]))
CELL_EQUIVALENT_UM = 8.0
LEVEL1 = CONFIG["context"]["level1_column"]
raw_path = ROOT / CONFIG["context"]["h5ad"]
raw_context = ad.read_h5ad(raw_path, backed="r")
full_coordinates = coordinates(raw_context); full_level1 = raw_context.obs[LEVEL1].astype(str).copy()
TISSUE_COORDINATES_UM = full_coordinates * CONFIG["spatial_region"]["microns_per_coordinate"]
TISSUE_XLIM = tuple(TISSUE_COORDINATES_UM["x"].agg(["min", "max"])); TISSUE_YLIM = tuple(TISSUE_COORDINATES_UM["y"].agg(["min", "max"])); ORIGIN_UM = (TISSUE_XLIM[0], TISSUE_YLIM[0])
ANATOMY = compute_anatomy_regions(full_coordinates=full_coordinates, full_level1_labels=full_level1, microns_per_coordinate=CONFIG["spatial_region"]["microns_per_coordinate"], candidate_window_sides_um=CONFIG["spatial_region"]["candidate_window_sides_um"], min_parent_units=CONFIG["spatial_region"]["min_parent_units"], cell_equivalent_um=CELL_EQUIVALENT_UM, **CONFIG["spatial_region"]["anatomy_region"])
write_anatomy_artifacts(OUTPUT_DIR, ANATOMY)
RUNS = {{}}
for comparison_config in CONFIG["partition_change"]["comparisons"]:
    parent = comparison_config["name"]
    recon_parent = ad.read_h5ad(ROOT / comparison_config["reconstructed_spatial_h5ad"])
    parent_ids, cohort_audit = select_raw_level1_parent_cohort(raw_context.obs[LEVEL1], recon_parent.obs_names, parent_value=parent)
    recon_parent = recon_parent[parent_ids].copy()
    raw_parent = raw_context[parent_ids].to_memory()
    partition = run_partition_analysis(raw_parent, recon_parent, level1_col=LEVEL1, final_cluster_key=comparison_config["reconstructed_cluster_key"], route_kind="sc_svc", resolution_mode="fixed_within_level1", within_level1_resolution=CONFIG["partition_change"]["within_level1_resolution"], random_state=CONFIG["partition_change"]["random_state"], n_top_genes=CONFIG["partition_change"]["n_top_genes"])
    comparison = next(iter(partition.comparisons.values())); assignments = comparison.assignments.set_index("unit_id")
    level2_mapping = raw_level2_mapping(raw_parent, parent)
    impact = compute_spatial_impact(full_coordinates=full_coordinates, full_level1_labels=full_level1, paired_coordinates=coordinates(raw_parent), raw_labels=assignments["raw_cluster"], raw_level2_labels=level2_mapping.labels, reconstructed_labels=assignments["recon_cluster"], unit_changed=assignments["unit_changed"], microns_per_coordinate=CONFIG["spatial_region"]["microns_per_coordinate"], candidate_window_sides_um=CONFIG["spatial_region"]["candidate_window_sides_um"], min_parent_units=CONFIG["spatial_region"]["min_parent_units"], rarefaction_draws=CONFIG["spatial_region"]["rarefaction_draws"], threshold_bootstraps=CONFIG["spatial_region"]["threshold_bootstraps"], cell_equivalent_um=CELL_EQUIVALENT_UM, anatomy_analysis=ANATOMY, **CONFIG["spatial_region"]["anatomy_region"])
    RUNS[parent] = {{"partition": partition, "comparison": comparison, "impact": impact, "raw_level2": level2_mapping, "cohort_audit": cohort_audit, "available": cohort_audit["carrier_units"], "used": raw_parent.n_obs, "expression_h5ad": comparison_config["expression_h5ad"]}}
for name, run in RUNS.items(): write_partition_artifacts(OUTPUT_DIR / name, run["partition"]); write_raw_level2_artifacts(OUTPUT_DIR / name, run["raw_level2"]); write_spatial_artifacts(OUTPUT_DIR / name, run["impact"])
manifest = {{"analysis_contract_version": 4, "route": "sc_svc", "anatomy_scale": ANATOMY.scale_audit, "inputs": {{"raw": {{"path": CONFIG["context"]["h5ad"], "sha256": file_sha256(raw_path)}}, "level2_reference": {{"path": CONFIG["raw_level2_mapping"]["reference_h5ad"], "sha256": file_sha256(ROOT / CONFIG["raw_level2_mapping"]["reference_h5ad"])}}}}, "parents": {{name: {{"expression_h5ad": run["expression_h5ad"], "spatial_expression_identical": run["partition"].representation_audit.get("spatial_expression_identical"), "raw_level2_mapping": run["raw_level2"].audit, "cohort_audit": run["cohort_audit"]}} for name, run in RUNS.items()}}}}
input_audit = pd.DataFrame([{{"role": "Raw full Level1 context", "units": raw_context.n_obs, "paired": True}}, *[{{"role": f"{{name}} spatial carrier after Raw Level1 audit", "units": run["used"], "excluded_units": run["cohort_audit"]["excluded_raw_level1_mismatch"], "paired": True}} for name, run in RUNS.items()], {{"role": "Raw Level2 single-cell reference", "units": sum(run["raw_level2"].audit["n_reference_cells"] for run in RUNS.values()), "paired": False}}])
write_analysis_artifacts(OUTPUT_DIR, config=CONFIG, manifest=manifest, input_audit=input_audit)
raw_context.file.close()''')


def with_cached_state(body: str) -> str:
    """Allow the constrained execution runner to restore a real workflow checkpoint."""
    indented = "\n".join(f"    {line}" if line else "" for line in body.splitlines())
    return '''CACHE_PATH = os.environ.get("RECONSTRUCTION_IMPACT_CACHE_PATH")
if CACHE_PATH:
    import pickle
    with Path(CACHE_PATH).open("rb") as handle:
        globals().update(pickle.load(handle))
else:
''' + indented


def parent_cells(number: int, parent: str):
    return [
        ("markdown", f"## {number}. {parent} internal diversity and Regions\n\nThis parent is analysed after the fixed macro anatomy. Raw Level2 is mapped from the original Raw expression after Level1 parent selection; reconstructed-carrier Level2 labels are not used as this baseline."),
        ("markdown", "### Cohort and window decision\n\nThe cell-equivalent side is fixed at 8 µm. The selected side comes only from parent occupancy; support is four units, and each of 200 paired draws samples four units."),
        ("code", f'''RUN = RUNS["{parent}"]; IMPACT = RUN["impact"]
fig = plot_support_curve(IMPACT, "{parent}: occupancy-only window decision")
save_figure(fig, "{parent.lower()}_window_decision"); plt.show()
occupancy = IMPACT.unit_assignments.groupby("window_id").size()
sampled = RUN.get("sampled", RUN["used"])
decision = pd.DataFrame([{{"available units": RUN["available"], "sampled units": sampled, "Raw-QC retained units": RUN["used"], "sampling mode": "full" if RUN["available"] == sampled else "deterministic 30k", "seed": CONFIG["partition_change"]["random_state"], "window side (cells per side)": IMPACT.scale_audit["main_window_cells_per_side"], "window side (µm)": IMPACT.scale_audit["main_window_side_um"], "occupancy median [Q1,Q3]": compact_iqr(occupancy), "minimum support": IMPACT.scale_audit["min_parent_units"], "rarefaction units": IMPACT.scale_audit["min_parent_units"], "paired draws": IMPACT.scale_audit["rarefaction_draws"]}}])
display(decision)'''),
        ("markdown", "### Internal-state baselines and matched-K assignments\n\nThe uniform Raw Level1 label is only the unexpanded baseline (`Kobs = Neff = evenness = 1`). The two informative Raw baselines are (i) expression-derived matched-K Leiden and (ii) reference-derived Raw Level2. Both are compared with the same reconstructed assignment on identical windows and identical rarefaction draws."),
        ("code", f'''METRICS = RUN["impact"].window_metrics.copy(); METRICS.attrs["side_um"] = IMPACT.scale_audit["main_window_side_um"]
fig = plot_reconstructed_clusters(RUN, "{parent}: reconstructed internal state")
save_figure(fig, "{parent.lower()}_reconstructed_clusters"); plt.show()
display(pd.DataFrame([{{"baseline": "Raw Level1 (uniform parent)", "Kobs": 1.0, "Neff": 1.0, "evenness": 1.0}}, {{"baseline": "Raw Level2 mapping", "Kobs": RUN["raw_level2"].audit["n_mapped_level2"], "Neff": "window-specific", "evenness": "window-specific"}}]))
display(pd.DataFrame([RUN["raw_level2"].audit]).loc[:, ["source", "method", "n_raw_units", "n_reference_cells", "n_reference_level2", "n_mapped_level2", "n_missing"]])'''),
        ("markdown", "### Matched-K reconstruction-associated internal structure\n\nThis is a **parent-internal reassignment** analysis, not a Level1 identity change. Raw Leiden is tuned to approximately the same cluster count as the reconstructed/final assignment. Hungarian matching controls label permutation for unit-change metrics; the diversity matrices below use the original cluster compositions within each window."),
        ("code", f'''fig = plot_cluster_pair(RUN, "{parent}: matched-K cluster assignment")
save_figure(fig, "{parent.lower()}_matched_clusters"); plt.show()
row = RUN["comparison"].summary.iloc[0]
display(pd.DataFrame([{{"paired units": row.n_units, "Raw K": row.n_raw_clusters, "Recon K": row.n_recon_clusters, "parent-internal reassignment": row.st_unit_change_fraction, "balanced change": row.balanced_cluster_change, "ARI": row["ARI"]}}]).round(4))'''),
        ("markdown", "### Kobs: local subtype richness\n\n`Kobs` is the number of cluster labels observed in a rarefied window. Here it asks how many distinct internal states coexist locally, without accounting for whether one state dominates."),
        ("code", f'''fig = plot_metric_comparison(RUN, "k_obs", "{parent}: Kobs evidence matrix")
save_figure(fig, "{parent.lower()}_kobs_matrix"); plt.show()
valid = METRICS.loc[METRICS.valid_window]
display(pd.DataFrame([{{"baseline": "Raw Leiden", "state median [Q1,Q3]": compact_iqr(valid["k_obs_raw"]), "Recon minus baseline": compact_iqr(valid["delta_k_obs_vs_raw_leiden"])}}, {{"baseline": "Raw Level2", "state median [Q1,Q3]": compact_iqr(valid["k_obs_level2"]), "Recon minus baseline": compact_iqr(valid["delta_k_obs_vs_raw_level2"])}}]))
display(Markdown(f"Across valid windows, reconstructed Kobs is **{{compact_iqr(valid['k_obs_recon'])}}**. Its median difference is **{{valid['delta_k_obs_vs_raw_leiden'].median():.3g}}** versus Raw Leiden and **{{valid['delta_k_obs_vs_raw_level2'].median():.3g}}** versus Raw Level2."))'''),
        ("markdown", "### Neff: abundance-aware effective subtype count\n\n`Neff = exp(Shannon entropy)` converts subtype composition into an effective number of equally abundant clusters. It is the primary local-diversity state because it discounts rare labels and single-cluster dominance."),
        ("code", f'''fig = plot_metric_comparison(RUN, "neff", "{parent}: Neff evidence matrix")
save_figure(fig, "{parent.lower()}_neff_matrix"); plt.show()
valid = METRICS.loc[METRICS.valid_window]
display(pd.DataFrame([{{"baseline": "Raw Leiden", "state median [Q1,Q3]": compact_iqr(valid["neff_raw"]), "Recon minus baseline": compact_iqr(valid["delta_neff_vs_raw_leiden"])}}, {{"baseline": "Raw Level2", "state median [Q1,Q3]": compact_iqr(valid["neff_level2"]), "Recon minus baseline": compact_iqr(valid["delta_neff_vs_raw_level2"])}}]))
display(Markdown(f"Reconstructed Neff is **{{compact_iqr(valid['neff_recon'])}}**. The two delta maps separate reconstruction-associated differences from dependence on the chosen Raw baseline."))'''),
        ("markdown", "### Evenness: balance conditional on observed richness\n\n`evenness = Neff / Kobs` ranges from dominance toward balanced coexistence. It distinguishes windows with many labels but one dominant subtype from windows where those labels have comparable abundance."),
        ("code", f'''fig = plot_metric_comparison(RUN, "evenness", "{parent}: evenness evidence matrix")
save_figure(fig, "{parent.lower()}_evenness_matrix"); plt.show()
valid = METRICS.loc[METRICS.valid_window]
display(pd.DataFrame([{{"baseline": "Raw Leiden", "state median [Q1,Q3]": compact_iqr(valid["evenness_raw"]), "Recon minus baseline": compact_iqr(valid["delta_evenness_vs_raw_leiden"])}}, {{"baseline": "Raw Level2", "state median [Q1,Q3]": compact_iqr(valid["evenness_level2"]), "Recon minus baseline": compact_iqr(valid["delta_evenness_vs_raw_level2"])}}]))
display(Markdown(f"Reconstructed evenness is **{{compact_iqr(valid['evenness_recon'])}}**; interpret it together with Kobs, because a pure one-cluster window also has evenness 1."))'''),
        ("markdown", "### High-diversity Region\n\nThe Region is a parent-specific mask derived only from the data-driven breakpoint of reconstructed Neff. It is displayed beside the continuous Neff field and is not overlaid on anatomy."),
        ("code", f'''fig = plot_high_diversity(RUN, "{parent}: reconstructed diversity state and Region")
save_figure(fig, "{parent.lower()}_high_diversity_region"); plt.show()
extent = extent_table(IMPACT, "in_state_region").loc[:, ["level1_region", "region_windows", "valid_windows", "region_area_um2", "area_fraction", "unit_fraction"]]
display(extent.loc[extent.level1_region.isin(["Overall", "Tumor", "Normal", "Interface"])])
display(Markdown(f"The reconstructed-Neff breakpoint is **{{IMPACT.state_threshold.get('threshold')}}** (status: **{{IMPACT.state_threshold.get('status')}}**). The anatomy rows describe context after Region definition; anatomy does not enter the threshold."))'''),
    ]


def cells_for(kind: str):
    route = "VisiumHD sp-SVC" if kind == "visium" else "Xenium sc-SVC"
    route_detail = "Raw and reconstructed expression are paired; only Raw→reconstructed-expression is analysed." if kind == "visium" else "Raw Leiden is compared with final SVC clusters; the spatial carrier expression-identity audit prevents an artificial expression edge."
    cells = [
        ("markdown", f"# {route} Reconstruction Impact\n\nA route-specific, reproducible description of reconstruction-associated partition and local-composition changes."),
        ("markdown", f"## 1. Route semantics and carrier audit\n\n{route_detail} **Evidence boundary:** these are paired representation and spatial-pattern observations; they do not establish a mechanism, biological truth, or clinical meaning."),
        ("code", HELPERS),
        ("code", setup_cell(kind)),
        ("code", "display(input_audit)\ndisplay(pd.DataFrame([{\"8 µm cell-equivalent\": CELL_EQUIVALENT_UM, \"coordinate to µm\": CONFIG[\"spatial_region\"][\"microns_per_coordinate\"], \"window candidates (µm)\": str(CONFIG[\"spatial_region\"][\"candidate_window_sides_um\"])}]))"),
        ("markdown", "### Partition preprocessing contract\n\nFor sp-SVC, observation and gene QC are defined on Raw counts and applied to both paired carriers. Raw-derived Seurat-v3 HVGs are then shared by Raw and reconstructed expression. Leiden uses the explicitly recorded igraph backend and seed. This prevents reconstructed expression from defining the Raw feature space."),
        ("code", "if CONFIG[\"route_kind\"] == \"sp_svc\":\n    audit = GLOBAL_PARTITION.audit\n    display(pd.DataFrame([{\"sampled units\": audit[\"input_units\"], \"Raw-QC retained units\": audit[\"n_units\"], \"excluded units\": audit[\"excluded_raw_qc_units\"], \"retained shared genes\": audit[\"n_shared_genes\"], \"features\": audit[\"n_features\"], \"feature selection\": audit[\"feature_selection\"], \"Leiden backend\": audit[\"leiden_backend\"], \"seed\": audit[\"leiden_random_state\"]}]))\nelse:\n    first = next(iter(RUNS.values()))[\"partition\"]\n    display(pd.DataFrame([{\"carrier expression identity\": all(run[\"partition\"].representation_audit.get(\"spatial_expression_identical\", False) for run in RUNS.values()), \"feature selection\": first.audit[\"feature_selection\"], \"Leiden backend\": first.audit[\"leiden_backend\"], \"seed\": first.audit[\"leiden_random_state\"]}]))"),
        ("markdown", "### Local diversity metric definitions\n\nFor every valid parent window and every paired rarefaction draw: **Kobs** counts observed subtype/cluster labels; **Neff** is the abundance-aware effective cluster count; **evenness = Neff / Kobs** measures balance conditional on richness. The same sampled units are used for Raw Leiden, Raw Level2 and reconstructed assignments. These definitions apply identically to all three parent sections."),
        ("markdown", "## 2. Partition complexity diagnostic\n\nThis diagnostic retains the Raw resolution. It asks whether reconstruction yields a finer or reorganized partition; it is not the controlled assignment-change headline."),
        ("code", "complexity_rows = []\nfor name, run in RUNS.items():\n    row = next(iter(run[\"partition\"].complexity_comparisons.values())).summary.iloc[0]\n    complexity_rows.append({\"scope\": name, \"Raw K\": row.n_raw_clusters, \"Reconstructed K at Raw resolution\": row.n_recon_clusters, \"ARI\": row[\"ARI\"]})\nif 'GLOBAL_PARTITION' in globals():\n    row = next(iter(GLOBAL_PARTITION.complexity_comparisons.values())).summary.iloc[0]\n    complexity_rows.insert(0, {\"scope\": \"Global\", \"Raw K\": row.n_raw_clusters, \"Reconstructed K at Raw resolution\": row.n_recon_clusters, \"ARI\": row[\"ARI\"]})\ncomplexity_table = pd.DataFrame(complexity_rows); display(complexity_table.round(3))\ndisplay(Markdown(\"A larger reconstructed K supports the focused structural observation: coarse expression structure can remain, while smaller or boundary-ambiguous Raw clusters are split or reorganized.\"))"),
        ("markdown", "## 3. Matched-complexity change and Level1 localization\n\nFor each comparison, the reconstructed resolution is independently selected to match Raw K as closely as possible. Hungarian matching is then global and is reused for every Level1 proportion and spatial map."),
        ("code", "COMP = ({name: run[\"comparison\"] for name, run in RUNS.items()})\nCOMP_STATUS = {name: run[\"partition\"].matched_cluster_status for name, run in RUNS.items()}\nif 'GLOBAL_PARTITION' in globals():\n    COMP = {\"Global\": next(iter(GLOBAL_PARTITION.comparisons.values())), **COMP}\n    COMP_STATUS = {\"Global\": GLOBAL_PARTITION.matched_cluster_status, **COMP_STATUS}\nplot_comp = {(name if COMP_STATUS[name] != \"unmatched_cluster_complexity\" else f\"{name} [audit: unmatched K]\"): comparison for name, comparison in COMP.items()}\nfig = plot_contingency(plot_comp, \"Cluster-count-controlled contingency (unmatched K is audit only)\")\nsave_figure(fig, \"matched_k_contingency\"); plt.show()"),
        ("code", "change_rows=[]\nfor name, comparison in COMP.items():\n    row=comparison.summary.iloc[0]; status=COMP_STATUS[name]; change_rows.append({\"scope\": name, \"paired units\": row.n_units, \"Raw K\": row.n_raw_clusters, \"Recon K\": row.n_recon_clusters, \"status\": status, \"headline eligible\": status != \"unmatched_cluster_complexity\", \"ST-unit change\": row.st_unit_change_fraction, \"balanced change\": row.balanced_cluster_change, \"ARI\": row[\"ARI\"]})\nchange_table=pd.DataFrame(change_rows); display(change_table.round(4))\nif (change_table[\"headline eligible\"] == False).any(): display(Markdown(\"Rows marked `unmatched_cluster_complexity` are retained as audits; their change metrics are not interpreted as matched-K headlines.\"))"),
        ("code", "level1_tables=[]\nfor name, run in RUNS.items():\n    table=run[\"partition\"].change_by_level1.copy(); table.insert(0, \"scope\", name); level1_tables.append(table)\nif 'GLOBAL_PARTITION' in globals():\n    table=GLOBAL_PARTITION.change_by_level1.copy(); table.insert(0, \"scope\", \"Global\"); level1_tables.insert(0, table)\nlevel1_change=pd.concat(level1_tables, ignore_index=True)\nfig, ax=plt.subplots(figsize=(7, 5.8))\nif CONFIG[\"route_kind\"] == \"sc_svc\":\n    parent_order=[\"Fibroblast\", \"Mono_Macro\", \"T\"]\n    plot_table=(level1_change.loc[level1_change[\"scope\"].isin(parent_order)].sort_values(\"scope\", key=lambda x: x.map({name:i for i,name in enumerate(parent_order)})).drop_duplicates(\"scope\"))\n    lower=plot_table[\"change_fraction\"]-plot_table[\"wilson_ci_lower\"]; upper=plot_table[\"wilson_ci_upper\"]-plot_table[\"change_fraction\"]\n    ax.bar(plot_table[\"scope\"], plot_table[\"change_fraction\"], color=\"#4c78a8\", yerr=np.vstack([lower, upper]), capsize=4)\n    ax.set(title=\"Parent-internal matched-K reassignment\", ylabel=\"reassigned fraction\")\n    display_table=plot_table.loc[:, [\"scope\",\"level1\",\"total_units\",\"changed_units\",\"change_fraction\",\"wilson_ci_lower\",\"wilson_ci_upper\"]]\nelse:\n    plot_table=level1_change.query(\"scope == 'Global'\")\n    if GLOBAL_PARTITION.matched_cluster_status == \"unmatched_cluster_complexity\":\n        row=next(iter(GLOBAL_PARTITION.comparisons.values())).summary.iloc[0]\n        ax.axis(\"off\"); ax.text(.5,.55,f\"No headline Level1 change estimate\\nRaw K={{int(row.n_raw_clusters)}}, closest Recon K={{int(row.n_recon_clusters)}}\",ha=\"center\",va=\"center\",fontsize=13)\n        display_table=pd.DataFrame([{\"scope\":\"Global\",\"status\":GLOBAL_PARTITION.matched_cluster_status,\"Raw K\":int(row.n_raw_clusters),\"closest Recon K\":int(row.n_recon_clusters),\"Level1 audit artifact\":\"global/st_unit/change_by_level1.csv\"}])\n    else:\n        plot_table=plot_table.loc[plot_table[\"level1\"] != \"Overall\"].sort_values(\"change_fraction\")\n        lower=plot_table[\"change_fraction\"]-plot_table[\"wilson_ci_lower\"]; upper=plot_table[\"wilson_ci_upper\"]-plot_table[\"change_fraction\"]\n        ax.barh(plot_table[\"level1\"], plot_table[\"change_fraction\"], color=\"#4c78a8\", xerr=np.vstack([lower, upper]), capsize=3)\n        ax.set(title=\"Level1 matched-K change fraction\", xlabel=\"changed fraction\", ylabel=\"\")\n        focus_level1=[\"Overall\",\"Fibroblast\",\"Mono/Macro\",\"T\"]\n        focus=level1_change.loc[(level1_change[\"scope\"] == \"Global\") & level1_change[\"level1\"].isin(focus_level1)]\n        top=plot_table.nlargest(3, \"change_fraction\")\n        display_table=pd.concat([focus,top],ignore_index=True).drop_duplicates(\"level1\").loc[:, [\"scope\",\"level1\",\"total_units\",\"changed_units\",\"change_fraction\",\"wilson_ci_lower\",\"wilson_ci_upper\"]]\nsave_figure(fig, \"level1_change_fraction\"); plt.show()\ndisplay(display_table.round(4))"),
        ("code", "fig=plot_changed_units(RUNS, \"Changed-unit spatial localization\")\nsave_figure(fig, \"changed_units\"); plt.show()\nkey_mapping=[]\nfor name, run in RUNS.items():\n    table=run[\"comparison\"].mapping.copy(); table.insert(0, \"scope\", name); key_mapping.append(table.nlargest(2, \"overlap_n\"))\ndisplay(pd.concat(key_mapping, ignore_index=True).head(6))"),
        ("markdown", "## 4. Level1 anatomy Regions\n\nThis is an independent, full-tissue macro context. It is neither sampled with the VisiumHD partition switch nor used to tune cluster complexity, parent windows, or diversity thresholds."),
        ("code", "fig=plot_support_curve(ANATOMY, \"Route-level anatomy regions: occupancy decision\")\nsave_figure(fig, \"anatomy_window_decision\"); plt.show()\ndisplay(pd.DataFrame([{\"selected cells per side\": ANATOMY.scale_audit[\"main_window_cells_per_side\"], \"selected side (µm)\": ANATOMY.scale_audit[\"main_window_side_um\"], \"valid windows\": int(ANATOMY.anatomy_windows.shape[0]), \"retained units\": ANATOMY.support_selection.get(\"retained_parent_unit_fraction\")}]))"),
        ("code", "fig, axes = plt.subplots(2, 2, figsize=(10, 9))\nanatomy_map(axes[0,0]); anatomy_map(axes[0,1], \"Tumor\"); anatomy_map(axes[1,0], \"Normal\"); anatomy_map(axes[1,1], \"Interface\")\nfig.suptitle(\"Level1 anatomy Regions\"); fig.tight_layout()\nsave_figure(fig, \"anatomy_regions\"); plt.show()\ndisplay(ANATOMY.anatomy_context_summary.loc[ANATOMY.anatomy_context_summary[\"level1_region\"].isin([\"Tumor\",\"Normal\",\"Interface\"]), [\"level1_region\",\"full_level1_units\",\"tissue_windows\",\"area_um2\",\"area_fraction\"]])"),
    ]
    for number, parent in enumerate(("Fibroblast", "Mono_Macro", "T"), start=5):
        cells.extend(parent_cells(number, parent))
    cells.extend([
        ("markdown", "## 8. Cross-parent summary and evidence boundary\n\nThe summary compares observations within each route and parent. Region masks use route- and parent-specific data-driven thresholds, so their fractions are not direct cross-platform biological comparisons."),
        ("code", "summary=[]\nfor parent, run in RUNS.items():\n    row=run[\"comparison\"].summary.iloc[0]; metrics=run[\"impact\"].window_metrics.loc[lambda x: x.valid_window]\n    region=extent_table(run[\"impact\"], \"in_state_region\").query(\"level1_region == 'Overall'\")\n    summary.append({\"parent\": parent, \"matched-K status\": run[\"partition\"].matched_cluster_status, \"parent-internal reassignment\": row.st_unit_change_fraction, \"ARI\": row[\"ARI\"], \"median ΔNeff vs Raw Leiden\": metrics[\"delta_neff_vs_raw_leiden\"].median(), \"median ΔNeff vs Raw Level2\": metrics[\"delta_neff_vs_raw_level2\"].median(), \"High-diversity unit fraction\": region[\"unit_fraction\"].iloc[0]})\nsummary=pd.DataFrame(summary); display(summary.round(4))\ndisplay(Markdown(\"The notebook separates global Level1-stratified change from parent-internal matched-K reassignment and diversity under two Raw baselines. High-diversity Regions are reconstructed-state descriptions, not anatomy, mechanism, pathology, or biological validation.\"))"),
    ])
    return cells


def build(kind: str, filename: str) -> None:
    notebook = nbf.v4.new_notebook()
    notebook["metadata"] = {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}, "language_info": {"name": "python", "version": "3.10"}}
    converted = []
    for index, (cell_type, source) in enumerate(cells_for(kind), start=1):
        source = source.replace("{{int(row.n_raw_clusters)}}", "{int(row.n_raw_clusters)}")
        source = source.replace("{{int(row.n_recon_clusters)}}", "{int(row.n_recon_clusters)}")
        cell = nbf.v4.new_markdown_cell(source) if cell_type == "markdown" else nbf.v4.new_code_cell(source)
        cell["id"] = f"impact-{kind}-{index:02d}"
        converted.append(cell)
    notebook["cells"] = converted
    nbf.write(notebook, NOTEBOOK_DIR / filename)


if __name__ == "__main__":
    build("visium", "VisiumHD_sp_SVC_Reconstruction_Impact.ipynb")
    build("xenium", "Xenium_sc_SVC_Fibroblast_Reconstruction_Impact.ipynb")
