"""Generate the two scientific notebooks in their four-layer reading order.

Only cell rendering and presentation helpers live here. Scientific decisions and
calculations are emitted as explicit notebook cells, using existing algorithms.
"""
from pathlib import Path
from textwrap import dedent
from inspect import cleandoc
from string import Template
import nbformat as nbf

NOTEBOOK_DIR = Path(__file__).resolve().parent
try:
    from reproduce.case.reconstruction_impact.content_contract import NODE_SPECS
except ModuleNotFoundError:
    from content_contract import NODE_SPECS
NODE_SPECS_BY_ID = {spec["node_id"]: spec for spec in NODE_SPECS}
NOTEBOOK_BLOCK_TEMPLATE_PATH = NOTEBOOK_DIR.parents[2] / "docs" / "design" / "reconstruction-impact" / "templates" / "notebook-block.md"
NOTEBOOK_BLOCK_TEMPLATE = Template(NOTEBOOK_BLOCK_TEMPLATE_PATH.read_text(encoding="utf-8"))


def notebook_block(heading, purpose, method_links, observations=""):
    """Render a short English purpose/method block for source notebooks."""

    return NOTEBOOK_BLOCK_TEMPLATE.safe_substitute(
        heading=heading,
        purpose=purpose,
        method_links=method_links,
        observations=observations,
    ).strip()

HELPERS = r'''import os
import sys
import json
from types import SimpleNamespace
from pathlib import Path

os.environ["TQDM_DISABLE"] = "1"  # Disable progress at the producer; retain warnings.
os.environ.setdefault("KMP_WARNINGS", "0")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("KMP_USE_SHM", "0")
os.environ.setdefault("KMP_AFFINITY", "disabled")

PROJECT_ROOT = Path(os.environ.get("REVISE_REPOSITORY_ROOT", Path.cwd())).resolve()
if not (PROJECT_ROOT / "configs").is_dir():
    PROJECT_ROOT = Path.cwd().resolve()
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from IPython.display import Markdown, display

from revise.analysis.reconstruction_impact import load_reconstruction_impact_config

ROOT = PROJECT_ROOT

from reproduce.case.reconstruction_impact.notebook_helpers import (
    NotebookFigures, coordinates, compact_iqr, extent_table, raw_level2_mapping,
    plot_threshold_reliability, prepare_aucell_overall, prepare_aucell_parent,
    prepare_route_inputs, prepare_expression_views, save_partition_evidence,
    save_anatomy_evidence, save_spatial_feature_values, prepare_local_state,
    save_local_state_evidence, save_final_evidence,
)
'''

def parent_cells(number: int, parent: str):
    return [
        ("markdown", f"## {number}. {parent} internal diversity and Regions\n\nThis parent is analysed after the fixed macro anatomy. Raw Level2 is mapped from the original Raw expression after Level1 parent selection; reconstructed-carrier Level2 labels are not used as this baseline."),
        ("markdown", "### Cohort and window decision\n\nThe cell-equivalent side is fixed at 8 µm. The selected side comes only from parent occupancy; support is four units, and each of 200 paired draws samples four units."),
        ("code", f'''RUN = RUNS["{parent}"]; IMPACT = RUN["impact"]
fig = figures.plot_support_curve(IMPACT, "{parent}: occupancy-only window decision")
figures.save_figure(fig, "{parent.lower()}_window_decision"); plt.show()
occupancy = IMPACT.unit_assignments.groupby("window_id").size()
sampled = RUN.get("sampled", RUN["used"])
decision = pd.DataFrame([{{"available units": RUN["available"], "sampled units": sampled, "Raw-QC retained units": RUN["used"], "sampling mode": "full" if RUN["available"] == sampled else "deterministic 30k", "seed": CONFIG["partition_change"]["random_state"], "window side (cells per side)": IMPACT.scale_audit["main_window_cells_per_side"], "window side (µm)": IMPACT.scale_audit["main_window_side_um"], "occupancy median [Q1,Q3]": compact_iqr(occupancy), "minimum support": IMPACT.scale_audit["min_parent_units"], "rarefaction units": IMPACT.scale_audit["min_parent_units"], "paired draws": IMPACT.scale_audit["rarefaction_draws"]}}])
display(decision)'''),
        ("markdown", "### Internal-state baselines and matched-K assignments\n\nThe uniform Raw Level1 label is only the unexpanded baseline (`Kobs = Neff = evenness = 1`). The two informative Raw baselines are (i) expression-derived matched-K Leiden and (ii) reference-derived Raw Level2. Both are compared with the same reconstructed assignment on identical windows and identical rarefaction draws."),
        ("code", f'''METRICS = RUN["impact"].window_metrics.copy(); METRICS.attrs["side_um"] = IMPACT.scale_audit["main_window_side_um"]
fig = figures.plot_reconstructed_clusters(RUN, "{parent}: reconstructed internal state")
figures.save_figure(fig, "{parent.lower()}_reconstructed_clusters"); plt.show()
display(pd.DataFrame([{{"baseline": "Raw Level1 (uniform parent)", "Kobs": 1.0, "Neff": 1.0, "evenness": 1.0}}, {{"baseline": "Raw Level2 mapping", "Kobs": RUN["raw_level2"].audit["n_mapped_level2"], "Neff": "window-specific", "evenness": "window-specific"}}]))
display(pd.DataFrame([RUN["raw_level2"].audit]).loc[:, ["source", "method", "n_raw_units", "n_reference_cells", "n_reference_level2", "n_mapped_level2", "n_missing"]])'''),
        ("markdown", "### Matched-K reconstruction-associated internal structure\n\nThis is a **parent-internal reassignment** analysis, not a Level1 identity change. Raw Leiden is tuned to approximately the same cluster count as the reconstructed/final assignment. Hungarian matching controls label permutation for unit-change metrics; the diversity matrices below use the original cluster compositions within each window."),
        ("code", f'''fig = figures.plot_cluster_pair(RUN, "{parent}: matched-K cluster assignment")
figures.save_figure(fig, "{parent.lower()}_matched_clusters"); plt.show()
row = RUN["comparison"].summary.iloc[0]
display(pd.DataFrame([{{"paired units": row.n_units, "Raw K": row.n_raw_clusters, "Recon K": row.n_recon_clusters, "parent-internal reassignment": row.st_unit_change_fraction, "balanced change": row.balanced_cluster_change, "ARI": row["ARI"]}}]).round(4))'''),
        ("markdown", "### Kobs: local subtype richness\n\n`Kobs` is the number of cluster labels observed in a rarefied window. Here it asks how many distinct internal states coexist locally, without accounting for whether one state dominates."),
        ("code", f'''fig = figures.plot_metric_comparison(RUN, "k_obs", "{parent}: Kobs evidence matrix")
figures.save_figure(fig, "{parent.lower()}_kobs_matrix"); plt.show()
valid = METRICS.loc[METRICS.valid_window]
display(pd.DataFrame([{{"baseline": "Raw Leiden", "state median [Q1,Q3]": compact_iqr(valid["k_obs_raw"]), "Recon minus baseline": compact_iqr(valid["delta_k_obs_vs_raw_leiden"])}}, {{"baseline": "Raw Level2", "state median [Q1,Q3]": compact_iqr(valid["k_obs_level2"]), "Recon minus baseline": compact_iqr(valid["delta_k_obs_vs_raw_level2"])}}]))
display(valid[["k_obs_recon", "delta_k_obs_vs_raw_leiden", "delta_k_obs_vs_raw_level2"]].agg(["count", "median", "mean"]).T)'''),
        ("markdown", "### Neff: abundance-aware effective subtype count\n\n`Neff = exp(Shannon entropy)` converts subtype composition into an effective number of equally abundant clusters. It is the primary local-diversity state because it discounts rare labels and single-cluster dominance."),
        ("code", f'''fig = figures.plot_metric_comparison(RUN, "neff", "{parent}: Neff evidence matrix")
figures.save_figure(fig, "{parent.lower()}_neff_matrix"); plt.show()
valid = METRICS.loc[METRICS.valid_window]
display(pd.DataFrame([{{"baseline": "Raw Leiden", "state median [Q1,Q3]": compact_iqr(valid["neff_raw"]), "Recon minus baseline": compact_iqr(valid["delta_neff_vs_raw_leiden"])}}, {{"baseline": "Raw Level2", "state median [Q1,Q3]": compact_iqr(valid["neff_level2"]), "Recon minus baseline": compact_iqr(valid["delta_neff_vs_raw_level2"])}}]))
display(valid[["neff_recon", "delta_neff_vs_raw_leiden", "delta_neff_vs_raw_level2"]].agg(["count", "median", "mean"]).T)'''),
        ("markdown", "### Evenness: balance conditional on observed richness\n\n`evenness = Neff / Kobs` ranges from dominance toward balanced coexistence. It distinguishes windows with many labels but one dominant subtype from windows where those labels have comparable abundance."),
        ("code", f'''fig = figures.plot_metric_comparison(RUN, "evenness", "{parent}: evenness evidence matrix")
figures.save_figure(fig, "{parent.lower()}_evenness_matrix"); plt.show()
valid = METRICS.loc[METRICS.valid_window]
display(pd.DataFrame([{{"baseline": "Raw Leiden", "state median [Q1,Q3]": compact_iqr(valid["evenness_raw"]), "Recon minus baseline": compact_iqr(valid["delta_evenness_vs_raw_leiden"])}}, {{"baseline": "Raw Level2", "state median [Q1,Q3]": compact_iqr(valid["evenness_level2"]), "Recon minus baseline": compact_iqr(valid["delta_evenness_vs_raw_level2"])}}]))
display(valid[["evenness_recon", "delta_evenness_vs_raw_leiden", "delta_evenness_vs_raw_level2"]].agg(["count", "median", "mean"]).T)'''),
        ("markdown", "### High-diversity Region\n\nThe Region is a parent-specific mask derived only from the data-driven breakpoint of reconstructed Neff. It is displayed beside the continuous Neff field and is not overlaid on anatomy."),
        ("code", f'''fig = figures.plot_high_diversity(RUN, "{parent}: reconstructed diversity state and Region")
figures.save_figure(fig, "{parent.lower()}_high_diversity_region"); plt.show()
extent = extent_table(IMPACT, "in_state_region").loc[:, ["level1_region", "region_windows", "valid_windows", "region_area_um2", "area_fraction", "unit_fraction"]]
display(extent.loc[extent.level1_region.isin(["Overall", "Tumor", "Normal", "Interface"])])
display(pd.DataFrame([IMPACT.state_threshold]))'''),
    ]

def md(source):
    return ("markdown", cleandoc(source))


def code(source):
    return ("code", dedent(source).strip())


def foundation_cells(kind):
    config_name = "reconstruction_impact_visiumhd_p1crc.yaml" if kind == "visium" else "reconstruction_impact_xenium_p2crc_fibroblast.yaml"
    return [
        md(
            notebook_block(
                "## 1. Comparison foundation — what are we comparing?",
                "Raw and reconstruction views are paired on explicit observation, gene, label, and coordinate axes; each side retains its full gene space before analysis-specific preprocessing.",
                "[Comparison foundation](../../../docs/design/reconstruction-impact/input-views.md#comparison-foundation)",
                "Batch output: carrier, observation, and gene-availability audits.",
            )
            + """
### 1.1 Route semantics and carrier audit
Raw expression, fixed Level1 labels, spatial coordinates, reconstruction assignments, and the separate Raw Level2 reference remain distinct evidence inputs.
"""
        ),
        code(HELPERS),
        code(f"""
        CONFIG = load_reconstruction_impact_config(ROOT / "configs/analysis/{config_name}")
        OUTPUT_DIR = Path(os.environ.get("REVISE_ANALYSIS_OUTPUT_ROOT", ROOT / CONFIG["output"]["dir"]))
        SAMPLE_ID = CONFIG["sample"]["id"]
        LEVEL1 = CONFIG["context"]["level1_column"]
        SEED = int(CONFIG["partition_change"]["random_state"])
        VISIUMHD_SAMPLE_N_UNITS = 30_000
        USE_FULL_VISIUMHD_COHORT = False
        CELL_EQUIVALENT_UM = 8.0
        PARENT_SOURCE = {{"Fibroblast": "Fibroblast", "Mono_Macro": "Mono/Macro", "T": "T"}}
        SOURCE_FILES = {{"raw": ROOT / CONFIG["context"]["h5ad"], "raw_level2_reference": ROOT / CONFIG["raw_level2_mapping"]["reference_h5ad"]}}
        raw_context = ad.read_h5ad(SOURCE_FILES["raw"], backed="r")
        full_coordinates = coordinates(raw_context)
        full_level1 = raw_context.obs[LEVEL1].astype(str).copy()
        if not raw_context.obs_names.is_unique or not raw_context.var_names.is_unique:
            raise ValueError("Raw observation and gene axes must be unique")
        if not np.isfinite(full_coordinates.to_numpy()).all():
            raise ValueError("Raw spatial coordinates must be finite")
        TISSUE_COORDINATES_UM = full_coordinates * CONFIG["spatial_region"]["microns_per_coordinate"]
        TISSUE_XLIM = tuple(TISSUE_COORDINATES_UM["x"].agg(["min", "max"]))
        TISSUE_YLIM = tuple(TISSUE_COORDINATES_UM["y"].agg(["min", "max"]))
        ORIGIN_UM = (TISSUE_XLIM[0], TISSUE_YLIM[0])
        figures = NotebookFigures(OUTPUT_DIR, TISSUE_XLIM, TISSUE_YLIM, ORIGIN_UM, CELL_EQUIVALENT_UM)
        RUNS = {{}}
        input_rows = [{{"role": "Raw full context", "n_units": raw_context.n_obs, "n_genes": raw_context.n_vars, "spatial_axis": True}}]
        # global_ids and parent IDs are sampled inside prepare_route_inputs.
        """),
        md(
            notebook_block(
                "### 1.2 Cohort selection and coordinate audit",
                "Freeze Raw-defined cohorts and sampling decisions before partition QC; cohort denominators remain in the saved input tables.",
                "[Comparison foundation](../../../docs/design/reconstruction-impact/input-views.md#comparison-foundation)",
                "Batch output: included and excluded IDs with sampling reasons.",
            )
        ),
        code("""
        def deterministic_same_id_sample(ids, limit):
            ids = pd.Index(ids); return ids if USE_FULL_VISIUMHD_COHORT or len(ids) <= limit else pd.Index(np.random.default_rng(SEED).choice(ids.to_numpy(), size=limit, replace=False))
        prepared = prepare_route_inputs(
            CONFIG, ROOT, raw_context, full_coordinates, full_level1, figures,
            SOURCE_FILES, input_rows, seed=SEED,
            sample_limit=VISIUMHD_SAMPLE_N_UNITS,
            use_full_visiumhd_cohort=USE_FULL_VISIUMHD_COHORT,
            parent_source=PARENT_SOURCE, sample_fn=deterministic_same_id_sample,
        )
        SOURCE_FILES = prepared["source_files"]; RUNS = prepared["runs"]
        GLOBAL_INPUTS = prepared["global_inputs"]; input_rows = prepared["input_rows"]
        display(pd.DataFrame(input_rows))
        """),
        md(
            notebook_block(
                "### 1.3 Pairing and spatial coordinates",
                "Pair by exact spatial observation ID and record the existing coordinate equality range, unit, and tolerance before analysis-specific filtering.",
                "[Comparison foundation](../../../docs/design/reconstruction-impact/input-views.md#comparison-foundation)",
                "Batch output: pairing_audit.csv with full carrier ranges and boolean checks.",
            )
        ),
        md(
            notebook_block(
                "### 1.4 Expression views and gene availability",
                "Retain the complete gene space on both sides and record whether each gene is measured, projected, or unmeasured.",
                "[Comparison foundation](../../../docs/design/reconstruction-impact/input-views.md#comparison-foundation)",
                "Batch output: per-scope carrier mapping and full-gene availability table.",
            )
        ),
        code("""
        # Xenium expression_h5ad is projected inside this complete-gene helper.
        GENE_AVAILABILITY, input_rows, EXPRESSION_AUDITS = prepare_expression_views(
            RUNS, CONFIG, ROOT, figures, SOURCE_FILES, input_rows,
        )
        for table in EXPRESSION_AUDITS:
            display(table)
        display(GENE_AVAILABILITY.groupby(
            ["scope", "side", "status"], sort=False
        ).size().rename("genes").reset_index())
        """),
        md(
            notebook_block(
                "### 1.5 Analysis-specific preprocessing",
                "Keep Raw-derived partition QC and feature selection separate from full-gene Moran and original-view AUCell calculations.",
                "[Partition complexity](../../../docs/design/reconstruction-impact/impact-analysis.md#partition-complexity); [Moran method](../../../docs/design/reconstruction-impact/gene-and-function.md#moran); [EMT score](../../../docs/design/reconstruction-impact/gene-and-function.md#emt-score)",
                "Batch output: frozen parameter and preprocessing audit.",
            )
        ),
        code("""
        display(pd.DataFrame([{"analysis": "partition", **CONFIG["partition_change"]}, {"analysis": "Moran", "normalization": "normalize_total(1e4), log1p", "min_units": 51}, {"analysis": "AUCell", "resource": "Hallmark 2025.1 EMT", "AUC_threshold": 0.01, "seed": SEED}]))
        figures.save_json({"sample": SAMPLE_ID, "route": CONFIG["route_kind"], "seed": SEED, "partition": CONFIG["partition_change"], "spatial": CONFIG["spatial_region"]}, "parameters.json")
        """),
    ]


def partition_cells(kind):
    cells = [
        md(
            notebook_block(
                "## 2. Overall changes — what changed?",
                "Separate partition complexity from cluster-count-controlled assignment change before Moran, EMT, or spatial localization.",
                "[Partition complexity](../../../docs/design/reconstruction-impact/impact-analysis.md#partition-complexity); [matched-K](../../../docs/design/reconstruction-impact/impact-analysis.md#matched-k)",
                "Batch output: complexity sweeps, matched assignments, and change summaries.",
            )
            + """
### 2.1 Partition complexity diagnostic
The same-resolution diagnostic describes complexity; matched-K controls cluster count before interpreting assignment change.
"""
        ),
        code("""
        from revise.analysis.reconstruction_impact import run_partition_analysis
        PARTITIONS = {}
        if CONFIG["route_kind"] == "sp_svc":
            raw_global, recon_global = GLOBAL_INPUTS
            GLOBAL_PARTITION = run_partition_analysis(
                raw_global, recon_global, level1_col=LEVEL1, route_kind="sp_svc",
                resolution_mode="level1_ari", resolution_candidates=CONFIG["partition_change"]["level1_resolution_candidates"],
                within_level1_resolution=CONFIG["partition_change"]["within_level1_resolution"], random_state=SEED,
                n_top_genes=CONFIG["partition_change"]["n_top_genes"],
                raw_qc_min_genes=CONFIG["partition_change"]["raw_qc_min_genes"], raw_qc_min_cells=CONFIG["partition_change"]["raw_qc_min_cells"],
            )
            PARTITIONS["All"] = GLOBAL_PARTITION
        for parent, run in RUNS.items():
            raw = run["raw_expression"]
            spatial = run.get("spatial_carrier", run["recon_expression"])
            partition = run_partition_analysis(
                raw, spatial, level1_col=LEVEL1, route_kind=CONFIG["route_kind"],
                final_cluster_key="SVC_cluster" if CONFIG["route_kind"] == "sc_svc" else None,
                resolution_mode="fixed_within_level1", within_level1_resolution=CONFIG["partition_change"]["within_level1_resolution"],
                random_state=SEED, n_top_genes=CONFIG["partition_change"]["n_top_genes"],
                raw_qc_min_genes=CONFIG["partition_change"].get("raw_qc_min_genes", 50),
                raw_qc_min_cells=CONFIG["partition_change"].get("raw_qc_min_cells", 3),
            )
            comparison = next(iter(partition.comparisons.values()))
            ids = comparison.assignments.index
            run.update(partition=partition, comparison=comparison, used=len(ids))
            run["partition_coordinates"] = coordinates(raw).loc[ids] * CONFIG["spatial_region"]["microns_per_coordinate"]
            PARTITIONS[parent] = partition
        COMP = {scope: next(iter(partition.comparisons.values())) for scope, partition in PARTITIONS.items()}
        COMP_STATUS = {scope: partition.matched_cluster_status for scope, partition in PARTITIONS.items()}
        display(pd.DataFrame([{"scope": scope, "Raw-QC retained units": len(COMP[scope].assignments), "n_features": len(part.feature_names), "representation_audit": part.representation_audit, **part.audit} for scope, part in PARTITIONS.items()]))
        """),
        code("""
        complexity_table = pd.DataFrame([
            {"scope": scope, "Raw K": c.summary.iloc[0].n_raw_clusters, "Recon K at Raw resolution": c.summary.iloc[0].n_recon_clusters, "ARI": c.summary.iloc[0]["ARI"]}
            for scope, part in PARTITIONS.items() for c in part.complexity_comparisons.values()
        ])
        display(complexity_table.round(4))
        figures.save_table(complexity_table, "reconstruction_impact/partition/complexity_summary.csv")
        display(Markdown("Cluster counts describe representation complexity; they are not assignment-change fractions."))
        for scope, part in PARTITIONS.items():
            display(part.complexity_sweep)
            display(part.sweep)
        """),
        md(
            notebook_block(
                "### 2.2 Matched-complexity change",
                "Use one Hungarian mapping and Raw-row-normalized contingency to quantify assignment change at controlled cluster complexity.",
                "[Matched-K method](../../../docs/design/reconstruction-impact/impact-analysis.md#matched-k)",
                "Batch output: contingency, matched assignment, and headline-eligibility tables.",
            )
        ),
        code("""
        plot_comp = {(scope if COMP_STATUS[scope] != "unmatched_cluster_complexity" else scope + " [audit: unmatched K]"): c for scope, c in COMP.items()}
        fig = figures.plot_contingency(plot_comp, "Cluster-count-controlled contingency")
        figures.save_figure(fig, "matched_k_contingency")
        plt.show()
        change_table = pd.DataFrame([
            {"scope": scope, "paired units": c.summary.iloc[0].n_units, "Raw K": c.summary.iloc[0].n_raw_clusters,
             "Recon K": c.summary.iloc[0].n_recon_clusters, "ST-unit change": c.summary.iloc[0].st_unit_change_fraction,
             "balanced change": c.summary.iloc[0].balanced_cluster_change, "ARI": c.summary.iloc[0]["ARI"],
             "status": COMP_STATUS[scope], "headline_eligible": COMP_STATUS[scope] != "unmatched_cluster_complexity"}
            for scope, c in COMP.items()
        ])
        display(change_table.round(4))
        figures.save_table(change_table, "reconstruction_impact/partition/matched_k_summary.csv")
        """),
        md(
            notebook_block(
                "### 2.3 Change by cell type",
                "Report global Level1-stratified change for HD and parent-internal reassignment for Xenium using the corresponding Raw-defined cohorts.",
                "[Changed units and assignment basis](../../../docs/design/reconstruction-impact/impact-analysis.md#changed-units)",
                "Batch output: per-scope change fractions, Wilson intervals, and audit status.",
            )
        ),
        code("""
        level1_change = pd.concat([part.change_by_level1.assign(scope=scope) for scope, part in PARTITIONS.items()], ignore_index=True)
        fig, ax = plt.subplots(figsize=(7, 5.8))
        if "All" in PARTITIONS:
            plot_table = level1_change.loc[level1_change.scope.eq("All") & ~level1_change.level1.eq("Overall")].sort_values("change_fraction")
            eligible = COMP_STATUS["All"] != "unmatched_cluster_complexity"
        else:
            plot_table = level1_change.drop_duplicates("scope").copy()
            eligible = True
        if eligible:
            low = plot_table.change_fraction - plot_table.wilson_ci_lower
            high = plot_table.wilson_ci_upper - plot_table.change_fraction
            labels = plot_table.level1 if "All" in PARTITIONS else plot_table.scope
            ax.barh(labels, plot_table.change_fraction, xerr=np.vstack([low, high]), color="#4c78a8", capsize=3)
            ax.set(xlabel="changed fraction", title="Global Level1 change" if "All" in PARTITIONS else "Parent-internal reassignment")
            if "All" not in PARTITIONS:
                ax.set_yticks(range(len(labels)), [label + (" [audit: unmatched K]" if COMP_STATUS[label] == "unmatched_cluster_complexity" else "") for label in labels])
        else:
            ax.axis("off")
            ax.text(.5, .5, "No matched-K headline estimate: unmatched cluster complexity", ha="center", wrap=True)
        fig.tight_layout()
        figures.save_figure(fig, "level1_change_fraction")
        plt.show()
        display(level1_change.round(4))
        figures.save_table(level1_change, "reconstruction_impact/partition/level1_summary.csv")
        """),
        code("""
        COMPARISONS = save_partition_evidence(
            PARTITIONS, GLOBAL_INPUTS, RUNS, SAMPLE_ID, SEED, figures,
            complexity_kind=("raw_reference_vs_final"
                             if CONFIG["route_kind"] == "sc_svc"
                             else "same_resolution"),
        )
        """),
    ]

    if kind == "xenium":
        cells[0] = md(
            notebook_block(
                "## 2. Overall changes — what changed?",
                "Separate Xenium Raw reference-resolution complexity from matched comparison with the fixed reconstructed final-cluster count.",
                "[Partition complexity](../../../docs/design/reconstruction-impact/impact-analysis.md#partition-complexity); [matched-K](../../../docs/design/reconstruction-impact/impact-analysis.md#matched-k)",
                "Batch output: Xenium complexity sweep, matched assignments, and parent-internal change audit.",
            )
            + """
### 2.1 Partition complexity diagnostic
Xenium uses the Raw reference resolution and fixed final-cluster count; matched-K tunes Raw complexity toward that fixed count.
"""
        )
        cells[2] = ("code", cells[2][1].replace('"Recon K at Raw resolution"', '"Fixed final-cluster K"').replace('"Raw K"', '"Raw K at reference resolution"'))
        cells[-1] = ("code", cells[-1][1].replace('("same_resolution", part.complexity_comparisons)', '("raw_reference_vs_final", part.complexity_comparisons)'))
    return cells


AU_CELL_MARKDOWN = notebook_block(
    "### 2.5 Pre-specified EMT sentinel: AUCell",
    "This single EMT sentinel uses the local Hallmark 2025.1 resource with the fixed provider parameters. Raw and reconstruction retain their complete expression spaces, and coverage/status rows remain visible.",
    "[AUCell method](../../../docs/design/reconstruction-impact/gene-and-function.md#emt-coverage) and [score contract](../../../docs/design/reconstruction-impact/gene-and-function.md#emt-score)",
    "A missing resource gene or zero provider rank cutoff is recorded as unavailable; provider failures remain execution errors.",
)

AU_CELL_CODE = r'''# Score complete expression views with the fixed EMT resource; publish each audit.
from reproduce.case.reconstruction_impact.notebook_analysis import compute_emt
AUCELL_AUC_THRESHOLD = 0.01  # Quantile of detected genes, not a fixed rank fraction.
AUCELL_SEED = 42
AUCELL, AVAIL, aucell_summary, EMT_NAME, GMT_PATH = compute_emt(
    RUNS, ROOT, OUTPUT_DIR, CONFIG,
    AUCELL_AUC_THRESHOLD=AUCELL_AUC_THRESHOLD, AUCELL_SEED=AUCELL_SEED,
)
AUCELL_AVAIL = AVAIL
display(aucell_summary.round(4))
'''

AUCELL_CELLS = [
    ("markdown", AU_CELL_MARKDOWN),
    ("code", AU_CELL_CODE),
]

AUCELL_COVERAGE_MARKDOWN = r"""#### AUCell resource coverage and computability

Coverage is shown before score distributions.  It records the complete EMT
resource, measured resource genes on each side, the provider-derived dynamic
rank cutoff, and the explicit status/reason.  Coverage is descriptive and is
not used as an additional score eligibility threshold.
"""

AUCELL_COVERAGE_CODE = r'''# Display the already-published coverage audit; do not rerun AUCell.
_AU_SCOPE_LABELS = {
    "Fibroblast": "Fibroblast",
    "Mono_Macro": "Mono/Macro",
    "Mono/Macro": "Mono/Macro",
    "T": "T",
}
_AU_SCOPE_ORDER = [
    scope for scope in ("Fibroblast", "Mono_Macro", "Mono/Macro", "T")
    if scope in set(AVAIL.get("scope", pd.Series(dtype=object)).astype(str))
]
_AU_COVERAGE = AVAIL.copy()
if not _AU_COVERAGE.empty:
    _AU_COVERAGE["scope"] = _AU_COVERAGE["scope"].astype(str)
    _AU_COVERAGE["scope_label"] = _AU_COVERAGE["scope"].map(_AU_SCOPE_LABELS).fillna(_AU_COVERAGE["scope"])
    _AU_COVERAGE["scope_order"] = _AU_COVERAGE["scope"].map({value: index for index, value in enumerate(_AU_SCOPE_ORDER)})
    _AU_COVERAGE = _AU_COVERAGE.sort_values(["scope_order", "side"], kind="stable")
    display(_AU_COVERAGE.loc[:, [
        "scope_label", "side", "pathway", "resource_gene_count",
        "available_gene_count", "coverage", "detected_signature_gene_count",
        "n_observations", "n_genes", "detected_count_quantile",
        "provider_auc_threshold", "effective_rank_length",
        "rank_cutoff_zero_based", "status", "reason",
    ]].reset_index(drop=True).round(4))
else:
    display(Markdown("No AUCell coverage rows were published."))
'''

AUCELL_OVERALL_MARKDOWN = r"""#### Overall AUCell distributions and paired change

The first panel preserves the existing Raw-versus-reconstruction distribution
format by parent.  The second panel shows the paired unit-level
`reconstruction − Raw` score, using only units with both scores.  These are
overall EMT sentinel score changes; their spatial fields are displayed in the later
localization layer.
"""

AUCELL_OVERALL_CODE = r'''# Display published AUCell scores; the scorer is called in the prior cell.
_AU_LONG, _AU_DELTA, _AU_SUMMARY, _AU_COMPUTED = prepare_aucell_overall(
    AUCELL, aucell_summary, _AU_SCOPE_LABELS, _AU_SCOPE_ORDER
)
_AU_FIG = figures.plot_aucell_overall(_AU_LONG, _AU_DELTA, _AU_SCOPE_ORDER, _AU_SCOPE_LABELS, EMT_NAME)
figures.save_figure(_AU_FIG, "aucell_overall_distribution_and_delta")
plt.show()
_AU_SUMMARY_COLUMNS = [
    "scope", "pathway", "raw_status", "reconstruction_status", "comparison_status",
    "raw_n_valid", "raw_median", "raw_q1", "raw_q3", "reconstruction_n_valid",
    "reconstruction_median", "reconstruction_q1", "reconstruction_q3", "paired_n",
    "paired_delta_median", "paired_delta_q1", "paired_delta_q3",
]
if _AU_SUMMARY.empty:
    display(Markdown("No AUCell summary rows were published."))
else:
    display(_AU_SUMMARY.loc[:, [c for c in _AU_SUMMARY_COLUMNS if c in _AU_SUMMARY.columns]].reset_index(drop=True).round(4))
'''

def _parent_display_cells(scope_key: str, scope_label: str, ordinal: int):
    markdown = f"""#### {ordinal}. {scope_label}: EMT sentinel score comparison

This parent-level view reuses the published EMT scores and summary.  It keeps
Raw, reconstruction, and paired `reconstruction − Raw` values together for
one parent before the later spatial-field section.
"""
    code = f'''# Parent display only: all values come from AUCELL/aucell_summary above.
_AU_PARENT_SCOPE = {scope_key!r}
_AU_PARENT_LABEL = {scope_label!r}
_AU_PARENT, _AU_PARENT_LONG, _AU_PARENT_DELTA, _AU_PARENT_SUMMARY, _ = prepare_aucell_parent(
    AUCELL, aucell_summary, _AU_PARENT_SCOPE, _AU_PARENT_LABEL
)
_AU_PARENT_FIG = figures.plot_aucell_parent(_AU_PARENT_LONG, _AU_PARENT_DELTA, _AU_PARENT_LABEL, EMT_NAME)
figures.save_figure(_AU_PARENT_FIG, f"aucell_{{_AU_PARENT_SCOPE.lower().replace('/', '_')}}_distribution_and_delta")
plt.show()
if _AU_PARENT_SUMMARY.empty:
    display(Markdown("No AUCell summary row was published."))
else:
    display(_AU_PARENT_SUMMARY.round(4))
'''
    return [("markdown", markdown), ("code", code)]

AUCELL_DISPLAY_CELLS = [
    ("markdown", AUCELL_COVERAGE_MARKDOWN),
    ("code", AUCELL_COVERAGE_CODE),
    ("markdown", AUCELL_OVERALL_MARKDOWN),
    ("code", AUCELL_OVERALL_CODE),
]

for _ordinal, (_scope_key, _scope_label) in enumerate(
    (("Fibroblast", "Fibroblast"), ("Mono_Macro", "Mono/Macro"), ("T", "T")),
    start=1,
):
    AUCELL_DISPLAY_CELLS.extend(_parent_display_cells(_scope_key, _scope_label, _ordinal))


MORAN_MARKDOWN = notebook_block(
    "### 2.4 Gene spatial autocorrelation: Moran I",
    "This is a descriptive full-gene comparison. Each scope uses one Raw-coordinate graph for independently normalized Raw and reconstruction views; missing and non-computable genes remain explicit.",
    "[Moran method](../../../docs/design/reconstruction-impact/gene-and-function.md#moran)",
    "The shared-valid subset is a paired display; it does not replace the full-gene status table or define a Region.",
)

MORAN_SOURCE_CELL = r'''# Prepare and audit exact paired IDs before normalization.
from reproduce.case.reconstruction_impact.notebook_analysis import prepare_moran
MORAN_MIN_UNITS, MORAN_N_NEIGHS = 51, 6
MORAN_SCOPE_ORDER = []  # Keep the preparation stage visible in the notebook.
MORAN_SOURCES = {}
MORAN_COORDINATE_SOURCES = {}
(
    MORAN_SCOPE_ORDER, MORAN_SOURCES, MORAN_COORDINATE_SOURCES,
    _coordinate_unit, _coordinate_to_microns, MORAN_INPUT_AUDIT,
) = prepare_moran(
    RUNS, CONFIG,
    global_partition=globals().get("GLOBAL_PARTITION"),
    raw_global=globals().get("raw_global"), recon_global=globals().get("recon_global"),
    moran_min_units=MORAN_MIN_UNITS, moran_n_neighs=MORAN_N_NEIGHS,
)
display(MORAN_INPUT_AUDIT)'''

MORAN_COMPUTE_CELL = r'''# Normalize full-gene views and compute Moran on the same graph within each scope.
from reproduce.case.reconstruction_impact.notebook_analysis import compute_moran
MORAN, MORAN_GRAPH_AUDIT = compute_moran(
    MORAN_SOURCES, MORAN_SCOPE_ORDER, RUNS, CONFIG, OUTPUT_DIR,
    MORAN_COORDINATE_SOURCES, _coordinate_unit, _coordinate_to_microns,
    MORAN_MIN_UNITS=MORAN_MIN_UNITS, MORAN_N_NEIGHS=MORAN_N_NEIGHS,
)
display(MORAN_GRAPH_AUDIT)
display(MORAN.head())
'''

MORAN_AUDIT_CELL = r'''# Summarize the complete gene union; preserve every status and denominator.
from reproduce.case.reconstruction_impact.notebook_analysis import summarize_moran
MORAN_GENE_AVAILABILITY, MORAN_SUMMARY, MORAN_DISPLAY_SUMMARY = summarize_moran(
    MORAN, MORAN_SOURCES, MORAN_SCOPE_ORDER
)
figures.save_table(MORAN_GENE_AVAILABILITY, "gene_availability.csv", scope="moran")
figures.save_table(MORAN_SUMMARY, "moran_summary.csv", scope="moran")
figures.save_table(MORAN_DISPLAY_SUMMARY, "moran_distribution_summary.csv", scope="moran")
display(MORAN_SUMMARY)
display(MORAN_DISPLAY_SUMMARY)
display(MORAN_GENE_AVAILABILITY.groupby(
    ["scope", "raw_status", "reconstruction_status"], dropna=False
).size().rename("n_genes").reset_index())'''

MORAN_Q75_PLOT_CELL = r'''# Compare all-valid and shared-valid Q75 summaries using one color scale.
fig = figures.plot_moran_q75(MORAN_SUMMARY, MORAN_SCOPE_ORDER)
figures.save_figure(fig, "moran_q75_heatmap")
plt.show()
display(MORAN_DISPLAY_SUMMARY.loc[:, ["scope", "gene_set", "side", "n_valid", "median", "q1", "q75"]])'''

MORAN_DISTRIBUTION_PLOT_CELL = r'''# Display all computed genes and the paired shared-gene subset.
from reproduce.case.reconstruction_impact.notebook_analysis import build_moran_distribution
MORAN_DISTRIBUTION = build_moran_distribution(MORAN, MORAN_SCOPE_ORDER)
fig = figures.plot_moran_distribution(MORAN_DISTRIBUTION, MORAN_SCOPE_ORDER)
figures.save_figure(fig, "moran_all_gene_distribution")
plt.show()
display(MORAN_DISPLAY_SUMMARY.loc[:, ["scope", "gene_set", "side", "n_valid", "median", "q1", "q75"]])'''

MORAN_SCATTER_PLOT_CELL = r'''# Shared-gene comparison remains a descriptive paired view.
fig = figures.plot_moran_shared_scatter(MORAN, MORAN_SCOPE_ORDER)
figures.save_figure(fig, "moran_shared_gene_scatter")
plt.show()'''

MORAN_PAIRED_DELTA_PLOT_CELL = r'''# Paired deltas complement the side-specific full-gene distributions.
from reproduce.case.reconstruction_impact.notebook_analysis import build_moran_paired_delta
MORAN_PAIRED_DELTA = build_moran_paired_delta(MORAN, MORAN_SCOPE_ORDER)
fig = figures.plot_moran_paired_delta(MORAN_PAIRED_DELTA, MORAN_SCOPE_ORDER)
figures.save_figure(fig, "moran_shared_gene_delta_distribution")
plt.show()
display(MORAN_SUMMARY.loc[:, ["scope", "paired_delta_n", "paired_delta_median", "paired_delta_q1", "paired_delta_q3"]])'''

def moran_cells() -> list[tuple[str, str]]:
    """Return cells in the intended second-layer order."""

    return [
        ("markdown", MORAN_MARKDOWN),
        ("code", MORAN_SOURCE_CELL),
        ("code", MORAN_COMPUTE_CELL),
        ("code", MORAN_AUDIT_CELL),
        ("code", MORAN_Q75_PLOT_CELL),
        ("code", MORAN_DISTRIBUTION_PLOT_CELL),
        ("code", MORAN_SCATTER_PLOT_CELL),
        ("code", MORAN_PAIRED_DELTA_PLOT_CELL),
    ]

MORAN_CELL_SOURCES = {
    "markdown": MORAN_MARKDOWN,
    "source": MORAN_SOURCE_CELL,
    "compute": MORAN_COMPUTE_CELL,
    "audit": MORAN_AUDIT_CELL,
    "q75_plot": MORAN_Q75_PLOT_CELL,
    "distribution_plot": MORAN_DISTRIBUTION_PLOT_CELL,
    "shared_gene_scatter_plot": MORAN_SCATTER_PLOT_CELL,
    "shared_gene_delta_plot": MORAN_PAIRED_DELTA_PLOT_CELL,
}

def spatial_cells():
    cells = [
        md(
            notebook_block(
                "## 3. Spatial localization — where are the observed changes?",
                "Localize assignment changes, EMT scores, and internal diversity on fixed tissue coordinates after the overall comparisons are complete.",
                "[Anatomy method](../../../docs/design/reconstruction-impact/impact-analysis.md#anatomy); [changed units](../../../docs/design/reconstruction-impact/impact-analysis.md#changed-units)",
                "Batch output: anatomy context, changed-unit maps, spatial fields, and local-state evidence.",
            )
            + """
### 3.1 Level1 anatomy context
Anatomy is a fixed descriptive background and does not select complexity, cohort, scale, or threshold.
"""
        ),
        code("""
        from revise.analysis.reconstruction_impact import compute_anatomy_regions, _map_to_anatomy_regions, _summarize_parent_window_anatomy
        from revise.analysis.basic.spatial_region import (
            assign_square_windows,
            compute_rarefied_window_diversity,
            convert_coordinates_to_microns,
            select_region_threshold,
            select_window_scale,
            summarize_region_extent_by_anatomy,
        )
        ANALYSIS_ROOT = OUTPUT_DIR / "analysis"
        ANALYSIS_ROOT.mkdir(parents=True, exist_ok=True)
        ANATOMY = compute_anatomy_regions(
            full_coordinates=full_coordinates, full_level1_labels=full_level1,
            microns_per_coordinate=CONFIG["spatial_region"]["microns_per_coordinate"],
            candidate_window_sides_um=CONFIG["spatial_region"]["candidate_window_sides_um"],
            min_parent_units=4, cell_equivalent_um=CELL_EQUIVALENT_UM,
            **CONFIG["spatial_region"]["anatomy_region"],
        )
        anatomy_units = save_anatomy_evidence(
            ANATOMY, full_level1, CONFIG, ANALYSIS_ROOT, figures,
        )
        fig, axes = plt.subplots(2, 2, figsize=(10, 9))
        figures.anatomy_map(axes[0, 0], ANATOMY); figures.anatomy_map(axes[0, 1], ANATOMY, "Tumor")
        figures.anatomy_map(axes[1, 0], ANATOMY, "Normal"); figures.anatomy_map(axes[1, 1], ANATOMY, "Interface")
        fig.suptitle("Level1 anatomy Regions"); fig.tight_layout()
        figures.save_figure(fig, "anatomy_regions"); plt.show()
        display(ANATOMY.anatomy_context_summary.loc[lambda x: x.level1_region.isin(["Tumor", "Normal", "Interface"])])

        """),
        md(
            notebook_block(
                "#### Anatomy window support",
                "Show the selected anatomy window scale and retained support before reconstruction-dependent localization.",
                "[Window support method](../../../docs/design/reconstruction-impact/impact-analysis.md#window-support)",
                "Batch output: support sensitivity table and anatomy window decision figure.",
            )
        ),
        code("""
        fig = figures.plot_support_curve(ANATOMY, "Anatomy window support selection")
        figures.save_figure(fig, "anatomy_window_decision")
        plt.show()
        display(ANATOMY.support_sensitivity)
        """),
        md(
            notebook_block(
                "### 3.2 Changed-unit spatial localization",
                "Map matched-K assignment changes using the mapping already established in the overall layer; partitioning is not rerun here.",
                "[Changed units](../../../docs/design/reconstruction-impact/impact-analysis.md#changed-units)",
                "Batch output: per-parent matched assignment coordinates and changed-unit map.",
            )
        ),
        code("""
        for parent, run in RUNS.items():
            assignments = run["comparison"].assignments.set_index("unit_id")
            raw = run["raw_expression"][assignments.index].copy()
            xy = coordinates(raw) * CONFIG["spatial_region"]["microns_per_coordinate"]
            changed_dir = ANALYSIS_ROOT / "changed_units" / parent
            changed_dir.mkdir(parents=True, exist_ok=True)
            assignments.assign(x=xy["x"], y=xy["y"]).reset_index(names="unit_id").to_csv(changed_dir / "matched_k_assignments.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
        fig = figures.plot_changed_units(RUNS, "Changed-unit spatial localization")
        figures.save_figure(fig, "changed_units"); plt.show()
        display(change_table[["scope", "status", "headline_eligible"]])
        display(Markdown("Unmatched-K maps are diagnostic reassignment locations, not cluster-count-controlled change estimates."))
        """),
        md(
            notebook_block(
                "### 3.3 EMT sentinel spatial field",
                "Localize the pre-specified EMT AUCell score from the overall layer on the same tissue coordinates, retaining coverage and computability status.",
                "[EMT spatial fields](../../../docs/design/reconstruction-impact/gene-and-function.md#emt-spatial-fields)",
                "Batch output: Raw, reconstruction, delta, and status values for each parent field.",
            )
        ),
        code("""
        FEATURES = save_spatial_feature_values(
            RUNS, AUCELL, EMT_NAME, CONFIG, ANALYSIS_ROOT,
        )
        """),
        *[("code", f'''parent = {parent!r}
feature = {feature!r}
field = FEATURES.loc[FEATURES.scope.eq(parent) & FEATURES.feature.eq(feature)]
fig = figures.plot_feature_field(field, feature, f"{{parent}}: {{feature}} spatial field")
figures.save_figure(fig, f"{{parent.lower()}}_{{feature.lower()}}_spatial")
plt.show()
display(field[["raw_value", "reconstruction_value", "delta"]].agg(["count", "median", "mean"]).T)
display(field[["raw_status", "reconstruction_status", "expression_layer"]].drop_duplicates())
display(Markdown("Inspect the EMT score against the same tissue coordinates; unavailable scores remain unavailable."))''') for parent in ("Fibroblast", "Mono_Macro", "T") for feature in ("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",)],
        md(
            notebook_block(
                "### 3.4 Internal-state baselines and window decision",
                "Map Raw Level2 from the selected Raw parent and independent reference, then select the primary window scale from parent occupancy.",
                "[Local diversity and window definition](../../../docs/design/reconstruction-impact/impact-analysis.md#local-diversity); [window support](../../../docs/design/reconstruction-impact/impact-analysis.md#window-support)",
                "Batch output: Raw Level2 mapping, paired window assignments, and scale decision.",
            )
        ),
        code("""
        for parent, run in RUNS.items():
            prepare_local_state(
                parent, run, CONFIG, ROOT, ANATOMY, ANALYSIS_ROOT, figures,
                CELL_EQUIVALENT_UM,
            )
        """),
        *[entry for parent in ("Fibroblast", "Mono_Macro", "T") for entry in
          [parent_cells(0, parent)[2], parent_cells(0, parent)[4], parent_cells(0, parent)[6]]],
        md(
            notebook_block(
                "### 3.5 Local richness: Kobs",
                "Count observed internal states in the same paired rarefaction windows for Raw Leiden, reconstruction, and Raw Level2.",
                "[Local diversity method](../../../docs/design/reconstruction-impact/impact-analysis.md#local-diversity)",
                "Batch output: 2×3 Kobs state and delta evidence matrix.",
            )
        ),
        code("""
        for parent, run in RUNS.items():
            spatial = run["spatial"]
            metrics = compute_rarefied_window_diversity(spatial["windows"], spatial["assignments"]["raw_cluster"], spatial["assignments"]["recon_cluster"], raw_level2_labels=spatial["level2"].labels, min_parent_units=4, n_draws=200, random_state=SEED)
            metrics["scope"] = parent
            metrics["scale_um"] = spatial["side_um"]
            run["spatial"]["metrics"] = metrics
            run["impact"].window_metrics = metrics
            metrics.to_csv(ANALYSIS_ROOT / "local_state" / parent / "window_metrics.csv", index=False)
        """),
        *[("code", f'RUN = RUNS[{parent!r}]; IMPACT = RUN["impact"]; METRICS = IMPACT.window_metrics\n' + parent_cells(0, parent)[8][1]) for parent in ("Fibroblast", "Mono_Macro", "T")],
        md(
            notebook_block(
                "### 3.6 Effective diversity: Neff",
                "Compare abundance-aware effective subtype counts under identical paired rarefaction support.",
                "[Local diversity method](../../../docs/design/reconstruction-impact/impact-analysis.md#local-diversity)",
                "Batch output: 2×3 Neff state and delta evidence matrix.",
            )
        ),
        *[("code", f'RUN = RUNS[{parent!r}]; IMPACT = RUN["impact"]; METRICS = IMPACT.window_metrics\n' + parent_cells(0, parent)[10][1]) for parent in ("Fibroblast", "Mono_Macro", "T")],
        md(
            notebook_block(
                "### 3.7 Evenness",
                "Report balance conditional on local richness and interpret it with the paired Kobs and Neff baselines.",
                "[Local diversity method](../../../docs/design/reconstruction-impact/impact-analysis.md#local-diversity)",
                "Batch output: 2×3 evenness state and delta evidence matrix.",
            )
        ),
        *[("code", f'RUN = RUNS[{parent!r}]; IMPACT = RUN["impact"]; METRICS = IMPACT.window_metrics\n' + parent_cells(0, parent)[12][1]) for parent in ("Fibroblast", "Mono_Macro", "T")],
        md(
            notebook_block(
                "### 3.8 Anatomy-stratified summaries",
                "Stratify local diversity and assignment change by the fixed anatomy context after window metrics are computed.",
                "[Anatomy method](../../../docs/design/reconstruction-impact/impact-analysis.md#anatomy)",
                "Batch output: anatomy-stratified diversity and change tables for each parent.",
            )
        ),
        code("""
        for parent, run in RUNS.items():
            change_summary, anatomy_summary = save_local_state_evidence(
                parent, run, ANATOMY, ANALYSIS_ROOT,
            )
            display(change_summary)
            display(anatomy_summary)
        """),
    ]
    return cells

def region_cells():
    return [
        md(
            notebook_block(
                "## 4. Regions: local-state coverage and reliability",
                "Define and audit the reconstructed-Neff State Region only after local diversity, anatomy, and window support are complete.",
                "[Region reliability](../../../docs/design/reconstruction-impact/impact-analysis.md#region-reliability); [Region extent](../../../docs/design/reconstruction-impact/impact-analysis.md#region-extent)",
                "Batch output: threshold reliability, Region masks, anatomy extent, gain audit, and scale sensitivity.",
            )
            + """
### 4.1 Threshold reliability — State Region
State Region is a reconstructed-Neff state mask, not a maximum Raw-to-reconstruction difference; invalid windows remain NA.
"""
        ),
        code("""
        for parent, run in RUNS.items():
            spatial = run["spatial"]; metrics = spatial["metrics"].copy()
            valid = metrics.loc[metrics.valid_window]
            state_threshold, bootstrap = select_region_threshold(valid.neff_recon.to_numpy(), n_bootstrap=CONFIG["spatial_region"]["threshold_bootstraps"], random_state=SEED)
            gain_values = valid.loc[valid.delta_neff_vs_raw_leiden.gt(0), "delta_neff_vs_raw_leiden"].to_numpy()
            gain_threshold, gain_bootstrap = select_region_threshold(gain_values, n_bootstrap=500, random_state=SEED)
            metrics["in_state_region"] = pd.Series(pd.NA, index=metrics.index, dtype="boolean")
            metrics["in_gain_region"] = pd.Series(pd.NA, index=metrics.index, dtype="boolean")
            if gain_threshold["status"] == "ok":
                metrics.loc[metrics.valid_window, "in_gain_region"] = metrics.loc[metrics.valid_window, "delta_neff_vs_raw_leiden"].ge(float(gain_threshold["threshold"]))
            if state_threshold["status"] == "ok":
                metrics.loc[metrics.valid_window, "in_state_region"] = metrics.loc[metrics.valid_window, "neff_recon"].ge(float(state_threshold["threshold"]))
            impact = SimpleNamespace(window_metrics=metrics, scale_audit={"main_window_side_um": spatial["side_um"]}, state_threshold=state_threshold)
            run["impact"] = impact
            spatial.update(metrics=metrics, state_threshold=state_threshold, state_bootstrap=bootstrap, gain_threshold=gain_threshold, gain_bootstrap=gain_bootstrap)
        """),
        code("""
        # Diagnose all three parents before showing any Region masks.
        THRESHOLD_RELIABILITY, fig = plot_threshold_reliability(
            RUNS, CONFIG["spatial_region"]["threshold_bootstraps"])
        figures.save_table(THRESHOLD_RELIABILITY, "threshold_reliability.csv")
        figures.save_figure(fig, "threshold_reliability")
        plt.show()
        display(THRESHOLD_RELIABILITY)
        """),
        md(
            notebook_block(
                "### 4.2 Continuous reconstructed-Neff state and mask",
                "Display the continuous reconstructed-Neff field beside its thresholded State Region mask, preserving unavailable windows.",
                "[Region extent](../../../docs/design/reconstruction-impact/impact-analysis.md#region-extent)",
                "Batch output: per-parent continuous state and State Region figures.",
            )
        ),
        *[("code", f'''run = RUNS[{parent!r}]
fig = figures.plot_high_diversity(run, "{parent}: reconstructed Neff and State Region")
figures.save_figure(fig, "{parent.lower()}_state_region")
plt.show()
display(pd.DataFrame([run["impact"].state_threshold]))
display(Markdown("The continuous state and Region mask answer different questions; an unavailable threshold yields no coverage estimate."))''') for parent in ("Fibroblast", "Mono_Macro", "T")],
        md(
            notebook_block(
                "### 4.3 Region extent by anatomy",
                "Summarize State Region and gain-region extent within the fixed anatomy context while retaining threshold status.",
                "[Region extent](../../../docs/design/reconstruction-impact/impact-analysis.md#region-extent)",
                "Batch output: per-parent overall, Tumor, Normal, and Interface extent rows.",
            )
        ),
        code("""
        REGION_ROWS = []
        for parent, run in RUNS.items():
            spatial = run["spatial"]; metrics = spatial["metrics"].copy(); metrics["in_region"] = metrics["in_state_region"]
            extent = summarize_region_extent_by_anatomy(metrics, window_side_length=float(spatial["side_um"])).assign(scope=parent, scale_um=spatial["side_um"], threshold_status=spatial["state_threshold"]["status"])
            metrics = metrics.drop(columns="in_region")
            spatial.update(metrics=metrics, extent=extent)
            run["impact"].region_extent_by_anatomy = extent
            local_dir = ANALYSIS_ROOT / "local_state" / parent
            metrics.to_csv(local_dir / "window_metrics_with_state_region.csv", index=False)
            extent.to_csv(local_dir / "state_region_extent_by_anatomy.csv", index=False)
            gain_extent = summarize_region_extent_by_anatomy(metrics.assign(in_region=metrics.in_gain_region), window_side_length=float(spatial["side_um"]))
            gain_extent.to_csv(local_dir / "gain_region_extent_audit.csv", index=False)
            display(extent.loc[extent.level1_region.isin(["Overall", "Tumor", "Normal", "Interface"])])
            REGION_ROWS.append(extent.assign(parent=parent))
        REGION_EXTENT = pd.concat(REGION_ROWS, ignore_index=True)
        """),
        md(
            notebook_block(
                "#### Threshold audit",
                "Publish bootstrap support and the audit-only gain threshold before interpreting Region coverage.",
                "[Region reliability](../../../docs/design/reconstruction-impact/impact-analysis.md#region-reliability)",
                "Batch output: state and gain bootstrap tables plus threshold JSON audits.",
            )
        ),
        code("""
        for parent, run in RUNS.items():
            spatial = run["spatial"]; local_dir = ANALYSIS_ROOT / "local_state" / parent
            spatial["state_bootstrap"].to_csv(local_dir / "state_threshold_bootstrap.csv", index=False)
            spatial["gain_bootstrap"].to_csv(local_dir / "gain_threshold_bootstrap_audit.csv", index=False)
            pd.DataFrame([spatial["state_threshold"]]).to_json(local_dir / "state_threshold.json", orient="records", indent=2)
            pd.DataFrame([spatial["gain_threshold"] | {"role": "audit_only"}]).to_json(local_dir / "gain_threshold_audit.json", orient="records", indent=2)
            display(pd.DataFrame([spatial["state_threshold"]])); display(spatial["state_bootstrap"].head())
        """),
        md(
            notebook_block(
                "### 4.4 Scale sensitivity",
                "Repeat the paired rarefaction definition across candidate window scales as a sensitivity audit of local diversity.",
                "[Scale sensitivity](../../../docs/design/reconstruction-impact/impact-analysis.md#scale-sensitivity)",
                "Batch output: per-parent Kobs, Neff, evenness, and delta summaries across scales.",
            )
        ),
        code("""
        SCALE_ROWS = []
        for parent, run in RUNS.items():
            spatial = run["spatial"]; rows = []
            for side_um in CONFIG["spatial_region"]["candidate_window_sides_um"]:
                windows = assign_square_windows(spatial["paired_um"], window_side_length=float(side_um), origin=spatial["origin"])
                values = compute_rarefied_window_diversity(windows, spatial["assignments"].raw_cluster, spatial["assignments"].recon_cluster, raw_level2_labels=spatial["level2"].labels, min_parent_units=4, n_draws=200, random_state=SEED)
                valid = values.loc[values.valid_window]
                row = {"parent": parent, "window_side_length": float(side_um), "n_valid_windows": len(valid), "min_parent_units": 4, "rarefaction_draws": 200, "random_state": SEED}
                for metric in ("k_obs", "neff", "evenness"):
                    for suffix in ("raw", "recon", "level2"):
                        row[f"median_{metric}_{suffix}"] = valid[f"{metric}_{suffix}"].median()
                    for baseline in ("raw_leiden", "raw_level2"):
                        row[f"median_delta_{metric}_vs_{baseline}"] = valid[f"delta_{metric}_vs_{baseline}"].median()
                rows.append(row)
            sensitivity = pd.DataFrame(rows); sensitivity.to_csv(ANALYSIS_ROOT / "local_state" / parent / "scale_sensitivity.csv", index=False)
            SCALE_ROWS.append(sensitivity); display(sensitivity)
        SCALE_SENSITIVITY = pd.concat(SCALE_ROWS, ignore_index=True)
        SCALE_SENSITIVITY.to_csv(ANALYSIS_ROOT / "local_state" / "scale_sensitivity_all_parents.csv", index=False)
        """),
        md(
            notebook_block(
                "### 4.5 Region conclusion",
                "Leave the Region interpretation to the post-run evidence record after reliability, continuous state, extent, and scale checks.",
                "[Region reliability](../../../docs/design/reconstruction-impact/impact-analysis.md#region-reliability)",
                "Batch output: reviewed Region conclusion placeholder.",
            )
        ),
    ]

def overall_summary_cells():
    return [
        md(
            notebook_block(
                "### 2.6 Summary of observed changes",
                "Keep partition, full-gene Moran, and the pre-specified EMT sentinel as separate evidence streams before spatial follow-up.",
                "[Changed units](../../../docs/design/reconstruction-impact/impact-analysis.md#changed-units); [Moran method](../../../docs/design/reconstruction-impact/gene-and-function.md#moran); [EMT score](../../../docs/design/reconstruction-impact/gene-and-function.md#emt-score)",
                "Batch output: compact overall tables with missing Raw measurements still visible.",
            )
        ),
        code('''display(change_table.round(4))
display(MORAN_SUMMARY.round(4))
display(aucell_summary.round(4))
display(Markdown("Next inspect assignment changes and the EMT sentinel field, then local diversity under both Raw baselines."))''')]


def final_cells():
    return [
        md(
            notebook_block(
                "## 5. Cross-parent summary and evidence boundary",
                "Summarize completed cross-parent evidence without combining partition, gene, pathway, local-diversity, or Region results into one score.",
                "[Comparison foundation](../../../docs/design/reconstruction-impact/input-views.md#comparison-foundation); [Region reliability](../../../docs/design/reconstruction-impact/impact-analysis.md#region-reliability)",
                "Batch output: cross-parent local-diversity, Region-extent, summary, source, and artifact audits.",
            )
            + """
The interpretation remains descriptive and is constrained by Raw coverage, projected Xenium expression, sampling, graph construction, and window support.
"""
        ),
        code('''SOURCE_FILES["content_contract.py"] = ROOT / "reproduce/case/reconstruction_impact/content_contract.py"
FINAL_EVIDENCE = save_final_evidence(
    RUNS, CONFIG, OUTPUT_DIR, COMPARISONS, change_table, MORAN_SUMMARY,
    aucell_summary, SAMPLE_ID, EMT_NAME, GMT_PATH, SOURCE_FILES, figures,
)
LOCAL_DIVERSITY_CROSS_PARENT = FINAL_EVIDENCE["local_diversity"]
REGION_EXTENT_CROSS_PARENT = FINAL_EVIDENCE["region_extent"]
LOCAL_WIDE = FINAL_EVIDENCE["local_wide"]
CROSS_PARENT = FINAL_EVIDENCE["cross_parent"]
ARTIFACTS = FINAL_EVIDENCE["artifacts"]
display(CROSS_PARENT.round(4))
display(LOCAL_DIVERSITY_CROSS_PARENT.round(4))
display(REGION_EXTENT_CROSS_PARENT.loc[
    REGION_EXTENT_CROSS_PARENT["level1_region"].isin(["Overall", "Tumor", "Normal", "Interface"])
].round(4))'''),
        code('''display(ARTIFACTS)
raw_context.file.close()
if "recon_context" in globals():
    recon_context.file.close()''')]



def impact_node(node_id):
    """Return a post-evidence placeholder consumed by the report runner."""

    spec = NODE_SPECS_BY_ID[node_id]
    method_file, method_anchor = spec["method"]
    method_target = method_file + (f"#{method_anchor}" if method_anchor else "")
    return (
        "markdown",
        f"<!-- impact-node:{node_id} -->\n"
        f"**{spec['title_en']} — batch interpretation**\n\n"
        f"Method: [{method_target}](../../../docs/design/reconstruction-impact/{method_target})",
    )


def _insert_before_heading(cells, heading, node_id):
    for index, (cell_type, source) in enumerate(cells):
        if cell_type == "markdown" and heading in source:
            cells.insert(index, impact_node(node_id))
            return
    raise ValueError(f"Cannot place impact node {node_id}: missing {heading!r}")


def _add_impact_nodes(cells):
    placements = (
        ("### 1.2 Cohort selection", "1.1"),
        ("### 1.3 Pairing and spatial coordinates", "1.2"),
        ("### 1.4 Expression views and gene availability", "1.3"),
        ("### 1.5 Analysis-specific preprocessing", "1.4"),
        ("## 2. Overall changes", "1.5"),
        ("## 2. Overall changes", "1.6"),
        ("### 2.2 Matched-complexity change", "2.1"),
        ("### 2.3 Change by cell type", "2.2"),
        ("### 2.4 Gene spatial autocorrelation", "2.3"),
        ("### 2.5 Pre-specified EMT sentinel", "2.4"),
        ("### 2.6 Summary of observed changes", "2.5"),
        ("## 3. Spatial localization", "2.6"),
        ("### 3.2 Changed-unit spatial localization", "3.1"),
        ("### 3.3 EMT sentinel spatial field", "3.2"),
        ("### 3.4 Internal-state baselines and window decision", "3.3"),
        ("### 3.5 Local richness", "3.4"),
        ("### 3.6 Effective diversity", "3.5"),
        ("### 3.7 Evenness", "3.6"),
        ("### 3.8 Anatomy-stratified summaries", "3.7"),
        ("## 4. Regions", "3.8"),
        ("## 4. Regions", "3.9"),
        ("### 4.2 Continuous reconstructed-Neff state and mask", "4.1"),
        ("### 4.3 Region extent by anatomy", "4.2"),
        ("#### Threshold audit", "4.3"),
        ("### 4.5 Region conclusion", "4.4"),
        ("## 5. Cross-parent summary", "4.5"),
    )
    for heading, node_id in placements:
        _insert_before_heading(cells, heading, node_id)
    cells.append(impact_node("summary"))
    return cells


def cells_for(kind):
    route = "VisiumHD" if kind == "visium" else "Xenium"
    cells = [md(f'''# {route} reconstruction impact: Fibroblast, Mono/Macro and T

Read the first half for **what changed**, and the second half for **where those
changes occur**. Figures are followed by summaries; analysis tables and figures
are saved as each block completes. The source files and exact cohorts are audited
below. Each route keeps its own spatial and expression semantics.

1. [Comparison foundation](#1.-Comparison-foundation-—-what-are-we-comparing?)
2. [Overall changes](#2.-Overall-changes-—-what-changed?)
3. Spatial localization: anatomy, changed units, expression and diversity
4. Regions: state, extent, reliability and scale sensitivity
5. Cross-parent summary and evidence boundary

Design reference: `docs/design/reconstruction-impact/README.md`.
''')] + foundation_cells(kind) + partition_cells(kind) + moran_cells() + AUCELL_CELLS + AUCELL_DISPLAY_CELLS + overall_summary_cells() + spatial_cells() + region_cells() + final_cells()
    return _add_impact_nodes(cells)


def build(kind, filename):
    cells = []
    for index, (cell_type, source) in enumerate(cells_for(kind)):
        if cell_type == "markdown" and source.startswith("#") and not source.startswith("# "):
            # Keep navigation and method access; reviewed batch findings follow the evidence.
            lines = source.splitlines()
            source = "\n\n".join(
                line for line in lines
                if line.startswith("#") or line.startswith("**Method:")
            )
        cell = nbf.v4.new_markdown_cell(source) if cell_type == "markdown" else nbf.v4.new_code_cell(source)
        cell["id"] = f"{kind}-{index:03d}"
        if source.startswith("<!-- impact-node:"):
            cell["metadata"]["impact_node"] = source.split(":", 1)[1].split("-->", 1)[0].strip()
        cells.append(cell)
    notebook = nbf.v4.new_notebook(cells=cells, metadata={"kernelspec": {"display_name": "Python 3.10", "language": "python", "name": "python3.10"}, "language_info": {"name": "python", "version": "3.10"}})
    nbf.write(notebook, NOTEBOOK_DIR / filename)


if __name__ == "__main__":
    build("visium", "VisiumHD_sp_SVC_Reconstruction_Impact.ipynb")
    build("xenium", "Xenium_sc_SVC_Fibroblast_Reconstruction_Impact.ipynb")
