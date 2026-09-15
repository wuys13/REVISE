"""Small real-AnnData contracts for the notebook Moran and AUCell cells.

These tests execute the source emitted by ``build_notebooks.py``.  They do not
execute either full route and do not copy the temporary AUCell fragment into
the repository.  The input deliberately contains a partial gene union, a
constant gene, and a parent with no resource-gene hit so that missing biology
cannot be hidden by a summary-only check.
"""

from __future__ import annotations

import importlib
import json
from pathlib import Path
import re
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
BUILDER_MODULE = "reproduce.case.reconstruction_impact.build_notebooks"


def _builder_source_cells() -> list[str]:
    """Collect actual cell strings exposed by the notebook builder."""

    builder = importlib.import_module(BUILDER_MODULE)
    return [source for kind, source in builder.cells_for("visium") if kind == "code"]


def _scientific_cell(marker: str, label: str) -> str:
    stage_marker = r"\bcompute_moran\s*\(" if label == "Moran" else r"\bcompute_emt\s*\("
    candidates = [source for source in _builder_source_cells() if re.search(stage_marker, source)]
    assert candidates, f"builder does not emit a {label} scientific cell"
    # Prefer the executable calculation cell over a short markdown/helper
    # string containing the same term.
    return max(candidates, key=len)


def _execute_stage(source, namespace):
    from reproduce.case.reconstruction_impact import notebook_analysis
    with pytest.MonkeyPatch.context() as patch:
        for name in ('compare_moran', 'score_gene_set_aucell', 'get_aucell_provider_metadata'):
            if name in namespace:
                patch.setattr(notebook_analysis, name, namespace[name])
        exec(compile(source, '<notebook-stage>', 'exec'), namespace)


def _expression_pair(
    *,
    raw_gene_ids: list[str],
    reconstruction_gene_ids: list[str],
    prefix: str,
    n_units: int = 60,
) -> tuple[ad.AnnData, ad.AnnData]:
    """Build deterministic sparse Raw/reconstruction carriers with exact IDs."""

    unit_ids = pd.Index([f"{prefix}-{index:03d}" for index in range(n_units)])
    obs = pd.DataFrame(index=unit_ids)
    raw_columns = []
    for column, gene in enumerate(raw_gene_ids):
        if gene.startswith("constant"):
            raw_columns.append(np.zeros(n_units, dtype=float))
        else:
            raw_columns.append((np.arange(n_units, dtype=float) + column + 1) % 17 + 1)
    recon_columns = []
    for column, gene in enumerate(reconstruction_gene_ids):
        if gene.startswith("constant"):
            recon_columns.append(np.zeros(n_units, dtype=float))
        else:
            recon_columns.append((np.arange(n_units, dtype=float) + 2 * column + 2) % 19 + 1)
    raw = ad.AnnData(
        sparse.csr_matrix(np.column_stack(raw_columns)),
        obs=obs.copy(),
        var=pd.DataFrame(index=raw_gene_ids),
    )
    reconstruction = ad.AnnData(
        sparse.csr_matrix(np.column_stack(recon_columns)),
        obs=obs.copy(),
        var=pd.DataFrame(index=reconstruction_gene_ids),
    )
    coordinates = np.column_stack(
        [np.arange(n_units, dtype=float), np.mod(np.arange(n_units, dtype=float) ** 2, 31)]
    )
    raw.obsm["spatial"] = coordinates.copy()
    reconstruction.obsm["spatial"] = coordinates.copy()
    return raw, reconstruction


def _save_table_factory(output_dir: Path):
    def save_table(frame: pd.DataFrame, name: str, *, scope: str | None = None):
        destination = output_dir / "analysis"
        if scope is not None:
            destination /= str(scope)
        destination /= name
        destination.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(destination, index=False)

    return save_table


def test_builder_moran_cell_uses_one_shared_graph_and_preserves_full_gene_statuses(tmp_path):
    pytest.importorskip("scanpy")
    pytest.importorskip("squidpy")
    import scanpy as sc
    import squidpy as sq

    source = _scientific_cell(r"\bcompare_moran\s*\(", "Moran")
    raw, reconstruction = _expression_pair(
        raw_gene_ids=["shared", "constant_raw", "raw_only"],
        reconstruction_gene_ids=["shared", "constant_recon", "recon_only"],
        prefix="moran",
    )
    parent_pairs = {
        "Fibroblast": (raw, reconstruction),
        "Mono_Macro": _expression_pair(
            raw_gene_ids=["shared", "constant_raw", "raw_only"],
            reconstruction_gene_ids=["shared", "constant_recon", "recon_only"],
            prefix="mono",
        ),
        "T": _expression_pair(
            raw_gene_ids=["shared", "constant_raw", "raw_only"],
            reconstruction_gene_ids=["shared", "constant_recon", "recon_only"],
            prefix="t",
        ),
    }
    output_dir = Path(tmp_path)
    calls = []
    from revise.analysis.paired_moran import compare_moran as real_compare_moran

    def compare_moran_spy(raw_view, reconstruction_view, weights, **kwargs):
        calls.append(sparse.csr_matrix(weights, copy=True))
        return real_compare_moran(raw_view, reconstruction_view, weights, **kwargs)

    namespace = {
        "RUNS": {
            parent: {"raw_expression": pair[0], "recon_expression": pair[1]}
            for parent, pair in parent_pairs.items()
        },
        "GLOBAL_PARTITION": None,
        "raw_global": raw,
        "recon_global": reconstruction,
        "np": np,
        "pd": pd,
        "sc": sc,
        "sq": sq,
        "sparse": sparse,
        "compare_moran": compare_moran_spy,
        "save_table": _save_table_factory(output_dir),
        "display": lambda *_args, **_kwargs: None,
        "plot_moran_summary": lambda *_args, **_kwargs: None,
        "save_figure": lambda *_args, **_kwargs: None,
        "OUTPUT_DIR": output_dir,
        "CONFIG": {
            "route_kind": "test",
            "sample": {"id": "notebook-cell-fixture"},
            "spatial_region": {"coordinate_unit": "fixture", "microns_per_coordinate": 1.0},
        },
    }
    # Avoid taking the optional global branch in cells that support it.
    namespace.pop("GLOBAL_PARTITION")
    setup_source = next(
        source_setup
        for source_setup in _builder_source_cells()
        if "MORAN_SCOPE_ORDER = []" in source_setup
    )
    exec(compile(setup_source, "<builder-moran-source-cell>", "exec"), namespace)
    _execute_stage(source, namespace)

    assert len(calls) == 3
    assert all(graph.shape == (60, 60) and graph.nnz > 0 for graph in calls)
    moran = namespace["MORAN"]
    assert set(moran["scope"]) == {"Fibroblast", "Mono_Macro", "T"}
    moran = moran.loc[moran["scope"].eq("Fibroblast")]
    by_gene = moran.set_index("gene_id")
    assert by_gene.loc["constant_raw", "raw_status"] == "not_computable"
    assert by_gene.loc["constant_recon", "reconstruction_status"] == "not_computable"
    assert by_gene.loc["raw_only", "reconstruction_status"] == "unmeasured"
    assert by_gene.loc["recon_only", "raw_status"] == "unmeasured"
    assert by_gene["n_units"].eq(60).all()
    assert by_gene["n_edges_raw"].eq(calls[0].nnz).all()
    assert by_gene.loc["raw_only", "delta"] != by_gene.loc["raw_only", "delta"]

    moran_csv = output_dir / "analysis" / "moran" / "moran_by_cell_type.csv"
    assert moran_csv.exists(), "Moran cell did not publish its result CSV"
    reloaded = pd.read_csv(moran_csv)
    reloaded = reloaded.loc[reloaded["scope"].eq("Fibroblast")].set_index("gene_id")
    assert reloaded.loc["constant_raw", "raw_status"] == "not_computable"
    assert reloaded.loc["raw_only", "reconstruction_status"] == "unmeasured"


def test_builder_aucell_cell_runs_actual_helper_and_keeps_unmeasured_parent(tmp_path):
    pytest.importorskip("omicverse")
    pytest.importorskip("ctxcore")

    source = _scientific_cell(r"\bscore_gene_set_aucell\s*\(", "AUCell")
    from revise.analysis.advanced.aucell import (
        get_aucell_provider_metadata,
        score_gene_set_aucell,
    )
    from revise.analysis.basic.gene_set_scoring import read_gmt

    gmt_path = ROOT / "reproduce" / "case" / "pathway" / "h.all.v2025.1.Hs.symbols.gmt"
    emt_gene = str(read_gmt(gmt_path)["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"][0])
    raw, reconstruction = _expression_pair(
        raw_gene_ids=[emt_gene, "shared", "constant_raw", "raw_only"],
        reconstruction_gene_ids=[emt_gene, "shared", "constant_recon", "recon_only", "new_gene"],
        prefix="paired",
    )
    nohit_raw, nohit_reconstruction = _expression_pair(
        raw_gene_ids=["nohit_a", "constant_raw", "nohit_b"],
        reconstruction_gene_ids=["nohit_a", "constant_recon", "nohit_c"],
        prefix="nohit",
    )
    output_dir = Path(tmp_path)
    real_scorer_calls = []

    def actual_helper_with_audit(adata, genes, *, score_name, AUC_threshold, seed):
        real_scorer_calls.append({
            "unit_ids": tuple(str(value) for value in adata.obs_names),
            "gene_axis": tuple(str(value) for value in adata.var_names),
            "genes": tuple(genes),
            "AUC_threshold": AUC_threshold,
            "seed": seed,
            "X": adata.X.copy(),
        })
        return score_gene_set_aucell(
            adata,
            genes,
            score_name=score_name,
            AUC_threshold=AUC_threshold,
            seed=seed,
        )

    namespace = {
        "RUNS": {
            "Fibroblast": {"raw_expression": raw, "recon_expression": reconstruction},
            "NoHit": {"raw_expression": nohit_raw, "recon_expression": nohit_reconstruction},
        },
        "ROOT": ROOT,
        "OUTPUT_DIR": output_dir,
        "SEED": 42,
        "CONFIG": {"sample": {"id": "notebook-cell-fixture"}},
        "np": np,
        "pd": pd,
        "sc": sc if "sc" in locals() else None,
        "sparse": sparse,
        "score_gene_set_aucell": actual_helper_with_audit,
        "get_aucell_provider_metadata": get_aucell_provider_metadata,
        "display": lambda *_args, **_kwargs: None,
    }
    _execute_stage(source, namespace)

    aucell = namespace["AUCELL"]
    assert {"scope", "unit_id", "raw_score", "reconstruction_score", "delta", "status", "reason"}.issubset(aucell.columns)
    paired = aucell.loc[aucell["scope"].eq("Fibroblast")]
    assert len(paired) == 60
    assert paired["status"].eq("computed").all()
    np.testing.assert_allclose(
        paired["delta"].to_numpy(),
        (paired["reconstruction_score"] - paired["raw_score"]).to_numpy(),
    )
    nohit = aucell.loc[aucell["scope"].eq("NoHit")]
    assert nohit["raw_score"].isna().all()
    assert nohit["reconstruction_score"].isna().all()
    assert nohit["status"].eq("unavailable").all()
    assert nohit["reason"].str.contains("resource_genes_unmeasured").all()

    # The real helper receives each full gene axis and the fixed provider
    # parameters.  No scorer call should be made for the no-hit parent.
    assert len(real_scorer_calls) == 2
    assert all(call["unit_ids"][0].startswith("paired-") for call in real_scorer_calls)
    assert all(call["AUC_threshold"] == pytest.approx(0.01) for call in real_scorer_calls)
    assert all(call["seed"] == 42 for call in real_scorer_calls)
    assert tuple(raw.var_names.astype(str)) in {call["gene_axis"] for call in real_scorer_calls}
    assert tuple(reconstruction.var_names.astype(str)) in {call["gene_axis"] for call in real_scorer_calls}

    summaries = namespace["aucell_summary"].set_index("scope")
    assert summaries.loc["Fibroblast", "paired_n"] == 60
    assert summaries.loc["NoHit", "paired_n"] == 0
    assert summaries.loc["NoHit", "comparison_status"] == "unavailable"

    score_paths = list(output_dir.rglob("scores.csv.gz"))
    resource_paths = list(output_dir.rglob("resource.json"))
    availability_paths = list(output_dir.rglob("availability.csv"))
    assert score_paths and resource_paths and availability_paths
    scores_csv = pd.read_csv(score_paths[0])
    assert len(scores_csv) == len(aucell)
    availability = pd.read_csv(availability_paths[0])
    assert availability.loc[availability["scope"].eq("NoHit"), "available_gene_count"].eq(0).all()
    resource = json.loads(resource_paths[0].read_text(encoding="utf-8"))
    assert resource["auc_threshold_quantile"] == pytest.approx(0.01)
    assert resource["seed"] == 42
    assert resource["resource_sha256"]
    side_cutoffs = resource["side_cutoffs"]
    if isinstance(side_cutoffs, dict):
        side_cutoffs = list(side_cutoffs.values())
    assert side_cutoffs
    assert all(int(row["effective_rank_length"]) > 0 for row in side_cutoffs if row["status"] == "computed")


def test_moran_summaries_distinguish_full_and_shared_genes(tmp_path):
    """Unequal gene spaces must not change the denominator of paired summaries."""
    builder = importlib.import_module(BUILDER_MODULE)
    rows = []
    for scope in ("Fibroblast", "T"):
        for gene, raw, recon in (("a", .1, .4), ("b", .5, .3), ("raw_only", .9, np.nan), ("recon_only", np.nan, -.2)):
            if scope == "T":
                recon = np.nan
            rows.append({"scope": scope, "comparison_id": scope, "gene_id": gene,
                         "raw_moran": raw, "reconstruction_moran": recon,
                         "raw_status": "computed" if np.isfinite(raw) else "unmeasured",
                         "reconstruction_status": "computed" if np.isfinite(recon) else "unmeasured",
                         "raw_reason": "", "reconstruction_reason": "",
                         "comparison_status": "computed" if np.isfinite(raw + recon) else "unavailable",
                         "delta": recon - raw})
    raw, recon = _expression_pair(raw_gene_ids=["a", "b", "raw_only"], reconstruction_gene_ids=["a", "b", "recon_only"], prefix="summary")
    namespace = {"np": np, "pd": pd, "MORAN": pd.DataFrame(rows),
                 "MORAN_SCOPE_ORDER": ["Fibroblast", "T"],
                 "MORAN_SOURCES": {s: (raw, recon, raw) for s in ("Fibroblast", "T")},
                 "save_table": _save_table_factory(tmp_path), "display": lambda *_: None}
    namespace["figures"] = SimpleNamespace(save_table=namespace["save_table"])
    exec(builder.MORAN_AUDIT_CELL, namespace)
    summary = namespace["MORAN_SUMMARY"].set_index("scope")
    assert summary.loc["Fibroblast", "raw_nvalid"] == 3
    assert summary.loc["Fibroblast", "reconstruction_nvalid"] == 3
    assert summary.loc["Fibroblast", "raw_q1"] == pytest.approx(.3)
    assert summary.loc["Fibroblast", "reconstruction_q1"] == pytest.approx(.05)
    assert summary.loc["Fibroblast", "matched_gene_n"] == 2
    assert summary.loc["Fibroblast", "shared_raw_median"] == pytest.approx(.3)
    assert summary.loc["Fibroblast", "shared_reconstruction_median"] == pytest.approx(.35)
    assert summary.loc["Fibroblast", "paired_delta_median"] == pytest.approx(.05)
    assert summary.loc["Fibroblast", "paired_delta_q1"] == pytest.approx(-.075)
    assert summary.loc["Fibroblast", "paired_delta_q3"] == pytest.approx(.175)
    assert summary.loc["T", "matched_gene_n"] == 0
    assert pd.isna(summary.loc["T", "paired_delta_median"])
    assert pd.isna(summary.loc["T", "shared_raw_median"])
    compact = pd.read_csv(tmp_path / "analysis/moran/moran_distribution_summary.csv")
    assert len(compact) == 8
    assert not compact.duplicated(["scope", "gene_set", "side"]).any()
    empty = compact.loc[compact.scope.eq("T") & compact.gene_set.eq("shared_valid")]
    assert empty.n_valid.eq(0).all()
    assert empty[["median", "q1", "q75"]].isna().all().all()
