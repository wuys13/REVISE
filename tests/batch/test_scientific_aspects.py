"""Real impact calculations across the batch handoff, using a solver stand-in.

This verifies publication/reuse and scientific table plumbing, not reconstruction
quality or biological validity.
"""

import json

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import sparse

from test_runner import batch_config, fake_solver as _fake_solver, make_sample, result_root

fake_solver = _fake_solver


def test_hst_impact_publishes_real_calculations_and_reuses(tmp_path, fake_solver):
    from revise.batch.analysis import run_analysis_batch
    from revise.batch.runner import run_batch

    root = tmp_path / "data" / "one"
    make_sample(root, modality="hST")
    rng = np.random.default_rng(19)
    counts = rng.poisson(3, size=(24, 6)).astype(float) + 1
    counts[:12, :3] += 15
    counts[12:, 3:] += 15
    data = ad.AnnData(
        sparse.csr_matrix(counts),
        obs=pd.DataFrame(
            {"Level1": ["T"] * 20 + ["Tumor"] * 4,
             "Level2": ["a"] * 12 + ["b"] * 12},
            index=[f"u{i}" for i in range(24)],
        ),
        var=pd.DataFrame(index=[f"g{i}" for i in range(6)]),
    )
    data.obsm["spatial"] = np.column_stack((np.arange(24) % 6, np.arange(24) // 6))
    data.write_h5ad(root / "spatial.h5ad")
    data.write_h5ad(root / "reference.h5ad")
    config = batch_config(tmp_path)
    document = yaml.safe_load(config.read_text())
    document["analysis"] = {
        "reconstruction_impact": {
            "entrypoint": "revise.analysis.reconstruction_impact_adapter:run",
            "version": "1",
            "requires_pairing": True,
            "resources": {"raw_level2_reference": str(root / "reference.h5ad")},
            "parameters": {
                "level1_column": "Level1",
                "parent_labels": {"T": "T"},
                "partition_change": {
                    "mode": "level1_ari",
                    "level1_resolution_candidates": [0.3, 0.5],
                    "within_level1_resolution": 0.5,
                    "random_state": 19,
                    "n_top_genes": 6,
                    "raw_qc_min_genes": 1,
                    "raw_qc_min_cells": 1,
                    "sample_n_units": 24,
                    "parent_sample_n_units": 20,
                },
                "spatial_region": {
                    "microns_per_coordinate": 1.0,
                    "candidate_window_sides_um": [3.0],
                    "min_parent_units": 2,
                    "rarefaction_draws": 2,
                    "threshold_bootstraps": 2,
                    "cell_equivalent_um": 1.0,
                },
                "raw_level2_mapping": {
                    "method": "pot", "level1_column": "Level1", "level2_column": "Level2",
                    "pot": {"reg": 0.1, "reg_m": 0.0, "reg_type": "entropy"},
                },
            },
        },
    }
    config.write_text(yaml.safe_dump(document))
    assert run_batch(config)["summary"]["succeeded"] == 1
    result = run_analysis_batch(config)
    target = result_root(root)
    log = target / ".revise" / "analysis" / "reconstruction_impact.log"
    assert result["summary"]["succeeded"] == 1, log.read_text()
    output = target / "analysis" / "reconstruction_impact"
    comparisons = pd.read_csv(output / "comparisons.csv")
    summary = pd.read_csv(output / "partition" / "summary.csv")
    assert set(summary["scope"]) == {"global", "T"}
    assert set(summary["comparison_id"]) <= set(comparisons["comparison_id"])
    assert comparisons["comparison_id"].is_unique
    assert "multiple" not in set(comparisons["comparison_edge"])
    windows = pd.read_csv(output / "spatial" / "window_metrics.csv")
    assert set(windows["scope"]) == {"T"}
    assert {"raw_leiden", "raw_level2"} <= set(windows["baseline"])
    raw_neff = windows.loc[(windows.baseline == "raw_leiden") & (windows.metric == "neff"), "raw"]
    assert raw_neff.gt(1).any(), "Raw baseline must retain Leiden diversity within the T parent"
    audit = json.loads((output / "audit.json").read_text())
    assert audit["scope_audits"]["global"]["spatial_status"] == "not_applicable_parent_baseline"
    assert run_analysis_batch(config)["summary"]["reused"] == 1
    assert run_batch(config)["summary"]["reused"] == 1
    from revise.analysis.impact_report import render_impact_report

    report = render_impact_report([target], target / "reports/reconstruction_impact")
    assert report["status"] == "succeeded", report.get("error")
    assert report["tasks"][0]["task_id"] == "one"
    assert report["execution"]["status"] == "succeeded"
    assert any(item["role"].startswith("figure:") for item in report["artifacts"])


def test_ist_pathway_projects_before_scoring_and_publishes(tmp_path, fake_solver, monkeypatch):
    from revise.analysis import pathway_activity_adapter as pathway
    from revise.batch import analysis
    from revise.batch.analysis import run_analysis_batch
    from revise.batch.runner import run_batch

    root = tmp_path / "data" / "one"
    make_sample(root, cell_types=["T"])
    resource = tmp_path / "gene_sets.json"
    resource.write_text(json.dumps({"measured": ["0"], "absent": ["missing"]}))
    config = batch_config(tmp_path)
    document = yaml.safe_load(config.read_text())
    document["analysis"] = {
        "pathway_activity": {
            "entrypoint": "revise.analysis.pathway_activity_adapter:run",
            "version": "1", "requires_pairing": True,
            "resources": {"gene_sets": str(resource)},
            "parameters": {
                "resource": {"name": "fixture-pathways", "species": "human", "gene_id_type": "symbol", "version": "fixture"},
                "level1_column": "Level1", "parent_labels": {"T": "T"},
                "min_coverage": 0.5, "min_genes": 1, "min_detected_genes": 1,
                "seed": 17, "cutoff_policy": "provider_detected_gene_quantile",
                "auc_threshold_quantile": 0.05,
            },
        },
    }
    config.write_text(yaml.safe_dump(document))
    calls = []

    def score(adata, genes, *, score_name, AUC_threshold, seed):
        calls.append((adata.obs_names.tolist(), adata.var_names.tolist()))
        result = adata.copy()
        result.obs[f"{score_name}_aucell"] = 0.5
        return result, f"{score_name}_aucell"

    # Attest production source before installing the deliberate provider stand-in.
    analysis._attest_framework_helpers()
    analysis._attest_adapter_module(pathway.__name__)
    monkeypatch.setattr(pathway, "score_gene_set_aucell", score)
    monkeypatch.setattr(pathway, "get_aucell_provider_metadata", lambda: {"name": "fixture", "version": "1"})
    assert run_batch(config)["summary"]["succeeded"] == 1
    result = run_analysis_batch(config)
    target = result_root(root) / "T"
    assert result["summary"]["succeeded"] == 1, (target / ".revise/analysis/pathway_activity.log").read_text()
    assert calls == [(["u1"], ["0", "1", "2"])] * 2
    output = target / "analysis/pathway_activity"
    scores = pd.read_csv(output / "scores.csv.gz")
    absent = scores.loc[scores.pathway == "absent"]
    assert absent.raw_status.eq("unmeasured").all()
    assert absent.raw_score.isna().all()
    assert absent.reconstruction_score.isna().all()
    registry = pd.read_csv(output / "comparisons.csv")
    assert registry.reconstruction_view.eq("cluster_mean_projection").all()
    assert run_analysis_batch(config)["summary"]["reused"] == 1
    assert len(calls) == 2
