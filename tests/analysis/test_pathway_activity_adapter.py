from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData


def _view(matrix, genes, labels, *, mapping="native_id", source="carrier.h5ad", provenance=None):
    n_units = len(labels)
    adata = AnnData(
        np.asarray(matrix, dtype=float),
        obs=pd.DataFrame({"Level1": labels}, index=[f"unit-{index + 1}" for index in range(n_units)]),
        var=pd.DataFrame(index=genes),
        obsm={"spatial": np.column_stack([np.arange(n_units), np.arange(n_units)])},
    )
    return SimpleNamespace(
        adata=adata,
        mapping=mapping,
        source=source,
        provenance=dict(provenance or {}),
        uncovered_observations=(),
    )


def _context(tmp_path: Path, raw, reconstruction, resource_mapping, **parameters):
    resource_path = tmp_path / "gene_sets.json"
    resource_path.write_text(json.dumps(resource_mapping), encoding="utf-8")
    calls = []

    class Inputs:
        def aligned_raw(self):
            calls.append("aligned_raw")
            return raw

        def reconstruction_expression(self, *, strategy=None):
            calls.append("reconstruction_expression")
            context.expression_strategies.append(strategy)
            return reconstruction

    context = SimpleNamespace(
        output_dir=tmp_path / "output",
        parameters=parameters,
        resources={"gene_sets": resource_path},
        inputs=Inputs(),
        reconstruction={
            "sample_id": "sample-1",
            "cell_type": "T",
            "modality": "hST",
            "expression_semantics": "fixture",
        },
    )
    context.inputs_calls = calls
    context.expression_strategies = []
    return context


def _parameters(**overrides):
    values = {
        "resource": {
            "name": "fixture-pathways",
            "species": "human",
            "gene_id_type": "symbol",
            "version": "fixture-v1",
        },
        "min_coverage": 0.5,
        "min_genes": 1,
        "min_detected_genes": 1,
        "seed": 17,
        "cutoff_policy": "provider_detected_gene_quantile",
        "auc_threshold_quantile": 0.05,
        "level1_column": "Level1",
    }
    values.update(overrides)
    return values


def test_pathway_adapter_preserves_all_units_and_pathways_with_independent_gene_spaces(
    tmp_path, monkeypatch
):
    from revise.analysis import pathway_activity_adapter as adapter

    raw = _view([[2, 0], [0, 0]], ["G1", "G2"], ["Raw-A", "Raw-B"])
    reconstruction = _view(
        [[2, 4, 0], [0, 0, 3]],
        ["G1", "G2", "G3"],
        ["Recon-A", "Recon-B"],
        mapping="cluster_mean_projection",
        provenance={"expression_view": "cluster_mean_projection"},
    )
    context = _context(
        tmp_path,
        raw,
        reconstruction,
        {"pathway-1": ["G1", "G3"], "unmeasured": ["MISSING"]},
        **_parameters(),
    )
    provider_calls = []

    def fake_provider(adata, genes, *, score_name, AUC_threshold=None, seed=None):
        provider_calls.append((score_name, list(genes), AUC_threshold, seed, list(adata.var_names)))
        scored = adata.copy()
        scored.obs[f"{score_name}_aucell"] = np.arange(scored.n_obs, dtype=float) + 1.0
        return scored, f"{score_name}_aucell"

    monkeypatch.setattr(adapter, "score_gene_set_aucell", fake_provider)
    monkeypatch.setattr(
        adapter,
        "get_aucell_provider_metadata",
        lambda: {"name": "omicverse", "version": "1.7.5", "scorer": "single.geneset_aucell"},
    )

    result = adapter.run(context)

    assert context.inputs_calls == ["aligned_raw", "reconstruction_expression"]
    assert context.expression_strategies == ["native"]
    assert set(result["artifacts"]) == {
        "audit",
        "comparisons",
        "observations",
        "genes",
        "mapping",
        "availability",
        "scores",
        "summary",
        "resource",
        "resource_genes",
    }
    scores = pd.read_csv(context.output_dir / "scores.csv.gz")
    assert len(scores) == 2 * 2
    assert set(scores["pathway"]) == {"pathway-1", "unmeasured"}
    assert set(scores["unit_id"]) == {"unit-1", "unit-2"}
    pathway_scores = scores[scores["pathway"] == "pathway-1"]
    assert set(pathway_scores["raw_status"]) == {"computed"}
    assert set(pathway_scores["reconstruction_status"]) == {"computed"}
    missing_scores = scores[scores["pathway"] == "unmeasured"]
    assert set(missing_scores["raw_status"]) == {"unmeasured"}
    assert set(missing_scores["reconstruction_status"]) == {"unmeasured"}
    assert missing_scores["raw_score"].isna().all()
    assert missing_scores["reconstruction_score"].isna().all()

    comparisons = pd.read_csv(context.output_dir / "comparisons.csv")
    assert comparisons["scope"].tolist() == ["Raw-A", "Raw-B"]
    assert comparisons["label_source"].eq("raw.Level1").all()
    assert comparisons["reconstruction_view"].eq("cluster_mean_projection").all()
    assert len(pd.read_csv(context.output_dir / "summary.csv")) == 4
    assert provider_calls == [
        ("pathway-1", ["G1"], 0.05, 17, ["G1", "G2"]),
        ("pathway-1", ["G1", "G3"], 0.05, 17, ["G1", "G2", "G3"]),
    ]

    resource = json.loads((context.output_dir / "resource.json").read_text())
    assert resource["species"] == "human"
    assert resource["provider"]["version"] == "1.7.5"
    assert resource["cutoff_policy"] == "provider_detected_gene_quantile"
    assert resource["auc_threshold_quantile"] == 0.05
    assert resource["side_cutoffs"]["raw"]["n_vars"] == 2
    assert resource["side_cutoffs"]["reconstruction"]["n_vars"] == 3
    assert resource["side_cutoffs"]["raw"]["provider_auc_threshold"] == pytest.approx(0.025)
    assert resource["side_cutoffs"]["reconstruction"]["provider_auc_threshold"] == pytest.approx(0.35)
    assert "NaN" not in (context.output_dir / "audit.json").read_text()


def test_pathway_adapter_marks_zero_expression_as_not_computable(tmp_path, monkeypatch):
    from revise.analysis import pathway_activity_adapter as adapter

    raw = _view([[0, 0], [0, 0]], ["G1", "G2"], ["Raw-A", "Raw-B"])
    reconstruction = _view([[1, 0], [0, 1]], ["G1", "G2"], ["Recon-A", "Recon-B"])
    context = _context(
        tmp_path,
        raw,
        reconstruction,
        {"pathway": ["G1", "G2"]},
        **_parameters(),
    )
    def fake_provider(adata, genes, *, score_name, AUC_threshold=None, seed=None):
        scored = adata.copy()
        scored.obs[f"{score_name}_aucell"] = np.ones(scored.n_obs)
        return scored, f"{score_name}_aucell"

    monkeypatch.setattr(adapter, "score_gene_set_aucell", fake_provider)
    monkeypatch.setattr(
        adapter,
        "get_aucell_provider_metadata",
        lambda: {"name": "omicverse", "version": "1.7.5", "scorer": "single.geneset_aucell"},
    )

    adapter.run(context)

    availability = pd.read_csv(context.output_dir / "availability.csv")
    raw_row = availability.query("side == 'raw'").iloc[0]
    assert raw_row["status"] == "not_computable"
    assert raw_row["reason"] == "zero_expression_no_detectable_signature"
    scores = pd.read_csv(context.output_dir / "scores.csv.gz")
    assert scores["raw_score"].isna().all()
    assert scores["reconstruction_score"].notna().all()


def test_pathway_adapter_marks_zero_detected_gene_rank_cutoff_not_computable(
    tmp_path, monkeypatch
):
    from revise.analysis import pathway_activity_adapter as adapter

    n_units = 20
    raw = _view(
        np.vstack([np.zeros((2, 2)), np.tile([1.0, 0.0], (n_units - 2, 1))]),
        ["G1", "G2"],
        ["Raw-A"] * n_units,
    )
    reconstruction = _view(
        np.tile([1.0, 0.0], (n_units, 1)), ["G1", "G2"], ["Recon-A"] * n_units
    )
    context = _context(
        tmp_path,
        raw,
        reconstruction,
        {"pathway": ["G1"]},
        **_parameters(),
    )
    provider_calls = []

    def fake_provider(adata, genes, *, score_name, AUC_threshold=None, seed=None):
        provider_calls.append(list(adata.obs_names))
        return _scored_copy(adata, score_name)

    monkeypatch.setattr(adapter, "score_gene_set_aucell", fake_provider)
    monkeypatch.setattr(
        adapter,
        "get_aucell_provider_metadata",
        lambda: {"name": "omicverse", "version": "1.7.5", "scorer": "fixture"},
    )

    adapter.run(context)

    availability = pd.read_csv(context.output_dir / "availability.csv")
    raw_row = availability.query("side == 'raw'").iloc[0]
    reconstruction_row = availability.query("side == 'reconstruction'").iloc[0]
    assert raw_row["status"] == "not_computable"
    assert raw_row["reason"] == "zero_detected_gene_rank_cutoff"
    assert reconstruction_row["status"] == "computed"
    scores = pd.read_csv(context.output_dir / "scores.csv.gz")
    assert scores["raw_score"].isna().all()
    assert scores["reconstruction_score"].notna().all()
    assert provider_calls == [[f"unit-{index + 1}" for index in range(n_units)]]


def test_pathway_adapter_records_explicit_zero_signature_policy(tmp_path, monkeypatch):
    from revise.analysis import pathway_activity_adapter as adapter

    raw = _view([[0, 5], [0, 4]], ["G1", "G2"], ["Raw-A", "Raw-B"])
    reconstruction = _view([[0, 5], [0, 4]], ["G1", "G2"], ["Recon-A", "Recon-B"])
    context = _context(
        tmp_path,
        raw,
        reconstruction,
        {"pathway": ["G1"]},
        **_parameters(min_detected_genes=0),
    )
    calls = []

    def fake_provider(adata, genes, *, score_name, AUC_threshold=None, seed=None):
        calls.append((score_name, list(genes), AUC_threshold, seed))
        return _scored_copy(adata, score_name)

    monkeypatch.setattr(adapter, "score_gene_set_aucell", fake_provider)
    monkeypatch.setattr(
        adapter,
        "get_aucell_provider_metadata",
        lambda: {"name": "omicverse", "version": "1.7.5", "scorer": "fixture"},
    )

    adapter.run(context)

    availability = pd.read_csv(context.output_dir / "availability.csv")
    assert set(availability["status"]) == {"computed"}
    assert set(availability["detected_signature_gene_count"]) == {0}
    assert calls == [
        ("pathway", ["G1"], 0.05, 17),
        ("pathway", ["G1"], 0.05, 17),
    ]


def test_pathway_adapter_marks_below_minimum_detected_signature_support(tmp_path, monkeypatch):
    from revise.analysis import pathway_activity_adapter as adapter

    raw = _view([[0, 5], [0, 4]], ["G1", "G2"], ["Raw-A", "Raw-B"])
    reconstruction = _view([[0, 5], [0, 4]], ["G1", "G2"], ["Recon-A", "Recon-B"])
    context = _context(
        tmp_path,
        raw,
        reconstruction,
        {"pathway": ["G1"]},
        **_parameters(min_detected_genes=1),
    )
    monkeypatch.setattr(
        adapter,
        "get_aucell_provider_metadata",
        lambda: {"name": "omicverse", "version": "1.7.5", "scorer": "fixture"},
    )

    def provider_should_not_run(*_args, **_kwargs):
        raise AssertionError("provider must not run below min_detected_genes")

    monkeypatch.setattr(adapter, "score_gene_set_aucell", provider_should_not_run)
    adapter.run(context)

    availability = pd.read_csv(context.output_dir / "availability.csv")
    assert set(availability["status"]) == {"insufficient_support"}
    assert set(availability["reason"]) == {"insufficient_detected_signature_genes"}
    scores = pd.read_csv(context.output_dir / "scores.csv.gz")
    assert scores["raw_score"].isna().all()
    assert scores["reconstruction_score"].isna().all()


def test_pathway_adapter_does_not_hide_provider_failures(tmp_path, monkeypatch):
    from revise.analysis import pathway_activity_adapter as adapter

    raw = _view([[1, 0], [0, 1]], ["G1", "G2"], ["Raw-A", "Raw-B"])
    reconstruction = _view([[1, 0], [0, 1]], ["G1", "G2"], ["Recon-A", "Recon-B"])
    context = _context(
        tmp_path,
        raw,
        reconstruction,
        {"pathway": ["G1"]},
        **_parameters(),
    )
    monkeypatch.setattr(
        adapter,
        "score_gene_set_aucell",
        lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("provider exploded")),
    )
    monkeypatch.setattr(
        adapter,
        "get_aucell_provider_metadata",
        lambda: {"name": "omicverse", "version": "1.7.5", "scorer": "single.geneset_aucell"},
    )

    with pytest.raises(RuntimeError, match="provider exploded"):
        adapter.run(context)


def test_pathway_adapter_requires_explicit_scientific_parameters(tmp_path):
    from revise.analysis import pathway_activity_adapter as adapter

    raw = _view([[1, 0], [0, 1]], ["G1", "G2"], ["Raw-A", "Raw-B"])
    reconstruction = _view([[1, 0], [0, 1]], ["G1", "G2"], ["Recon-A", "Recon-B"])
    context = _context(
        tmp_path,
        raw,
        reconstruction,
        {"pathway": ["G1"]},
        **{key: value for key, value in _parameters().items() if key != "auc_threshold_quantile"},
    )

    with pytest.raises(ValueError, match="auc_threshold_quantile"):
        adapter.run(context)


def test_pathway_adapter_requires_a_named_resource_identity(tmp_path, monkeypatch):
    from revise.analysis import pathway_activity_adapter as adapter

    raw = _view([[1, 0], [0, 1]], ["G1", "G2"], ["Raw-A", "Raw-B"])
    reconstruction = _view(
        [[1, 0], [0, 1]], ["G1", "G2"], ["Recon-A", "Recon-B"]
    )
    resource = dict(_parameters()["resource"])
    del resource["name"]
    context = _context(
        tmp_path,
        raw,
        reconstruction,
        {"pathway": ["G1"]},
        **_parameters(resource=resource),
    )
    monkeypatch.setattr(
        adapter,
        "get_aucell_provider_metadata",
        lambda: {"name": "omicverse", "version": "1.7.5", "scorer": "fixture"},
    )
    monkeypatch.setattr(
        adapter,
        "score_gene_set_aucell",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("resource identity must be validated before scoring")
        ),
    )

    with pytest.raises(ValueError, match="resource.name"):
        adapter.run(context)


def test_pathway_comparison_ids_include_exact_sample_task_and_raw_scope():
    from revise.analysis import pathway_activity_adapter as adapter

    raw = _view([[1, 0], [0, 1]], ["G1", "G2"], ["A::B", "A::B"])
    reconstruction = _view([[1, 0], [0, 1]], ["G1", "G2"], ["ignored", "ignored"])
    context = SimpleNamespace(
        reconstruction={"sample_id": "sample::1", "cell_type": "T::cell"}
    )

    comparisons = adapter._comparison_rows(
        context, raw, reconstruction, ["A::B"], "Level1"
    )

    assert comparisons[0]["comparison_id"] == (
        '{"raw_level1":"A::B","sample_id":"sample::1","task_cell_type":"T::cell"}'
    )


def test_pathway_adapter_accepts_only_provider_quantiles(tmp_path):
    from revise.analysis import pathway_activity_adapter as adapter

    raw = _view([[1, 0], [0, 1]], ["G1", "G2"], ["Raw-A", "Raw-B"])
    reconstruction = _view([[1, 0], [0, 1]], ["G1", "G2"], ["Recon-A", "Recon-B"])
    context = _context(
        tmp_path,
        raw,
        reconstruction,
        {"pathway": ["G1"]},
        **_parameters(auc_threshold_quantile=0.02),
    )

    with pytest.raises(ValueError, match="provider quantiles"):
        adapter.run(context)


def test_pathway_adapter_requires_exact_parent_label_for_ist(tmp_path, monkeypatch):
    from revise.analysis import pathway_activity_adapter as adapter

    raw = _view([[1, 0], [0, 1]], ["G1", "G2"], ["Mono/Macro", "Mono/Macro"])
    reconstruction = _view([[1, 0], [0, 1]], ["G1", "G2"], ["ignored", "ignored"])
    context = _context(
        tmp_path,
        raw,
        reconstruction,
        {"pathway": ["G1"]},
        **_parameters(parent_labels={"Mono_Macro": "Mono/Macro"}),
    )
    context.reconstruction["modality"] = "iST"
    context.reconstruction["cell_type"] = "Mono_Macro"
    monkeypatch.setattr(
        adapter,
        "get_aucell_provider_metadata",
        lambda: {"name": "omicverse", "version": "1.7.5", "scorer": "fixture"},
    )
    monkeypatch.setattr(
        adapter,
        "score_gene_set_aucell",
        lambda adata, genes, *, score_name, AUC_threshold, seed: _scored_copy(
            adata, score_name
        ),
    )

    # The exact mapped Raw label is accepted, while a normalized alias is not.
    adapter.run(context)
    context.parameters["parent_labels"] = {"Mono_Macro": "Mono_Macro"}
    with pytest.raises(ValueError, match="exact Raw Level1"):
        adapter.run(context)


def _scored_copy(adata, score_name):
    scored = adata.copy()
    scored.obs[f"{score_name}_aucell"] = np.ones(scored.n_obs)
    return scored, f"{score_name}_aucell"
