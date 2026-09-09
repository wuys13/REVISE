from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

from revise.analysis import bio
from revise.backend.runners import sc_svc_application


def _adata() -> AnnData:
    obs = pd.DataFrame(
        {
            "SVC_cluster": pd.Categorical(["0", "1", "0", "1"]),
            "Level2": pd.Categorical(["A", "B", "A", "B"]),
        },
        index=["cell-0", "cell-1", "cell-2", "cell-3"],
    )
    return AnnData(
        X=np.array(
            [
                [1.0, 0.0],
                [0.0, 1.0],
                [2.0, 0.0],
                [0.0, 2.0],
            ]
        ),
        obs=obs,
        var=pd.DataFrame(index=["G1", "G2"]),
    )


def test_constructor_keeps_input_objects_and_eagerly_computes_degs(monkeypatch):
    spatial = _adata()
    expression = _adata()
    expected = pd.DataFrame({"group": ["0"], "gene": ["G1"]})
    calls = []

    def fake_get_degs(adata, **kwargs):
        calls.append((adata, kwargs))
        return expected

    monkeypatch.setattr(sc_svc_application, "get_degs", fake_get_degs)

    analysis = sc_svc_application.ScSVCAnalysis(
        spatial, expression, "SVC_cluster"
    )

    assert analysis.sc_SVC_adata_spatial is spatial
    assert analysis.sc_SVC_adata_expr is expression
    assert analysis.sc_SVC_degs is expected
    assert calls == [
        (
            expression,
            {
                "groupby": "SVC_cluster",
                "method": "t-test",
                "fc_threshold": None,
            },
        )
    ]


def test_empty_target_cluster_keeps_legacy_full_data_behavior(monkeypatch):
    expression = _adata()
    calls = []

    def fake_get_degs(adata, **kwargs):
        calls.append((adata, kwargs))
        return pd.DataFrame()

    monkeypatch.setattr(sc_svc_application, "get_degs", fake_get_degs)
    analysis = sc_svc_application.ScSVCAnalysis(_adata(), expression, "SVC_cluster")
    calls.clear()

    analysis.get_svc_degs([], fc_threshold=0)

    assert calls == [
        (
            expression,
            {
                "groupby": "SVC_cluster",
                "method": "t-test",
                "fc_threshold": 0,
            },
        )
    ]


def test_pathway_normalize_keeps_legacy_unused_subset_call_sequence(monkeypatch):
    expression = _adata()
    monkeypatch.setattr(
        sc_svc_application,
        "get_degs",
        lambda *args, **kwargs: pd.DataFrame(),
    )
    analysis = sc_svc_application.ScSVCAnalysis(_adata(), expression, "SVC_cluster")
    events = []
    expected_degs = pd.DataFrame({"group": ["0"], "gene": ["G1"]})
    expected_pathways = pd.DataFrame({"Term": ["P1"]})

    monkeypatch.setattr(
        sc_svc_application.sc.pp,
        "normalize_total",
        lambda adata, **kwargs: events.append(("normalize", adata, kwargs)),
    )
    monkeypatch.setattr(
        sc_svc_application.sc.pp,
        "log1p",
        lambda adata: events.append(("log1p", adata)),
    )

    def fake_get_svc_degs(cluster_nums, fc_threshold):
        events.append(("degs", cluster_nums, fc_threshold))
        return expected_degs

    def fake_conclusions(deg_df, geneset_file, **kwargs):
        events.append(("conclusions", deg_df, geneset_file, kwargs))
        return expected_pathways

    monkeypatch.setattr(analysis, "get_svc_degs", fake_get_svc_degs)
    monkeypatch.setattr(sc_svc_application, "conclusions_write", fake_conclusions)

    result = analysis.get_pathway_conclusion(
        ["0"],
        fc_threshold=1,
        pathway_num=5,
        gene_num=60,
        geneset_file=["hallmark"],
        normalize=True,
    )

    assert result is expected_pathways
    assert [event[0] for event in events] == [
        "normalize",
        "log1p",
        "degs",
        "conclusions",
    ]
    normalized_subset = events[0][1]
    assert normalized_subset is events[1][1]
    assert normalized_subset is not expression
    assert events[2] == ("degs", ["0"], 1)
    assert events[3][1] is expected_degs


def test_conclusions_write_resolves_runtime_bio_enrichment_patch(monkeypatch):
    degs = pd.DataFrame(
        {
            "group": ["0", "0", "1"],
            "gene": ["G1", "G2", "G3"],
        }
    )
    calls = []

    def patched_enrichment(genes, geneset_file, cutoff=0.05):
        calls.append((genes, geneset_file, cutoff))
        return pd.DataFrame(
            {"Term": [f"term-{genes[0]}"], "Adjusted P-value": [0.01]}
        )

    monkeypatch.setattr(bio, "get_enrichment", patched_enrichment)

    result = bio.conclusions_write(
        degs,
        "genes.gmt",
        gene_num=1,
        pathway_num=1,
        print_flag=False,
    )

    assert calls == [(["G1"], "genes.gmt", 0.05), (["G3"], "genes.gmt", 0.05)]
    assert result[["Term", "group"]].to_dict("records") == [
        {"Term": "term-G1", "group": "0"},
        {"Term": "term-G3", "group": "1"},
    ]
