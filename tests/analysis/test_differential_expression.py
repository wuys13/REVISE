from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

from revise.analysis import bio
from revise.analysis.basic import differential_expression


def _adata() -> AnnData:
    return AnnData(
        X=np.array(
            [
                [1.0, 0.0, 2.0],
                [0.0, 1.0, 2.0],
                [1.0, 1.0, 0.0],
            ]
        ),
        obs=pd.DataFrame(
            {"cluster": ["A", "B", "reference"]},
            index=["cell-a", "cell-b", "cell-reference"],
        ),
        var=pd.DataFrame(index=["G1", "G2", "G3"]),
    )


def _install_rank_result(monkeypatch, events):
    def normalize_total(adata, **kwargs):
        events.append(("normalize_total", adata, kwargs))
        adata.X = np.asarray(adata.X) + 10

    def log1p(adata):
        events.append(("log1p", adata))

    def rank_genes_groups(adata, **kwargs):
        events.append(("rank_genes_groups", adata, kwargs))
        adata.uns["rank_genes_groups"] = {
            "names": {
                "A": np.array(["G1", "G2", "G3"]),
                "B": np.array(["G2", "G1", "G3"]),
                "reference": np.array(["G3", "G1", "G2"]),
            },
            "logfoldchanges": {
                "A": np.array([2.0, -2.0, 0.0]),
                "B": np.array([1.5, -1.5, 0.0]),
                "reference": np.array([3.0, 2.0, 1.0]),
            },
            "pvals": {
                "A": np.array([0.01, 0.02, 0.5]),
                "B": np.array([0.03, 0.04, 0.6]),
                "reference": np.array([0.01, 0.01, 0.01]),
            },
            "pvals_adj": {
                "A": np.array([0.02, 0.04, 0.8]),
                "B": np.array([0.03, 0.05, 0.9]),
                "reference": np.array([0.02, 0.02, 0.02]),
            },
        }

    monkeypatch.setattr(
        differential_expression.sc.pp, "normalize_total", normalize_total
    )
    monkeypatch.setattr(differential_expression.sc.pp, "log1p", log1p)
    monkeypatch.setattr(
        differential_expression.sc.tl, "rank_genes_groups", rank_genes_groups
    )


def test_get_degs_uses_a_normalized_copy_and_preserves_result_contract(monkeypatch):
    adata = _adata()
    original_x = adata.X.copy()
    original_cluster = adata.obs["cluster"].copy()
    events = []
    _install_rank_result(monkeypatch, events)

    result = differential_expression.get_degs(
        adata,
        groupby="cluster",
        method="wilcoxon",
        reference="reference",
    )

    assert [event[0] for event in events] == [
        "normalize_total",
        "log1p",
        "rank_genes_groups",
    ]
    working_copy = events[0][1]
    assert working_copy is events[1][1] is events[2][1]
    assert working_copy is not adata
    assert events[0][2] == {"target_sum": 1e4}
    assert events[2][2] == {
        "reference": "reference",
        "groupby": "cluster",
        "method": "wilcoxon",
        "use_raw": False,
    }
    np.testing.assert_array_equal(adata.X, original_x)
    pd.testing.assert_series_equal(adata.obs["cluster"], original_cluster)
    assert result.columns.tolist() == [
        "group",
        "gene",
        "logfoldchanges",
        "pvals",
        "pvals_adj",
        "log_q",
    ]
    assert set(result["group"]) == {"A", "B"}
    assert "reference" not in result["group"].tolist()
    assert result["log_q"].is_monotonic_decreasing
    assert result.index.tolist() == list(range(len(result)))


def test_get_degs_preserves_positive_and_negative_fold_change_filters(monkeypatch):
    events = []
    _install_rank_result(monkeypatch, events)

    positive = differential_expression.get_degs(
        _adata(), "cluster", fc_threshold=1, reference="reference"
    )
    negative = differential_expression.get_degs(
        _adata(), "cluster", fc_threshold=-1, reference="reference"
    )

    assert positive["logfoldchanges"].tolist() == [2.0, 1.5]
    assert negative["logfoldchanges"].tolist() == [-2.0, -1.5]


def test_bio_get_degs_forwards_to_the_basic_authority(monkeypatch):
    expected = pd.DataFrame({"group": ["A"], "gene": ["G1"]})
    calls = []

    def authority(*args, **kwargs):
        calls.append((args, kwargs))
        return expected

    monkeypatch.setattr(differential_expression, "get_degs", authority)
    adata = object()

    result = bio.get_degs(
        adata,
        "cluster",
        method="wilcoxon",
        fc_threshold=1,
        reference="reference",
    )

    assert result is expected
    assert calls == [
        (
            (adata,),
            {
                "groupby": "cluster",
                "method": "wilcoxon",
                "fc_threshold": 1,
                "reference": "reference",
            },
        )
    ]
