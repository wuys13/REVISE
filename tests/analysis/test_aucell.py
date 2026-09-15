from __future__ import annotations

import importlib
import sys
from types import ModuleType, SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

from revise.analysis.advanced.aucell import score_gene_set_aucell


def _expression_fixture() -> AnnData:
    return AnnData(
        np.arange(24, dtype=float).reshape(6, 4),
        obs=pd.DataFrame(index=[f"cell-{i}" for i in range(6)]),
        var=pd.DataFrame(index=["G1", "G2", "G3", "G4"]),
    )


def test_aucell_import_is_lazy(monkeypatch):
    import revise.analysis.advanced.aucell as module

    monkeypatch.delitem(sys.modules, "omicverse", raising=False)
    importlib.reload(module)
    assert "omicverse" not in sys.modules


def test_score_aucell_passes_overlap_and_returns_validated_copy(monkeypatch):
    calls = []

    def geneset_aucell(*, adata, geneset_name, geneset):
        calls.append((adata, geneset_name, geneset))
        adata.obs[f"{geneset_name}_aucell"] = np.arange(adata.n_obs, dtype=float)

    fake = ModuleType("omicverse")
    fake.single = SimpleNamespace(geneset_aucell=geneset_aucell)
    monkeypatch.setitem(sys.modules, "omicverse", fake)
    adata = _expression_fixture()
    original_x = adata.X.copy()

    scored, score_key = score_gene_set_aucell(
        adata,
        ["G3", "MISSING", "G1"],
        score_name="SIGNATURE",
    )

    assert score_key == "SIGNATURE_aucell"
    assert calls == [(scored, "SIGNATURE", ["G3", "G1"])]
    assert score_key in scored.obs
    assert score_key not in adata.obs
    assert np.array_equal(adata.X, original_x)
    assert np.array_equal(scored.X.toarray(), original_x)


def test_score_aucell_missing_dependency_names_pathway_extra(monkeypatch):
    import revise.analysis.advanced.aucell as module

    real_import = module.importlib.import_module

    def missing(name):
        if name == "omicverse":
            raise ModuleNotFoundError("missing omicverse", name="omicverse")
        return real_import(name)

    monkeypatch.setattr(module.importlib, "import_module", missing)

    with pytest.raises(ImportError, match=r"revise-svc\[pathway\]"):
        score_gene_set_aucell(_expression_fixture(), ["G1"], score_name="SIGNATURE")


def test_score_aucell_preserves_nested_provider_import_failure(monkeypatch):
    import revise.analysis.advanced.aucell as module

    def missing_torch(name):
        if name == "omicverse":
            raise ModuleNotFoundError("No module named 'torch'", name="torch")
        raise AssertionError(name)

    monkeypatch.setattr(module.importlib, "import_module", missing_torch)

    with pytest.raises(ImportError, match=r"OmicVerse import failed.*torch"):
        score_gene_set_aucell(_expression_fixture(), ["G1"], score_name="SIGNATURE")


def test_score_aucell_distinguishes_failure_missing_result_and_zero_overlap(monkeypatch):
    fake = ModuleType("omicverse")
    fake.single = SimpleNamespace(
        geneset_aucell=lambda **_kwargs: (_ for _ in ()).throw(ValueError("third party failed"))
    )
    monkeypatch.setitem(sys.modules, "omicverse", fake)
    with pytest.raises(ValueError, match="third party failed"):
        score_gene_set_aucell(_expression_fixture(), ["G1"], score_name="SIGNATURE")

    fake.single = SimpleNamespace(geneset_aucell=lambda **_kwargs: None)
    with pytest.raises(RuntimeError, match="did not create 'SIGNATURE_aucell'"):
        score_gene_set_aucell(_expression_fixture(), ["G1"], score_name="SIGNATURE")

    with pytest.raises(ValueError, match="overlap"):
        score_gene_set_aucell(_expression_fixture(), ["MISSING"], score_name="SIGNATURE")


def test_score_aucell_forwards_explicit_provider_parameters(monkeypatch):
    calls = []

    def geneset_aucell(*, adata, geneset_name, geneset, AUC_threshold, seed):
        calls.append({
            "geneset_name": geneset_name,
            "geneset": geneset,
            "AUC_threshold": AUC_threshold,
            "seed": seed,
        })
        adata.obs[f"{geneset_name}_aucell"] = np.ones(adata.n_obs)

    fake = ModuleType("omicverse")
    fake.single = SimpleNamespace(geneset_aucell=geneset_aucell)
    monkeypatch.setitem(sys.modules, "omicverse", fake)

    score_gene_set_aucell(
        _expression_fixture(),
        ["G1"],
        score_name="SIGNATURE",
        AUC_threshold=0.05,
        seed=17,
    )

    assert calls == [{
        "geneset_name": "SIGNATURE",
        "geneset": ["G1"],
        "AUC_threshold": 0.05,
        "seed": 17,
    }]


def test_score_aucell_rejects_nonfinite_provider_scores(monkeypatch):
    fake = ModuleType("omicverse")

    def geneset_aucell(*, adata, geneset_name, geneset):
        adata.obs[f"{geneset_name}_aucell"] = [0.1, np.nan, 0.3, 0.4, 0.5, 0.6]

    fake.single = SimpleNamespace(geneset_aucell=geneset_aucell)
    monkeypatch.setitem(sys.modules, "omicverse", fake)

    with pytest.raises(RuntimeError, match="non-finite"):
        score_gene_set_aucell(_expression_fixture(), ["G1"], score_name="SIGNATURE")


def test_score_aucell_converts_dense_expression_to_the_provider_sparse_carrier(monkeypatch):
    fake = ModuleType("omicverse")
    received = []

    def geneset_aucell(*, adata, geneset_name, geneset):
        received.append(adata.X)
        adata.obs[f"{geneset_name}_aucell"] = np.ones(adata.n_obs)

    fake.single = SimpleNamespace(geneset_aucell=geneset_aucell)
    monkeypatch.setitem(sys.modules, "omicverse", fake)

    score_gene_set_aucell(_expression_fixture(), ["G1"], score_name="SIGNATURE")

    assert len(received) == 1
    assert sparse.isspmatrix_csr(received[0])
