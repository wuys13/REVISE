from __future__ import annotations

from types import SimpleNamespace

import pandas as pd
import pytest

from revise.analysis import bio
from revise.analysis.advanced import enrichment


def test_empty_gene_list_returns_canonical_result_without_loading_gseapy(
    monkeypatch,
):
    monkeypatch.setattr(
        enrichment,
        "_require_gseapy",
        lambda: (_ for _ in ()).throw(
            AssertionError("empty input attempted to load gseapy")
        ),
    )

    result = enrichment.get_enrichment([], "genes.gmt")

    assert result.empty
    assert result.columns.tolist() == list(enrichment.ENRICHMENT_RESULT_COLUMNS)


def test_missing_gseapy_names_the_pathway_extra(monkeypatch):
    def missing():
        raise ImportError(
            'gseapy is required; install it with "revise-svc[pathway]"'
        )

    monkeypatch.setattr(enrichment, "_require_gseapy", missing)

    with pytest.raises(ImportError, match=r"revise-svc\[pathway\]"):
        enrichment.get_enrichment(["G1"], "genes.gmt")


def test_strict_authority_returns_provider_results(monkeypatch):
    expected = pd.DataFrame(
        {"Term": ["P1"], "Adjusted P-value": [0.01], "Genes": ["G1"]}
    )
    calls = []

    def enrichr(**kwargs):
        calls.append(kwargs)
        return SimpleNamespace(results=expected)

    monkeypatch.setattr(
        enrichment, "_require_gseapy", lambda: SimpleNamespace(enrichr=enrichr)
    )

    result = enrichment.get_enrichment(["G1"], "genes.gmt", cutoff=0.1)

    assert result is expected
    assert calls == [
        {
            "gene_list": ["G1"],
            "gene_sets": "genes.gmt",
            "organism": "human",
            "cutoff": 0.1,
        }
    ]


def test_provider_empty_result_is_returned_with_canonical_schema(monkeypatch):
    monkeypatch.setattr(
        enrichment,
        "_require_gseapy",
        lambda: SimpleNamespace(
            enrichr=lambda **kwargs: SimpleNamespace(results=pd.DataFrame())
        ),
    )

    result = enrichment.get_enrichment(["G1"], "genes.gmt")

    assert result.empty
    assert result.columns.tolist() == list(enrichment.ENRICHMENT_RESULT_COLUMNS)


def test_strict_authority_distinguishes_third_party_execution_failure(monkeypatch):
    def fail(**kwargs):
        raise RuntimeError("provider unavailable")

    monkeypatch.setattr(
        enrichment, "_require_gseapy", lambda: SimpleNamespace(enrichr=fail)
    )

    with pytest.raises(enrichment.EnrichmentExecutionError) as exc_info:
        enrichment.get_enrichment(["G1"], "genes.gmt")

    assert isinstance(exc_info.value.__cause__, RuntimeError)
    assert str(exc_info.value.__cause__) == "provider unavailable"


def test_local_enrichment_uses_gseapy_enrich_with_explicit_background(monkeypatch):
    expected = pd.DataFrame(
        {"Term": ["P1"], "Adjusted P-value": [0.01], "Genes": ["G1"]}
    )
    calls = []

    def enrich(**kwargs):
        calls.append(kwargs)
        return SimpleNamespace(results=expected)

    monkeypatch.setattr(
        enrichment, "_require_gseapy", lambda: SimpleNamespace(enrich=enrich)
    )

    result = enrichment.get_enrichment_local(
        ["G1"], "genes.gmt", cutoff=0.1, background=None
    )

    assert result is expected
    assert calls == [
        {
            "gene_list": ["G1"],
            "gene_sets": "genes.gmt",
            "background": None,
            "cutoff": 0.1,
            "no_plot": True,
        }
    ]


def test_local_enrichment_empty_input_and_provider_empty_use_canonical_schema(monkeypatch):
    monkeypatch.setattr(
        enrichment,
        "_require_gseapy",
        lambda: SimpleNamespace(enrich=lambda **kwargs: SimpleNamespace(results=pd.DataFrame())),
    )

    for genes in ([], ["G1"]):
        result = enrichment.get_enrichment_local(genes, "genes.gmt")
        assert result.empty
        assert result.columns.tolist() == list(enrichment.ENRICHMENT_RESULT_COLUMNS)


def test_local_enrichment_surfaces_missing_dependency(monkeypatch):
    def missing():
        raise ImportError("install revise-svc[pathway]")

    monkeypatch.setattr(enrichment, "_require_gseapy", missing)

    with pytest.raises(ImportError, match=r"revise-svc\[pathway\]"):
        enrichment.get_enrichment_local(["G1"], "genes.gmt")


def test_local_enrichment_preserves_provider_failure(monkeypatch):
    def fail(**kwargs):
        raise RuntimeError("local provider failed")

    monkeypatch.setattr(
        enrichment, "_require_gseapy", lambda: SimpleNamespace(enrich=fail)
    )

    with pytest.raises(
        enrichment.EnrichmentExecutionError, match="local enrichment execution failed"
    ) as exc_info:
        enrichment.get_enrichment_local(["G1"], "genes.gmt")

    assert isinstance(exc_info.value.__cause__, RuntimeError)
    assert str(exc_info.value.__cause__) == "local provider failed"


def test_bio_get_enrichment_converts_execution_failure_to_legacy_empty_result(
    monkeypatch,
):
    def fail(*args, **kwargs):
        try:
            raise RuntimeError("provider unavailable")
        except RuntimeError as exc:
            raise enrichment.EnrichmentExecutionError(
                "Enrichr execution failed"
            ) from exc

    monkeypatch.setattr(enrichment, "get_enrichment", fail)
    messages = []
    monkeypatch.setattr(
        "builtins.print",
        lambda *args, **kwargs: messages.append(" ".join(map(str, args))),
    )

    result = bio.get_enrichment(["G1"], "genes.gmt")

    assert result.empty
    assert result.columns.tolist() == list(enrichment.ENRICHMENT_RESULT_COLUMNS)
    assert messages == [
        "Skipping enrichment analysis: RuntimeError: provider unavailable"
    ]


def test_bio_get_enrichment_does_not_hide_missing_dependency(monkeypatch):
    def missing(*args, **kwargs):
        raise ImportError("install revise-svc[pathway]")

    monkeypatch.setattr(enrichment, "get_enrichment", missing)

    with pytest.raises(ImportError, match=r"revise-svc\[pathway\]"):
        bio.get_enrichment(["G1"], "genes.gmt")
