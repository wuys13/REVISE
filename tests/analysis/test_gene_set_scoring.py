from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise.analysis.basic.gene_set_scoring import read_gmt, score_genes


def _expression_fixture() -> AnnData:
    rng = np.random.default_rng(20260812)
    x = rng.normal(size=(30, 100))
    return AnnData(
        x,
        obs=pd.DataFrame(index=[f"cell-{i}" for i in range(30)]),
        var=pd.DataFrame(index=[f"G{i}" for i in range(100)]),
    )


def test_read_gmt_preserves_set_and_gene_order(tmp_path):
    path = tmp_path / "sets.gmt"
    path.write_text(
        "SET_B\tdescription\tG3\tG1\nSET_A\tna\tG2\tG4\tG6\n",
        encoding="utf-8",
    )

    result = read_gmt(path)

    assert list(result) == ["SET_B", "SET_A"]
    assert result == {"SET_B": ["G3", "G1"], "SET_A": ["G2", "G4", "G6"]}


def test_read_gmt_distinguishes_missing_malformed_and_duplicate_resources(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_gmt(tmp_path / "missing.gmt")

    malformed = tmp_path / "malformed.gmt"
    malformed.write_text("SET_ONLY\tdescription\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Malformed GMT row 1"):
        read_gmt(malformed)

    duplicate = tmp_path / "duplicate.gmt"
    duplicate.write_text(
        "SET\tdescription\tG1\nSET\tdescription\tG2\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="Duplicate gene set 'SET'"):
        read_gmt(duplicate)


def test_score_genes_scores_only_overlap_on_a_copy():
    adata = _expression_fixture()
    original_x = adata.X.copy()

    scored, score_key = score_genes(
        adata,
        ["G1", "MISSING", "G3"],
        score_name="HYPOXIA",
    )

    assert score_key == "HYPOXIA"
    assert score_key in scored.obs
    assert score_key not in adata.obs
    np.testing.assert_array_equal(adata.X, original_x)
    np.testing.assert_array_equal(scored.X, original_x)


def test_score_genes_rejects_zero_overlap():
    with pytest.raises(ValueError, match="overlap"):
        score_genes(_expression_fixture(), ["MISSING"], score_name="EMPTY")
