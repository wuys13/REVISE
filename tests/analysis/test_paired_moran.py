"""Paired descriptive Moran preserves unavailable genes and shared weights."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse


def test_full_gene_moran_preserves_missing_constant_and_valid_results():
    from revise.analysis.paired_moran import compare_moran

    obs = pd.DataFrame(index=["a", "b", "c"])
    raw = AnnData(sparse.csr_matrix([[0, 0], [1, 0], [2, 0]]), obs=obs,
                  var=pd.DataFrame(index=["shared", "zero"]))
    recon = AnnData(sparse.csr_matrix([[0, 0, 0], [0, 0, 1], [3, 0, 0]]), obs=obs,
                    var=pd.DataFrame(index=["shared", "zero", "new"]))
    weights = sparse.csr_matrix([[0, 1, 0], [0.5, 0, 0.5], [0, 1, 0]])

    result = compare_moran(raw, recon, weights, min_units=3, min_edges=2).set_index("gene_id")

    assert result.index.tolist() == ["shared", "zero", "new"]
    assert result.loc["shared", "raw_moran"] == pytest.approx(0)
    assert result.loc["shared", "reconstruction_moran"] == pytest.approx(-0.25)
    assert result.loc["shared", "delta"] == pytest.approx(-0.25)
    assert result.loc["zero", "raw_status"] == "not_computable"
    assert result.loc["zero", "raw_reason"] == "constant_expression"
    assert result.loc["new", "raw_status"] == "unmeasured"
    assert result.loc["new", "reconstruction_status"] == "computed"
    assert pd.isna(result.loc["new", "raw_moran"])
    assert pd.isna(result.loc["new", "delta"])
    assert result.loc["new", "comparison_status"] == "unavailable"


def test_moran_without_edges_keeps_missing_distinct_from_low_support():
    from revise.analysis.paired_moran import compare_moran

    obs = pd.DataFrame(index=["a", "b"])
    raw = AnnData(np.array([[0.0], [1.0]]), obs=obs, var=pd.DataFrame(index=["raw"]))
    recon = AnnData(np.array([[1.0], [0.0]]), obs=obs, var=pd.DataFrame(index=["new"]))
    result = compare_moran(raw, recon, sparse.csr_matrix((2, 2)), min_units=2, min_edges=1)

    assert result.loc[0, "raw_status"] == "insufficient_support"
    assert result.loc[0, "raw_reason"] == "no_edges"
    assert result.loc[0, "reconstruction_status"] == "unmeasured"
    assert result["delta"].isna().all()


def test_moran_refuses_row_order_mismatch_instead_of_using_wrong_graph():
    from revise.analysis.paired_moran import compare_moran

    raw = AnnData(np.array([[0.0], [1.0]]), obs=pd.DataFrame(index=["a", "b"]))
    recon = raw[::-1].copy()
    with pytest.raises(ValueError, match="observation"):
        compare_moran(raw, recon, sparse.csr_matrix([[0, 1], [1, 0]]), min_units=2, min_edges=1)


def test_moran_matches_scanpy_with_the_same_final_weights():
    import scanpy as sc
    from revise.analysis.paired_moran import compare_moran

    values = sparse.csr_matrix([[0.0, 4], [2, 1], [4, 2], [5, 0]])
    adata = AnnData(values, obs=pd.DataFrame(index=list("abcd")),
                    var=pd.DataFrame(index=["g1", "g2"]))
    weights = sparse.csr_matrix([[0, 1, 0, 0], [0.5, 0, 0.5, 0],
                                 [0, 0.25, 0, 0.75], [0, 0, 1, 0]])
    expected = sc.metrics.morans_i(weights, values.T.tocsr())
    result = compare_moran(adata, adata.copy(), weights, min_units=2, min_edges=1)
    np.testing.assert_allclose(result["raw_moran"], expected, atol=1e-12)
    np.testing.assert_array_equal(result["delta"], 0)
