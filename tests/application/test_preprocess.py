from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse


def _adata(*, with_transcript_counts: bool = True) -> AnnData:
    obs = pd.DataFrame(
        {"Patient": ["P1", "P2", "P2"]},
        index=["o1", "o2", "o3"],
    )
    if with_transcript_counts:
        obs["transcript_counts"] = [59, 60, 70]
    return AnnData(
        X=np.array(
            [
                [1.0, 1.0, 0.0],
                [1.0, 0.0, 1.0],
                [1.0, 0.0, 1.0],
            ]
        ),
        obs=obs,
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
    )


def test_filter_reference_returns_an_exact_filtered_copy():
    from revise.application.preprocess import filter_reference

    source = _adata()
    result = filter_reference(source, filter_column="Patient", filter_value="P2")

    assert result.obs_names.tolist() == ["o2", "o3"]
    assert result is not source
    result.obs.loc["o2", "Patient"] = "changed"
    assert source.obs.loc["o2", "Patient"] == "P2"


def test_preprocess_spatial_filters_existing_transcript_counts_then_genes():
    from revise.application.preprocess import preprocess_spatial

    source = _adata()
    result = preprocess_spatial(
        source,
        min_transcript_counts=60,
        min_cell_counts=2,
    )

    assert result.obs_names.tolist() == ["o2", "o3"]
    assert result.var_names.tolist() == ["g1", "g3"]
    assert source.shape == (3, 3)


def test_preprocess_reference_defaults_to_no_transcript_filtering():
    from revise.application.preprocess import preprocess_reference

    source = _adata(with_transcript_counts=False)
    result = preprocess_reference(source, min_cell_counts=2)

    assert result.obs_names.tolist() == ["o1", "o2", "o3"]
    assert result.var_names.tolist() == ["g1", "g3"]


def test_preprocess_spatial_can_filter_cells_by_x_counts_for_sp_svc():
    from revise.application.preprocess import preprocess_spatial

    result = preprocess_spatial(
        _adata(with_transcript_counts=False),
        min_transcript_counts=None,
        min_counts=2,
        min_cell_counts=2,
    )

    assert result.obs_names.tolist() == ["o1", "o2", "o3"]
    assert result.var_names.tolist() == ["g1", "g3"]


def test_preprocess_reference_can_filter_cells_by_detected_genes_for_sp_svc():
    from revise.application.preprocess import preprocess_reference

    source = _adata(with_transcript_counts=False)
    source.X[0, :] = 0.0
    result = preprocess_reference(
        source,
        min_genes=1,
        min_cell_counts=2,
    )

    assert result.obs_names.tolist() == ["o2", "o3"]
    assert result.var_names.tolist() == ["g1", "g3"]


def test_transcript_filter_does_not_infer_missing_counts_from_x():
    from revise.application.preprocess import preprocess_spatial

    with pytest.raises(KeyError, match="transcript_counts"):
        preprocess_spatial(_adata(with_transcript_counts=False))


def test_prepare_sc_svc_pair_normalizes_labels_and_limits_only_spatial_genes():
    from revise.application.preprocess import prepare_sc_svc_pair

    spatial = AnnData(
        X=np.ones((2, 3)),
        obs=pd.DataFrame(index=["spot-1", "spot-2"]),
        var=pd.DataFrame(index=["g1", "g2", "spatial-only"]),
    )
    reference = AnnData(
        X=np.ones((2, 3)),
        obs=pd.DataFrame(
            {"Level1": [" Mono/Macro ", "T"], "Level2": ["M/1", "T/1"]},
            index=["cell-1", "cell-2"],
        ),
        var=pd.DataFrame(index=["g1", "g2", "reference-only"]),
    )

    prepared_spatial, prepared_reference = prepare_sc_svc_pair(
        spatial,
        reference,
        broad_column="Level1",
        subtype_column="Level2",
    )

    assert prepared_spatial.var_names.tolist() == ["g1", "g2"]
    assert prepared_reference.var_names.tolist() == ["g1", "g2", "reference-only"]
    assert prepared_reference.obs.columns.tolist() == ["Level1", "Level2"]
    assert prepared_reference.obs["Level1"].tolist() == ["Mono_Macro", "T"]
    assert prepared_reference.obs["Level2"].tolist() == ["M_1", "T_1"]


def test_sp_sr_reference_label_normalization_preserves_surrounding_whitespace_by_default():
    from revise.application.preprocess import normalize_reference_labels

    reference = AnnData(
        X=np.ones((1, 1)),
        obs=pd.DataFrame({"Level1": [" Mono/Macro "]}, index=["cell-1"]),
        var=pd.DataFrame(index=["g1"]),
    )

    normalized = normalize_reference_labels(reference, ["Level1"])

    assert normalized.obs["Level1"].tolist() == [" Mono_Macro "]


@pytest.mark.parametrize("matrix_format", ["dense", "sparse"])
def test_cyto_linear_v1_expression_score_matches_formula_without_sparse_densification(
    monkeypatch,
    matrix_format,
):
    from revise.utils.spot_sr_input import _cyto_linear_v1_expression_score

    values = np.array([[1.0, 3.0, 0.0], [0.0, 0.0, 0.0]])
    expected = np.log2(1.0 + 1e6 * np.array([0.25, 0.75, 0.0])).sum()
    matrix = values
    if matrix_format == "sparse":
        matrix = sparse.csr_matrix(values)
        monkeypatch.setattr(
            sparse.csr_matrix,
            "toarray",
            lambda *_args, **_kwargs: pytest.fail("sparse path densified X"),
        )
    adata = AnnData(X=matrix)

    score = _cyto_linear_v1_expression_score(adata)

    np.testing.assert_allclose(score, [expected, 0.0])


def test_cyto_linear_v1_mapping_has_minimum_one_cell_and_no_twelve_cell_cap():
    from revise.utils.spot_sr_input import ensure_all_cells_in_spot

    adata = AnnData(
        X=sparse.csr_matrix(
            np.vstack(
                [
                    np.zeros(6000, dtype=np.float64),
                    np.ones(6000, dtype=np.float64),
                ]
            )
        ),
        obs=pd.DataFrame(index=["empty", "large"]),
    )

    ensure_all_cells_in_spot(adata, cell_count_method="cyto_linear_v1")

    mapping = adata.uns["all_cells_in_spot"]
    assert len(mapping["empty"]) == 1
    assert len(mapping["large"]) > 12
    assert adata.obs["estimated_cell_count"].to_dict() == {
        spot: len(cells) for spot, cells in mapping.items()
    }


def test_cyto_linear_v1_preserves_supplied_spot_cell_mapping():
    from revise.utils.spot_sr_input import ensure_all_cells_in_spot

    adata = AnnData(
        X=np.ones((1, 2)),
        obs=pd.DataFrame(index=["spot-1"]),
    )
    supplied = {"spot-1": ["cell-a", "cell-b"]}
    adata.uns["all_cells_in_spot"] = supplied

    ensure_all_cells_in_spot(adata, cell_count_method="cyto_linear_v1")

    assert adata.uns["all_cells_in_spot"] == supplied
    assert adata.obs["estimated_cell_count"].tolist() == [2]
