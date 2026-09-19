import numpy as np
import pandas as pd
from anndata import AnnData, read_h5ad

from revise.backend.kernels.ot import _publish_annotation
from revise.backend.ops.assignment import GlobalAssignment
from revise.reference_preparation.scoring import row_metrics
from revise.utils.confidence import entropy_certainty
from revise.utils.confidence import max_confidence
from revise.utils.confidence import top2_margin


def test_shared_metrics_match_reference_row_metrics_exactly():
    matrix = np.array(
        [[1.0, 0.0, 0.0], [0.7, 0.2, 0.1], [0.5, 0.5, 0.0]],
        dtype=float,
    )

    with np.errstate(divide="raise", invalid="raise"):
        metrics = row_metrics(matrix)
        np.testing.assert_array_equal(max_confidence(matrix), metrics["max_confidence"])
        np.testing.assert_array_equal(
            entropy_certainty(matrix), metrics["normalized_certainty"]
        )
        np.testing.assert_array_equal(top2_margin(matrix), metrics["top2_margin"])


def test_annotation_publication_keeps_legacy_pandas_max_and_records_stage_metadata(
    tmp_path,
):
    target = AnnData(
        X=np.ones((2, 1)),
        obs=pd.DataFrame(index=["spot1", "spot2"]),
        var=pd.DataFrame(index=["g1"]),
    )
    posterior = pd.DataFrame(
        [[0.2, np.nan], [np.nan, np.nan]],
        index=target.obs_names,
        columns=["A", "B"],
    )
    assignment = GlobalAssignment(
        labels=pd.Series(["A", "Unknown"], index=target.obs_names),
        posterior=posterior,
    )

    result = _publish_annotation(
        target,
        assignment,
        annotation_key="Level1",
        confidence_key="Confidence",
        method="tacco",
        reference_origin="reference.h5ad",
    )

    expected = posterior.max(axis=1).to_numpy(copy=True)
    np.testing.assert_allclose(result.obs["Confidence"].to_numpy(), expected, equal_nan=True)
    metadata = result.uns["revise_confidence"]
    assert metadata["stage"] == "annotation:Level1"
    assert metadata["method"] == "tacco"
    assert metadata["annotation_key"] == "Level1"
    assert metadata["confidence_key"] == "Confidence"
    assert metadata["candidate_categories"] == ["A", "B"]
    assert metadata["candidate_categories_meaning"] == (
        "posterior columns are candidate annotation categories"
    )
    assert metadata["reference_origin"] == "reference.h5ad"

    output_path = tmp_path / "annotated.h5ad"
    result.write_h5ad(output_path)
    loaded = read_h5ad(output_path)
    loaded_metadata = loaded.uns["revise_confidence"]
    for key, value in metadata.items():
        if key == "candidate_categories":
            assert loaded_metadata[key].tolist() == value
        else:
            assert loaded_metadata[key] == value


def test_assembly_preserves_published_confidence_column():
    spatial = AnnData(
        X=np.zeros((2, 1)),
        obs=pd.DataFrame(
            {"Confidence": [0.8, 0.6], "SVC_cluster": ["a", "a"]},
            index=["spot1", "spot2"],
        ),
        var=pd.DataFrame(index=["g1"]),
    )
    expression = AnnData(
        X=np.array([[3.0], [5.0]]),
        obs=pd.DataFrame({"SVC_cluster": ["a", "a"]}, index=["cell1", "cell2"]),
        var=pd.DataFrame(index=["g1"]),
    )

    from revise.application.ist_assembly import assemble_ist

    result = assemble_ist(spatial, expression, mapping="mean", seed=42)
    np.testing.assert_allclose(result.obs["Confidence"].to_numpy(), [0.8, 0.6])
