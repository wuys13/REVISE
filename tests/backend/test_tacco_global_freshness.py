from __future__ import annotations

import importlib.metadata
import logging
import sys
import types
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData


pytestmark = pytest.mark.filterwarnings("error")


def _inputs():
    target = AnnData(
        X=np.array([[2.0, 1.0], [1.0, 2.0]]),
        obs=pd.DataFrame(index=["spot1", "spot2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    reference = AnnData(
        X=np.array([[2.0, 0.0], [0.0, 2.0]]),
        obs=pd.DataFrame({"Level1": ["A", "B"]}, index=["cell1", "cell2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    return target, reference


def _config(mode):
    return SimpleNamespace(
        annotate_mode=mode,
        cell_type_col="Level1",
        confidence_col="Confidence",
        unknown_key="Unknown",
    )


def _install_tacco(monkeypatch, annotate, version="0.5.0"):
    module = types.ModuleType("tacco")
    module.__version__ = version
    module.tl = SimpleNamespace(annotate=annotate)
    monkeypatch.setitem(sys.modules, "tacco", module)
    monkeypatch.setattr(importlib.metadata, "version", lambda package: version)
    return module


def _kernel_class(monkeypatch):
    ot_module = types.ModuleType("ot")
    ot_module.unbalanced = SimpleNamespace(
        sinkhorn_unbalanced=lambda *args, **kwargs: np.eye(2) * 0.5
    )
    monkeypatch.setitem(sys.modules, "ot", ot_module)
    distance = types.ModuleType("revise.backend.ops.distance")
    distance.bhattacharyya_distance = lambda profiles, expression: np.array(
        [[0.0, 1.0], [1.0, 0.0]]
    )
    monkeypatch.setitem(sys.modules, "revise.backend.ops.distance", distance)
    from revise.backend.kernels import global_anchoring, ot as ot_kernel

    monkeypatch.setattr(
        ot_kernel,
        "bhattacharyya_distance",
        distance.bhattacharyya_distance,
    )

    return global_anchoring.GlobalAnchoringKernel


def _valid_result(adata):
    return pd.DataFrame(
        [[0.8, 0.2], [0.1, 0.9]],
        index=adata.obs_names,
        columns=["A", "B"],
    )


def test_application_sc_pot_does_not_require_tacco_annotation_parameters(
    monkeypatch,
):
    from revise.backend.kernels import global_anchoring
    from revise.config.runner_conf import ApplicationScConf

    target, reference = _inputs()
    config = ApplicationScConf(
        sample_name="sample",
        raw_data_path="data",
        result_root_path="output",
        st_file="sp.h5ad",
        sc_ref_file="sc.h5ad",
        annotate_mode="pot",
        rec_ot_method="pot",
        cell_type_col="Level1",
        confidence_col="Confidence",
        unknown_key="Unknown",
        tacco_annotate_multi_center=None,
        tacco_annotate_lamb=None,
        annotate_pot_reg=0.1,
        annotate_pot_reg_m=0.0,
        annotate_pot_reg_type="entropy",
        rec_graph_n_neighbors=10,
        rec_graph_exp_neighbor_num=15,
        rec_graph_spatial_neighbor_num=6,
        rec_graph_method="joint",
        rec_graph_alpha=0.2,
        rec_random_state=0,
        rec_pot_reg=0.1,
        rec_pot_reg_m=0.0,
        rec_pot_reg_type="entropy",
        rec_alpha=0.5,
        rec_match_spot_sum=False,
    )
    calls = []
    monkeypatch.setattr(
        global_anchoring,
        "OTKernel",
        SimpleNamespace(
            annotate=lambda target, reference, **kwargs: calls.append(kwargs)
            or target.copy()
        ),
    )

    global_anchoring.GlobalAnchoringKernel(
        config,
        logging.getLogger("test"),
    ).run(
        target,
        reference,
        annotate_pot_reg=0.1,
        annotate_pot_reg_m=0.0,
        annotate_pot_reg_type="entropy",
    )

    assert calls[0]["method"] == "pot"


def test_tacco_rejects_zero_mass_posterior_row(monkeypatch):
    GlobalAnchoringKernel = _kernel_class(monkeypatch)
    target, reference = _inputs()

    def annotate(adata, ref, annotation_key, *, result_key, return_reference):
        adata.obsm[result_key] = pd.DataFrame(
            [[0.0, 0.0], [3.0, 7.0]],
            index=adata.obs_names,
            columns=["A", "B"],
        )
        return adata, ref

    _install_tacco(monkeypatch, annotate)

    with pytest.raises(ValueError, match="positive mass"):
        GlobalAnchoringKernel(
            _config("tacco"),
            logging.getLogger("test"),
        ).run(target, reference)


def test_stale_canonical_key_cannot_satisfy_tacco_noop(monkeypatch):
    GlobalAnchoringKernel = _kernel_class(monkeypatch)
    target, reference = _inputs()
    target.obsm["Level1"] = _valid_result(target)

    def annotate(adata, ref, annotation_key, *, result_key, return_reference):
        assert return_reference is True
        return adata, ref

    _install_tacco(monkeypatch, annotate)
    kernel = GlobalAnchoringKernel(_config("tacco"), logging.getLogger("test"))

    with pytest.raises(KeyError, match="fresh"):
        kernel.run(target, reference)


def test_tacco_uses_unique_result_key_then_promotes_only_fresh_output(monkeypatch):
    GlobalAnchoringKernel = _kernel_class(monkeypatch)
    target, reference = _inputs()
    target.obsm["Level1"] = pd.DataFrame(
        [[1.0, 0.0], [1.0, 0.0]],
        index=target.obs_names,
        columns=["A", "B"],
    )
    seen = []

    def annotate(adata, ref, annotation_key, *, result_key, return_reference):
        seen.append(result_key)
        assert result_key != "Level1"
        assert return_reference is True
        adata.obsm[result_key] = _valid_result(adata)
        return adata, ref

    _install_tacco(monkeypatch, annotate)
    kernel = GlobalAnchoringKernel(
        _config("tacco"),
        logging.getLogger("test"),
    )

    result = kernel.run(target, reference)

    assert len(seen) == 1
    assert seen[0] not in result.obsm
    pd.testing.assert_frame_equal(result.obsm["Level1"], _valid_result(result))
    assert result.obs["Level1"].tolist() == ["A", "B"]


def test_tacco_rejects_all_nan_rows_without_reference_prior_repair(
    monkeypatch,
):
    GlobalAnchoringKernel = _kernel_class(monkeypatch)
    target = AnnData(
        X=np.array(
            [
                [2.0, 1.0, 0.0, 0.0],
                [1.0, 2.0, 0.0, 0.0],
                [0.0, 0.0, 3.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
                [0.0, 0.0, 0.0, 0.0],
            ]
        ),
        obs=pd.DataFrame(
            index=[
                "spot1",
                "spot2",
                "spot_filtered",
                "spot_unshared",
                "spot_zero",
            ]
        ),
        var=pd.DataFrame(
            index=["g1", "g2", "filtered_shared", "target_only"]
        ),
    )
    reference = AnnData(
        X=np.array(
            [
                [1.0, 0.0, 1.0],
                [0.0, 9.0, 1.0],
                [1.0, 0.0, 1.0],
            ]
        ),
        obs=pd.DataFrame(
            {"Level1": ["A", "B", "A"]},
            index=["cell1", "cell2", "cell3"],
        ),
        var=pd.DataFrame(index=["g1", "g2", "filtered_shared"]),
    )
    processed_reference = reference[:, ["g1", "g2"]].copy()

    def annotate(adata, ref, annotation_key, *, result_key, return_reference):
        assert return_reference is True
        assert adata.obs_names.tolist() == target.obs_names.tolist()
        ref.uns["tacco_side_effect"] = True
        adata.obsm[result_key] = pd.DataFrame(
            [
                [0.8, 0.2],
                [0.1, 0.9],
                [np.nan, np.nan],
                [np.nan, np.nan],
                [np.nan, np.nan],
            ],
            index=adata.obs_names,
            columns=["A", "B"],
        )
        return adata, processed_reference

    _install_tacco(monkeypatch, annotate)
    with pytest.raises(ValueError, match="finite"):
        GlobalAnchoringKernel(
            _config("tacco"),
            logging.getLogger("test"),
        ).run(target, reference)

    assert "tacco_side_effect" not in reference.uns


def test_tacco_does_not_hide_partial_nonfinite_final_zero_row(monkeypatch):
    GlobalAnchoringKernel = _kernel_class(monkeypatch)
    target, reference = _inputs()
    target = AnnData(
        X=np.vstack([target.X, np.zeros(target.n_vars)]),
        obs=pd.DataFrame(index=["spot1", "spot2", "spot_zero"]),
        var=target.var.copy(),
    )

    def annotate(adata, ref, annotation_key, *, result_key, return_reference):
        assert return_reference is True
        adata.obsm[result_key] = pd.DataFrame(
            [[0.8, 0.2], [0.1, 0.9], [np.nan, 0.0]],
            index=adata.obs_names,
            columns=["A", "B"],
        )
        return adata, ref

    _install_tacco(monkeypatch, annotate)

    with pytest.raises(ValueError, match="finite"):
        GlobalAnchoringKernel(
            _config("tacco"),
            logging.getLogger("test"),
        ).run(target, reference)


def test_tacco_rejects_target_with_no_informative_shared_gene_row(monkeypatch):
    GlobalAnchoringKernel = _kernel_class(monkeypatch)
    _, reference = _inputs()
    target = AnnData(
        X=np.array([[1.0], [2.0]]),
        obs=pd.DataFrame(index=["spot1", "spot2"]),
        var=pd.DataFrame(index=["target_only"]),
    )

    def annotate(*args, **kwargs):
        raise AssertionError("TACCO must not run without an informative row")

    _install_tacco(monkeypatch, annotate)
    with pytest.raises(ValueError, match="at least one target row"):
        GlobalAnchoringKernel(
            _config("tacco"),
            logging.getLogger("test"),
        ).run(target, reference)


def test_tacco_wrong_result_key_fails_closed(monkeypatch):
    GlobalAnchoringKernel = _kernel_class(monkeypatch)
    target, reference = _inputs()

    def annotate(adata, ref, annotation_key, *, result_key, return_reference):
        assert return_reference is True
        adata.obsm["some_other_key"] = _valid_result(adata)
        return adata, ref

    _install_tacco(monkeypatch, annotate)

    with pytest.raises(KeyError, match="fresh"):
        GlobalAnchoringKernel(_config("tacco"), logging.getLogger("test")).run(
            target, reference
        )


@pytest.mark.parametrize(
    "values, message",
    [
        ([[np.nan, 0.2], [0.1, 0.9]], "finite"),
        ([[-0.1, 1.1], [0.1, 0.9]], "non-negative"),
    ],
)
def test_malformed_fresh_tacco_result_fails_before_completed(
    monkeypatch, values, message
):
    GlobalAnchoringKernel = _kernel_class(monkeypatch)
    target, reference = _inputs()

    def annotate(adata, ref, annotation_key, *, result_key, return_reference):
        assert return_reference is True
        adata.obsm[result_key] = pd.DataFrame(
            values, index=adata.obs_names, columns=["A", "B"]
        )
        return adata, ref

    _install_tacco(monkeypatch, annotate)
    with pytest.raises(ValueError, match=message):
        GlobalAnchoringKernel(
            _config("tacco"),
            logging.getLogger("test"),
        ).run(target, reference)


@pytest.mark.parametrize(
    "coupling, message",
    [
        (np.ones((1, 2)), "shape"),
        (np.array([[np.nan, 0.0], [0.0, 1.0]]), "finite"),
        (np.array([[-0.1, 0.6], [0.0, 0.5]]), "non-negative"),
        (np.array([[0.0, 0.0], [0.0, 0.5]]), "positive mass"),
    ],
)
def test_malformed_pot_coupling_fails_before_completed(
    monkeypatch, coupling, message
):
    GlobalAnchoringKernel = _kernel_class(monkeypatch)
    target, reference = _inputs()
    monkeypatch.setattr(
        sys.modules["ot"].unbalanced,
        "sinkhorn_unbalanced",
        lambda *args, **kwargs: coupling,
    )
    kernel = GlobalAnchoringKernel(
        _config("pot"),
        logging.getLogger("test"),
    )

    with pytest.raises(ValueError, match=message):
        kernel.run(
            target,
            reference,
            annotate_pot_reg=0.1,
            annotate_pot_reg_m=0.0,
            annotate_pot_reg_type="entropy",
        )


def test_pot_global_records_completed_only_after_validated_result(monkeypatch):
    GlobalAnchoringKernel = _kernel_class(monkeypatch)
    target, reference = _inputs()
    monkeypatch.setattr(
        sys.modules["ot"].unbalanced,
        "sinkhorn_unbalanced",
        lambda *args, **kwargs: np.array([[0.4, 0.1], [0.2, 0.3]]),
    )
    result = GlobalAnchoringKernel(
        _config("pot"),
        logging.getLogger("test"),
    ).run(
        target,
        reference,
        annotate_pot_reg=0.1,
        annotate_pot_reg_m=0.0,
        annotate_pot_reg_type="entropy",
    )

    assert result.obsm["Level1"].shape == (2, 2)


def test_pot_publishes_ordered_row_normalized_posterior_and_argmax_labels(
    monkeypatch,
):
    GlobalAnchoringKernel = _kernel_class(monkeypatch)
    target, reference = _inputs()
    reference.obs["Level1"] = ["B", "A"]
    monkeypatch.setattr(
        sys.modules["ot"].unbalanced,
        "sinkhorn_unbalanced",
        lambda *args, **kwargs: np.array([[0.4, 0.1], [0.2, 0.3]]),
    )

    result = GlobalAnchoringKernel(
        _config("pot"),
        logging.getLogger("test"),
    ).run(
        target,
        reference,
        annotate_pot_reg=0.1,
        annotate_pot_reg_m=0.0,
        annotate_pot_reg_type="entropy",
    )

    expected = pd.DataFrame(
        [[0.8, 0.2], [0.4, 0.6]],
        index=target.obs_names,
        columns=["B", "A"],
    )
    pd.testing.assert_frame_equal(result.obsm["Level1"], expected)
    assert result.obs["Level1"].tolist() == expected.idxmax(axis=1).tolist()


def test_tacco_publishes_ordered_values_without_reordering_or_normalizing(
    monkeypatch,
):
    GlobalAnchoringKernel = _kernel_class(monkeypatch)
    target, reference = _inputs()
    reference.obs["Level1"] = ["B", "A"]
    expected = pd.DataFrame(
        [[0.7, 0.3], [0.4000005, 0.6]],
        index=target.obs_names,
        columns=["B", "A"],
    )

    def annotate(adata, ref, annotation_key, *, result_key, return_reference):
        adata.obsm[result_key] = expected.copy()
        return adata, ref

    _install_tacco(monkeypatch, annotate)

    result = GlobalAnchoringKernel(
        _config("tacco"),
        logging.getLogger("test"),
    ).run(target, reference)

    pd.testing.assert_frame_equal(result.obsm["Level1"], expected)
    assert result.obs["Level1"].tolist() == expected.idxmax(axis=1).tolist()


def test_tacco_accepts_returned_category_order_when_reference_set_matches(
    monkeypatch,
):
    GlobalAnchoringKernel = _kernel_class(monkeypatch)
    target, reference = _inputs()
    reference.obs["Level1"] = ["B", "A"]
    expected = pd.DataFrame(
        [[0.3, 0.7], [0.6, 0.4]],
        index=target.obs_names,
        columns=["A", "B"],
    )

    def annotate(adata, ref, annotation_key, *, result_key, return_reference):
        adata.obsm[result_key] = expected.copy()
        return adata, ref

    _install_tacco(monkeypatch, annotate)
    result = GlobalAnchoringKernel(
        _config("tacco"),
        logging.getLogger("test"),
    ).run(target, reference)

    pd.testing.assert_frame_equal(result.obsm["Level1"], expected)
    assert result.obs["Level1"].tolist() == expected.idxmax(axis=1).tolist()


@pytest.mark.parametrize(
    ("index", "columns", "message"),
    [
        (["spot2", "spot1"], ["A", "B"], "observation.*order"),
        (["spot1", "spot2"], ["A", "C"], "category.*mismatch"),
    ],
)
def test_tacco_rejects_permuted_axes_instead_of_reordering(
    monkeypatch,
    index,
    columns,
    message,
):
    GlobalAnchoringKernel = _kernel_class(monkeypatch)
    target, reference = _inputs()
    candidate = pd.DataFrame(
        [[0.8, 0.2], [0.1, 0.9]],
        index=index,
        columns=columns,
    )

    def annotate(adata, ref, annotation_key, *, result_key, return_reference):
        annotated = adata[candidate.index].copy()
        annotated.obsm[result_key] = candidate.copy()
        return annotated, ref

    _install_tacco(monkeypatch, annotate)
    with pytest.raises(ValueError, match=message):
        GlobalAnchoringKernel(
            _config("tacco"),
            logging.getLogger("test"),
        ).run(target, reference)


def test_tacco_publishes_owned_posterior(monkeypatch):
    GlobalAnchoringKernel = _kernel_class(monkeypatch)
    target, reference = _inputs()
    candidate = _valid_result(target)

    def annotate(adata, ref, annotation_key, *, result_key, return_reference):
        adata.obsm[result_key] = candidate
        return adata, ref

    _install_tacco(monkeypatch, annotate)
    result = GlobalAnchoringKernel(
        _config("tacco"),
        logging.getLogger("test"),
    ).run(target, reference)

    assert result.obsm["Level1"] is not candidate
    pd.testing.assert_frame_equal(result.obsm["Level1"], candidate)


def test_tacco_invalid_values_use_global_assignment_contract_error(monkeypatch):
    from revise.backend.ops.assignment import GlobalAssignmentContractError

    GlobalAnchoringKernel = _kernel_class(monkeypatch)
    target, reference = _inputs()
    candidate = _valid_result(target)
    candidate.iloc[0, 0] = np.nan

    def annotate(adata, ref, annotation_key, *, result_key, return_reference):
        adata.obsm[result_key] = candidate
        return adata, ref

    _install_tacco(monkeypatch, annotate)
    with pytest.raises(GlobalAssignmentContractError, match="finite"):
        GlobalAnchoringKernel(
            _config("tacco"),
            logging.getLogger("test"),
        ).run(target, reference)
