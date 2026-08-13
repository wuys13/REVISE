from __future__ import annotations

import ast
from dataclasses import fields
import logging
from pathlib import Path
import sys
import types
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse
from scipy.stats import pearsonr

from revise.backend.kernels.graph_aggregate import GraphAggregateKernel
from revise.backend.ops.posterior_conditioning import posterior_reference_allocation
from revise.recon.pipeline import UnifiedReconstructionPipeline
from revise.svc import SVC


REPO_ROOT = Path(__file__).resolve().parents[2]


def _global_assignment(labels, posterior):
    from revise.backend.ops import assignment

    return assignment.GlobalAssignment(labels=labels, posterior=posterior)


def _validate_global_assignment(
    labels,
    posterior,
    *,
    expected_observations=("spot-1", "spot-2"),
    expected_categories=("A", "B"),
):
    from revise.backend.ops import assignment

    return assignment.validate_global_assignment(
        _global_assignment(labels, posterior),
        expected_observations=expected_observations,
        expected_categories=expected_categories,
    )


def _ordered_global_assignment():
    posterior = pd.DataFrame(
        [[0.75, 0.25], [0.2, 0.8]],
        index=["spot-1", "spot-2"],
        columns=["A", "B"],
    )
    labels = pd.Series(["A", "B"], index=posterior.index)
    return labels, posterior


def test_global_assignment_carrier_contains_only_labels_and_posterior():
    from revise.backend.ops import assignment

    assert [field.name for field in fields(assignment.GlobalAssignment)] == [
        "labels",
        "posterior",
    ]


def test_global_assignment_is_a_mutable_identity_carrier():
    labels, posterior = _ordered_global_assignment()
    carrier = _global_assignment(labels, posterior)
    equivalent = _global_assignment(labels.copy(), posterior.copy())
    replacement = labels.copy()

    carrier.labels = replacement

    assert carrier.labels is replacement
    assert carrier != equivalent


def test_global_assignment_accepts_ordered_posterior_without_changing_values():
    labels, posterior = _ordered_global_assignment()

    validated = _validate_global_assignment(labels, posterior)

    pd.testing.assert_series_equal(validated.labels, labels)
    pd.testing.assert_frame_equal(validated.posterior, posterior)


@pytest.mark.parametrize(
    ("axis", "actual", "message"),
    [
        ("observation", ["spot-1"], "observation.*missing"),
        ("observation", ["spot-1", "spot-2", "spot-3"], "observation.*extra"),
        ("observation", ["spot-2", "spot-1"], "observation.*order"),
        ("category", ["A"], "category.*missing"),
        ("category", ["A", "B", "C"], "category.*extra"),
        ("category", ["B", "A"], "category.*order"),
    ],
)
def test_global_assignment_rejects_nonexact_axes(axis, actual, message):
    from revise.backend.ops import assignment

    labels, posterior = _ordered_global_assignment()
    if axis == "observation":
        posterior = pd.DataFrame(
            np.tile([[0.75, 0.25]], (len(actual), 1)),
            index=actual,
            columns=posterior.columns,
        )
        labels = pd.Series(["A"] * len(actual), index=actual)
    else:
        posterior = pd.DataFrame(
            np.full((2, len(actual)), 1.0 / len(actual)),
            index=posterior.index,
            columns=actual,
        )
        labels = posterior.idxmax(axis=1)

    with pytest.raises(assignment.GlobalAssignmentContractError, match=message):
        _validate_global_assignment(labels, posterior)


@pytest.mark.parametrize(
    ("axis", "actual", "message"),
    [
        ("observation", ["spot-1", "spot-1"], "observation.*duplicate"),
        ("observation", ["spot/1", "spot_1"], "observation.*collid"),
        ("category", ["A", "A"], "category.*duplicate"),
        ("category", ["A/B", "A_B"], "category.*collid"),
        ("observation", ["spot-1", None], "observation.*null"),
        ("category", ["A", ""], "category.*empty"),
    ],
)
def test_global_assignment_rejects_ambiguous_or_invalid_axes(
    axis,
    actual,
    message,
):
    from revise.backend.ops import assignment

    labels, posterior = _ordered_global_assignment()
    if axis == "observation":
        posterior.index = actual
        labels.index = actual
        expected_observations = actual
        expected_categories = posterior.columns
    else:
        posterior.columns = actual
        labels = pd.Series([actual[0], actual[0]], index=posterior.index)
        expected_observations = posterior.index
        expected_categories = actual

    with pytest.raises(assignment.GlobalAssignmentContractError, match=message):
        _validate_global_assignment(
            labels,
            posterior,
            expected_observations=expected_observations,
            expected_categories=expected_categories,
        )


@pytest.mark.parametrize(
    ("values", "message"),
    [
        ([[1.1, -0.1], [0.2, 0.8]], "non-negative"),
        ([[np.nan, 0.0], [0.2, 0.8]], "finite"),
        ([[np.inf, 0.0], [0.2, 0.8]], "finite"),
        ([[0.6, 0.400002], [0.2, 0.8]], "row-normalized"),
    ],
)
def test_global_assignment_rejects_invalid_posterior_values(values, message):
    from revise.backend.ops import assignment

    posterior = pd.DataFrame(
        values,
        index=["spot-1", "spot-2"],
        columns=["A", "B"],
    )
    labels = posterior.idxmax(axis=1)

    with pytest.raises(assignment.GlobalAssignmentContractError, match=message):
        _validate_global_assignment(labels, posterior)


def test_global_assignment_accepts_in_tolerance_row_mass_without_normalizing():
    posterior = pd.DataFrame(
        [[0.6, 0.4000005], [0.2, 0.8]],
        index=["spot-1", "spot-2"],
        columns=["A", "B"],
    )
    labels = posterior.idxmax(axis=1)

    validated = _validate_global_assignment(labels, posterior)

    pd.testing.assert_frame_equal(validated.posterior, posterior)
    assert validated.posterior.iloc[0].sum() == pytest.approx(1.0000005)


@pytest.mark.parametrize(
    ("values", "labels"),
    [
        ([[0.8, 0.2], [0.1, 0.9]], ["B", "B"]),
        ([[0.5, 0.5], [0.2, 0.8]], ["B", "B"]),
    ],
)
def test_global_assignment_rejects_labels_that_differ_from_pandas_idxmax(
    values,
    labels,
):
    from revise.backend.ops import assignment

    posterior = pd.DataFrame(
        values,
        index=["spot-1", "spot-2"],
        columns=["A", "B"],
    )
    hard_labels = pd.Series(labels, index=posterior.index)

    with pytest.raises(
        assignment.GlobalAssignmentContractError,
        match="labels.*argmax",
    ):
        _validate_global_assignment(hard_labels, posterior)


def test_global_assignment_rejects_hard_labels_without_posterior():
    from revise.backend.ops import assignment

    labels = pd.Series(["A", "B"], index=["spot-1", "spot-2"])

    with pytest.raises(
        assignment.GlobalAssignmentContractError,
        match="posterior.*DataFrame",
    ):
        _validate_global_assignment(labels, None)


def test_global_assignment_validator_does_not_mutate_inputs_or_share_outputs():
    labels, posterior = _ordered_global_assignment()
    original_labels = labels.copy(deep=True)
    original_posterior = posterior.copy(deep=True)

    validated = _validate_global_assignment(labels, posterior)
    validated.labels.iloc[0] = "changed"
    validated.posterior.iloc[0, 0] = 999.0

    pd.testing.assert_series_equal(labels, original_labels)
    pd.testing.assert_frame_equal(posterior, original_posterior)


def _load_functions(relative_path, names, namespace):
    path = REPO_ROOT / relative_path
    tree = ast.parse(path.read_text())
    selected = []
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name in names:
            selected.append(node)
        elif isinstance(node, ast.ClassDef):
            selected.extend(
                child
                for child in node.body
                if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef))
                and child.name in names
            )
    exec(compile(ast.Module(body=selected, type_ignores=[]), path, "exec"), namespace)
    return SimpleNamespace(**{name: namespace[name] for name in names})


def _load_metrics(ssim=None):
    def normalize_total(adata, target_sum):
        values = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)
        row_sums = values.sum(axis=1, keepdims=True)
        adata.X = values / row_sums * float(target_sum)

    if ssim is None:
        def ssim(left, right, data_range):
            return 1.0 - float(
                np.mean(np.abs(np.asarray(left) - np.asarray(right)))
            ) / data_range

    return _load_functions(
        "revise/analysis/metrics.py",
        ["_to_numpy_matrix", "normalize_data", "compute_metric"],
        {
            "issparse": sparse.issparse,
            "np": np,
            "pd": pd,
            "pearsonr": pearsonr,
            "sc": SimpleNamespace(pp=SimpleNamespace(normalize_total=normalize_total)),
            "ssim": ssim,
        },
    )


def test_pipeline_evaluation_passes_aligned_ground_truth_before_prediction(
    monkeypatch,
    tmp_path,
):
    captured = {}

    def compute_metric(gt, pred, *args, **kwargs):
        captured["gt"] = np.asarray(gt.X).copy()
        captured["pred"] = np.asarray(pred.X).copy()
        captured["obs"] = (list(gt.obs_names), list(pred.obs_names))
        return pd.DataFrame({"NRMSE": [0.0]}, index=["gene"])

    metrics = types.ModuleType("revise.analysis.metrics")
    metrics.compute_metric = compute_metric
    monkeypatch.setitem(sys.modules, "revise.analysis.metrics", metrics)

    pred = AnnData(
        X=np.array([[200.0], [100.0]]),
        obs=pd.DataFrame(index=["cell-b", "cell-a"]),
        var=pd.DataFrame(index=["gene"]),
    )
    gt = AnnData(
        X=np.array([[10.0], [20.0]]),
        obs=pd.DataFrame(index=["cell-a", "cell-b"]),
        var=pd.DataFrame(index=["gene"]),
    )
    ctx = SimpleNamespace(
        svc=SimpleNamespace(
            artifacts={"outputs": {"prediction": pred}},
            provenance={},
            quality_metrics={},
        ),
        real_st_adata=gt,
        run_dir=tmp_path,
        logger=logging.getLogger("u15-pipeline-evaluation"),
        quality_metrics={},
        compatibility_mode=False,
        stage_trace=[],
        record_artifact=lambda *_args: None,
    )
    pipeline = UnifiedReconstructionPipeline(
        strategy=object(),
        validation_policy=object(),
        evaluation_policy=SimpleNamespace(should_evaluate=lambda _ctx: True),
    )

    pipeline.evaluate(ctx)

    assert captured["obs"] == (["cell-b", "cell-a"], ["cell-b", "cell-a"])
    np.testing.assert_array_equal(captured["gt"].ravel(), [20.0, 10.0])
    np.testing.assert_array_equal(captured["pred"].ravel(), [200.0, 100.0])


def _metric_adata(values):
    values = np.asarray(values, dtype=np.float64)
    return AnnData(
        X=values,
        obs=pd.DataFrame(index=[f"cell-{i}" for i in range(values.shape[0])]),
        var=pd.DataFrame(index=[f"gene-{i}" for i in range(values.shape[1])]),
    )


def test_nrmse_orientation_is_non_symmetric():
    metrics = _load_metrics()
    gt = _metric_adata(
        [
            [10.0, 90.0],
            [20.0, 80.0],
            [30.0, 70.0],
            [40.0, 60.0],
            [50.0, 50.0],
            [60.0, 40.0],
            [70.0, 30.0],
            [80.0, 20.0],
        ]
    )
    pred = _metric_adata(
        [
            [10.0, 90.0],
            [11.0, 89.0],
            [12.0, 88.0],
            [13.0, 87.0],
            [14.0, 86.0],
            [15.0, 85.0],
            [16.0, 84.0],
            [80.0, 20.0],
        ]
    )

    correct = metrics.compute_metric(gt, pred, logging.getLogger("u15-metric"), normalize=True)
    swapped = metrics.compute_metric(pred, gt, logging.getLogger("u15-metric"), normalize=True)

    assert correct.loc[0, "MSE"] == pytest.approx(0.1880357142857143)
    assert correct.loc[0, "NRMSE"] == pytest.approx(0.8672617004934883)
    assert swapped.loc[0, "NRMSE"] == pytest.approx(2.668497539979964)
    assert correct.loc[0, "NRMSE"] != pytest.approx(swapped.loc[0, "NRMSE"])


def test_ssim_receives_observation_order_while_nrmse_is_order_invariant():
    def ordered_similarity(left, right, data_range):
        weights = np.arange(1, len(left) + 1, dtype=np.float64)
        squared_error = (np.asarray(left) - np.asarray(right)) ** 2
        return 1.0 - float(np.average(squared_error, weights=weights)) / data_range**2

    metrics = _load_metrics(ssim=ordered_similarity)
    gt_values = np.array(
        [[10.0, 90.0], [20.0, 80.0], [30.0, 70.0], [40.0, 60.0],
         [50.0, 50.0], [60.0, 40.0], [70.0, 30.0], [80.0, 20.0]]
    )
    pred_values = np.array(
        [[10.0, 90.0], [11.0, 89.0], [12.0, 88.0], [13.0, 87.0],
         [14.0, 86.0], [15.0, 85.0], [16.0, 84.0], [80.0, 20.0]]
    )
    order = np.array([0, 7, 1, 6, 2, 5, 3, 4])

    original = metrics.compute_metric(
        _metric_adata(gt_values),
        _metric_adata(pred_values),
        logging.getLogger("u15-metric-order"),
        normalize=True,
    )
    permuted = metrics.compute_metric(
        _metric_adata(gt_values[order]),
        _metric_adata(pred_values[order]),
        logging.getLogger("u15-metric-order"),
        normalize=True,
    )

    assert original.loc[0, "NRMSE"] == pytest.approx(permuted.loc[0, "NRMSE"])
    assert original.loc[0, "SSIM"] != pytest.approx(permuted.loc[0, "SSIM"])


def test_metric_undefined_constant_and_zero_mean_boundaries():
    metrics = _load_metrics()
    constant = _metric_adata(np.ones((8, 2)))

    constant_result = metrics.compute_metric(
        constant,
        constant,
        logging.getLogger("u15-metric"),
        normalize=True,
    )

    assert np.isnan(constant_result.loc[0, "PCC"])
    assert np.isnan(constant_result.loc[0, "NRMSE"])

    varying = _metric_adata([[value, 10.0 - value] for value in range(1, 9)])
    zero_mean_result = metrics.compute_metric(
        constant,
        varying,
        logging.getLogger("u15-metric"),
        normalize=True,
    )

    assert np.isposinf(zero_mean_result.loc[0, "NRMSE"])


def test_reference_allocation_uses_category_labels_not_positions():
    allocated = posterior_reference_allocation(
        np.array([[5.0]]),
        pd.DataFrame([[0.0, 1.0]], columns=["B", "A"]),
        pd.DataFrame([[1.0], [10.0]], index=["A", "B"]),
    )

    np.testing.assert_allclose(allocated, [[[5.0, 0.0]]])


@pytest.mark.parametrize(
    ("posterior", "reference", "message"),
    [
        (np.ones((1, 2)), pd.DataFrame([[1.0], [2.0]], index=["A", "B"]), "posterior must"),
        (pd.DataFrame([[1.0, 0.0]], columns=["A", "B"]), np.ones((2, 1)), "reference_profiles must"),
        (pd.DataFrame([[1.0, 0.0]], columns=["A", None]), pd.DataFrame([[1.0], [2.0]], index=["A", "B"]), "null"),
        (pd.DataFrame([[1.0, 0.0]], columns=["A", "A"]), pd.DataFrame([[1.0], [2.0]], index=["A", "B"]), "duplicate"),
        (pd.DataFrame([[1.0, 0.0]], columns=["A/B", "A_B"]), pd.DataFrame([[1.0], [2.0]], index=["A", "B"]), "collide"),
        (pd.DataFrame([[1.0]], columns=["A"]), pd.DataFrame([[1.0], [2.0]], index=["A", "B"]), "missing=['B'], extra=[]"),
        (pd.DataFrame([[1.0, 0.0, 0.0]], columns=["A", "B", "C"]), pd.DataFrame([[1.0], [2.0]], index=["A", "B"]), "missing=[], extra=['C']"),
    ],
)
def test_reference_allocation_rejects_ambiguous_category_axes(
    posterior,
    reference,
    message,
):
    with pytest.raises(ValueError) as exc_info:
        posterior_reference_allocation(np.array([[5.0]]), posterior, reference)

    assert message in str(exc_info.value)


@pytest.mark.parametrize("storage", ["dense", "sparse"])
def test_graph_aggregate_honors_mask_and_preserves_zero_support(storage):
    raw = np.array([[1.0, 0.0], [0.0, 10.0], [100.0, 100.0]])
    X = sparse.csr_matrix(raw) if storage == "sparse" else raw
    adata = AnnData(X=X)
    kernel = GraphAggregateKernel(SimpleNamespace(rec_alpha=1.0), None)

    result = kernel.run(
        adata,
        neighbor_idx_matrix=np.array([[1, 2], [0, 2], [0, 1]]),
        coupling_matrix=np.array([[2.0, 0.0, 3.0], [0.0, 0.0, 0.0]]),
        valid_neighbor_mask=np.array([[True, False], [False, False], [True, False]]),
    )

    values = result.X.toarray() if sparse.issparse(result.X) else np.asarray(result.X)
    np.testing.assert_allclose(values, [[0.0, 10.0], [0.0, 10.0], [1.0, 0.0]])


def test_graph_aggregate_rejects_mass_on_padded_slot():
    kernel = GraphAggregateKernel(SimpleNamespace(rec_alpha=1.0), None)
    with pytest.raises(ValueError, match="mass on invalid neighbor slots"):
        kernel.run(
            AnnData(X=np.eye(2)),
            neighbor_idx_matrix=np.array([[1, 0], [0, 1]]),
            coupling_matrix=np.array([[1.0, 1.0], [0.25, 0.0]]),
            valid_neighbor_mask=np.array([[True, False], [True, False]]),
        )


@pytest.mark.parametrize(
    ("coupling", "message"),
    [
        (np.ones((3, 2)), "must have shape"),
        (np.array([[np.inf, 1.0], [0.0, 0.0]]), "finite"),
        (np.array([[-1.0, 1.0], [0.0, 0.0]]), "non-negative"),
    ],
)
def test_graph_aggregate_rejects_invalid_coupling(coupling, message):
    kernel = GraphAggregateKernel(SimpleNamespace(rec_alpha=1.0), None)
    with pytest.raises(ValueError, match=message):
        kernel.run(
            AnnData(X=np.eye(2)),
            neighbor_idx_matrix=np.array([[1, 0], [0, 1]]),
            coupling_matrix=coupling,
        )


def test_graph_aggregate_with_no_neighbor_slots_is_identity():
    raw = np.array([[1.0, 2.0], [3.0, 4.0]])
    result = GraphAggregateKernel(SimpleNamespace(rec_alpha=1.0), None).run(
        AnnData(X=raw.copy()),
        neighbor_idx_matrix=np.empty((2, 0), dtype=np.int64),
        coupling_matrix=np.empty((0, 2), dtype=np.float64),
        valid_neighbor_mask=np.empty((2, 0), dtype=bool),
    )

    np.testing.assert_array_equal(result.X, raw)


def test_similarity_to_distance_preserves_minimization_order_and_padding():
    distance_op = _load_functions(
        "revise/backend/ops/distance.py",
        ["similarity_to_distance"],
        {"Any": object, "np": np},
    )
    distance = distance_op.similarity_to_distance(
        np.array([[0.2, 999.0], [0.8, 0.4]]),
        np.array([[True, False], [True, True]]),
    )

    assert distance[1, 0] < distance[1, 1] < distance[0, 0]
    assert np.isinf(distance[0, 1])


def test_application_svc_fields_preserve_their_labeled_row_identity():
    adapter = _load_functions(
        "revise/backend/adapters.py",
        ["_extract_probs", "_build_svc"],
        {"Any": object, "Dict": dict, "SVC": SVC, "np": np, "pd": pd},
    )

    expr = AnnData(
        X=np.array([[20.0], [10.0]]),
        obs=pd.DataFrame(
            {"Level1": ["B", "A"], "Confidence": [0.8, 0.9]},
            index=["cell-b", "cell-a"],
        ),
    )
    spatial = AnnData(
        X=np.array([[2.0], [1.0]]),
        obs=pd.DataFrame(index=["spot-b", "spot-a"]),
    )
    spatial.obsm["Level1"] = pd.DataFrame(
        [[0.2, 0.8], [0.9, 0.1]],
        index=spatial.obs_names,
        columns=["A", "B"],
    )
    ctx = SimpleNamespace(
        columns={"cell_type_col": "Level1", "confidence_col": "Confidence"},
        runtime={"strategy": "ScSvcApplicationStrategy", "svc_kind": "sc"},
        route="application:sc",
        stage_trace=[],
        quality_metrics={},
    )

    svc = adapter._build_svc(
        ctx,
        {"sc_svc_expr": expr, "sc_svc_spatial": spatial},
        default_key="sc_svc_expr",
        expr=expr,
        spatial=spatial,
    )

    assert list(svc.cell_type_label.index) == ["cell-b", "cell-a"]
    assert list(svc.cell_type_label) == ["B", "A"]
    assert list(svc.confidence.index) == ["cell-b", "cell-a"]
    assert list(svc.cell_type_probs.index) == ["spot-b", "spot-a"]
    np.testing.assert_array_equal(svc.expr.X.ravel(), [20.0, 10.0])


def test_sr_benchmark_output_preserves_svc_obs_order():
    runner = _load_functions(
        "revise/backend/runners/sc_svc_sr_benchmark.py",
        ["_build_svc_adata"],
        {"sc": SimpleNamespace(AnnData=AnnData)},
    )
    target = SimpleNamespace(
        svc_obs=pd.DataFrame(
            {
                "cell_id": ["cell-b", "cell-a"],
                "cell_type": ["B", "A"],
            }
        )
    )

    output = runner._build_svc_adata(
        target,
        np.array([[20.0], [10.0]]),
        pd.Index(["gene"]),
    )

    assert list(output.obs_names) == ["cell-b", "cell-a"]
    assert list(output.obs["cell_type"]) == ["B", "A"]
    np.testing.assert_array_equal(output.X.ravel(), [20.0, 10.0])
