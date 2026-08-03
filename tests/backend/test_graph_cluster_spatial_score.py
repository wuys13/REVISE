import importlib
import sys
import types
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

from revise.backend.ops.assignment_guidance import (
    AssignmentGuidanceCollector,
    NotApplicableReason,
)


ISOLATED_GRAPH_CLUSTER_MODULE_NAMES = (
    "scanpy",
    "squidpy",
    "revise.backend.kernels.graph_cluster",
)
_MISSING = object()


def _snapshot_modules(module_names):
    modules = {
        module_name: sys.modules.get(module_name, _MISSING)
        for module_name in module_names
    }
    parent_attributes = {}
    for module_name in module_names:
        parent_name, separator, attribute = module_name.rpartition(".")
        if not separator:
            continue
        parent = sys.modules.get(parent_name)
        parent_attributes[(parent_name, attribute)] = (
            getattr(parent, attribute, _MISSING)
            if parent is not None
            else _MISSING
        )
    return modules, parent_attributes


def _restore_modules(snapshot) -> None:
    modules, parent_attributes = snapshot
    for module_name in modules:
        sys.modules.pop(module_name, None)
    for module_name, module in modules.items():
        if module is not _MISSING:
            sys.modules[module_name] = module
    for (parent_name, attribute), value in parent_attributes.items():
        parent = sys.modules.get(parent_name)
        if parent is None:
            continue
        if value is _MISSING:
            if hasattr(parent, attribute):
                delattr(parent, attribute)
        else:
            setattr(parent, attribute, value)


@pytest.fixture
def graph_cluster_module():
    snapshot = _snapshot_modules(ISOLATED_GRAPH_CLUSTER_MODULE_NAMES)
    for module_name in ISOLATED_GRAPH_CLUSTER_MODULE_NAMES:
        sys.modules.pop(module_name, None)
    scanpy = types.ModuleType("scanpy")
    scanpy.pp = SimpleNamespace()
    scanpy.tl = SimpleNamespace()
    scanpy.pl = SimpleNamespace()
    squidpy = types.ModuleType("squidpy")
    squidpy.gr = SimpleNamespace()
    sys.modules["scanpy"] = scanpy
    sys.modules["squidpy"] = squidpy
    try:
        yield importlib.import_module("revise.backend.kernels.graph_cluster")
    finally:
        _restore_modules(snapshot)


@pytest.fixture
def spatial_score(graph_cluster_module):
    return graph_cluster_module.get_spatial_score


def _config(collector=None):
    return SimpleNamespace(
        rec_random_state=11,
        rec_graph_alpha=0.25,
        rec_graph_method="joint",
        plot_flag=False,
        # Temporary compatibility fields are used only by
        # record_not_applicable until U7 removes the old collector.
        assignment_guidance_policy="off",
        assignment_guidance_callback=(
            None if collector is None else collector.callback
        ),
        assignment_guidance_route="sc_svc:segmentation",
        posterior_conditioning_beta=1.0,
        posterior_conditioning_min_affinity=0.05,
        posterior_conditioning_cost_strength=0.2,
    )


def _adata(level1_q):
    names = ["cell-1", "cell-2", "cell-3", "cell-4"]
    adata = AnnData(
        X=np.ones((4, 3), dtype=np.float64),
        obs=pd.DataFrame(
            {"Level1": ["A", "A", "B", "B"]}, index=names
        ),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
    )
    adata.obsm["Level1"] = pd.DataFrame(
        level1_q, index=names, columns=["A", "B"]
    )
    adata.obsm["Level2"] = pd.DataFrame(
        [[0.8, 0.2], [0.7, 0.3], [0.2, 0.8], [0.1, 0.9]],
        index=names,
        columns=["a1", "a2"],
    )
    return adata


def _patch_graph_runtime(module, monkeypatch):
    gene_graph = sparse.csr_matrix(
        [
            [0.0, 1.0, 0.5, 0.0],
            [1.0, 0.0, 0.0, 0.5],
            [0.5, 0.0, 0.0, 1.0],
            [0.0, 0.5, 1.0, 0.0],
        ]
    )
    spatial_graph = sparse.csr_matrix(
        [
            [0.0, 1.0, 0.0, 0.0],
            [1.0, 0.0, 1.0, 0.0],
            [0.0, 1.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, 0.0],
        ]
    )
    captured = {"leiden": []}

    def highly_variable_genes(adata, **_kwargs):
        adata.var["highly_variable"] = True

    def neighbors(adata, **_kwargs):
        adata.obsp["connectivities"] = gene_graph.copy()

    def spatial_neighbors(adata):
        adata.obsp["spatial_connectivities"] = spatial_graph.copy()

    def leiden(adata, *, adjacency, resolution, key_added, random_state):
        captured["leiden"].append(
            {
                "resolution": resolution,
                "key_added": key_added,
                "random_state": random_state,
                "adjacency": adjacency.copy(),
            }
        )
        labels = (
            ["0", "1", "0", "1"]
            if float(resolution) == 0.5
            else ["0", "0", "1", "1"]
        )
        adata.obs[key_added] = pd.Categorical(labels)

    module.sc.pp = SimpleNamespace(
        normalize_total=lambda *_args, **_kwargs: None,
        log1p=lambda *_args, **_kwargs: None,
        highly_variable_genes=highly_variable_genes,
        pca=lambda *_args, **_kwargs: None,
        neighbors=neighbors,
    )
    module.sc.tl = SimpleNamespace(leiden=leiden)
    module.sc.pl = SimpleNamespace(scatter=lambda *_args, **_kwargs: None)
    module.sq.gr = SimpleNamespace(spatial_neighbors=spatial_neighbors)

    from revise.backend.ops import coefficients

    monkeypatch.setattr(
        coefficients,
        "get_weighted_align_score",
        lambda _adata, *, res, label: 0.3 if float(res) == 0.5 else 0.9,
    )
    return gene_graph, spatial_graph, captured


def _run_kernel(module, monkeypatch, level1_q, *, collector=None):
    gene_graph, spatial_graph, captured = _patch_graph_runtime(
        module, monkeypatch
    )
    kernel = module.GraphClusterKernel(
        _config(collector),
        SimpleNamespace(
            info=lambda *_args, **_kwargs: None,
            warning=lambda *_args, **_kwargs: None,
        ),
    )
    output, metrics, best_res = kernel.run(
        _adata(level1_q),
        resolution=[0.5, 1.0],
        label="Level2",
    )
    return output, metrics, best_res, captured, gene_graph, spatial_graph


def test_same_argmax_different_ga_q_preserves_graph_edges_and_leiden(
    graph_cluster_module, monkeypatch
):
    first_q = np.array(
        [[0.9, 0.1], [0.6, 0.4], [0.2, 0.8], [0.1, 0.9]]
    )
    second_q = np.array(
        [[0.51, 0.49], [0.99, 0.01], [0.49, 0.51], [0.3, 0.7]]
    )

    first = _run_kernel(graph_cluster_module, monkeypatch, first_q)
    second = _run_kernel(graph_cluster_module, monkeypatch, second_q)

    (
        first_output,
        first_metrics,
        first_best,
        first_calls,
        gene_graph,
        spatial_graph,
    ) = first
    second_output, second_metrics, second_best, second_calls, *_ = second
    assert first_best == second_best == 1.0
    pd.testing.assert_frame_equal(first_metrics, second_metrics)
    assert len(first_calls["leiden"]) == len(second_calls["leiden"]) == 2
    for left, right in zip(first_calls["leiden"], second_calls["leiden"]):
        assert left["resolution"] == right["resolution"]
        assert left["key_added"] == right["key_added"]
        assert left["random_state"] == right["random_state"] == 11
        np.testing.assert_allclose(
            left["adjacency"].toarray(), right["adjacency"].toarray()
        )
        np.testing.assert_allclose(
            left["adjacency"].toarray(),
            (0.75 * gene_graph + 0.25 * spatial_graph).toarray(),
        )
    assert first_output.obs["leiden_0.5"].tolist() == second_output.obs[
        "leiden_0.5"
    ].tolist()
    assert first_output.obs["leiden_1.0"].tolist() == second_output.obs[
        "leiden_1.0"
    ].tolist()
    assert "assignment_guided_connectivities" not in first_output.obsp
    assert "assignment_guided_connectivities" not in second_output.obsp


def test_temporary_not_applicable_event_does_not_change_clustering(
    graph_cluster_module, monkeypatch
):
    collector = AssignmentGuidanceCollector()
    q = np.array(
        [[0.9, 0.1], [0.6, 0.4], [0.2, 0.8], [0.1, 0.9]]
    )
    kernel = graph_cluster_module.GraphClusterKernel(
        _config(collector),
        SimpleNamespace(
            info=lambda *_args, **_kwargs: None,
            warning=lambda *_args, **_kwargs: None,
        ),
    )
    kernel.record_not_applicable(
        problem_key="standard-sc:skipped",
        reason=NotApplicableReason.INSUFFICIENT_UNITS,
        reason_details={"observed": 1, "required": 2},
    )
    _gene, _space, calls = _patch_graph_runtime(
        graph_cluster_module, monkeypatch
    )
    output, _metrics, _best = kernel.run(
        _adata(q), resolution=[0.5, 1.0], label="Level2"
    )

    assert len(calls["leiden"]) == 2
    assert output.obs["leiden_1.0"].tolist() == ["0", "0", "1", "1"]
    [event] = collector.events
    assert event["outcome"] == "not_applicable"


def test_graph_clustering_needs_no_assignment_policy_config(
    graph_cluster_module, monkeypatch
):
    q = np.array(
        [[0.9, 0.1], [0.6, 0.4], [0.2, 0.8], [0.1, 0.9]]
    )
    _gene, _space, calls = _patch_graph_runtime(
        graph_cluster_module, monkeypatch
    )
    kernel = graph_cluster_module.GraphClusterKernel(
        SimpleNamespace(
            rec_random_state=11,
            rec_graph_alpha=0.25,
            rec_graph_method="joint",
            plot_flag=False,
        ),
        SimpleNamespace(
            info=lambda *_args, **_kwargs: None,
            warning=lambda *_args, **_kwargs: None,
        ),
    )

    output, _metrics, best_res = kernel.run(
        _adata(q), resolution=[0.5, 1.0], label="Level2"
    )

    assert best_res == 1.0
    assert len(calls["leiden"]) == 2
    assert output.obs["leiden_1.0"].tolist() == ["0", "0", "1", "1"]


def test_first_resolution_leiden_failure_propagates_without_continuing(
    graph_cluster_module, monkeypatch
):
    q = np.array(
        [[0.9, 0.1], [0.6, 0.4], [0.2, 0.8], [0.1, 0.9]]
    )
    _gene, _space, _captured = _patch_graph_runtime(
        graph_cluster_module, monkeypatch
    )
    calls = []

    def fail_first_leiden(*_args, **kwargs):
        calls.append(float(kwargs["resolution"]))
        raise RuntimeError("base Leiden failed")

    graph_cluster_module.sc.tl.leiden = fail_first_leiden
    kernel = graph_cluster_module.GraphClusterKernel(
        SimpleNamespace(
            rec_random_state=11,
            rec_graph_alpha=0.25,
            rec_graph_method="joint",
            plot_flag=False,
        ),
        SimpleNamespace(
            info=lambda *_args, **_kwargs: None,
            warning=lambda *_args, **_kwargs: None,
        ),
    )

    with pytest.raises(RuntimeError, match="base Leiden failed"):
        kernel.run(_adata(q), resolution=[0.5, 1.0], label="Level2")

    assert calls == [0.5]


def _dense_source_row_oracle(graph, labels):
    same_label = labels[:, None] == labels[None, :]
    return (graph.toarray() * same_label).sum(axis=1)


def _graph_snapshot(graph):
    if sparse.isspmatrix_coo(graph):
        return graph.row.copy(), graph.col.copy(), graph.data.copy()
    return graph.indptr.copy(), graph.indices.copy(), graph.data.copy()


@pytest.mark.parametrize("graph_format", ["coo", "csr"])
def test_spatial_score_matches_dense_source_row_oracle_without_mutating_graph(
    graph_format, spatial_score
):
    labels = np.array(["a", "a", "b", "b", "c", "a"])
    graph = sparse.coo_matrix(
        (
            np.array(
                [1.5, 2.0, -0.5, 4.0, 3.0, 1.25, -2.0, 0.0, 1.0, 2.0, -1.0]
            ),
            (
                np.array([0, 0, 1, 1, 2, 3, 3, 5, 5, 5, 2]),
                np.array([0, 1, 0, 2, 3, 2, 5, 0, 1, 1, 2]),
            ),
        ),
        shape=(6, 6),
    )
    graph = graph if graph_format == "coo" else graph.tocsr()
    expected = _dense_source_row_oracle(graph, labels)
    before = _graph_snapshot(graph)
    adata = SimpleNamespace(
        obs=pd.DataFrame({"leiden_0.5": labels}),
        obsp={"spatial_connectivities": graph},
    )

    result = spatial_score(adata, res=0.5)

    assert result is adata
    np.testing.assert_allclose(
        adata.obs["spatial_score_0.5"].to_numpy(), expected
    )
    after = _graph_snapshot(graph)
    for actual, original in zip(after, before):
        np.testing.assert_array_equal(actual, original)


class _BroadcastGuardLabels(np.ndarray):
    def __new__(cls, values):
        return np.asarray(values).view(cls)

    def __eq__(self, other):
        other_shape = getattr(other, "shape", ())
        if self.ndim == 2 and len(other_shape) == 2:
            if np.broadcast_shapes(self.shape, other_shape) == (
                self.size,
                self.size,
            ):
                raise AssertionError(
                    "dense n_obs by n_obs label comparison is forbidden"
                )
        return super().__eq__(other)


class _FakeObs(dict):
    def __init__(self, key, labels):
        super().__init__({key: SimpleNamespace(to_numpy=lambda: labels)})


def test_spatial_score_compares_only_sparse_edges(spatial_score):
    labels = _BroadcastGuardLabels(["a", "a", "b", "b", "a"])
    graph = sparse.csr_matrix(
        (
            np.array([2.0, 4.0, 3.0, -1.0]),
            (np.array([0, 1, 2, 4]), np.array([1, 2, 3, 0])),
        ),
        shape=(5, 5),
    )
    adata = SimpleNamespace(
        obsp={"spatial_connectivities": graph},
        obs=_FakeObs("leiden_1", labels),
    )

    spatial_score(adata, res=1)

    np.testing.assert_allclose(
        adata.obs["spatial_score_1"], [2.0, 0.0, 3.0, 0.0, -1.0]
    )


def test_spatial_score_preserves_sparse_duplicate_reduction_semantics(
    spatial_score,
):
    graph = sparse.coo_matrix(
        (
            np.array([1e10, 1.0, -1e10], dtype=np.float32),
            (np.zeros(3, dtype=int), np.zeros(3, dtype=int)),
        ),
        shape=(1, 1),
    )
    labels = np.array(["x"])
    expected = _dense_source_row_oracle(graph, labels)
    before = _graph_snapshot(graph)
    adata = SimpleNamespace(
        obs=pd.DataFrame({"leiden_1.5": labels}),
        obsp={"spatial_connectivities": graph},
    )

    spatial_score(adata, res=1.5)

    np.testing.assert_array_equal(expected, [0.0])
    np.testing.assert_array_equal(
        adata.obs["spatial_score_1.5"].to_numpy(), expected
    )
    after = _graph_snapshot(graph)
    for actual, original in zip(after, before):
        np.testing.assert_array_equal(actual, original)


def test_spatial_score_handles_large_sparse_ring_without_dense_label_broadcast(
    spatial_score,
):
    n_obs = 50_000
    rows = np.arange(n_obs)
    cols = (rows + 1) % n_obs
    graph = sparse.csr_matrix(
        (np.ones(n_obs), (rows, cols)), shape=(n_obs, n_obs)
    )
    labels = _BroadcastGuardLabels(rows % 2)
    adata = SimpleNamespace(
        obsp={"spatial_connectivities": graph},
        obs=_FakeObs("leiden_2", labels),
    )

    result = spatial_score(adata, res=2)

    assert result is adata
    scores = adata.obs["spatial_score_2"]
    assert scores.shape == (n_obs,)
    np.testing.assert_array_equal(scores, np.zeros(n_obs))
