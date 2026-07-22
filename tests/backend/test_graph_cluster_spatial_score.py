import importlib
import sys
import types
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy import sparse


@pytest.fixture
def spatial_score(monkeypatch):
    monkeypatch.setitem(sys.modules, "scanpy", types.ModuleType("scanpy"))
    monkeypatch.setitem(sys.modules, "squidpy", types.ModuleType("squidpy"))
    monkeypatch.delitem(sys.modules, "revise.backend.kernels.graph_cluster", raising=False)
    module = importlib.import_module("revise.backend.kernels.graph_cluster")
    return module.get_spatial_score


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
            np.array([1.5, 2.0, -0.5, 4.0, 3.0, 1.25, -2.0, 0.0, 1.0, 2.0, -1.0]),
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
    np.testing.assert_allclose(adata.obs["spatial_score_0.5"].to_numpy(), expected)
    after = _graph_snapshot(graph)
    for actual, original in zip(after, before):
        np.testing.assert_array_equal(actual, original)


class _BroadcastGuardLabels(np.ndarray):
    def __new__(cls, values):
        return np.asarray(values).view(cls)

    def __eq__(self, other):
        other_shape = getattr(other, "shape", ())
        if self.ndim == 2 and len(other_shape) == 2:
            if np.broadcast_shapes(self.shape, other_shape) == (self.size, self.size):
                raise AssertionError("dense n_obs by n_obs label comparison is forbidden")
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

    np.testing.assert_allclose(adata.obs["spatial_score_1"], [2.0, 0.0, 3.0, 0.0, -1.0])


def test_spatial_score_preserves_sparse_duplicate_reduction_semantics(spatial_score):
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
    np.testing.assert_array_equal(adata.obs["spatial_score_1.5"].to_numpy(), expected)
    after = _graph_snapshot(graph)
    for actual, original in zip(after, before):
        np.testing.assert_array_equal(actual, original)


def test_spatial_score_handles_large_sparse_ring_without_dense_label_broadcast(
    spatial_score,
):
    n_obs = 50_000
    rows = np.arange(n_obs)
    cols = (rows + 1) % n_obs
    graph = sparse.csr_matrix((np.ones(n_obs), (rows, cols)), shape=(n_obs, n_obs))
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
