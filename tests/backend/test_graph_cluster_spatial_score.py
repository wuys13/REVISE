import importlib
import sys
import types
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

from revise.backend.ops.assignment import AssignmentState
from revise.backend.ops.assignment_guidance import AssignmentGuidanceCollector


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
            getattr(parent, attribute, _MISSING) if parent is not None else _MISSING
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
    unrestored = [
        module_name
        for module_name, module in modules.items()
        if (
            (module is _MISSING and module_name in sys.modules)
            or (module is not _MISSING and sys.modules.get(module_name) is not module)
        )
    ]
    if unrestored:
        raise AssertionError(f"failed to restore isolated modules: {unrestored}")


@pytest.fixture
def graph_cluster_module():
    snapshot = _snapshot_modules(ISOLATED_GRAPH_CLUSTER_MODULE_NAMES)
    for module_name in ISOLATED_GRAPH_CLUSTER_MODULE_NAMES:
        sys.modules.pop(module_name, None)
    sys.modules["scanpy"] = types.ModuleType("scanpy")
    sys.modules["squidpy"] = types.ModuleType("squidpy")
    try:
        yield importlib.import_module("revise.backend.kernels.graph_cluster")
    finally:
        _restore_modules(snapshot)


@pytest.fixture
def spatial_score(graph_cluster_module):
    return graph_cluster_module.get_spatial_score


def _patch_graph_cluster_runtime(module, monkeypatch, captured):
    gene_graph = sparse.csr_matrix([[0.0, 1.0], [1.0, 0.0]])
    spatial_graph = sparse.csr_matrix([[1.0, 0.0], [0.0, 1.0]])

    def highly_variable_genes(adata, **_kwargs):
        adata.var["highly_variable"] = True

    def neighbors(adata, **_kwargs):
        adata.obsp["connectivities"] = gene_graph.copy()

    def spatial_neighbors(adata):
        adata.obsp["spatial_connectivities"] = spatial_graph.copy()

    def leiden(adata, *, adjacency, key_added, **_kwargs):
        captured["leiden_adjacency"] = adjacency.copy()
        captured.setdefault("leiden_adjacencies", []).append(adjacency.copy())
        adata.obs[key_added] = ["0", "1"]

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
        lambda *_args, **_kwargs: 1.0,
    )
    return gene_graph


def _run_graph_cluster(
    module,
    monkeypatch,
    *,
    obsm,
    guidance="prefer",
    guidance_state=None,
):
    captured = {"warnings": []}
    base_graph = _patch_graph_cluster_runtime(module, monkeypatch, captured)
    obs_names = ["cell-1", "cell-2"]
    matrices = {
        key: pd.DataFrame(values, index=obs_names)
        for key, values in obsm.items()
    }
    adata = AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame(index=obs_names),
        var=pd.DataFrame(index=["g1", "g2"]),
        obsm=matrices,
    )

    def capture_compatibility(left, right, *, support, **_kwargs):
        captured["left_state"] = left
        captured["right_state"] = right
        captured["support"] = support
        return np.ones(len(support[0]), dtype=np.float64)

    monkeypatch.setattr(module, "assignment_compatibility", capture_compatibility)
    monkeypatch.setattr(
        module,
        "graph_guidance",
        lambda weights, _affinity, _strength: weights * 3.0,
    )
    collector = AssignmentGuidanceCollector()
    logger = SimpleNamespace(
        info=lambda *_args, **_kwargs: None,
        warning=lambda *args, **_kwargs: captured["warnings"].append(args),
    )
    config = SimpleNamespace(
        rec_random_state=0,
        rec_graph_alpha=0.5,
        rec_graph_method="pca",
        posterior_conditioning_enabled=guidance != "off",
        posterior_conditioning_mode="cost",
        posterior_conditioning_key="Level1",
        posterior_conditioning_beta=1.0,
        posterior_conditioning_min_affinity=0.0,
        posterior_conditioning_strict=guidance == "require",
        posterior_conditioning_cost_strength=1.0,
        assignment_guidance_callback=collector.callback,
        assignment_guidance_route="real2real:cellular",
        plot_flag=False,
    )
    module.GraphClusterKernel(config, logger).run(
        adata,
        resolution=[0.5],
        label="Level2",
        guidance_state=guidance_state,
        problem_key="standard-sc:A",
    )
    captured["events"] = collector.events
    return captured, base_graph


def test_graph_cluster_uses_explicit_level2_state_not_global_level1_key(
    graph_cluster_module,
    monkeypatch,
):
    module = graph_cluster_module
    level1 = np.array([[0.9, 0.1], [0.2, 0.8]])
    level2 = np.array([[0.1, 0.9], [0.8, 0.2]])
    state = AssignmentState(
        values=level2,
        observation_labels=("cell-1", "cell-2"),
        category_labels=("sub-a", "sub-b"),
        source="local_anchoring:obsm[Level2]",
        level="Level2",
        value_semantics="soft",
        lineage=[],
    )
    captured, _ = _run_graph_cluster(
        module,
        monkeypatch,
        obsm={
            "Level1": level1,
            "Level2": level2,
        },
        guidance_state=state,
    )

    np.testing.assert_allclose(captured["left_state"].values, level2)
    assert captured["left_state"].level == "Level2"
    assert captured["right_state"].source == "local_anchoring:obsm[Level2]"
    np.testing.assert_allclose(
        captured["leiden_adjacency"].toarray(),
        [[0.0, 3.0], [3.0, 0.0]],
    )
    assert captured["events"][0]["outcome"] == "applied"


@pytest.mark.parametrize(
    ("guidance", "expected_outcome", "raises"),
    [
        ("prefer", "fallback", False),
        ("require", "failed", True),
    ],
)
def test_graph_cluster_missing_level2_obeys_prefer_require_policy(
    graph_cluster_module,
    monkeypatch,
    guidance,
    expected_outcome,
    raises,
):
    module = graph_cluster_module
    if raises:
        with pytest.raises(ValueError, match="assignment guidance unavailable"):
            _run_graph_cluster(
                module,
                monkeypatch,
                obsm={},
                guidance=guidance,
                guidance_state=None,
            )
        return

    captured, base_graph = _run_graph_cluster(
        module,
        monkeypatch,
        obsm={},
        guidance=guidance,
        guidance_state=None,
    )

    assert "left_state" not in captured
    assert captured["events"][0]["outcome"] == expected_outcome
    assert captured["events"][0]["reason"] == "assignment_missing"
    np.testing.assert_allclose(
        captured["leiden_adjacency"].toarray(),
        base_graph.toarray(),
    )


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
