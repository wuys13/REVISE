from __future__ import annotations

import sys
import types

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, read_h5ad
from scipy import sparse

from revise.application.ist_assembly import assemble_ist
from revise.application.ist_ot import smooth_spatial_overlap, typed_labels


def _carriers(*, spatial_rows=4):
    spatial = AnnData(
        X=np.array(
            [
                [0.0, 10.0],
                [10.0, 0.0],
                [5.0, 5.0],
                [1.0, 9.0],
            ][:spatial_rows]
        ),
        obs=pd.DataFrame(
            {
                "SVC_cluster": ["a", "a", "b", "c"][:spatial_rows],
                "parent": ["T", "T", "T", "B"][:spatial_rows],
                "Confidence": [0.91, 0.82, 0.73, 0.64][:spatial_rows],
            },
            index=["s0", "s1", "s2", "s3"][:spatial_rows],
        ),
        var=pd.DataFrame(index=["g2", "g1"]),
    )
    spatial.obsm["spatial"] = np.c_[np.arange(spatial_rows), np.zeros(spatial_rows)]
    spatial.obsp["spatial_connectivities"] = sparse.csr_matrix(
        np.ones((spatial_rows, spatial_rows))
    )
    expression = AnnData(
        X=sparse.csr_matrix(
            [
                [10.0, 0.0, 100.0],
                [0.0, 10.0, 200.0],
                [5.0, 5.0, 300.0],
                [9.0, 1.0, 400.0],
            ]
        ),
        obs=pd.DataFrame(
            {
                "SVC_cluster": ["a", "a", "b", "c"],
                "parent": ["T", "T", "T", "B"],
            },
            index=["d0", "d1", "d2", "d3"],
        ),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
    )
    return spatial, expression


def _install_spatial_graph(monkeypatch, adjacency, calls=None):
    module = types.ModuleType("squidpy")

    def spatial_neighbors(work, **kwargs):
        if calls is not None:
            calls.append((work, kwargs))
        assert not list(work.obsp.keys())
        assert kwargs["coord_type"] == "generic"
        assert kwargs["copy"] is True
        graph = sparse.csr_matrix(adjacency, dtype=np.float64)
        return graph, graph.copy()

    module.gr = types.SimpleNamespace(spatial_neighbors=spatial_neighbors)
    monkeypatch.setitem(sys.modules, "squidpy", module)


def _fixture_graph():
    # s0-s1 is same-cluster; s1-s2 is cross-cluster and must be discarded;
    # s2 and s3 are isolated after filtering.
    return np.array(
        [
            [1, 1, 0, 0],
            [1, 1, 1, 0],
            [0, 1, 1, 0],
            [0, 0, 0, 1],
        ],
        dtype=float,
    )


def test_smoothing_rebuilds_graph_drops_self_and_cross_cluster_edges(monkeypatch):
    spatial, _ = _carriers()
    calls = []
    _install_spatial_graph(monkeypatch, _fixture_graph(), calls)

    observed, metadata = smooth_spatial_overlap(
        spatial,
        np.array([0, 1]),
        typed_labels(spatial, "SVC_cluster", "spatial carrier"),
        spatial_weight=0.2,
    )

    np.testing.assert_allclose(
        observed,
        [
            [2.0, 8.0],
            [8.0, 2.0],
            [5.0, 5.0],
            [1.0, 9.0],
        ],
    )
    assert metadata == {"same_cluster_edges": 2, "isolated_rows": 2}
    assert len(calls) == 1
    np.testing.assert_array_equal(
        spatial.obsp["spatial_connectivities"].toarray(), np.ones((4, 4))
    )


def test_within_cluster_preserves_axes_projects_full_genes_and_reports_groups(
    monkeypatch,
):
    spatial, expression = _carriers()
    _install_spatial_graph(monkeypatch, _fixture_graph())
    observed_calls = []

    def couple(source, target, cost, *, method):
        observed_calls.append((source.copy(), target.copy(), cost.copy(), method))
        return np.array([[0.75, 0.25], [0.2, 0.8]])

    monkeypatch.setattr("revise.application.ist_ot.OTKernel.couple", couple)
    result = assemble_ist(
        spatial,
        expression,
        mapping="within_cluster",
        seed=42,
        broad_column="parent",
        ot_options={"gene_block_size": 1},
    )

    np.testing.assert_allclose(
        result.X,
        [
            [7.5, 2.5, 125.0],
            [2.0, 8.0, 180.0],
            [5.0, 5.0, 300.0],
            [9.0, 1.0, 400.0],
        ],
    )
    assert result.obs_names.tolist() == spatial.obs_names.tolist()
    assert result.var_names.tolist() == expression.var_names.tolist()
    np.testing.assert_array_equal(result.obs["Confidence"], spatial.obs["Confidence"])
    assert len(observed_calls) == 1
    source, target, cost, method = observed_calls[0]
    np.testing.assert_allclose(source, [10, 10])
    np.testing.assert_allclose(target, [10, 10])
    np.testing.assert_allclose(cost, [[0.1909830056, 1], [1, 0.1909830056]])
    assert method == "tacco"
    metadata = result.uns["revise_reconstruction"]
    assert metadata["ot_effective_options"] == {
        "method": "tacco",
        "spatial_weight": 0.2,
        "max_cost_entries": 2_000_000,
        "gene_block_size": 1,
    }
    assert metadata["ot_group_labels"] == ["str:a", "str:b", "str:c"]
    assert metadata["ot_group_spatial_sizes"] == [2, 1, 1]
    assert metadata["ot_group_candidate_sizes"] == [2, 1, 1]
    assert len(metadata["ot_group_timings_seconds"]) == 3
    assert "coupling" not in repr(metadata).lower()


def test_outside_cluster_uses_cluster_means_actual_broad_column_and_unequal_axes(
    monkeypatch,
):
    spatial, expression = _carriers()
    _install_spatial_graph(monkeypatch, _fixture_graph())
    observed_calls = []

    def couple(source, target, cost, *, method):
        observed_calls.append((source.copy(), target.copy(), cost.copy(), method))
        return np.array([[0.8, 0.2], [0.5, 0.5], [0.1, 0.9]])

    monkeypatch.setattr("revise.application.ist_ot.OTKernel.couple", couple)
    result = assemble_ist(
        spatial,
        expression,
        mapping="outside_cluster",
        seed=42,
        broad_column="parent",
    )

    # Broad type T receives cluster-mean profiles a=[5,5,150], b=[5,5,300].
    np.testing.assert_allclose(
        result.X,
        [
            [5, 5, 180],
            [5, 5, 225],
            [5, 5, 285],
            [9, 1, 400],
        ],
    )
    assert len(observed_calls) == 1
    source, target, _, _ = observed_calls[0]
    np.testing.assert_allclose(source, [10, 10, 10])
    np.testing.assert_allclose(target, [20, 10])
    metadata = result.uns["revise_reconstruction"]
    assert metadata["ot_broad_column"] == "parent"
    assert metadata["ot_group_labels"] == ["str:T", "str:B"]
    assert metadata["ot_group_spatial_sizes"] == [3, 1]
    assert metadata["ot_group_candidate_sizes"] == [2, 1]


def test_single_candidate_path_still_validates_masses_without_calling_solver(monkeypatch):
    spatial, expression = _carriers(spatial_rows=1)
    expression = expression[[0]].copy()
    monkeypatch.setattr(
        "revise.application.ist_ot.OTKernel.couple",
        lambda *args, **kwargs: pytest.fail("single donor must not invoke TACCO"),
    )

    result = assemble_ist(
        spatial, expression, mapping="within_cluster", seed=42
    )

    np.testing.assert_allclose(result.X, [[10, 0, 100]])


def test_ot_metadata_round_trips_through_h5ad(monkeypatch, tmp_path):
    spatial, expression = _carriers(spatial_rows=1)
    expression = expression[[0]].copy()
    result = assemble_ist(
        spatial, expression, mapping="within_cluster", seed=42
    )
    path = tmp_path / "assembled.h5ad"

    result.write_h5ad(path)
    restored = read_h5ad(path)

    metadata = restored.uns["revise_reconstruction"]
    assert metadata["ist_mapping"] == "within_cluster"
    assert metadata["ot_effective_options"]["method"] == "tacco"
    np.testing.assert_array_equal(metadata["ot_group_spatial_sizes"], [1])


def test_cost_cap_is_checked_before_distance_matrix(monkeypatch):
    spatial, expression = _carriers(spatial_rows=2)
    expression = expression[:2].copy()
    _install_spatial_graph(monkeypatch, np.ones((2, 2)))
    monkeypatch.setattr(
        "revise.application.ist_ot.bhattacharyya_distance",
        lambda *args: pytest.fail("cost function must not run beyond allocation cap"),
    )

    with pytest.raises(
        ValueError,
        match=(
            "requires 4 entries.*max_cost_entries=3.*outside_cluster.*raise "
            "max_cost_entries"
        ),
    ):
        assemble_ist(
            spatial,
            expression,
            mapping="within_cluster",
            seed=42,
            ot_options={"max_cost_entries": 3},
        )


@pytest.mark.parametrize(
    ("mutation", "mapping", "message"),
    [
        (
            lambda sp, ex: setattr(sp, "var_names", ["x", "y"]),
            "within_cluster",
            "no overlapping genes",
        ),
        (
            lambda sp, ex: ex.obs.__setitem__("SVC_cluster", ["z"] * 4),
            "within_cluster",
            "no within-cluster candidates",
        ),
        (
            lambda sp, ex: sp.obs.__setitem__(
                "SVC_cluster", [None, "a", "b", "c"]
            ),
            "within_cluster",
            "null SVC_cluster labels",
        ),
        (
            lambda sp, ex: sp.X.__setitem__((0, 0), -1),
            "within_cluster",
            "non-negative",
        ),
        (
            lambda sp, ex: sp.X.__setitem__((2, slice(None)), 0),
            "within_cluster",
            "strictly positive",
        ),
        (
            lambda sp, ex: ex.obs.__setitem__(
                "parent", ["T", "B", "T", "B"]
            ),
            "outside_cluster",
            "inconsistent parent",
        ),
    ],
)
def test_invalid_inputs_fail_clearly(monkeypatch, mutation, mapping, message):
    spatial, expression = _carriers()
    _install_spatial_graph(monkeypatch, _fixture_graph())
    mutation(spatial, expression)

    with pytest.raises(ValueError, match=message):
        assemble_ist(
            spatial,
            expression,
            mapping=mapping,
            seed=42,
            broad_column="parent",
        )


def test_invalid_solver_weights_are_not_repaired(monkeypatch):
    spatial, expression = _carriers(spatial_rows=2)
    expression = expression[:2].copy()
    _install_spatial_graph(monkeypatch, np.ones((2, 2)))
    monkeypatch.setattr(
        "revise.application.ist_ot.OTKernel.couple",
        lambda *args, **kwargs: np.array([[np.nan, 0.0], [0.0, 1.0]]),
    )

    with pytest.raises(ValueError, match="weights must be finite and non-negative"):
        assemble_ist(spatial, expression, mapping="within_cluster", seed=42)


def test_real_tacco_solves_small_nonmock_assembly(monkeypatch):
    spatial, expression = _carriers(spatial_rows=2)
    expression = expression[:2].copy()
    _install_spatial_graph(monkeypatch, np.ones((2, 2)))

    result = assemble_ist(
        spatial,
        expression,
        mapping="within_cluster",
        seed=42,
        ot_options={"gene_block_size": 2},
    )

    assert result.shape == (2, 3)
    assert np.isfinite(result.X).all()
    assert (result.X >= 0).all()
    np.testing.assert_allclose(result.X[:, :2].sum(axis=1), [10, 10])
