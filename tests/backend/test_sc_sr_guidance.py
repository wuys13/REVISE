from __future__ import annotations

import importlib
import logging
import sys
import types
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise.backend.kernels.spot_sr import SpotSrKernel
from revise.backend.ops.assignment import (
    GlobalAssignment,
    GlobalAssignmentContractError,
)


def _spot_and_virtual_inputs():
    spatial = AnnData(
        X=np.ones((2, 2), dtype=np.float64),
        obs=pd.DataFrame(index=["spot-b", "spot-a"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    spatial.obsm["major_type"] = pd.DataFrame(
        [[0.2, 0.8], [0.9, 0.1]],
        index=spatial.obs_names,
        columns=["B", "A"],
    )
    spatial.obs["major_type"] = spatial.obsm["major_type"].idxmax(axis=1)
    virtual = pd.DataFrame(
        {
            "cell_id": ["v3", "v1", "v2"],
            "spot_name": ["spot-a", "spot-b", "spot-a"],
            "cell_type": ["A", "B", "A"],
        }
    )
    return spatial, virtual


def test_strict_spot_posterior_projects_to_virtual_cells_in_mapping_order():
    from revise.backend.ops.sr_allocation import (
        project_spot_assignment_to_virtual_cells,
        spot_global_assignment,
    )

    spatial, virtual = _spot_and_virtual_inputs()
    spot_assignment = spot_global_assignment(
        spatial,
        broad_key="major_type",
        expected_categories=pd.Index(["B", "A"]),
    )
    projected = project_spot_assignment_to_virtual_cells(
        spot_assignment,
        virtual,
    )

    assert isinstance(projected, GlobalAssignment)
    assert projected.posterior.index.tolist() == ["v3", "v1", "v2"]
    assert projected.posterior.columns.tolist() == ["B", "A"]
    np.testing.assert_allclose(
        projected.posterior.to_numpy(),
        [[0.9, 0.1], [0.2, 0.8], [0.9, 0.1]],
    )
    assert projected.labels.tolist() == ["B", "A", "B"]


def test_strict_spot_posterior_rejects_missing_labels_without_fallback():
    from revise.backend.ops.sr_allocation import spot_global_assignment

    spatial, virtual = _spot_and_virtual_inputs()
    del spatial.obs["major_type"]
    virtual["pm_A"] = 1.0

    with pytest.raises(GlobalAssignmentContractError, match=r"missing obs\[major_type\]"):
        spot_global_assignment(
            spatial,
            broad_key="major_type",
            expected_categories=pd.Index(["B", "A"]),
        )


def test_strict_spot_posterior_rejects_permuted_category_axis():
    from revise.backend.ops.sr_allocation import spot_global_assignment

    spatial, _virtual = _spot_and_virtual_inputs()

    with pytest.raises(GlobalAssignmentContractError, match="category axis order"):
        spot_global_assignment(
            spatial,
            broad_key="major_type",
            expected_categories=pd.Index(["A", "B"]),
        )


def test_virtual_projection_rejects_missing_spot_before_local_ot():
    from revise.backend.ops.sr_allocation import (
        project_spot_assignment_to_virtual_cells,
        spot_global_assignment,
    )

    spatial, virtual = _spot_and_virtual_inputs()
    virtual.loc[0, "spot_name"] = "missing-spot"
    spot_assignment = spot_global_assignment(
        spatial,
        broad_key="major_type",
        expected_categories=pd.Index(["B", "A"]),
    )

    with pytest.raises(GlobalAssignmentContractError, match="projection spot mapping"):
        project_spot_assignment_to_virtual_cells(spot_assignment, virtual)


def test_virtual_subset_preserves_requested_order_and_rejects_missing_cells():
    from revise.backend.ops.sr_allocation import (
        project_spot_assignment_to_virtual_cells,
        spot_global_assignment,
        subset_virtual_assignment,
    )

    spatial, virtual = _spot_and_virtual_inputs()
    projected = project_spot_assignment_to_virtual_cells(
        spot_global_assignment(
            spatial,
            broad_key="major_type",
            expected_categories=pd.Index(["B", "A"]),
        ),
        virtual,
    )

    subset = subset_virtual_assignment(projected, ["v2", "v3"])
    assert subset.posterior.index.tolist() == ["v2", "v3"]
    np.testing.assert_allclose(subset.posterior, [[0.9, 0.1], [0.9, 0.1]])

    with pytest.raises(GlobalAssignmentContractError, match="missing"):
        subset_virtual_assignment(projected, ["absent"])


def test_route_validates_each_complete_assignment_once_before_group_subsetting(
    monkeypatch,
):
    from revise.backend.ops import sr_allocation

    real_validate = sr_allocation.validate_global_assignment
    validated_observation_counts = []

    def validate_once(assignment, *, expected_observations, expected_categories):
        validated_observation_counts.append(len(expected_observations))
        return real_validate(
            assignment,
            expected_observations=expected_observations,
            expected_categories=expected_categories,
        )

    monkeypatch.setattr(
        sr_allocation,
        "validate_global_assignment",
        validate_once,
    )
    spatial, virtual = _spot_and_virtual_inputs()
    spot_assignment = sr_allocation.spot_global_assignment(
        spatial,
        broad_key="major_type",
        expected_categories=pd.Index(["B", "A"]),
    )
    projected = sr_allocation.project_spot_assignment_to_virtual_cells(
        spot_assignment,
        virtual,
    )

    first = sr_allocation.subset_virtual_assignment(projected, ["v2", "v3"])
    second = sr_allocation.subset_virtual_assignment(projected, ["v1"])

    assert first.posterior.index.tolist() == ["v2", "v3"]
    assert second.posterior.index.tolist() == ["v1"]
    assert validated_observation_counts == [2, 3]


def test_virtual_cell_conditioning_uses_projected_soft_q_on_existing_support():
    from revise.backend.ops.sr_allocation import condition_virtual_cell_ot_cost

    posterior = pd.DataFrame(
        [[0.9, 0.1], [0.2, 0.8], [0.9, 0.1]],
        index=["v3", "v1", "v2"],
        columns=["B", "A"],
    )
    assignment = GlobalAssignment(
        labels=posterior.idxmax(axis=1),
        posterior=posterior,
    )
    support = np.array([[0, 1], [1, 2], [2, 0]], dtype=np.int64)
    cost = np.zeros_like(support, dtype=np.float64)

    conditioned = condition_virtual_cell_ot_cost(
        cost,
        assignment=assignment,
        neighbor_indices=support,
        strength=2.0,
    )

    affinity = np.einsum(
        "id,ikd->ik",
        posterior.to_numpy(),
        posterior.to_numpy()[support],
    )
    np.testing.assert_allclose(
        conditioned,
        2.0 * -np.log(np.maximum(affinity, 1e-12)),
    )


def test_zero_strength_preserves_the_unconditioned_cost_object():
    from revise.backend.ops.sr_allocation import condition_virtual_cell_ot_cost

    posterior = pd.DataFrame(
        [[0.9, 0.1]],
        index=["v1"],
        columns=["A", "B"],
    )
    assignment = GlobalAssignment(
        labels=posterior.idxmax(axis=1),
        posterior=posterior,
    )
    cost = np.array([[0.25]], dtype=np.float64)

    result = condition_virtual_cell_ot_cost(
        cost,
        assignment=assignment,
        neighbor_indices=np.array([[0]], dtype=np.int64),
        strength=0.0,
    )

    assert result is cost


def test_mandatory_reference_allocation_remains_closed_form_and_mass_preserving():
    from revise.backend.ops.sr_allocation import mandatory_reference_allocation

    expression = np.array([[2.0, 8.0]], dtype=np.float64)
    composition = pd.DataFrame(
        [[0.75, 0.25]],
        index=["spot-a"],
        columns=["A", "B"],
    )
    reference = pd.DataFrame(
        [[9.0, 1.0], [1.0, 9.0]],
        index=["A", "B"],
        columns=["g1", "g2"],
    )

    first = mandatory_reference_allocation(expression, composition, reference)
    second = mandatory_reference_allocation(expression, composition, reference)

    np.testing.assert_allclose(first, second)
    np.testing.assert_allclose(first.sum(axis=2), expression)


def _sr_runner_inputs(*, cells_per_spot=50):
    spatial = AnnData(
        X=np.array([[2.0, 8.0], [8.0, 2.0]], dtype=np.float64),
        obs=pd.DataFrame(index=["spot-1", "spot-2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    spatial.obsm["major_type"] = pd.DataFrame(
        [[0.8, 0.2], [0.2, 0.8]],
        index=spatial.obs_names,
        columns=["A", "B"],
    )
    spatial.obs["major_type"] = spatial.obsm["major_type"].idxmax(axis=1)
    reference = AnnData(
        X=np.array([[9.0, 1.0], [1.0, 9.0]], dtype=np.float64),
        obs=pd.DataFrame(
            {"major_type": ["A", "B"]},
            index=["ref-a", "ref-b"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    rows = []
    for spot_index, spot_name in enumerate(spatial.obs_names):
        for cell_index in range(cells_per_spot):
            rows.append(
                {
                    "spot_name": str(spot_name),
                    "cell_id": f"cell-{spot_index}-{cell_index}",
                    "x": float(spot_index),
                    "y": float(cell_index),
                    "true_cell_type": "Unknown",
                }
            )
    return spatial, reference, pd.DataFrame(rows)


def _runner_config(*, graph_enabled=True, strength=0.0, seed=17):
    return SimpleNamespace(
        pm_on_cell_file="/path/that/does/not/exist.csv",
        svc_completeness=True,
        sr_assignment_seed=seed,
        cell_type_col="major_type",
        local_refinement_strength=strength,
        rec_graph_agg_enabled=graph_enabled,
        rec_graph_agg_low_conf_only=False,
        rec_graph_agg_anchor_only=False,
        rec_graph_agg_conf_weighted_alpha=False,
        rec_graph_n_neighbors=2,
        rec_graph_method="joint",
        rec_graph_alpha=0.2,
        rec_graph_exp_neighbor_num=2,
        rec_graph_spatial_neighbor_num=2,
        rec_ot_method="pot",
        rec_pot_reg=0.05,
        rec_pot_reg_m=1.0,
        rec_pot_reg_type="kl",
        rec_alpha=1.0,
        rec_match_spot_sum=True,
    )


@pytest.fixture
def load_sr_runner(monkeypatch):
    def normalize_total(adata, target_sum=1e4):
        values = (
            adata.X.toarray()
            if hasattr(adata.X, "toarray")
            else np.asarray(adata.X)
        )
        row_sums = values.sum(axis=1, keepdims=True)
        adata.X = values / np.maximum(row_sums, 1e-12) * float(target_sum)

    scanpy = types.ModuleType("scanpy")
    scanpy.pp = SimpleNamespace(normalize_total=normalize_total)
    scanpy.pl = SimpleNamespace()
    scanpy.tl = SimpleNamespace()
    scanpy.AnnData = AnnData
    monkeypatch.setitem(sys.modules, "scanpy", scanpy)
    distance = types.ModuleType("revise.backend.ops.distance")
    distance.similarity_to_distance = lambda values, support: np.where(
        support,
        np.asarray(values, dtype=np.float64),
        np.inf,
    )
    monkeypatch.setitem(sys.modules, "revise.backend.ops.distance", distance)
    kernels = importlib.import_module("revise.backend.kernels")

    class _UnusedKernel:
        def __init__(self, *_args, **_kwargs):
            pass

    monkeypatch.setitem(kernels.__dict__, "GlobalAnchoringKernel", _UnusedKernel)
    monkeypatch.setitem(kernels.__dict__, "GraphAggregateKernel", _UnusedKernel)
    monkeypatch.setitem(kernels.__dict__, "SpotSrKernel", SpotSrKernel)

    def load(name):
        module_name = f"revise.backend.runners.{name}"
        sys.modules.pop(module_name, None)
        return importlib.import_module(module_name)

    return load


def _build_runner(module, *, benchmark, strength, graph_enabled=True, seed=17):
    spatial, reference, svc_obs = _sr_runner_inputs()
    config = _runner_config(
        graph_enabled=graph_enabled,
        strength=strength,
        seed=seed,
    )
    runner = module.ScSVCSr.__new__(module.ScSVCSr)
    runner.st_adata = spatial
    runner.sc_ref_adata = reference
    runner.config = config
    runner.svc_obs = svc_obs
    runner.spot_sr = SpotSrKernel(
        config,
        logging.getLogger(f"test-sr-{benchmark}-{strength}"),
    )
    runner.graph_aggregate = SimpleNamespace(
        run=lambda *, adata, **_kwargs: adata
    )
    runner.logger = logging.getLogger(f"test-sr-{benchmark}-{strength}")
    runner.svc = {}
    if benchmark:
        runner._graphagg_confidence_cache = None
        runner._graphagg_confidence_source = None
        runner._graphagg_alpha_weight_cache = None
        runner._graphagg_posterior_source = None
    return runner


@pytest.mark.parametrize("solver", ["pot", "tacco"])
def test_application_always_conditions_executed_local_ot_and_preserves_allocation(
    solver,
    load_sr_runner,
    monkeypatch,
):
    from scipy import sparse

    snapshots = []
    for strength in (0.0, 3.0):
        module = load_sr_runner("sc_svc_sr_application")
        runner = _build_runner(
            module,
            benchmark=False,
            strength=strength,
        )
        runner.config.rec_ot_method = solver
        calls = []
        monkeypatch.setattr(
            module,
            "get_adjacency_graph",
            lambda adata, **_kwargs: sparse.eye(adata.n_obs, format="csr"),
        )

        def solve(source, target, cost, **kwargs):
            calls.append(
                {
                    "cost": np.asarray(cost).copy(),
                    "reference_measure": kwargs["reference_measure"],
                    "method": kwargs["method"],
                    "valid_support_mask": np.asarray(
                        kwargs["valid_support_mask"]
                    ).copy(),
                }
            )
            return np.outer(
                np.asarray(source) / np.sum(source),
                np.asarray(target) / np.sum(target),
            )

        monkeypatch.setattr(module, "solve_local_ot", solve)
        applied = runner.local_refinement()
        snapshots.append(
            {
                "applied": applied,
                "calls": calls,
                "allocation": runner.svc_obs[
                    ["cell_id", "spot_name", "cell_type"]
                ].reset_index(drop=True),
            }
        )

    baseline, conditioned = snapshots
    assert baseline["applied"] is True
    assert conditioned["applied"] is True
    assert len(baseline["calls"]) == len(conditioned["calls"]) == 2
    assert {call["method"] for call in conditioned["calls"]} == {solver}
    assert all(
        call["reference_measure"] is None
        for call in conditioned["calls"]
    )
    assert all(
        np.isposinf(call["cost"][~call["valid_support_mask"]]).all()
        for call in conditioned["calls"]
    )
    assert any(
        np.any(
            guided["cost"][np.isfinite(guided["cost"])]
            > base["cost"][np.isfinite(base["cost"])]
        )
        for base, guided in zip(
            baseline["calls"],
            conditioned["calls"],
        )
    )
    pd.testing.assert_frame_equal(
        conditioned["allocation"],
        baseline["allocation"],
    )


@pytest.mark.parametrize("solver", ["pot", "tacco"])
def test_benchmark_enabled_graph_conditions_every_solver_call(
    solver,
    load_sr_runner,
    monkeypatch,
):
    from scipy import sparse

    module = load_sr_runner("sc_svc_sr_benchmark")
    runner = _build_runner(
        module,
        benchmark=True,
        strength=2.0,
    )
    runner.config.rec_ot_method = solver
    calls = []
    monkeypatch.setattr(
        module,
        "get_adjacency_graph",
        lambda adata, **_kwargs: sparse.eye(adata.n_obs, format="csr"),
    )

    def solve(source, target, cost, **kwargs):
        calls.append(
            {
                "cost": np.asarray(cost).copy(),
                "reference_measure": kwargs["reference_measure"],
                "method": kwargs["method"],
                "valid_support_mask": np.asarray(
                    kwargs["valid_support_mask"]
                ).copy(),
            }
        )
        return np.outer(
            np.asarray(source) / np.sum(source),
            np.asarray(target) / np.sum(target),
        )

    monkeypatch.setattr(module, "solve_local_ot", solve)

    assert runner.local_refinement() is True
    assert len(calls) == 2
    assert {call["method"] for call in calls} == {solver}
    assert all(call["reference_measure"] is None for call in calls)
    assert all(
        np.isposinf(call["cost"][~call["valid_support_mask"]]).all()
        for call in calls
    )
    assert all(
        np.any(call["cost"][np.isfinite(call["cost"])] > 1.0)
        for call in calls
    )


def test_benchmark_disabled_graph_control_reports_no_local_ot(
    load_sr_runner,
):
    module = load_sr_runner("sc_svc_sr_benchmark")
    runner = _build_runner(
        module,
        benchmark=True,
        strength=4.0,
        graph_enabled=False,
    )

    assert runner.local_refinement() is False
    assert "sc_svc_dec_graphagg" not in runner.svc


def test_allocation_completes_before_invalid_virtual_projection_fails(
    load_sr_runner,
):
    module = load_sr_runner("sc_svc_sr_application")
    runner = _build_runner(
        module,
        benchmark=False,
        strength=1.0,
    )
    runner.svc_obs.loc[0, "cell_id"] = "cell/collision"
    runner.svc_obs.loc[1, "cell_id"] = "cell_collision"
    allocation_records = []
    runner.config.sr_allocation_callback = allocation_records.append

    with pytest.raises(GlobalAssignmentContractError, match="collide"):
        runner.local_refinement()

    assert allocation_records == [
        {
            "status": "completed",
            "broad_key": "major_type",
            "n_spots": 2,
            "n_virtual_cells": 100,
            "allocation_method": "posterior_reference_allocation",
        }
    ]
