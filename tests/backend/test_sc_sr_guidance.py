from __future__ import annotations

import importlib
import json
import logging
import sys
import types
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise.backend.kernels.spot_sr import SpotSrKernel
from revise.backend.ops.assignment import AssignmentStateError
from revise.backend.ops.assignment_guidance import AssignmentGuidanceCollector


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
    virtual = pd.DataFrame(
        {
            "cell_id": ["v3", "v1", "v2"],
            "spot_name": ["spot-a", "spot-b", "spot-a"],
            "cell_type": ["A", "B", "A"],
        }
    )
    return spatial, virtual


def _guidance_config(*, guidance="prefer", compatibility_mode="cost"):
    collector = AssignmentGuidanceCollector()
    return SimpleNamespace(
        cell_type_col="major_type",
        posterior_conditioning_enabled=guidance != "off",
        posterior_conditioning_mode=compatibility_mode,
        posterior_conditioning_strict=guidance == "require",
        posterior_conditioning_beta=1.0,
        posterior_conditioning_min_affinity=0.05,
        posterior_conditioning_cost_strength=2.0,
        rec_ot_method="pot",
        assignment_guidance_route="sc_svc_sr:spot_size",
        assignment_guidance_callback=collector.callback,
    ), collector


@pytest.mark.parametrize(
    ("canonical", "legacy", "expected"),
    [("require", "off", "require"), ("off", "require", "off")],
)
def test_sr_route_prefers_canonical_guidance_over_conflicting_legacy_fields(
    canonical,
    legacy,
    expected,
):
    from revise.backend.ops import sr_allocation

    config, _collector = _guidance_config(guidance=legacy)
    config.assignment_guidance_policy = canonical

    assert sr_allocation.guidance_mode(config) == expected


def test_projected_spot_assignment_is_the_canonical_virtual_cell_state():
    from revise.backend.ops.sr_allocation import projected_virtual_assignment

    spatial, virtual = _spot_and_virtual_inputs()

    state = projected_virtual_assignment(
        spatial,
        virtual,
        broad_key="major_type",
    )

    assert state.observation_labels == ("v3", "v1", "v2")
    assert state.category_labels == ("B", "A")
    np.testing.assert_allclose(
        state.values,
        [[0.9, 0.1], [0.2, 0.8], [0.9, 0.1]],
    )
    assert state.value_semantics == "soft"
    assert state.source == "project(obsm[major_type])"
    assert state.level == "virtual_cell"
    assert state.lineage[-1] == {
        "operation": "spot_to_virtual_projection",
        "mapping": "svc_obs[cell_id->spot_name]",
        "source_axis": "spot",
        "target_axis": "virtual_cell",
        "category_axis": "major_type",
    }


def test_projected_spot_assignment_does_not_fallback_to_pm_or_final_type():
    from revise.backend.ops.sr_allocation import projected_virtual_assignment

    spatial, virtual = _spot_and_virtual_inputs()
    del spatial.obsm["major_type"]
    virtual["pm_A"] = [1.0, 1.0, 1.0]

    with pytest.raises(AssignmentStateError, match="assignment_state_unavailable"):
        projected_virtual_assignment(
            spatial,
            virtual,
            broad_key="major_type",
        )


def test_mandatory_reference_allocation_is_independent_of_guidance_numerics():
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


def test_off_guidance_does_not_load_or_validate_virtual_assignment():
    from revise.backend.ops.sr_allocation import prepare_virtual_cell_guidance

    config, collector = _guidance_config(guidance="off")
    cost = np.array([[0.0, 1.0], [1.0, 0.0]])

    result, reference, attempted = prepare_virtual_cell_guidance(
        config,
        problem_key="sr:A",
        state_loader=lambda: pytest.fail("off must not load assignment"),
        neighbor_support=np.array([[0, 1], [1, 0]], dtype=np.int64),
        distance_matrix=cost,
        source_mass=np.ones(2),
        target_mass=np.ones(2),
    )

    np.testing.assert_array_equal(result, cost)
    assert reference is None
    assert attempted is False
    assert collector.events[0]["outcome"] == "off"


@pytest.mark.parametrize(
    ("guidance", "expected_outcome"),
    [("prefer", "fallback"), ("require", "failed")],
)
def test_missing_projected_assignment_follows_policy_after_allocation(
    guidance,
    expected_outcome,
):
    from revise.backend.ops.sr_allocation import prepare_virtual_cell_guidance

    config, collector = _guidance_config(guidance=guidance)

    if guidance == "require":
        with pytest.raises(ValueError, match="assignment_missing"):
            prepare_virtual_cell_guidance(
                config,
                problem_key="sr:A",
                state_loader=lambda: None,
                neighbor_support=np.array([[0]], dtype=np.int64),
                distance_matrix=np.zeros((1, 1)),
                source_mass=np.ones(1),
                target_mass=np.ones(1),
            )
    else:
        result, reference, attempted = prepare_virtual_cell_guidance(
            config,
            problem_key="sr:A",
            state_loader=lambda: None,
            neighbor_support=np.array([[0]], dtype=np.int64),
            distance_matrix=np.zeros((1, 1)),
            source_mass=np.ones(1),
            target_mass=np.ones(1),
        )
        np.testing.assert_array_equal(result, np.zeros((1, 1)))
        assert reference is None
        assert attempted is False

    assert collector.events[0]["outcome"] == expected_outcome
    assert collector.events[0]["reason"] == "assignment_missing"


def test_cost_guidance_uses_projected_soft_state_and_records_bilateral_lineage():
    from revise.backend.ops.sr_allocation import (
        prepare_virtual_cell_guidance,
        projected_virtual_assignment,
        record_virtual_cell_guidance_terminal,
    )

    spatial, virtual = _spot_and_virtual_inputs()
    state = projected_virtual_assignment(
        spatial,
        virtual,
        broad_key="major_type",
    )
    config, collector = _guidance_config(guidance="prefer")
    support = np.array([[0, 1], [1, 2], [2, 0]], dtype=np.int64)
    cost = np.zeros_like(support, dtype=np.float64)

    guided, reference, attempted = prepare_virtual_cell_guidance(
        config,
        problem_key="sr:A",
        state_loader=lambda: state,
        neighbor_support=support,
        distance_matrix=cost,
        source_mass=np.ones(2),
        target_mass=np.ones(3),
    )
    record_virtual_cell_guidance_terminal(
        config,
        problem_key="sr:A",
        attempted=attempted,
        outcome="applied",
    )

    assert attempted is True
    assert reference is None
    assert np.all(guided >= cost)
    [event] = collector.events
    assert event["outcome"] == "applied"
    assert event["left_assignment"]["source"] == "project(obsm[major_type])"
    assert event["right_assignment"]["source"] == "project(obsm[major_type])"
    assert event["left_assignment"]["value_semantics"] == "soft"
    assert event["left_assignment"]["lineage"][-1]["mapping"] == (
        "svc_obs[cell_id->spot_name]"
    )


def test_reference_guidance_builds_pot_measure_without_changing_base_cost():
    from revise.backend.ops.sr_allocation import (
        prepare_virtual_cell_guidance,
        projected_virtual_assignment,
    )

    spatial, virtual = _spot_and_virtual_inputs()
    state = projected_virtual_assignment(
        spatial,
        virtual,
        broad_key="major_type",
    )
    config, collector = _guidance_config(
        guidance="prefer",
        compatibility_mode="reference",
    )
    support = np.array([[0, 1], [1, 2], [2, 0]], dtype=np.int64)
    cost = np.zeros_like(support, dtype=np.float64)

    unchanged, reference, attempted = prepare_virtual_cell_guidance(
        config,
        problem_key="sr-reference:A",
        state_loader=lambda: state,
        neighbor_support=support,
        distance_matrix=cost,
        source_mass=np.ones(2),
        target_mass=np.ones(3),
    )

    np.testing.assert_array_equal(unchanged, cost)
    assert reference is not None
    assert reference.shape == (2, 3)
    assert attempted is True
    assert collector.events[0]["attempted"] is True


def test_not_applicable_block_does_not_load_assignment():
    from revise.backend.ops.sr_allocation import record_virtual_cell_not_applicable

    config, collector = _guidance_config(guidance="require")
    record_virtual_cell_not_applicable(
        config,
        problem_key="sr:small",
        reason="insufficient_virtual_cells",
    )

    [event] = collector.events
    assert event["applicability"] == "not_applicable"
    assert event["availability"] == "not_checked"
    assert event["outcome"] == "not_applicable"


def test_allocation_callback_records_completion_before_projection_failure():
    from revise.backend.ops.sr_allocation import (
        projected_virtual_assignment,
        record_mandatory_allocation,
    )

    records = []
    config = SimpleNamespace(sr_allocation_callback=records.append)
    record_mandatory_allocation(
        config,
        status="completed",
        broad_key="major_type",
        n_spots=2,
        n_virtual_cells=3,
        allocation_method="posterior_reference_allocation",
    )
    spatial, virtual = _spot_and_virtual_inputs()
    virtual.loc[0, "spot_name"] = "missing-spot"

    with pytest.raises(AssignmentStateError, match="projection_source_missing"):
        projected_virtual_assignment(
            spatial,
            virtual,
            broad_key="major_type",
        )

    assert records == [
        {
            "status": "completed",
            "broad_key": "major_type",
            "n_spots": 2,
            "n_virtual_cells": 3,
            "allocation_method": "posterior_reference_allocation",
        }
    ]


def _sr_runner_inputs():
    spatial = AnnData(
        X=np.array([[2.0, 8.0]], dtype=np.float64),
        obs=pd.DataFrame(index=["spot-1"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    spatial.obsm["major_type"] = pd.DataFrame(
        [[0.75, 0.25]],
        index=spatial.obs_names,
        columns=["A", "B"],
    )
    reference = AnnData(
        X=np.array([[9.0, 1.0], [1.0, 9.0]], dtype=np.float64),
        obs=pd.DataFrame(
            {"major_type": ["A", "B"], "clusters": ["legacy-1", "legacy-2"]},
            index=["ref-a", "ref-b"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    svc_obs = pd.DataFrame(
        {
            "spot_name": ["spot-1", "spot-1", "spot-1"],
            "cell_id": ["virtual-1", "virtual-2", "virtual-3"],
            "x": [0.0, 0.1, 0.2],
            "y": [0.0, 0.1, 0.2],
            "true_cell_type": ["Unknown"] * 3,
        }
    )
    return spatial, reference, svc_obs


def _runner_config(
    *,
    guidance,
    graph_enabled=True,
    compatibility_mode="cost",
    seed=17,
):
    return SimpleNamespace(
        pm_on_cell_file="/path/that/does/not/exist.csv",
        svc_completeness=True,
        sr_assignment_seed=seed,
        cell_type_col="major_type",
        posterior_conditioning_enabled=guidance != "off",
        posterior_conditioning_mode=compatibility_mode,
        posterior_conditioning_strict=guidance == "require",
        posterior_conditioning_beta=7.0,
        posterior_conditioning_min_affinity=0.2,
        posterior_conditioning_cost_strength=13.0,
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
    distance.similarity_to_distance = lambda values, _support: np.asarray(
        values,
        dtype=np.float64,
    )
    monkeypatch.setitem(sys.modules, "revise.backend.ops.distance", distance)
    kernels = importlib.import_module("revise.backend.kernels")

    class _UnusedKernel:
        def __init__(self, *_args, **_kwargs):
            pass

    monkeypatch.setitem(
        kernels.__dict__,
        "GlobalAnchoringKernel",
        _UnusedKernel,
    )
    monkeypatch.setitem(
        kernels.__dict__,
        "GraphAggregateKernel",
        _UnusedKernel,
    )
    monkeypatch.setitem(kernels.__dict__, "SpotSrKernel", SpotSrKernel)

    def load(name):
        module_name = f"revise.backend.runners.{name}"
        sys.modules.pop(module_name, None)
        return importlib.import_module(module_name)

    return load


@pytest.mark.parametrize("guidance", ["off", "prefer", "require"])
def test_application_mandatory_allocation_is_identical_across_guidance_modes(
    guidance,
    load_sr_runner,
):
    ScSVCSr = load_sr_runner("sc_svc_sr_application").ScSVCSr

    spatial, reference, svc_obs = _sr_runner_inputs()
    config = _runner_config(guidance=guidance)
    collector = AssignmentGuidanceCollector()
    allocation_records = []
    config.assignment_guidance_callback = collector.callback
    config.assignment_guidance_route = "sc_svc_sr:spot_size"
    config.sr_allocation_callback = allocation_records.append
    runner = ScSVCSr.__new__(ScSVCSr)
    runner.st_adata = spatial
    runner.sc_ref_adata = reference
    runner.config = config
    runner.svc_obs = svc_obs
    runner.spot_sr = SpotSrKernel(
        config,
        logging.getLogger(f"test-application-allocation-{guidance}"),
    )
    runner.graph_aggregate = SimpleNamespace()
    runner.logger = logging.getLogger(
        f"test-application-allocation-{guidance}"
    )
    runner.svc = {}

    runner.local_refinement()

    assert allocation_records[0]["status"] == "completed"
    assert allocation_records[0]["broad_key"] == "major_type"
    assert runner.svc_obs["cell_type"].tolist() == ["A", "B", "A"]
    np.testing.assert_allclose(
        np.asarray(runner.svc["sc_svc_dec"].X),
        np.array(
            [
                [981.81818182, 1600.0],
                [36.36363636, 4800.0],
                [981.81818182, 1600.0],
            ]
        ),
        rtol=1e-6,
    )
    assert collector.summary() == "not_applicable"


def test_benchmark_custom_broad_column_and_disabled_graph_prefer_fallback(
    load_sr_runner,
):
    ScSVCSr = load_sr_runner("sc_svc_sr_benchmark").ScSVCSr

    spatial, reference, svc_obs = _sr_runner_inputs()
    config = _runner_config(guidance="prefer", graph_enabled=False)
    collector = AssignmentGuidanceCollector()
    allocation_records = []
    config.assignment_guidance_callback = collector.callback
    config.assignment_guidance_route = "sc_svc_sr:batch_effect"
    config.sr_allocation_callback = allocation_records.append
    runner = ScSVCSr.__new__(ScSVCSr)
    runner.st_adata = spatial
    runner.sc_ref_adata = reference
    runner.config = config
    runner.svc_obs = svc_obs
    runner.spot_sr = SpotSrKernel(
        config,
        logging.getLogger("test-benchmark-custom-broad"),
    )
    runner.graph_aggregate = SimpleNamespace()
    runner.logger = logging.getLogger("test-benchmark-custom-broad")
    runner.svc = {}
    runner._graphagg_confidence_cache = None
    runner._graphagg_confidence_source = None
    runner._graphagg_alpha_weight_cache = None
    runner._graphagg_posterior_cache = None
    runner._graphagg_posterior_source = None

    runner.local_refinement()

    assert allocation_records[0]["status"] == "completed"
    assert "sc_svc_dec" in runner.svc
    assert "sc_svc_dec_graphagg" not in runner.svc
    assert runner.svc["sc_svc_dec"].uns["sr_allocation"]["broad_key"] == (
        "major_type"
    )
    assert runner.svc["sc_svc_dec"].uns["sr_allocation"]["posterior_key"] == (
        "major_type"
    )
    [event] = collector.events
    assert event["outcome"] == "fallback"
    assert event["reason"] == "graph_branch_disabled"


def test_benchmark_required_disabled_graph_fails_before_allocation(
    load_sr_runner,
):
    ScSVCSr = load_sr_runner("sc_svc_sr_benchmark").ScSVCSr

    spatial, reference, svc_obs = _sr_runner_inputs()
    config = _runner_config(guidance="require", graph_enabled=False)
    collector = AssignmentGuidanceCollector()
    allocation_records = []
    config.assignment_guidance_callback = collector.callback
    config.assignment_guidance_route = "sc_svc_sr:batch_effect"
    config.sr_allocation_callback = allocation_records.append
    runner = ScSVCSr.__new__(ScSVCSr)
    runner.st_adata = spatial
    runner.sc_ref_adata = reference
    runner.config = config
    runner.svc_obs = svc_obs
    runner.spot_sr = SpotSrKernel(
        config,
        logging.getLogger("test-benchmark-required-disabled"),
    )
    runner.graph_aggregate = SimpleNamespace()
    runner.logger = logging.getLogger("test-benchmark-required-disabled")
    runner.svc = {}

    with pytest.raises(ValueError, match="graph_branch_disabled"):
        runner.local_refinement()

    assert allocation_records == []
    [event] = collector.events
    assert event["outcome"] == "failed"
    assert event["attempted"] is False


def test_missing_configured_broad_assignment_records_allocation_failure(
    load_sr_runner,
):
    ScSVCSr = load_sr_runner("sc_svc_sr_application").ScSVCSr
    spatial, reference, svc_obs = _sr_runner_inputs()
    del spatial.obsm["major_type"]
    config = _runner_config(guidance="off")
    records = []
    config.sr_allocation_callback = records.append
    runner = ScSVCSr.__new__(ScSVCSr)
    runner.st_adata = spatial
    runner.sc_ref_adata = reference
    runner.config = config
    runner.svc_obs = svc_obs
    runner.spot_sr = SimpleNamespace(
        run=lambda _runner: pytest.fail(
            "missing broad assignment must fail before row assignment"
        )
    )
    runner.logger = logging.getLogger("test-missing-broad-allocation")

    with pytest.raises(KeyError, match="major_type"):
        runner.local_refinement()

    assert records == [
        {
            "status": "failed",
            "broad_key": "major_type",
            "n_spots": 1,
            "n_virtual_cells": 3,
            "allocation_method": "posterior_reference_allocation",
            "reason": "allocation_failed",
        }
    ]


def test_benchmark_runner_reads_configured_ground_truth_broad_key(
    load_sr_runner,
):
    module = load_sr_runner("sc_svc_sr_benchmark")
    spatial = AnnData(
        X=np.ones((2, 2), dtype=np.float64),
        obs=pd.DataFrame(index=["spot-1", "spot-2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    spatial.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 1.0]])
    spatial.uns["all_cells_in_spot"] = {
        "spot-1": ["cell-1"],
        "spot-2": ["cell-2"],
    }
    ground_truth = AnnData(
        X=np.ones((2, 2), dtype=np.float64),
        obs=pd.DataFrame(
            {
                "cell_id": ["cell-1", "cell-2"],
                "major_type": ["A", "B"],
                "x": [0.2, 1.2],
                "y": [0.3, 1.3],
            },
            index=["gt-1", "gt-2"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    runner = module.ScSVCSr.__new__(module.ScSVCSr)
    runner.st_adata = spatial
    runner.real_st_adata = ground_truth
    runner.config = SimpleNamespace(cell_type_col="major_type")

    svc_obs = runner._get_svc_obs()

    assert svc_obs["true_cell_type"].tolist() == ["A", "B"]
    np.testing.assert_allclose(
        svc_obs[["x", "y"]].to_numpy(),
        [[0.2, 0.3], [1.2, 1.3]],
    )


def test_benchmark_runner_default_level1_accepts_clusters_only_ground_truth(
    load_sr_runner,
):
    module = load_sr_runner("sc_svc_sr_benchmark")
    spatial = AnnData(
        X=np.ones((1, 1), dtype=np.float64),
        obs=pd.DataFrame(index=["spot-1"]),
        var=pd.DataFrame(index=["g1"]),
    )
    spatial.obsm["spatial"] = np.array([[0.0, 0.0]])
    spatial.uns["all_cells_in_spot"] = {"spot-1": ["cell-1"]}
    ground_truth = AnnData(
        X=np.ones((1, 1), dtype=np.float64),
        obs=pd.DataFrame(
            {
                "cell_id": ["cell-1"],
                "clusters": ["A"],
                "x": [0.2],
                "y": [0.3],
            },
            index=["gt-1"],
        ),
        var=pd.DataFrame(index=["g1"]),
    )
    runner = module.ScSVCSr.__new__(module.ScSVCSr)
    runner.st_adata = spatial
    runner.real_st_adata = ground_truth
    runner.config = SimpleNamespace(cell_type_col="Level1")

    svc_obs = runner._get_svc_obs()

    assert svc_obs["true_cell_type"].tolist() == ["A"]
    assert runner.ground_truth_label_source == "clusters"


def test_benchmark_runner_custom_ground_truth_key_is_strict(
    load_sr_runner,
):
    module = load_sr_runner("sc_svc_sr_benchmark")
    spatial = AnnData(
        X=np.ones((1, 1), dtype=np.float64),
        obs=pd.DataFrame(index=["spot-1"]),
        var=pd.DataFrame(index=["g1"]),
    )
    spatial.obsm["spatial"] = np.array([[0.0, 0.0]])
    spatial.uns["all_cells_in_spot"] = {"spot-1": ["cell-1"]}
    ground_truth = AnnData(
        X=np.ones((1, 1), dtype=np.float64),
        obs=pd.DataFrame(
            {
                "cell_id": ["cell-1"],
                "clusters": ["A"],
                "x": [0.2],
                "y": [0.3],
            },
            index=["gt-1"],
        ),
        var=pd.DataFrame(index=["g1"]),
    )
    runner = module.ScSVCSr.__new__(module.ScSVCSr)
    runner.st_adata = spatial
    runner.real_st_adata = ground_truth
    runner.config = SimpleNamespace(cell_type_col="major_type")

    with pytest.raises(KeyError, match="major_type"):
        runner._get_svc_obs()


@pytest.mark.parametrize(
    "strategy_name",
    ["ScSvcSrApplicationStrategy", "ScSvcSrBenchmarkStrategy"],
)
def test_sr_public_strategy_wires_allocation_and_guidance_callbacks(
    strategy_name,
    monkeypatch,
):
    scanpy = types.ModuleType("scanpy")
    scanpy.pp = SimpleNamespace()
    scanpy.pl = SimpleNamespace()
    scanpy.tl = SimpleNamespace()
    monkeypatch.setitem(sys.modules, "scanpy", scanpy)
    adapters = importlib.import_module("revise.backend.adapters")
    records = []
    guidance = []
    guidance_callback = guidance.append
    allocation_callback = records.append
    config = SimpleNamespace()

    class Runner:
        def local_refinement(self):
            assert config.assignment_guidance_callback is guidance_callback
            config.sr_allocation_callback({"status": "completed"})

    ctx = SimpleNamespace(
        runner_config=config,
        runner=Runner(),
        route_key="sc_svc_sr:spot_size",
        assignment_guidance_callback=guidance_callback,
        record_sr_allocation=allocation_callback,
    )

    getattr(adapters, strategy_name)().solve_ot(ctx)

    assert records == [{"status": "completed"}]


def _install_public_sr_route_stubs(monkeypatch):
    scanpy = types.ModuleType("scanpy")
    scanpy.pp = SimpleNamespace()
    scanpy.pl = SimpleNamespace()
    scanpy.tl = SimpleNamespace()
    scanpy.AnnData = AnnData
    monkeypatch.setitem(sys.modules, "scanpy", scanpy)
    adapters = importlib.import_module("revise.backend.adapters")
    backend_package = importlib.import_module("revise.backend")
    monkeypatch.setattr(
        backend_package,
        "adapters",
        adapters,
        raising=False,
    )
    from revise.backend.ops.sr_allocation import (
        record_virtual_cell_not_applicable,
    )
    from revise.backend.policies import ModeEvaluationPolicy, ModeValidationPolicy
    from revise.backend.registry import StrategyRegistry
    from revise.recon.pipeline import UnifiedReconstructionPipeline

    spatial = AnnData(
        X=np.ones((1, 2), dtype=np.float64),
        obs=pd.DataFrame(index=["spot-1"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    spatial.obsm["spatial"] = np.array([[0.0, 0.0]])
    spatial.uns["all_cells_in_spot"] = {"spot-1": ["cell-1"]}
    reference = AnnData(
        X=np.ones((1, 2), dtype=np.float64),
        obs=pd.DataFrame(
            {"Patient": ["sample"], "Level1": ["A"]},
            index=["ref-1"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    ground_truth = AnnData(
        X=np.ones((1, 2), dtype=np.float64),
        obs=pd.DataFrame(
            {
                "cell_id": ["cell-1"],
                "Level1": ["A"],
                "x": [0.0],
                "y": [0.0],
            },
            index=["cell-1"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    captured = {"strategies": [], "routes": []}

    monkeypatch.setattr(ModeValidationPolicy, "validate", lambda _self, _ctx: None)
    monkeypatch.setattr(
        ModeEvaluationPolicy,
        "should_evaluate",
        lambda _self, _ctx: False,
    )
    monkeypatch.setattr(
        UnifiedReconstructionPipeline,
        "_persist_outputs",
        lambda _self, _ctx: None,
    )
    monkeypatch.setattr(
        adapters,
        "_input_service",
        lambda _ctx: SimpleNamespace(
            read_st_adata=lambda _path: spatial.copy(),
            read_sc_ref_adata=lambda _path: reference.copy(),
            read_real_adata=lambda _path: ground_truth.copy(),
        ),
    )
    monkeypatch.setattr(
        adapters,
        "ensure_all_cells_in_spot",
        lambda *_args, **_kwargs: None,
    )

    class FakeSrRunner:
        def __init__(self, st_adata, sc_ref_adata, config, *rest):
            self.st_adata = st_adata
            self.sc_ref_adata = sc_ref_adata
            self.config = config
            self.svc = {}

        def local_refinement(self):
            broad_key = str(self.config.cell_type_col)
            self.config.sr_allocation_callback(
                {
                    "status": "completed",
                    "broad_key": broad_key,
                    "n_spots": 1,
                    "n_virtual_cells": 1,
                    "allocation_method": "posterior_reference_allocation",
                }
            )
            record_virtual_cell_not_applicable(
                self.config,
                problem_key=(
                    f"public-route:{self.config.assignment_guidance_route}"
                ),
                reason="insufficient_virtual_cells",
            )
            self.svc["sc_svc_dec"] = AnnData(
                X=np.ones((1, 2), dtype=np.float64),
                obs=pd.DataFrame(index=["cell-1"]),
                var=pd.DataFrame(index=["g1", "g2"]),
            )
            captured["routes"].append(
                self.config.assignment_guidance_route
            )

    for name in ("sc_svc_sr_application", "sc_svc_sr_benchmark"):
        module = types.ModuleType(f"revise.backend.runners.{name}")
        module.ScSVCSr = FakeSrRunner
        monkeypatch.setitem(sys.modules, module.__name__, module)
    monkeypatch.setattr(
        adapters.ScSvcSrApplicationStrategy,
        "global_anchoring",
        lambda self, ctx: self._attach_assignment_guidance_callback(ctx),
    )
    monkeypatch.setattr(
        adapters.ScSvcSrBenchmarkStrategy,
        "global_anchoring",
        lambda self, ctx: self._attach_assignment_guidance_callback(ctx),
    )
    real_get = StrategyRegistry.get

    def registry_get(registry, strategy_id):
        captured["strategies"].append(strategy_id)
        return real_get(registry, strategy_id)

    monkeypatch.setattr(StrategyRegistry, "get", registry_get)
    return captured


@pytest.mark.parametrize(
    ("surface", "profile", "expected_route"),
    [
        ("application", "application_sc_sr", "sc_svc_sr:spot_size"),
        ("batch", "benchmark_sr_batch", "sim2real:batch_effect"),
        ("spot_size", "benchmark_sr_spot_size", "sim2real:spot_size"),
    ],
)
def test_public_sr_routes_persist_allocation_and_guidance_manifest(
    monkeypatch,
    tmp_path,
    surface,
    profile,
    expected_route,
):
    from revise.application import service as application_service
    from revise.benchmark import cli as benchmark_cli

    captured = _install_public_sr_route_stubs(monkeypatch)
    output_root = tmp_path / f"output-{surface}"
    data_root = tmp_path / "data"
    data_root.mkdir()

    if surface == "application":
        args = SimpleNamespace(
            svc_type="sc-SVC-sr",
            config="revise/revise.yaml",
            seed=17,
            data_root=str(data_root),
            output_root=str(output_root),
            sample_name="sample",
            st_file="st.h5ad",
            sc_ref_file="sc.h5ad",
            patient_key="Patient",
            ot_method=None,
            cell_type_col=None,
            sub_cell_type_col=None,
            select_ct="all",
        )
        resolved_profile, _, _ = application_service._run_pipeline(args)
        assert resolved_profile == profile
    else:
        monkeypatch.setattr(benchmark_cli, "BATCH_NUMS", [1])
        monkeypatch.setattr(benchmark_cli, "SPOT_SIZES", [50])
        monkeypatch.setattr(
            benchmark_cli,
            "_discover_spot_sizes",
            lambda **_kwargs: [50],
        )
        args = SimpleNamespace(
            platform="sim2real",
            confounding="batch_effect" if surface == "batch" else "spot_size",
            data_root=str(data_root),
            dataset_task=None,
            sample_name="sample",
            st_file=None,
            gt_svc_file=None,
            sc_ref_file=None,
            output_root=str(output_root),
            sample_size=None,
            config="revise/revise.yaml",
            seed=17,
            seed_scope="run",
            local_refinement_guidance=None,
            local_refinement_compatibility_mode=None,
            posterior_mode=None,
            posterior_key=None,
            posterior_beta=None,
            posterior_min_affinity=None,
            posterior_cost_strength=None,
            posterior_strict=False,
            sr_refinement_preset=None,
        )
        benchmark_cli.main(args)

    paths = list(output_root.rglob("provenance.json"))
    assert len(paths) == 1
    manifest = json.loads(paths[0].read_text(encoding="utf-8"))
    assert manifest["run"]["status"] == "succeeded"
    assert manifest["profile"] == profile
    assert manifest["route_key"] == expected_route
    assert captured["strategies"] == [
        (
            "ScSvcSrApplicationStrategy"
            if surface == "application"
            else "ScSvcSrBenchmarkStrategy"
        )
    ]
    assert captured["routes"] == [expected_route]
    assert manifest["sr_allocation"] == [
        {
            "status": "completed",
            "broad_key": "Level1",
            "n_spots": 1,
            "n_virtual_cells": 1,
            "allocation_method": "posterior_reference_allocation",
        }
    ]
    [event] = manifest["assignment_guidance"]["events"]
    assert event["route"] == expected_route
    assert event["operator"] == "virtual_cell_ot"
    assert event["outcome"] == "not_applicable"


def _large_sr_runner_inputs():
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
        for cell_index in range(50):
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


@pytest.mark.parametrize("solver", ["pot", "tacco"])
def test_application_cost_guidance_reaches_actual_runner_solver_seam(
    solver,
    load_sr_runner,
    monkeypatch,
):
    from scipy import sparse

    module = load_sr_runner("sc_svc_sr_application")
    spatial, reference, svc_obs = _large_sr_runner_inputs()
    config = _runner_config(guidance="prefer", compatibility_mode="cost")
    config.rec_ot_method = solver
    collector = AssignmentGuidanceCollector()
    config.assignment_guidance_callback = collector.callback
    config.assignment_guidance_route = "sc_svc_sr:spot_size"
    calls = []

    monkeypatch.setattr(
        module,
        "get_adjacency_graph",
        lambda adata, **_kwargs: sparse.csr_matrix(
            np.ones((adata.n_obs, adata.n_obs), dtype=np.float64)
        ),
    )

    def solve(source, target, cost, **kwargs):
        calls.append(
            {
                "method": kwargs["method"],
                "cost": np.asarray(cost).copy(),
                "reference_measure": kwargs["reference_measure"],
            }
        )
        return np.outer(
            np.asarray(source) / np.sum(source),
            np.asarray(target) / np.sum(target),
        )

    monkeypatch.setattr(module, "solve_local_ot", solve)
    runner = module.ScSVCSr.__new__(module.ScSVCSr)
    runner.st_adata = spatial
    runner.sc_ref_adata = reference
    runner.config = config
    runner.svc_obs = svc_obs
    runner.spot_sr = SpotSrKernel(
        config,
        logging.getLogger(f"test-app-cost-{solver}"),
    )
    runner.graph_aggregate = SimpleNamespace(
        run=lambda *, adata, **_kwargs: adata
    )
    runner.logger = logging.getLogger(f"test-app-cost-{solver}")
    runner.svc = {}

    runner.local_refinement()

    assert len(calls) == 2
    assert {call["method"] for call in calls} == {solver}
    assert all(call["reference_measure"] is None for call in calls)
    assert any(
        np.any(call["cost"][np.isfinite(call["cost"])] > 1.0)
        for call in calls
    ), [event.get("reason") for event in collector.events]


def test_benchmark_reference_guidance_reaches_pot_reference_measure(
    load_sr_runner,
    monkeypatch,
):
    from scipy import sparse

    module = load_sr_runner("sc_svc_sr_benchmark")
    spatial, reference, svc_obs = _large_sr_runner_inputs()
    config = _runner_config(
        guidance="prefer",
        graph_enabled=True,
        compatibility_mode="reference",
    )
    config.rec_ot_method = "pot"
    collector = AssignmentGuidanceCollector()
    config.assignment_guidance_callback = collector.callback
    config.assignment_guidance_route = "sim2real:batch_effect"
    calls = []

    monkeypatch.setattr(
        module,
        "get_adjacency_graph",
        lambda adata, **_kwargs: sparse.csr_matrix(
            np.ones((adata.n_obs, adata.n_obs), dtype=np.float64)
        ),
    )

    def solve(source, target, cost, **kwargs):
        calls.append(
            {
                "method": kwargs["method"],
                "cost": np.asarray(cost).copy(),
                "reference_measure": np.asarray(
                    kwargs["reference_measure"]
                ).copy(),
            }
        )
        return np.outer(
            np.asarray(source) / np.sum(source),
            np.asarray(target) / np.sum(target),
        )

    monkeypatch.setattr(module, "solve_local_ot", solve)
    runner = module.ScSVCSr.__new__(module.ScSVCSr)
    runner.st_adata = spatial
    runner.sc_ref_adata = reference
    runner.config = config
    runner.svc_obs = svc_obs
    runner.spot_sr = SpotSrKernel(
        config,
        logging.getLogger("test-benchmark-reference-pot"),
    )
    runner.graph_aggregate = SimpleNamespace(
        run=lambda *, adata, **_kwargs: adata
    )
    runner.logger = logging.getLogger("test-benchmark-reference-pot")
    runner.svc = {}
    runner._graphagg_confidence_cache = None
    runner._graphagg_confidence_source = None
    runner._graphagg_alpha_weight_cache = None
    runner._graphagg_posterior_source = None

    runner.local_refinement()

    assert len(calls) == 2
    assert {call["method"] for call in calls} == {"pot"}
    assert all(
        call["reference_measure"].ndim == 2 for call in calls
    ), [event.get("reason") for event in collector.events]
    assert all(np.all(call["reference_measure"] >= 0) for call in calls)
    assert all(
        np.max(call["cost"][np.isfinite(call["cost"])]) <= 1.0
        for call in calls
    )


@pytest.mark.parametrize(
    ("mode", "solver", "message"),
    [
        ("application", "pot", "application local refinement"),
        ("benchmark", "tacco", "TACCO local refinement"),
    ],
)
def test_sr_reference_capability_rejects_before_strategy_allocation(
    tmp_path,
    mode,
    solver,
    message,
):
    from revise.backend.policies import (
        ModeEvaluationPolicy,
        ModeValidationPolicy,
    )
    from revise.recon.pipeline import UnifiedReconstructionPipeline

    allocation_attempts = []

    class Strategy:
        def prepare_context(self, _ctx):
            allocation_attempts.append("prepared")

    io = {
        "data_root": str(tmp_path),
        "output_root": str(tmp_path / "output"),
        "sample_name": "sample",
        "st_file": "st.h5ad",
        "sc_ref_file": "sc.h5ad",
    }
    if mode == "benchmark":
        io.update(gt_svc_file="gt.h5ad", spot_size=50)
    ctx = SimpleNamespace(
        runtime={"mode": mode, "task": "sc_svc_sr"},
        io=io,
        merged_config={
            "local_refinement": {
                "guidance": "prefer",
                "compatibility": {"mode": "reference"},
            },
            "ot": {
                "ga": {"solver": "pot"},
                "lr": {"solver": solver},
            },
            "sc": {"sr_graph_agg_enabled": True},
        },
    )
    pipeline = UnifiedReconstructionPipeline(
        strategy=Strategy(),
        validation_policy=ModeValidationPolicy(),
        evaluation_policy=ModeEvaluationPolicy(),
    )

    with pytest.raises(ValueError, match=message):
        pipeline.validate_inputs(ctx)

    assert allocation_attempts == []


def test_graph_enabled_benchmark_mandatory_snapshot_is_policy_invariant(
    load_sr_runner,
    monkeypatch,
):
    from scipy import sparse

    module = load_sr_runner("sc_svc_sr_benchmark")
    monkeypatch.setattr(
        module,
        "get_adjacency_graph",
        lambda adata, **_kwargs: sparse.csr_matrix(
            np.ones((adata.n_obs, adata.n_obs), dtype=np.float64)
        ),
    )
    active_mode = {"value": None}
    solver_calls = {"off": 0, "prefer": 0, "require": 0}

    def solve(source, target, _cost, **_kwargs):
        solver_calls[active_mode["value"]] += 1
        return np.outer(
            np.asarray(source) / np.sum(source),
            np.asarray(target) / np.sum(target),
        )

    monkeypatch.setattr(module, "solve_local_ot", solve)
    snapshots = {}

    for guidance in ("off", "prefer", "require"):
        active_mode["value"] = guidance
        spatial, reference, svc_obs = _large_sr_runner_inputs()
        config = _runner_config(
            guidance=guidance,
            graph_enabled=True,
            compatibility_mode="cost",
            seed=23,
        )
        collector = AssignmentGuidanceCollector()
        allocation_records = []
        config.assignment_guidance_callback = collector.callback
        config.assignment_guidance_route = "sim2real:spot_size"
        config.sr_allocation_callback = allocation_records.append
        runner = module.ScSVCSr.__new__(module.ScSVCSr)
        runner.st_adata = spatial
        runner.sc_ref_adata = reference
        runner.config = config
        runner.svc_obs = svc_obs
        runner.spot_sr = SpotSrKernel(
            config,
            logging.getLogger(f"test-policy-invariant-{guidance}"),
        )
        original_distribution = (
            runner.spot_sr.get_spot_cell_distribution
        )
        quota_holder = {}

        def capture_distribution(*args, **kwargs):
            quota = original_distribution(*args, **kwargs)
            quota_holder["value"] = quota.copy()
            return quota

        runner.spot_sr.get_spot_cell_distribution = capture_distribution
        runner.graph_aggregate = SimpleNamespace(
            run=lambda *, adata, **_kwargs: adata
        )
        runner.logger = logging.getLogger(
            f"test-policy-invariant-{guidance}"
        )
        runner.svc = {}
        runner._graphagg_confidence_cache = None
        runner._graphagg_confidence_source = None
        runner._graphagg_alpha_weight_cache = None
        runner._graphagg_posterior_source = None
        original_graph = runner._apply_graph_aggregation
        mandatory_holder = {}

        def capture_graph(self, values, **kwargs):
            mandatory_holder["value"] = np.asarray(values).copy()
            return original_graph(values, **kwargs)

        runner._apply_graph_aggregation = types.MethodType(
            capture_graph,
            runner,
        )

        runner.local_refinement()

        snapshots[guidance] = {
            "quota": quota_holder["value"],
            "rows": runner.svc_obs[
                ["cell_id", "spot_name", "cell_type"]
            ].reset_index(drop=True),
            "mandatory": mandatory_holder["value"],
            "raw_output": np.asarray(
                runner.svc["sc_svc_dec"].X
            ).copy(),
        }
        assert allocation_records[0]["status"] == "completed"
        assert len(collector.events) == 2
        assert {
            event["outcome"] for event in collector.events
        } == ({"off"} if guidance == "off" else {"applied"})

    baseline = snapshots["off"]
    for guidance in ("prefer", "require"):
        pd.testing.assert_frame_equal(
            snapshots[guidance]["quota"],
            baseline["quota"],
        )
        pd.testing.assert_frame_equal(
            snapshots[guidance]["rows"],
            baseline["rows"],
        )
        np.testing.assert_allclose(
            snapshots[guidance]["mandatory"],
            baseline["mandatory"],
        )
        np.testing.assert_allclose(
            snapshots[guidance]["raw_output"],
            baseline["raw_output"],
        )
    assert solver_calls == {"off": 2, "prefer": 2, "require": 2}
