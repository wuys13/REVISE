from __future__ import annotations

import importlib
import logging
import sys
import types
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, concat as anndata_concat
from scipy import sparse

from revise.backend.ops.assignment import AssignmentStateError
from revise.backend.ops.assignment_guidance import (
    AssignmentGuidanceCollector,
    FallbackReason,
    NotApplicableReason,
)


ISOLATED_MODULE_NAMES = (
    "scanpy",
    "squidpy",
    "numba",
    "revise.backend.kernels",
    "revise.backend.kernels.graph_aggregate",
    "revise.backend.kernels.seg_evaluate",
    "revise.backend.adapters",
    "revise.backend.ops.distance",
    "revise.backend.ops.meta",
    "revise.backend.ops.shaver",
    "revise.backend.ops.topology",
    "revise.backend.runners.application_svc",
    "revise.backend.runners.benchmark_svc",
    "revise.backend.runners.base_svc",
    "revise.backend.runners.base_svc_anchor",
    "revise.backend.runners.sp_svc_assignment_guidance",
    "revise.backend.runners.sp_svc_application",
    "revise.backend.runners.sp_svc_benchmark",
)
RELEVANT_MODULE_PREFIXES = (
    "scanpy",
    "squidpy",
    "numba",
    "revise.backend.adapters",
    "revise.backend.kernels",
    "revise.backend.ops",
    "revise.backend.runners",
)
SETUP_ISOLATION_PREFIXES = (
    "scanpy",
    "squidpy",
    "numba",
    "revise.backend.adapters",
    "revise.backend.kernels",
    "revise.backend.runners",
    "revise.backend.ops.distance",
    "revise.backend.ops.local_ot",
    "revise.backend.ops.meta",
    "revise.backend.ops.shaver",
    "revise.backend.ops.tacco_runtime",
    "revise.backend.ops.topology",
)
_MISSING = object()


def _module_names_for_prefixes(prefixes):
    return tuple(
        name
        for name in sys.modules
        if any(
            name == prefix or name.startswith(f"{prefix}.")
            for prefix in prefixes
        )
    )


def _snapshot_modules(module_names=ISOLATED_MODULE_NAMES):
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


def _remove_modules(module_names) -> None:
    for module_name in module_names:
        sys.modules.pop(module_name, None)
    for module_name in module_names:
        parent_name, separator, attribute = module_name.rpartition(".")
        if not separator:
            continue
        parent = sys.modules.get(parent_name)
        if parent is not None and hasattr(parent, attribute):
            delattr(parent, attribute)


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


@contextmanager
def _isolated_sp_module_environment():
    relevant_snapshot = _snapshot_modules(
        _module_names_for_prefixes(RELEVANT_MODULE_PREFIXES)
    )
    _remove_modules(
        _module_names_for_prefixes(SETUP_ISOLATION_PREFIXES)
    )

    scanpy = types.ModuleType("scanpy")
    scanpy.AnnData = AnnData
    scanpy.concat = anndata_concat
    scanpy.pp = SimpleNamespace()
    scanpy.pl = SimpleNamespace()
    scanpy.tl = SimpleNamespace()
    squidpy = types.ModuleType("squidpy")
    squidpy.gr = SimpleNamespace()
    numba = types.ModuleType("numba")
    numba.njit = lambda function=None, **_kwargs: (
        function if function is not None else lambda target: target
    )
    numba.prange = range
    sys.modules["scanpy"] = scanpy
    sys.modules["squidpy"] = squidpy
    sys.modules["numba"] = numba
    try:
        yield (
            importlib.import_module(
                "revise.backend.runners.sp_svc_application"
            ),
            importlib.import_module(
                "revise.backend.runners.sp_svc_benchmark"
            ),
        )
    finally:
        _remove_modules(
            _module_names_for_prefixes(RELEVANT_MODULE_PREFIXES)
        )
        _restore_modules(relevant_snapshot)


def test_module_isolation_restores_module_and_parent_attribute_together():
    module_name = "revise.backend.adapters"
    backend_package = importlib.import_module("revise.backend")
    original = _snapshot_modules()
    fake_adapters = types.ModuleType(module_name)
    try:
        sys.modules[module_name] = fake_adapters
        setattr(backend_package, "adapters", fake_adapters)
        fake_snapshot = _snapshot_modules()

        _remove_modules(ISOLATED_MODULE_NAMES)

        assert module_name not in sys.modules
        assert not hasattr(backend_package, "adapters")

        _restore_modules(fake_snapshot)

        assert sys.modules[module_name] is fake_adapters
        assert backend_package.adapters is fake_adapters
    finally:
        _restore_modules(original)


def test_sp_module_environment_restores_relevant_prefix_closure():
    prefixes = (
        "scanpy",
        "squidpy",
        "numba",
        "revise.backend.adapters",
        "revise.backend.kernels",
        "revise.backend.ops",
        "revise.backend.runners",
    )

    def relevant_modules():
        return {
            name: module
            for name, module in sys.modules.items()
            if any(
                name == prefix or name.startswith(f"{prefix}.")
                for prefix in prefixes
            )
        }

    before = relevant_modules()
    before_parent_attributes = {
        name: (
            getattr(sys.modules[parent_name], attribute, _MISSING)
            if separator and parent_name in sys.modules
            else _MISSING
        )
        for name in before
        for parent_name, separator, attribute in [name.rpartition(".")]
    }
    with _isolated_sp_module_environment():
        importlib.import_module("revise.backend.kernels.factory")
        importlib.import_module("revise.backend.ops.local_ot")
        importlib.import_module(
            "revise.backend.runners.sp_svc_application"
        )
        inside = relevant_modules()
        assert "revise.backend.kernels.factory" in inside
        assert "revise.backend.ops.local_ot" in inside
        assert "revise.backend.runners.sp_svc_application" in inside

    after = relevant_modules()
    assert after.keys() == before.keys()
    assert all(after[name] is module for name, module in before.items())
    for name, module in before.items():
        parent_name, separator, attribute = name.rpartition(".")
        if separator and parent_name in sys.modules:
            expected = before_parent_attributes[name]
            if expected is _MISSING:
                assert not hasattr(sys.modules[parent_name], attribute)
            else:
                assert getattr(sys.modules[parent_name], attribute) is expected


@pytest.fixture
def sp_modules():
    with _isolated_sp_module_environment() as modules:
        yield modules


def _config(
    collector: AssignmentGuidanceCollector,
    *,
    guidance: str = "prefer",
    compatibility_mode: str = "cost",
    solver: str = "pot",
):
    return SimpleNamespace(
        plot_flag=False,
        cell_type_col="major_type",
        assignment_guidance_policy=guidance,
        rec_graph_method="pca",
        rec_graph_alpha=0.0,
        rec_graph_exp_neighbor_num=1,
        rec_graph_spatial_neighbor_num=1,
        rec_graph_n_neighbors=1,
        rec_ot_method=solver,
        rec_pot_reg=0.1,
        rec_pot_reg_m=0.0,
        rec_pot_reg_type="kl",
        rec_alpha=1.0,
        posterior_conditioning_enabled=guidance != "off",
        posterior_conditioning_mode=compatibility_mode,
        posterior_conditioning_key="major_type",
        posterior_conditioning_strict=guidance == "require",
        posterior_conditioning_beta=1.0,
        posterior_conditioning_min_affinity=0.1,
        posterior_conditioning_cost_strength=2.0,
        assignment_guidance_callback=collector.callback,
        ot_event_callback=None,
    )


@pytest.mark.parametrize(
    ("canonical", "legacy", "expected"),
    [("require", "off", "require"), ("off", "require", "off")],
)
def test_sp_route_prefers_canonical_guidance_over_conflicting_legacy_fields(
    sp_modules,
    canonical,
    legacy,
    expected,
):
    _application, _benchmark = sp_modules
    collector = AssignmentGuidanceCollector()
    config = _config(collector, guidance=legacy)
    config.assignment_guidance_policy = canonical
    guidance = importlib.import_module(
        "revise.backend.runners.sp_svc_assignment_guidance"
    )

    assert guidance.guidance_mode(config) == expected


def _application_runner(
    module,
    collector,
    *,
    n_obs: int = 51,
    posterior: pd.DataFrame | np.ndarray | None = None,
    guidance: str = "prefer",
):
    obs_names = [f"spot-{index}" for index in range(n_obs)]
    st = AnnData(
        X=sparse.csr_matrix(np.ones((n_obs, 2))),
        obs=pd.DataFrame({"major_type": ["A"] * n_obs}, index=obs_names),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    if posterior is not None:
        st.obsm["major_type"] = posterior
    reference = AnnData(
        X=sparse.csr_matrix(np.ones((2, 2))),
        obs=pd.DataFrame(
            {"major_type": ["A", "B"]},
            index=["ref-a", "ref-b"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    runner = module.SpSVC.__new__(module.SpSVC)
    runner.st_adata = st
    runner.sc_ref_adata = reference
    runner.config = _config(collector, guidance=guidance)
    runner.logger = logging.getLogger("test-sp-application-guidance")
    runner.graph_aggregate = SimpleNamespace(
        run=lambda *, adata, **_kwargs: adata
    )
    runner.svc = {}
    return runner


def _patch_application_problem(module, monkeypatch, captured):
    monkeypatch.setattr(
        module,
        "trim_sp_adata",
        lambda adata, *_args, **_kwargs: (adata.copy(), {}),
    )
    monkeypatch.setattr(
        module,
        "get_adjacency_graph",
        lambda adata, **_kwargs: sparse.eye(adata.n_obs, format="csr"),
    )
    neighbor_idx = np.roll(np.arange(51), -1)[:, None].astype(np.int32)
    monkeypatch.setattr(
        module,
        "_compute_topk_expression",
        lambda **_kwargs: (
            np.ones((51, 1)),
            np.ones(1),
            neighbor_idx,
            np.ones((51, 1), dtype=bool),
            0,
        ),
    )
    monkeypatch.setattr(
        module,
        "stabilize_local_ot_support",
        lambda *_args, **_kwargs: (
            np.array([0]),
            np.arange(51),
            np.ones((1, 51), dtype=bool),
        ),
    )
    monkeypatch.setattr(
        module,
        "similarity_to_distance",
        lambda similarities, _mask: np.zeros_like(similarities),
    )

    def solve(_nu, _mu, cost, **kwargs):
        captured["cost"] = np.asarray(cost).copy()
        captured["reference_measure"] = kwargs["reference_measure"]
        return np.ones((1, 51), dtype=np.float64)

    monkeypatch.setattr(module, "solve_local_ot", solve)


def test_application_soft_level1_applies_cost_guidance_and_records_both_axes(
    sp_modules,
    monkeypatch,
):
    application, _benchmark = sp_modules
    collector = AssignmentGuidanceCollector()
    obs_names = [f"spot-{index}" for index in range(51)]
    posterior = pd.DataFrame(
        [[0.9, 0.1] if index % 2 == 0 else [0.1, 0.9] for index in range(51)],
        index=obs_names,
        columns=["A", "B"],
    )
    runner = _application_runner(application, collector, posterior=posterior)
    captured = {}
    _patch_application_problem(application, monkeypatch, captured)

    runner.local_refinement()

    assert np.all(captured["cost"] > 0.0)
    assert captured["reference_measure"] is None
    [event] = collector.events
    assert event["route"] == "sp_svc_application"
    assert event["operator"] == "neighbor_ot"
    assert event["attempted"] is True
    assert event["outcome"] == "applied"
    assert event["left_assignment"]["level"] == "major_type"
    assert event["right_assignment"]["level"] == "major_type"
    assert event["left_assignment"]["observation_axis"]["count"] == 51
    assert event["right_assignment"]["observation_axis"]["count"] == 51


def test_application_off_does_not_construct_compatibility_and_keeps_base_cost(
    sp_modules,
    monkeypatch,
):
    application, _benchmark = sp_modules
    collector = AssignmentGuidanceCollector()
    runner = _application_runner(application, collector, guidance="off")
    captured = {}
    _patch_application_problem(application, monkeypatch, captured)
    guidance = importlib.import_module(
        "revise.backend.runners.sp_svc_assignment_guidance"
    )
    monkeypatch.setattr(
        guidance,
        "assignment_compatibility",
        lambda *_args, **_kwargs: pytest.fail(
            "off guidance must not construct compatibility"
        ),
    )

    runner.local_refinement()

    np.testing.assert_array_equal(captured["cost"], np.zeros((1, 51)))
    [event] = collector.events
    assert event["availability"] == "not_checked"
    assert event["attempted"] is False
    assert event["outcome"] == "off"
    assert event["reason"] is None


def test_application_argmax_only_is_one_hot_guidance(
    sp_modules,
    monkeypatch,
):
    application, _benchmark = sp_modules
    collector = AssignmentGuidanceCollector()
    runner = _application_runner(application, collector)
    captured = {}
    _patch_application_problem(application, monkeypatch, captured)

    runner.local_refinement()

    [event] = collector.events
    assert event["attempted"] is True
    assert event["outcome"] == "applied"
    assert event["left_assignment"]["value_semantics"] == "one_hot"
    assert event["right_assignment"]["value_semantics"] == "one_hot"


def test_application_invalid_soft_assignment_prefer_falls_back_to_base_cost(
    sp_modules,
    monkeypatch,
):
    application, _benchmark = sp_modules
    collector = AssignmentGuidanceCollector()
    obs_names = [f"spot-{index}" for index in range(51)]
    posterior = pd.DataFrame(
        np.vstack(([-1.0, 2.0], np.ones((50, 2)))),
        index=obs_names,
        columns=["A", "B"],
    )
    runner = _application_runner(application, collector, posterior=posterior)
    captured = {}
    _patch_application_problem(application, monkeypatch, captured)

    runner.local_refinement()

    np.testing.assert_array_equal(captured["cost"], np.zeros((1, 51)))
    [event] = collector.events
    assert event["availability"] == "unavailable"
    assert event["attempted"] is False
    assert event["outcome"] == "fallback"
    assert event["reason"] == FallbackReason.ASSIGNMENT_INVALID.value
    assert event["reason_details"] == {"cause": "values_negative"}


def test_application_invalid_soft_assignment_require_fails_before_solver(
    sp_modules,
    monkeypatch,
):
    application, _benchmark = sp_modules
    collector = AssignmentGuidanceCollector()
    posterior = np.ones((51, 2))
    runner = _application_runner(
        application,
        collector,
        posterior=posterior,
        guidance="require",
    )
    runner.sc_ref_adata = runner.sc_ref_adata[["ref-a"], :].copy()
    captured = {}
    _patch_application_problem(application, monkeypatch, captured)
    monkeypatch.setattr(
        application,
        "solve_local_ot",
        lambda *_args, **_kwargs: pytest.fail(
            "require must fail before the solver"
        ),
    )

    with pytest.raises(AssignmentStateError):
        runner.local_refinement()

    [event] = collector.events
    assert event["availability"] == "unavailable"
    assert event["attempted"] is False
    assert event["outcome"] == "failed"
    assert event["reason"] is None


def test_application_insufficient_observations_is_not_applicable(
    sp_modules,
):
    application, _benchmark = sp_modules
    collector = AssignmentGuidanceCollector()
    runner = _application_runner(application, collector, n_obs=2)
    application.trim_sp_adata = lambda adata, *_args, **_kwargs: (
        adata.copy(),
        {},
    )

    runner.local_refinement()

    [event] = collector.events
    assert event["applicability"] == "not_applicable"
    assert event["availability"] == "not_checked"
    assert event["attempted"] is False
    assert event["outcome"] == "not_applicable"
    assert event["reason"] == NotApplicableReason.INSUFFICIENT_UNITS.value
    assert event["reason_details"] == {
        "unit": "observation",
        "observed": 2,
        "required": 51,
    }


def test_application_empty_active_support_is_not_applicable_before_solver(
    sp_modules,
    monkeypatch,
):
    application, _benchmark = sp_modules
    collector = AssignmentGuidanceCollector()
    runner = _application_runner(application, collector)
    captured = {}
    _patch_application_problem(application, monkeypatch, captured)
    monkeypatch.setattr(
        application,
        "stabilize_local_ot_support",
        lambda *_args, **_kwargs: (
            np.array([], dtype=np.int64),
            np.array([], dtype=np.int64),
            np.zeros((0, 0), dtype=bool),
        ),
    )
    monkeypatch.setattr(
        application,
        "solve_local_ot",
        lambda *_args, **_kwargs: pytest.fail(
            "empty support must not call the solver"
        ),
    )

    runner.local_refinement()

    [event] = collector.events
    assert event["applicability"] == "not_applicable"
    assert event["availability"] == "not_checked"
    assert event["attempted"] is False
    assert event["outcome"] == "not_applicable"
    assert event["reason"] == NotApplicableReason.EMPTY_SUPPORT.value
    assert event["reason_details"] == {"support": "active"}


def test_application_post_solver_update_failure_is_terminal_and_reraised(
    sp_modules,
    monkeypatch,
):
    application, _benchmark = sp_modules
    collector = AssignmentGuidanceCollector()
    runner = _application_runner(application, collector)
    captured = {}
    _patch_application_problem(application, monkeypatch, captured)
    runner.graph_aggregate = SimpleNamespace(
        run=lambda **_kwargs: (_ for _ in ()).throw(
            RuntimeError("application update exploded")
        )
    )

    with pytest.raises(RuntimeError, match="application update exploded"):
        runner.local_refinement()

    assert "cost" in captured
    [event] = collector.events
    assert event["availability"] == "available"
    assert event["attempted"] is True
    assert event["outcome"] == "failed"
    assert event["reason"] is None


def test_application_solver_failure_uses_exception_not_reason(
    sp_modules,
    monkeypatch,
):
    application, _benchmark = sp_modules
    collector = AssignmentGuidanceCollector()
    runner = _application_runner(application, collector)
    captured = {}
    _patch_application_problem(application, monkeypatch, captured)
    monkeypatch.setattr(
        application,
        "solve_local_ot",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            RuntimeError("application solver exploded")
        ),
    )

    with pytest.raises(RuntimeError, match="application solver exploded"):
        runner.local_refinement()

    [event] = collector.events
    assert event["attempted"] is True
    assert event["outcome"] == "failed"
    assert event["reason"] is None


def _benchmark_runner(
    module,
    collector,
    *,
    guidance: str = "prefer",
    compatibility_mode: str = "cost",
    solver: str = "pot",
    donor_count: int = 2,
    route: str | None = None,
):
    replace_count = 50
    total = replace_count + donor_count
    obs_names = [f"cell-{index}" for index in range(total)]
    st = AnnData(
        X=sparse.csr_matrix(np.ones((total, 2))),
        obs=pd.DataFrame(
            {
                "major_type": ["A"] * total,
                "no_effect": [False] * replace_count + [True] * donor_count,
            },
            index=obs_names,
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    st.obsm["major_type"] = pd.DataFrame(
        [[0.9, 0.1] if index < replace_count else [0.1, 0.9]
         for index in range(total)],
        index=obs_names,
        columns=["A", "B"],
    )
    runner = module.SpSVC.__new__(module.SpSVC)
    runner.st_adata = st
    runner.sc_ref_adata = AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame(
            {"major_type": ["A", "B"]},
            index=["ref-a", "ref-b"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    runner.config = _config(
        collector,
        guidance=guidance,
        compatibility_mode=compatibility_mode,
        solver=solver,
    )
    if route is not None:
        runner.config.assignment_guidance_route = route
    runner.logger = logging.getLogger("test-sp-benchmark-guidance")
    runner.seg_evaluate = SimpleNamespace(
        run=lambda adata, _logger: adata
    )
    runner.svc = {}
    return runner


def _patch_benchmark_problem(module, monkeypatch, captured):
    adjacency = sparse.lil_matrix((52, 52), dtype=np.float64)
    adjacency[:50, 50:] = 1.0
    monkeypatch.setattr(
        module,
        "get_adjacency_graph",
        lambda *_args, **_kwargs: adjacency.tocsr(),
    )
    monkeypatch.setattr(
        module,
        "stabilize_local_ot_support",
        lambda *_args, **_kwargs: (
            np.array([0]),
            np.arange(50),
            np.ones((1, 50), dtype=bool),
        ),
    )
    monkeypatch.setattr(
        module,
        "similarity_to_distance",
        lambda similarities, _mask: np.zeros_like(similarities),
    )

    def solve(_nu, _mu, cost, **kwargs):
        captured["cost"] = np.asarray(cost).copy()
        captured["reference_measure"] = kwargs["reference_measure"]
        return np.ones((1, 50), dtype=np.float64)

    monkeypatch.setattr(module, "solve_local_ot", solve)


def test_benchmark_off_keeps_base_cost_and_does_not_read_assignment(
    sp_modules,
    monkeypatch,
):
    _application, benchmark = sp_modules
    collector = AssignmentGuidanceCollector()
    runner = _benchmark_runner(benchmark, collector, guidance="off")
    captured = {}
    _patch_benchmark_problem(benchmark, monkeypatch, captured)
    guidance = importlib.import_module(
        "revise.backend.runners.sp_svc_assignment_guidance"
    )
    monkeypatch.setattr(
        guidance,
        "assignment_state_from_adata",
        lambda *_args, **_kwargs: pytest.fail(
            "off guidance must not read assignment state"
        ),
    )

    runner.local_refinement()

    np.testing.assert_array_equal(captured["cost"], np.zeros((1, 50)))
    assert captured["reference_measure"] is None
    [event] = collector.events
    assert event["availability"] == "not_checked"
    assert event["attempted"] is False
    assert event["outcome"] == "off"
    assert event["reason"] is None


@pytest.mark.parametrize(
    ("mode", "mismatch_side", "axis", "outcome", "solver_called"),
    [
        (
            "prefer",
            "replace",
            "observation",
            "fallback",
            True,
        ),
        (
            "require",
            "donor",
            "category",
            "failed",
            False,
        ),
    ],
)
def test_benchmark_axis_mismatch_obeys_pre_solver_policy(
    sp_modules,
    monkeypatch,
    mode,
    mismatch_side,
    axis,
    outcome,
    solver_called,
):
    _application, benchmark = sp_modules
    collector = AssignmentGuidanceCollector()
    runner = _benchmark_runner(benchmark, collector, guidance=mode)
    captured = {}
    _patch_benchmark_problem(benchmark, monkeypatch, captured)
    guidance = importlib.import_module(
        "revise.backend.runners.sp_svc_assignment_guidance"
    )
    real_loader = guidance.assignment_state_from_adata

    def mismatched_loader(adata, **kwargs):
        if mismatch_side == "replace" and adata.n_obs == 50:
            raise AssignmentStateError("observation_labels_mismatch")
        if mismatch_side == "donor" and adata.n_obs == 2:
            raise AssignmentStateError("category_labels_mismatch")
        return real_loader(adata, **kwargs)

    monkeypatch.setattr(
        guidance,
        "assignment_state_from_adata",
        mismatched_loader,
    )

    if mode == "require":
        with pytest.raises(AssignmentStateError):
            runner.local_refinement()
    else:
        runner.local_refinement()

    assert ("cost" in captured) is solver_called
    if solver_called:
        np.testing.assert_array_equal(captured["cost"], np.zeros((1, 50)))
        assert runner.svc["sp_svc"].n_obs == 52
    [event] = collector.events
    assert event["availability"] == "unavailable"
    assert event["attempted"] is False
    assert event["outcome"] == outcome
    if mode == "prefer":
        assert event["reason"] == FallbackReason.ASSIGNMENT_MISALIGNED.value
        assert event["reason_details"]["axis"] == axis
    else:
        assert event["reason"] is None
        assert event["reason_details"] == {}


def test_benchmark_post_solver_update_failure_is_terminal_and_reraised(
    sp_modules,
    monkeypatch,
):
    _application, benchmark = sp_modules
    collector = AssignmentGuidanceCollector()
    runner = _benchmark_runner(benchmark, collector)
    captured = {}
    _patch_benchmark_problem(benchmark, monkeypatch, captured)
    monkeypatch.setattr(
        benchmark.scipy.sparse,
        "lil_matrix",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            RuntimeError("benchmark update exploded")
        ),
    )

    with pytest.raises(RuntimeError, match="benchmark update exploded"):
        runner.local_refinement()

    assert "cost" in captured
    [event] = collector.events
    assert event["availability"] == "available"
    assert event["attempted"] is True
    assert event["outcome"] == "failed"
    assert event["reason"] is None


@pytest.mark.parametrize(
    "route",
    ["sp_svc:segmentation", "sp_svc:bin2cell"],
)
def test_benchmark_public_routes_use_replace_and_donor_assignment_axes(
    sp_modules,
    monkeypatch,
    route,
):
    _application, benchmark = sp_modules
    collector = AssignmentGuidanceCollector()
    runner = _benchmark_runner(
        benchmark,
        collector,
        compatibility_mode="reference",
        route=route,
    )
    captured = {}
    _patch_benchmark_problem(benchmark, monkeypatch, captured)

    runner.local_refinement()

    np.testing.assert_array_equal(captured["cost"], np.zeros((1, 50)))
    assert captured["reference_measure"] is not None
    assert captured["reference_measure"].shape == (1, 50)
    [event] = collector.events
    assert event["route"] == route
    assert event["operator"] == "replacement_ot"
    assert event["outcome"] == "applied"
    assert event["left_assignment"]["observation_axis"]["count"] == 50
    assert event["right_assignment"]["observation_axis"]["count"] == 2


def test_benchmark_cost_guidance_reaches_tacco_with_no_reference_measure(
    sp_modules,
    monkeypatch,
):
    _application, benchmark = sp_modules
    collector = AssignmentGuidanceCollector()
    runner = _benchmark_runner(
        benchmark,
        collector,
        compatibility_mode="cost",
        solver="tacco",
    )
    captured = {}
    _patch_benchmark_problem(benchmark, monkeypatch, captured)

    runner.local_refinement()

    assert np.all(captured["cost"] > 0.0)
    assert captured["reference_measure"] is None
    [event] = collector.events
    assert event["solver"] == "tacco"
    assert event["outcome"] == "applied"


@pytest.mark.parametrize(
    ("profile", "expected_route"),
    [
        ("benchmark_seg", "sim2real:segmentation"),
        ("benchmark_bin2cell", "sim2real:bin2cell"),
    ],
)
def test_public_benchmark_routes_resolve_default_and_emit_contextual_outcome(
    sp_modules,
    monkeypatch,
    tmp_path,
    profile,
    expected_route,
):
    _application, benchmark = sp_modules
    from revise.backend import adapters
    from revise.config import load_raw_config, merge_unified_config

    config_path = Path(__file__).resolve().parents[2] / "revise" / "revise.yaml"
    merged = merge_unified_config(
        raw_config=load_raw_config(config_path),
        profile=profile,
        runtime_overrides={},
        io_overrides={
            "data_root": str(tmp_path),
            "output_root": str(tmp_path),
            "sample_name": "sample",
        },
        algorithm_overrides={"graph": {"n_neighbors": 1}},
    )
    assert merged["runtime"]["strategy"] == "SpSvcBenchmarkSegStrategy"
    assert merged["local_refinement"] == {
        "guidance": "prefer",
        "compatibility": {
            "mode": "cost",
            "beta": 1.0,
            "min_affinity": 0.05,
            "strength": 0.2,
        },
    }
    route_key = (
        f"{merged['runtime']['platform']}:{merged['runtime']['confounding']}"
    )
    assert route_key == expected_route

    seed_collector = AssignmentGuidanceCollector()
    seed = _benchmark_runner(benchmark, seed_collector)
    st = seed.st_adata.copy()
    st.obs["Level1"] = st.obs["major_type"].copy()
    st.obsm["Level1"] = st.obsm["major_type"].copy()
    reference = seed.sc_ref_adata.copy()
    reference.obs["Level1"] = reference.obs["major_type"].copy()
    ground_truth = st.copy()

    class InputService:
        def read_st_adata(self, _path):
            return st.copy()

        def read_real_adata(self, _path):
            return ground_truth.copy()

        def read_sc_ref_adata(self, _path):
            return reference.copy()

    monkeypatch.setattr(
        adapters,
        "_input_service",
        lambda _ctx: InputService(),
    )
    collector = AssignmentGuidanceCollector()
    ctx = SimpleNamespace(
        merged_config=merged,
        runtime=merged["runtime"],
        io=merged["io"],
        columns=merged["columns"],
        route_key=route_key,
        run_dir=tmp_path,
        input_specs=[],
        compatibility_mode=True,
        assignment_guidance_callback=collector.callback,
        logger=logging.getLogger(f"test-public-{profile}"),
    )
    strategy = adapters.SpSvcBenchmarkSegStrategy()
    strategy.prepare_context(ctx)
    assert ctx.runner_config.posterior_conditioning_enabled is True
    assert ctx.runner_config.posterior_conditioning_strict is False
    assert ctx.runner_config.posterior_conditioning_mode == "cost"

    captured = {}
    _patch_benchmark_problem(benchmark, monkeypatch, captured)
    strategy.solve_ot(ctx)

    assert np.all(captured["cost"] > 0.0)
    [event] = collector.events
    assert event["route"] == expected_route
    assert event["mode"] == "prefer"
    assert event["operator"] == "replacement_ot"
    assert event["outcome"] == "applied"


@pytest.mark.parametrize(
    ("profile", "lr_solver", "expected_error"),
    [
        ("application_sp", "pot", "application"),
        ("benchmark_seg", "tacco", "TACCO"),
    ],
)
def test_public_sp_routes_reject_unsupported_reference_capability(
    profile,
    lr_solver,
    expected_error,
):
    from revise.backend.policies import ModeValidationPolicy
    from revise.config import load_raw_config, merge_unified_config

    config_path = Path(__file__).resolve().parents[2] / "revise" / "revise.yaml"
    merged = merge_unified_config(
        raw_config=load_raw_config(config_path),
        profile=profile,
        runtime_overrides={},
        io_overrides={},
        algorithm_overrides={
            "local_refinement": {
                "guidance": "prefer",
                "compatibility": {"mode": "reference"},
            },
            "ot": {"lr": {"solver": lr_solver}},
        },
    )
    ctx = SimpleNamespace(
        merged_config=merged,
        runtime=merged["runtime"],
    )

    with pytest.raises(ValueError, match=expected_error):
        ModeValidationPolicy._validate_solver_compatibility(ctx)


def test_public_benchmark_pot_reference_capability_is_allowed():
    from revise.backend.policies import ModeValidationPolicy
    from revise.config import load_raw_config, merge_unified_config

    config_path = Path(__file__).resolve().parents[2] / "revise" / "revise.yaml"
    merged = merge_unified_config(
        raw_config=load_raw_config(config_path),
        profile="benchmark_seg",
        runtime_overrides={},
        io_overrides={},
        algorithm_overrides={
            "local_refinement": {
                "guidance": "prefer",
                "compatibility": {"mode": "reference"},
            }
        },
    )

    ModeValidationPolicy._validate_solver_compatibility(
        SimpleNamespace(
            merged_config=merged,
            runtime=merged["runtime"],
        )
    )


def test_benchmark_missing_donors_is_not_applicable(sp_modules):
    _application, benchmark = sp_modules
    collector = AssignmentGuidanceCollector()
    runner = _benchmark_runner(benchmark, collector, donor_count=0)

    runner.local_refinement()

    [event] = collector.events
    assert event["applicability"] == "not_applicable"
    assert event["outcome"] == "not_applicable"
    assert event["reason"] == NotApplicableReason.REFERENCE_UNAVAILABLE.value
    assert event["reason_details"] == {"role": "donor"}
