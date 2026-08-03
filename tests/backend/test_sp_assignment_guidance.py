from __future__ import annotations

import importlib
import logging
import sys
import types
from contextlib import contextmanager
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, concat as anndata_concat
from scipy import sparse

from revise.backend.ops.assignment import GlobalAssignmentContractError


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
    "revise.backend.runners.sp_svc_assignment",
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


def _config(*, solver="pot", strength=0.2, callback=None):
    return SimpleNamespace(
        plot_flag=False,
        cell_type_col="major_type",
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
        posterior_conditioning_cost_strength=strength,
        assignment_guidance_callback=callback,
        ot_event_callback=None,
    )


def _posterior(obs_names, probabilities):
    return pd.DataFrame(
        probabilities,
        index=obs_names,
        columns=["A", "B"],
    )


def _application_runner(
    module,
    *,
    n_obs=51,
    probabilities=None,
    solver="pot",
    strength=0.2,
    callback=None,
):
    obs_names = [f"spot-{index}" for index in range(n_obs)]
    st = AnnData(
        X=sparse.csr_matrix(np.ones((n_obs, 2))),
        obs=pd.DataFrame({"major_type": ["A"] * n_obs}, index=obs_names),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    if probabilities is not None:
        st.obsm["major_type"] = _posterior(obs_names, probabilities)
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
    runner.config = _config(
        solver=solver,
        strength=strength,
        callback=callback,
    )
    runner.logger = logging.getLogger("test-sp-application-strict-ga")
    runner.graph_aggregate = SimpleNamespace(
        run=lambda *, adata, **_kwargs: adata
    )
    runner.svc = {}
    return runner


def _patch_application_problem(module, monkeypatch, captured, *, n_slots=1):
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
    neighbor_idx = np.column_stack(
        [np.roll(np.arange(51), -offset) for offset in range(1, n_slots + 1)]
    ).astype(np.int32)
    monkeypatch.setattr(
        module,
        "_compute_topk_expression",
        lambda **_kwargs: (
            np.ones((51, n_slots)),
            np.full(n_slots, 102.0 / n_slots),
            neighbor_idx,
            np.ones((51, n_slots), dtype=bool),
            0,
        ),
    )
    monkeypatch.setattr(
        module,
        "stabilize_local_ot_support",
        lambda *_args, **_kwargs: (
            np.arange(n_slots),
            np.arange(51),
            np.ones((n_slots, 51), dtype=bool),
        ),
    )
    monkeypatch.setattr(
        module,
        "similarity_to_distance",
        lambda similarities, _mask: np.zeros_like(similarities),
    )

    def solve(nu, mu, cost, **kwargs):
        captured.update(
            nu=np.asarray(nu).copy(),
            mu=np.asarray(mu).copy(),
            cost=np.asarray(cost).copy(),
            kwargs=kwargs.copy(),
        )
        return np.ones((n_slots, 51), dtype=np.float64)

    monkeypatch.setattr(module, "solve_local_ot", solve)


def test_strict_sp_assignment_requires_soft_q_without_one_hot_fallback(
    sp_modules,
):
    application, _benchmark = sp_modules
    runner = _application_runner(application)
    assignment = importlib.import_module(
        "revise.backend.runners.sp_svc_assignment"
    )

    with pytest.raises(
        GlobalAssignmentContractError,
        match=r"obsm\[major_type\]",
    ):
        assignment.global_assignment_from_adata(
            runner.st_adata,
            key="major_type",
            expected_categories=pd.Index(["A", "B"]),
        )


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("permuted", "observation axis"),
        ("negative", "non-negative"),
        ("argmax", "argmax"),
    ],
)
def test_strict_sp_assignment_rejects_invalid_q(
    sp_modules,
    mutation,
    message,
):
    application, _benchmark = sp_modules
    runner = _application_runner(
        application,
        probabilities=[[0.9, 0.1]] * 51,
    )
    if mutation == "permuted":
        runner.st_adata.obsm["major_type"].index = (
            runner.st_adata.obs_names[::-1]
        )
    elif mutation == "negative":
        runner.st_adata.obsm["major_type"].iloc[0] = [-0.1, 1.1]
    else:
        runner.st_adata.obs["major_type"] = "B"
    assignment = importlib.import_module(
        "revise.backend.runners.sp_svc_assignment"
    )

    with pytest.raises(GlobalAssignmentContractError, match=message):
        assignment.global_assignment_from_adata(
            runner.st_adata,
            key="major_type",
            expected_categories=pd.Index(["A", "B"]),
        )


@pytest.mark.parametrize(
    ("columns", "probabilities"),
    [
        (["A"], [1.0]),
        (["A", "B", "C"], [0.8, 0.1, 0.1]),
        (["B", "A"], [0.1, 0.9]),
    ],
)
def test_application_rejects_q_category_axis_that_differs_from_reference(
    sp_modules,
    monkeypatch,
    columns,
    probabilities,
):
    application, _benchmark = sp_modules
    runner = _application_runner(
        application,
        probabilities=[[0.9, 0.1]] * 51,
    )
    runner.st_adata.obsm["major_type"] = pd.DataFrame(
        [probabilities] * 51,
        index=runner.st_adata.obs_names,
        columns=columns,
    )
    monkeypatch.setattr(
        application,
        "trim_sp_adata",
        lambda adata, *_args, **_kwargs: (adata.copy(), {}),
    )
    monkeypatch.setattr(
        application,
        "solve_local_ot",
        lambda *_args, **_kwargs: pytest.fail(
            "category mismatch must fail before solve"
        ),
    )

    with pytest.raises(GlobalAssignmentContractError, match="category axis"):
        runner.local_refinement()


@pytest.mark.parametrize("solver", ["pot", "tacco"])
def test_application_conditions_cost_for_cost_capable_solver(
    sp_modules,
    monkeypatch,
    solver,
):
    application, _benchmark = sp_modules
    runner = _application_runner(
        application,
        probabilities=[[0.9, 0.1]] * 51,
        solver=solver,
        callback=lambda *_args, **_kwargs: pytest.fail(
            "strict sp route must not call legacy guidance callback"
        ),
    )
    captured = {}
    _patch_application_problem(application, monkeypatch, captured)

    assert runner.local_refinement() is True
    assert np.all(captured["cost"] > 0.0)
    assert captured["kwargs"]["method"] == solver
    assert captured["kwargs"]["reference_measure"] is None


def test_application_validates_full_ga_once_before_group_conditioning(
    sp_modules,
    monkeypatch,
):
    application, _benchmark = sp_modules
    assignment_module = importlib.import_module(
        "revise.backend.runners.sp_svc_assignment"
    )
    real_validate = assignment_module.validate_global_assignment
    validation_calls = []

    def count_validation(*args, **kwargs):
        validation_calls.append(args[0])
        return real_validate(*args, **kwargs)

    monkeypatch.setattr(
        assignment_module,
        "validate_global_assignment",
        count_validation,
    )
    runner = _application_runner(
        application,
        probabilities=[[0.9, 0.1]] * 51,
    )
    captured = {}
    _patch_application_problem(application, monkeypatch, captured)

    runner.local_refinement()

    assert len(validation_calls) == 1


def test_same_hard_labels_with_different_soft_q_change_solver_coupling(
    sp_modules,
    monkeypatch,
):
    application, _benchmark = sp_modules
    solutions = []
    posterior_variants = (
        [[0.95, 0.05] if i % 2 == 0 else [0.55, 0.45] for i in range(51)],
        [[0.95, 0.05] if i % 3 == 0 else [0.55, 0.45] for i in range(51)],
    )
    for probabilities in posterior_variants:
        runner = _application_runner(
            application,
            probabilities=probabilities,
        )
        runner.config.rec_graph_n_neighbors = 2
        captured = {}
        _patch_application_problem(
            application,
            monkeypatch,
            captured,
            n_slots=2,
        )

        def solve(nu, mu, cost, **_kwargs):
            nu = np.asarray(nu, dtype=np.float64)
            mu = np.asarray(mu, dtype=np.float64)
            kernel = np.exp(-np.asarray(cost, dtype=np.float64))
            column_scale = np.ones_like(mu)
            for _ in range(1000):
                row_scale = nu / (kernel @ column_scale)
                column_scale = mu / (kernel.T @ row_scale)
            coupling = (
                row_scale[:, None] * kernel * column_scale[None, :]
            )
            solutions.append(
                (coupling.copy(), nu, mu)
            )
            return coupling

        monkeypatch.setattr(application, "solve_local_ot", solve)
        runner.local_refinement()

    for coupling, nu, mu in solutions:
        np.testing.assert_allclose(
            coupling.sum(axis=1), nu, rtol=1e-6, atol=1e-8
        )
        np.testing.assert_allclose(
            coupling.sum(axis=0), mu, rtol=1e-6, atol=1e-8
        )
    assert not np.allclose(
        solutions[0][0], solutions[1][0], rtol=1e-6, atol=1e-8
    )


def test_application_zero_strength_keeps_baseline_solver_inputs(
    sp_modules,
    monkeypatch,
):
    application, _benchmark = sp_modules
    runner = _application_runner(
        application,
        probabilities=[[0.9, 0.1]] * 51,
        strength=0.0,
    )
    captured = {}
    _patch_application_problem(application, monkeypatch, captured)

    assert runner.local_refinement() is False
    np.testing.assert_allclose(
        captured["cost"],
        np.zeros((1, 51)),
        rtol=1e-6,
        atol=1e-8,
    )
    assert captured["kwargs"]["reference_measure"] is None


def test_application_all_skipped_reports_not_applied(
    sp_modules,
    monkeypatch,
):
    application, _benchmark = sp_modules
    runner = _application_runner(
        application,
        n_obs=2,
        probabilities=[[0.9, 0.1]] * 2,
        callback=lambda *_args, **_kwargs: pytest.fail(
            "strict sp route must not call legacy guidance callback"
        ),
    )
    monkeypatch.setattr(
        application,
        "trim_sp_adata",
        lambda adata, *_args, **_kwargs: (adata.copy(), {}),
    )

    assert runner.local_refinement() is False


def _benchmark_runner(module, *, solver="tacco", strength=0.2, callback=None):
    replace_count = 50
    donor_count = 2
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
    st.obsm["major_type"] = _posterior(
        obs_names,
        [[0.9, 0.1]] * replace_count + [[0.6, 0.4]] * donor_count,
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
        solver=solver,
        strength=strength,
        callback=callback,
    )
    runner.logger = logging.getLogger("test-sp-benchmark-strict-ga")
    runner.seg_evaluate = SimpleNamespace(run=lambda adata, _logger: adata)
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

    def solve(nu, mu, cost, **kwargs):
        captured.update(
            nu=np.asarray(nu).copy(),
            mu=np.asarray(mu).copy(),
            cost=np.asarray(cost).copy(),
            kwargs=kwargs.copy(),
        )
        return np.ones((1, 50), dtype=np.float64)

    monkeypatch.setattr(module, "solve_local_ot", solve)


def test_benchmark_conditions_replace_to_donor_q(sp_modules, monkeypatch):
    _application, benchmark = sp_modules
    runner = _benchmark_runner(
        benchmark,
        callback=lambda *_args, **_kwargs: pytest.fail(
            "strict sp route must not call legacy guidance callback"
        ),
    )
    captured = {}
    _patch_benchmark_problem(benchmark, monkeypatch, captured)

    assert runner.local_refinement() is True
    assert np.all(captured["cost"] > 0.0)
    assert captured["kwargs"]["reference_measure"] is None


def test_benchmark_zero_strength_matches_unconditioned_solver_call(
    sp_modules,
    monkeypatch,
):
    _application, benchmark = sp_modules
    real_condition = benchmark.condition_sp_local_ot_cost

    baseline = {}
    monkeypatch.setattr(
        benchmark,
        "condition_sp_local_ot_cost",
        lambda cost, **_kwargs: cost,
    )
    _patch_benchmark_problem(benchmark, monkeypatch, baseline)
    baseline_runner = _benchmark_runner(benchmark, strength=0.0)
    assert baseline_runner.local_refinement() is False

    zero_strength = {}
    monkeypatch.setattr(
        benchmark,
        "condition_sp_local_ot_cost",
        real_condition,
    )
    _patch_benchmark_problem(benchmark, monkeypatch, zero_strength)
    zero_runner = _benchmark_runner(benchmark, strength=0.0)
    assert zero_runner.local_refinement() is False

    for field in ("nu", "mu", "cost"):
        np.testing.assert_allclose(
            zero_strength[field],
            baseline[field],
            rtol=1e-6,
            atol=1e-8,
        )
    assert zero_strength["kwargs"].keys() == baseline["kwargs"].keys()
    for field, expected in baseline["kwargs"].items():
        actual = zero_strength["kwargs"][field]
        if isinstance(expected, np.ndarray):
            np.testing.assert_array_equal(actual, expected)
        else:
            assert actual == expected
