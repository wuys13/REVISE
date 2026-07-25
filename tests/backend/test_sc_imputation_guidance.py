from __future__ import annotations

import importlib
import json
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

from revise.backend.ops.assignment_guidance import AssignmentGuidanceCollector


_MISSING = object()
ISOLATION_PREFIXES = (
    "scanpy",
    "revise.backend.adapters",
    "revise.backend.kernels",
    "revise.backend.ops.distance",
    "revise.backend.ops.meta",
    "revise.backend.ops.shaver",
    "revise.backend.runners.benchmark_svc",
    "revise.backend.runners.sc_svc_impute_benchmark",
    "revise.benchmark.cli",
    "revise.framework",
)


def _module_names():
    return tuple(
        name
        for name in sys.modules
        if any(
            name == prefix or name.startswith(f"{prefix}.")
            for prefix in ISOLATION_PREFIXES
        )
    )


def _snapshot_modules(names):
    modules = {name: sys.modules.get(name, _MISSING) for name in names}
    parent_attributes = {}
    for name in names:
        parent_name, separator, attribute = name.rpartition(".")
        if separator and parent_name in sys.modules:
            parent_attributes[(parent_name, attribute)] = getattr(
                sys.modules[parent_name],
                attribute,
                _MISSING,
            )
    return modules, parent_attributes


def _remove_modules(names):
    for name in names:
        sys.modules.pop(name, None)
    for name in names:
        parent_name, separator, attribute = name.rpartition(".")
        parent = sys.modules.get(parent_name)
        if separator and parent is not None and hasattr(parent, attribute):
            delattr(parent, attribute)


def _restore_modules(snapshot):
    modules, parent_attributes = snapshot
    for name in modules:
        sys.modules.pop(name, None)
    for name, module in modules.items():
        if module is not _MISSING:
            sys.modules[name] = module
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
def _isolated_impute_module():
    names = _module_names()
    snapshot = _snapshot_modules(names)
    _remove_modules(names)

    scanpy = types.ModuleType("scanpy")
    scanpy.AnnData = AnnData
    scanpy.concat = anndata_concat
    scanpy.pp = SimpleNamespace()
    sys.modules["scanpy"] = scanpy

    kernels = types.ModuleType("revise.backend.kernels")
    kernels.GeneImputeKernel = object
    kernels.GeneUncertaintyKernel = object
    sys.modules[kernels.__name__] = kernels

    benchmark = types.ModuleType("revise.backend.runners.benchmark_svc")
    benchmark.BenchmarkSVC = object
    sys.modules[benchmark.__name__] = benchmark

    distance = types.ModuleType("revise.backend.ops.distance")
    distance.bhattacharyya_distance = lambda *_args, **_kwargs: None
    sys.modules[distance.__name__] = distance

    meta = types.ModuleType("revise.backend.ops.meta")
    meta.get_subcluster = lambda *_args, **_kwargs: None
    meta.merge_subcluster = lambda *_args, **_kwargs: None
    sys.modules[meta.__name__] = meta

    shaver = types.ModuleType("revise.backend.ops.shaver")
    shaver.get_prune_adata = lambda adata: adata
    sys.modules[shaver.__name__] = shaver
    try:
        yield importlib.import_module(
            "revise.backend.runners.sc_svc_impute_benchmark"
        )
    finally:
        _remove_modules(_module_names())
        _restore_modules(snapshot)


@pytest.fixture
def impute_module():
    with _isolated_impute_module() as module:
        yield module


def _config(
    collector,
    *,
    guidance="prefer",
    compatibility_mode="cost",
    solver="pot",
):
    return SimpleNamespace(
        cell_type_col="Level1",
        rec_ot_method=solver,
        rec_impute_pot_reg=0.1,
        rec_impute_pot_reg_m=0.0,
        rec_impute_pot_reg_type="kl",
        rec_merge_subcluster_method="mean",
        rec_impute_prune_flag=False,
        posterior_conditioning_enabled=guidance != "off",
        posterior_conditioning_mode=compatibility_mode,
        posterior_conditioning_key="Level1",
        posterior_conditioning_strict=guidance == "require",
        posterior_conditioning_beta=1.0,
        posterior_conditioning_min_affinity=0.1,
        posterior_conditioning_cost_strength=2.0,
        assignment_guidance_callback=collector.callback,
        assignment_guidance_route="sim2real:gene_panel",
        ot_event_callback=None,
    )


@pytest.mark.parametrize(
    ("canonical", "legacy", "expected"),
    [("require", "off", "require"), ("off", "require", "off")],
)
def test_imputation_route_prefers_canonical_guidance_over_conflicting_legacy_fields(
    impute_module,
    canonical,
    legacy,
    expected,
):
    collector = AssignmentGuidanceCollector()
    config = _config(collector, guidance=legacy)
    config.assignment_guidance_policy = canonical

    assert impute_module._guidance_mode(config) == expected


def _inputs(*, spot_posterior=True, reference_posterior=True):
    spot_names = ["spot-1", "spot-2"]
    reference_names = ["ref-1", "ref-2", "ref-3", "ref-4"]
    spatial = AnnData(
        X=sparse.csr_matrix([[2.0, 1.0], [1.0, 2.0]]),
        obs=pd.DataFrame({"Level1": ["A", "A"]}, index=spot_names),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    reference = AnnData(
        X=sparse.csr_matrix(
            [[2.0, 1.0], [1.0, 2.0], [2.0, 2.0], [3.0, 1.0]]
        ),
        obs=pd.DataFrame(
            {
                "Level1": ["A", "A", "A", "A"],
                "leiden_3": ["s1", "s1", "s2", "s2"],
            },
            index=reference_names,
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    if spot_posterior:
        spatial.obsm["Level1"] = pd.DataFrame(
            [[0.9, 0.1], [0.1, 0.9]],
            index=spot_names,
            columns=["A", "B"],
        )
    if reference_posterior:
        # The reversed category order is intentional: compatibility must align
        # by category label, not by column position.
        reference.obsm["Level1"] = pd.DataFrame(
            [
                [0.2, 0.8],
                [0.2, 0.8],
                [0.8, 0.2],
                [0.8, 0.2],
            ],
            index=reference_names,
            columns=["B", "A"],
        )
    return spatial, reference


def _runner(module, collector, *, guidance="prefer", compatibility_mode="cost"):
    spatial, reference = _inputs()
    runner = module.ScSVCImpute.__new__(module.ScSVCImpute)
    runner.st_adata = spatial
    runner.sc_ref_adata = reference
    runner.config = _config(
        collector,
        guidance=guidance,
        compatibility_mode=compatibility_mode,
    )
    runner.logger = logging.getLogger("test-sc-imputation-guidance")
    runner.svc = {}
    return runner, reference


def _patch_problem(module, monkeypatch, runner, *, solver_error=None, assembly_error=None):
    captured = {"solve_calls": [], "weights": []}

    monkeypatch.setattr(
        module,
        "bhattacharyya_distance",
        lambda profiles, spots: np.zeros((profiles.shape[0], spots.shape[0])),
    )

    def solve(source_mass, target_mass, cost, **kwargs):
        captured["solve_calls"].append(
            {
                "source_mass": np.asarray(source_mass).copy(),
                "target_mass": np.asarray(target_mass).copy(),
                "cost": np.asarray(cost).copy(),
                "reference_measure": kwargs.get("reference_measure"),
                "method": kwargs.get("method"),
            }
        )
        if solver_error is not None:
            raise solver_error
        return np.outer(source_mass, target_mass)

    monkeypatch.setattr(module, "solve_local_ot", solve)

    def merge(adata, *, subcluster, mode):
        labels = list(dict.fromkeys(adata.obs[subcluster].astype(str)))
        return AnnData(
            X=np.ones((len(labels), adata.n_vars)),
            obs=pd.DataFrame(
                {subcluster: labels},
                index=[f"merged-{label}" for label in labels],
            ),
            var=pd.DataFrame(index=adata.var_names.copy()),
        )

    monkeypatch.setattr(module, "merge_subcluster", merge)

    def impute(spatial, reference, *, genes_to_predict, neighbor_weights):
        captured["weights"].append(neighbor_weights.copy())
        if assembly_error is not None:
            raise assembly_error
        return spatial.copy()

    runner.gene_impute = SimpleNamespace(run=impute)
    return captured


def _event(collector):
    assert len(collector.events) == 1
    return collector.events[0]


def test_default_sc_guidance_off_preserves_base_problem_without_reading_assignment(
    impute_module,
    monkeypatch,
):
    collector = AssignmentGuidanceCollector()
    runner, reference = _runner(impute_module, collector, guidance="off")
    runner.st_adata.obsm["Level1"] = np.array([[np.nan], [np.nan]])
    reference.obsm["Level1"] = np.full((4, 1), np.nan)
    captured = _patch_problem(impute_module, monkeypatch, runner)
    monkeypatch.setattr(
        impute_module,
        "assignment_compatibility",
        lambda *_args, **_kwargs: pytest.fail(
            "off guidance must not construct compatibility"
        ),
        raising=False,
    )

    runner.local_impute(reference, "leiden_3", guidance_scope="panel")

    np.testing.assert_allclose(captured["solve_calls"][0]["cost"], 0.0)
    event = _event(collector)
    assert event["outcome"] == "off"
    assert event["availability"] == "not_checked"
    assert event["attempted"] is False
    assert event["left_assignment"] is None
    assert event["right_assignment"] is None


def test_soft_bilateral_assignments_align_reversed_categories_and_apply_cost(
    impute_module,
    monkeypatch,
):
    collector = AssignmentGuidanceCollector()
    runner, reference = _runner(impute_module, collector)
    captured = _patch_problem(impute_module, monkeypatch, runner)

    runner.local_impute(reference, "leiden_3", guidance_scope="panel")

    affinity = np.array([[0.74, 0.26], [0.26, 0.74]])
    conditioned = -2.0 * np.log(affinity)
    expected_solver_cost = conditioned / conditioned.max()
    np.testing.assert_allclose(
        captured["solve_calls"][0]["cost"],
        expected_solver_cost,
    )
    event = _event(collector)
    assert event["outcome"] == "applied"
    assert event["attempted"] is True
    assert event["operator"] == "imputation_ot"
    assert event["left_assignment"]["source"] == "obsm[Level1]"
    assert event["left_assignment"]["level"] == "Level1"
    assert event["left_assignment"]["value_semantics"] == "soft"
    assert event["right_assignment"]["source"] == "aggregate(obsm[Level1])"
    assert event["right_assignment"]["level"] == "leiden_3"
    assert event["right_assignment"]["value_semantics"] == "soft"
    assert [step["operation"] for step in event["right_assignment"]["lineage"]] == [
        "load",
        "aggregate",
    ]


def test_spot_only_broad_label_is_documented_as_one_hot_fallback(
    impute_module,
    monkeypatch,
):
    collector = AssignmentGuidanceCollector()
    runner, reference = _runner(impute_module, collector)
    del runner.st_adata.obsm["Level1"]
    _patch_problem(impute_module, monkeypatch, runner)

    runner.local_impute(reference, "leiden_3", guidance_scope="panel")

    event = _event(collector)
    assert event["outcome"] == "applied"
    assert event["left_assignment"]["source"] == "obs[Level1]"
    assert event["left_assignment"]["value_semantics"] == "one_hot"
    assert event["left_assignment"]["lineage"] == [
        {
            "operation": "fallback",
            "reason": "broad_label_only",
            "container": "obs",
            "key": "Level1",
        }
    ]
    assert event["right_assignment"]["source"] == "aggregate(obsm[Level1])"
    assert [step["operation"] for step in event["right_assignment"]["lineage"]] == [
        "load",
        "aggregate",
    ]


@pytest.mark.parametrize(
    ("guidance", "expected_outcome", "should_solve"),
    [
        ("prefer", "fallback", True),
        ("require", "failed", False),
    ],
)
def test_unlabeled_reference_posterior_is_ambiguous_not_positional(
    impute_module,
    monkeypatch,
    guidance,
    expected_outcome,
    should_solve,
):
    collector = AssignmentGuidanceCollector()
    runner, reference = _runner(impute_module, collector, guidance=guidance)
    reference.obsm["Level1"] = reference.obsm["Level1"].to_numpy()
    captured = _patch_problem(impute_module, monkeypatch, runner)

    if guidance == "require":
        with pytest.raises(ValueError, match="category_axis_mismatch"):
            runner.local_impute(reference, "leiden_3", guidance_scope="panel")
    else:
        runner.local_impute(reference, "leiden_3", guidance_scope="panel")

    assert bool(captured["solve_calls"]) is should_solve
    event = _event(collector)
    assert event["outcome"] == expected_outcome
    assert event["availability"] == "unavailable"
    assert event["reason"] == "category_axis_mismatch"
    assert event["attempted"] is False


@pytest.mark.parametrize("guidance", ["prefer", "require"])
def test_misaligned_categories_fallback_or_fail_before_solver(
    impute_module,
    monkeypatch,
    guidance,
):
    collector = AssignmentGuidanceCollector()
    runner, reference = _runner(impute_module, collector, guidance=guidance)
    reference.obsm["Level1"].columns = ["C", "A"]
    captured = _patch_problem(impute_module, monkeypatch, runner)

    if guidance == "require":
        with pytest.raises(ValueError, match="category_axis_mismatch"):
            runner.local_impute(reference, "leiden_3", guidance_scope="panel")
    else:
        runner.local_impute(reference, "leiden_3", guidance_scope="panel")

    assert bool(captured["solve_calls"]) is (guidance == "prefer")
    event = _event(collector)
    assert event["outcome"] == ("fallback" if guidance == "prefer" else "failed")
    assert event["reason"] == "category_axis_mismatch"


def test_reference_ablation_supplies_pot_reference_measure(
    impute_module,
    monkeypatch,
):
    collector = AssignmentGuidanceCollector()
    runner, reference = _runner(
        impute_module,
        collector,
        compatibility_mode="reference",
    )
    captured = _patch_problem(impute_module, monkeypatch, runner)

    runner.local_impute(reference, "leiden_3", guidance_scope="panel")

    reference_measure = captured["solve_calls"][0]["reference_measure"]
    assert reference_measure is not None
    assert np.asarray(reference_measure).shape == (2, 2)
    assert _event(collector)["outcome"] == "applied"


def test_tacco_cost_guidance_uses_conditioned_cost_without_reference_measure(
    impute_module,
    monkeypatch,
):
    collector = AssignmentGuidanceCollector()
    runner, reference = _runner(impute_module, collector)
    runner.config.rec_ot_method = "tacco"
    captured = _patch_problem(impute_module, monkeypatch, runner)

    runner.local_impute(reference, "leiden_3", guidance_scope="panel")

    call = captured["solve_calls"][0]
    assert call["method"] == "tacco"
    assert call["reference_measure"] is None
    assert np.count_nonzero(call["cost"]) == 4
    assert _event(collector)["outcome"] == "applied"


@pytest.mark.parametrize(
    ("drop_spatial", "drop_reference", "left_semantics", "right_source"),
    [
        (False, True, "soft", "aggregate(obs[Level1])"),
        (True, True, "one_hot", "aggregate(obs[Level1])"),
    ],
)
def test_reference_and_bilateral_broad_labels_use_documented_fallback_lineage(
    impute_module,
    monkeypatch,
    drop_spatial,
    drop_reference,
    left_semantics,
    right_source,
):
    collector = AssignmentGuidanceCollector()
    runner, reference = _runner(impute_module, collector)
    if drop_spatial:
        del runner.st_adata.obsm["Level1"]
    if drop_reference:
        del reference.obsm["Level1"]
    _patch_problem(impute_module, monkeypatch, runner)

    runner.local_impute(reference, "leiden_3", guidance_scope="panel")

    event = _event(collector)
    assert event["outcome"] == "applied"
    assert event["left_assignment"]["value_semantics"] == left_semantics
    assert event["right_assignment"]["source"] == right_source
    assert event["right_assignment"]["lineage"][-2]["operation"] == "fallback"
    assert event["right_assignment"]["lineage"][-1]["operation"] == "aggregate"


def _multi_type_runner(module, collector, spatial, reference, *, guidance="off"):
    runner = module.ScSVCImpute.__new__(module.ScSVCImpute)
    runner.st_adata = spatial
    runner.sc_ref_adata = reference
    runner.config = _config(collector, guidance=guidance)
    runner.logger = logging.getLogger("test-sc-imputation-multi-type")
    runner.svc = {}
    return runner


def test_problem_keys_are_unique_for_slash_and_underscore_cell_types(
    impute_module,
    monkeypatch,
):
    spatial = AnnData(
        X=sparse.csr_matrix(np.ones((4, 2))),
        obs=pd.DataFrame(
            {"Level1": ["A/B", "A/B", "A_B", "A_B"]},
            index=[f"spot-{index}" for index in range(4)],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    reference = AnnData(
        X=sparse.csr_matrix(np.ones((4, 2))),
        obs=pd.DataFrame(
            {
                "Level1": ["A/B", "A/B", "A_B", "A_B"],
                "leiden_3": ["s1", "s1", "s2", "s2"],
            },
            index=[f"ref-{index}" for index in range(4)],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    collector = AssignmentGuidanceCollector()
    runner = _multi_type_runner(
        impute_module,
        collector,
        spatial,
        reference,
    )
    captured = _patch_problem(impute_module, monkeypatch, runner)

    runner.local_impute(reference, "leiden_3", guidance_scope="panel")

    assert len(captured["solve_calls"]) == 2
    assert [event["ordinal"] for event in collector.events] == [1, 2]
    assert [event["outcome"] for event in collector.events] == ["off", "off"]
    problem_keys = [event["problem_key"] for event in collector.events]
    assert len(set(problem_keys)) == 2
    assert problem_keys[0].endswith('3:"A/B"')
    assert problem_keys[1].endswith('3:"A_B"')


def test_reference_only_cell_type_is_not_applicable_and_keeps_event_order(
    impute_module,
    monkeypatch,
):
    spatial = AnnData(
        X=sparse.csr_matrix([[2.0, 1.0], [1.0, 2.0]]),
        obs=pd.DataFrame(
            {"Level1": ["A", "A"]},
            index=["spot-a1", "spot-a2"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    spatial.obsm["Level1"] = pd.DataFrame(
        [[0.9, 0.1], [0.8, 0.2]],
        index=spatial.obs_names,
        columns=["A", "B"],
    )
    reference = AnnData(
        X=sparse.csr_matrix(
            [[2.0, 1.0], [1.0, 2.0], [1.0, 1.0], [2.0, 1.0]]
        ),
        obs=pd.DataFrame(
            {
                "Level1": ["A", "A", "B", "B"],
                "leiden_3": ["a1", "a1", "b1", "b1"],
            },
            index=["ref-a1", "ref-a2", "ref-b1", "ref-b2"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    reference.obsm["Level1"] = pd.DataFrame(
        [[0.9, 0.1], [0.8, 0.2], [0.1, 0.9], [0.2, 0.8]],
        index=reference.obs_names,
        columns=["A", "B"],
    )
    collector = AssignmentGuidanceCollector()
    runner = _multi_type_runner(
        impute_module,
        collector,
        spatial,
        reference,
        guidance="prefer",
    )
    captured = _patch_problem(impute_module, monkeypatch, runner)

    runner.local_impute(reference, "leiden_3", guidance_scope="panel")

    assert len(captured["solve_calls"]) == 1
    assert [event["ordinal"] for event in collector.events] == [1, 2]
    assert [event["outcome"] for event in collector.events] == [
        "applied",
        "not_applicable",
    ]
    assert collector.events[1]["reason"] == "empty_spatial_support"
    assert collector.summary() == "mixed"


@pytest.mark.parametrize(
    ("spatial_values", "reference_values", "categorical", "reason"),
    [
        (
            [[0.0, 0.0], [1.0, 1.0]],
            [[1.0, 1.0], [1.0, 1.0]],
            False,
            "zero_source_mass",
        ),
        (
            [[1.0, 1.0], [1.0, 1.0]],
            [[0.0, 0.0], [0.0, 0.0]],
            False,
            "zero_target_mass",
        ),
    ],
)
def test_zero_observed_marginals_are_not_applicable(
    impute_module,
    monkeypatch,
    spatial_values,
    reference_values,
    categorical,
    reason,
):
    spatial = AnnData(
        X=sparse.csr_matrix(spatial_values),
        obs=pd.DataFrame(
            {"Level1": ["A", "A"]},
            index=["spot-1", "spot-2"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    subclusters = (
        pd.Categorical(["s1", "s1"], categories=["s1", "unused"])
        if categorical
        else ["s1", "s1"]
    )
    reference = AnnData(
        X=sparse.csr_matrix(reference_values),
        obs=pd.DataFrame(
            {"Level1": ["A", "A"], "leiden_3": subclusters},
            index=["ref-1", "ref-2"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    collector = AssignmentGuidanceCollector()
    runner = _multi_type_runner(
        impute_module,
        collector,
        spatial,
        reference,
    )
    captured = _patch_problem(impute_module, monkeypatch, runner)

    result = runner.local_impute(
        reference,
        "leiden_3",
        guidance_scope="panel",
    )

    assert captured["solve_calls"] == []
    assert result.n_obs == 0
    event = _event(collector)
    assert event["outcome"] == "not_applicable"
    assert event["reason"] == reason
    assert event["availability"] == "not_checked"


def test_disjoint_categorical_subclusters_solve_each_observed_cell_type(
    impute_module,
    monkeypatch,
):
    spatial = AnnData(
        X=sparse.csr_matrix(np.ones((4, 2))),
        obs=pd.DataFrame(
            {"Level1": ["A", "A", "B", "B"]},
            index=[f"spot-{index}" for index in range(4)],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    reference = AnnData(
        X=sparse.csr_matrix(np.ones((4, 2))),
        obs=pd.DataFrame(
            {
                "Level1": ["A", "A", "B", "B"],
                "leiden_3": pd.Categorical(
                    ["a1", "a1", "b1", "b1"],
                    categories=["a1", "b1"],
                ),
            },
            index=[f"ref-{index}" for index in range(4)],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    collector = AssignmentGuidanceCollector()
    runner = _multi_type_runner(
        impute_module,
        collector,
        spatial,
        reference,
    )
    captured = _patch_problem(impute_module, monkeypatch, runner)

    result = runner.local_impute(
        reference,
        "leiden_3",
        guidance_scope="panel",
    )

    assert len(captured["solve_calls"]) == 2
    assert result.n_obs == 4
    assert [event["outcome"] for event in collector.events] == ["off", "off"]


def test_metadata_only_unused_subcluster_level_preserves_base_result(
    impute_module,
    monkeypatch,
):
    spatial = AnnData(
        X=sparse.csr_matrix([[2.0, 1.0], [1.0, 2.0]]),
        obs=pd.DataFrame(
            {"Level1": ["A", "A"]},
            index=["spot-1", "spot-2"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )

    def reference(subclusters):
        return AnnData(
            X=sparse.csr_matrix([[2.0, 1.0], [1.0, 2.0]]),
            obs=pd.DataFrame(
                {"Level1": ["A", "A"], "leiden_3": subclusters},
                index=["ref-1", "ref-2"],
            ),
            var=pd.DataFrame(index=["g1", "g2"]),
        )

    plain_reference = reference(["s1", "s1"])
    plain_collector = AssignmentGuidanceCollector()
    plain_runner = _multi_type_runner(
        impute_module,
        plain_collector,
        spatial,
        plain_reference,
    )
    plain_capture = _patch_problem(
        impute_module,
        monkeypatch,
        plain_runner,
    )
    plain_result = plain_runner.local_impute(
        plain_reference,
        "leiden_3",
        guidance_scope="plain",
    )

    categorical_reference = reference(
        pd.Categorical(
            ["s1", "s1"],
            categories=["s1", "unused"],
        )
    )
    categorical_collector = AssignmentGuidanceCollector()
    categorical_runner = _multi_type_runner(
        impute_module,
        categorical_collector,
        spatial,
        categorical_reference,
    )
    categorical_capture = _patch_problem(
        impute_module,
        monkeypatch,
        categorical_runner,
    )
    categorical_result = categorical_runner.local_impute(
        categorical_reference,
        "leiden_3",
        guidance_scope="categorical",
    )

    assert len(plain_capture["solve_calls"]) == 1
    assert len(categorical_capture["solve_calls"]) == 1
    np.testing.assert_allclose(
        categorical_capture["solve_calls"][0]["source_mass"],
        plain_capture["solve_calls"][0]["source_mass"],
    )
    np.testing.assert_allclose(
        categorical_capture["solve_calls"][0]["target_mass"],
        plain_capture["solve_calls"][0]["target_mass"],
    )
    np.testing.assert_allclose(
        categorical_result.X.toarray(),
        plain_result.X.toarray(),
    )
    assert _event(categorical_collector)["outcome"] == "off"


@pytest.mark.parametrize(
    ("failure_stage", "expected_outcome", "reason"),
    [
        ("solver", "failed", "solver_failed"),
        ("assembly", "failed", "result_assembly_failed"),
    ],
)
def test_failures_record_terminal_stage_and_reraise(
    impute_module,
    monkeypatch,
    failure_stage,
    expected_outcome,
    reason,
):
    collector = AssignmentGuidanceCollector()
    runner, reference = _runner(impute_module, collector)
    error = RuntimeError(f"{failure_stage} boom")
    _patch_problem(
        impute_module,
        monkeypatch,
        runner,
        solver_error=error if failure_stage == "solver" else None,
        assembly_error=error if failure_stage == "assembly" else None,
    )

    with pytest.raises(RuntimeError, match=f"{failure_stage} boom"):
        runner.local_impute(reference, "leiden_3", guidance_scope="panel")

    event = _event(collector)
    assert event["attempted"] is True
    assert event["outcome"] == expected_outcome
    assert event["reason"] == reason


def test_result_assembly_interrupt_records_interrupted_and_reraises(
    impute_module,
    monkeypatch,
):
    collector = AssignmentGuidanceCollector()
    runner, reference = _runner(impute_module, collector)
    _patch_problem(
        impute_module,
        monkeypatch,
        runner,
        assembly_error=KeyboardInterrupt("assembly interrupted"),
    )

    with pytest.raises(KeyboardInterrupt, match="assembly interrupted"):
        runner.local_impute(reference, "leiden_3", guidance_scope="panel")

    event = _event(collector)
    assert event["attempted"] is True
    assert event["outcome"] == "interrupted"
    assert event["reason"] == "result_assembly_interrupted"


def _run_public_cli_route(
    impute_module,
    monkeypatch,
    tmp_path,
    *,
    confounding,
    enable_guidance,
):
    from revise.backend import adapters
    from revise.backend.policies import (
        ModeEvaluationPolicy,
        ModeValidationPolicy,
    )
    from revise.backend.registry import StrategyRegistry
    from revise.benchmark import cli
    from revise.recon.pipeline import UnifiedReconstructionPipeline

    spatial, reference = _inputs()
    if not enable_guidance:
        spatial.obsm["Level1"] = np.full((spatial.n_obs, 1), np.nan)
        reference.obsm["Level1"] = np.full((reference.n_obs, 1), np.nan)

    captured = {"solve_calls": [], "assignment_reads": 0}
    observed = {
        "registry_strategy": None,
        "prepare_routes": [],
        "solve_routes": [],
        "runner_configs": [],
    }

    monkeypatch.setattr(adapters, "_install_safe_topology_patch", lambda: None)
    input_service = SimpleNamespace(
        read_st_adata=lambda _path: spatial.copy(),
        read_real_adata=lambda _path: spatial.copy(),
        read_sc_ref_adata=lambda _path: reference.copy(),
    )
    monkeypatch.setattr(adapters, "_input_service", lambda _ctx: input_service)
    monkeypatch.setattr(
        ModeValidationPolicy,
        "validate",
        lambda _self, _ctx: None,
    )
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

    real_registry_get = StrategyRegistry.get

    def registry_get(registry, strategy_id):
        observed["registry_strategy"] = strategy_id
        return real_registry_get(registry, strategy_id)

    monkeypatch.setattr(StrategyRegistry, "get", registry_get)

    real_prepare = adapters.ScSvcImputeBenchmarkStrategy.prepare_context
    real_solve = adapters.ScSvcImputeBenchmarkStrategy.solve_ot

    def prepare(strategy, ctx):
        observed["prepare_routes"].append(ctx.route_key)
        real_prepare(strategy, ctx)
        observed["runner_configs"].append(ctx.runner_config)

    def solve(strategy, ctx):
        observed["solve_routes"].append(ctx.route_key)
        return real_solve(strategy, ctx)

    monkeypatch.setattr(
        adapters.ScSvcImputeBenchmarkStrategy,
        "prepare_context",
        prepare,
    )
    monkeypatch.setattr(
        adapters.ScSvcImputeBenchmarkStrategy,
        "global_anchoring",
        lambda _self, _ctx: None,
    )
    monkeypatch.setattr(
        adapters.ScSvcImputeBenchmarkStrategy,
        "solve_ot",
        solve,
    )

    def runner_init(self, st_adata, sc_ref_adata, config, real_st_adata, logger):
        self.st_adata = st_adata.copy()
        self.sc_ref_adata = sc_ref_adata.copy()
        self.config = config
        self.real_st_adata = real_st_adata.copy()
        self.logger = logger
        self.svc = {}

        def impute(spatial_block, _reference_block, **_kwargs):
            return spatial_block.copy()

        self.gene_impute = SimpleNamespace(run=impute)

    def runner_local_refinement(self):
        self.svc["sc_svc_impute_in_panel"] = self.local_impute(
            self.sc_ref_adata,
            "leiden_3",
            guidance_scope=confounding,
        )

    monkeypatch.setattr(impute_module.ScSVCImpute, "__init__", runner_init)
    monkeypatch.setattr(
        impute_module.ScSVCImpute,
        "local_refinement",
        runner_local_refinement,
    )
    monkeypatch.setattr(
        impute_module,
        "bhattacharyya_distance",
        lambda profiles, spots: np.zeros(
            (profiles.shape[0], spots.shape[0])
        ),
    )

    def solve_local(source_mass, target_mass, cost, **kwargs):
        captured["solve_calls"].append(
            {
                "cost": np.asarray(cost).copy(),
                "reference_measure": kwargs.get("reference_measure"),
            }
        )
        return np.outer(source_mass, target_mass)

    monkeypatch.setattr(impute_module, "solve_local_ot", solve_local)

    def merge(adata, *, subcluster, mode):
        labels = list(dict.fromkeys(adata.obs[subcluster].astype(str)))
        return AnnData(
            X=np.ones((len(labels), adata.n_vars)),
            obs=pd.DataFrame(
                {subcluster: labels},
                index=[f"merged-{label}" for label in labels],
            ),
            var=pd.DataFrame(index=adata.var_names.copy()),
        )

    monkeypatch.setattr(impute_module, "merge_subcluster", merge)
    real_assignment_loader = impute_module._assignment_state_from_adata

    def assignment_loader(*args, **kwargs):
        captured["assignment_reads"] += 1
        if not enable_guidance:
            pytest.fail("route-default off must not read assignment state")
        return real_assignment_loader(*args, **kwargs)

    monkeypatch.setattr(
        impute_module,
        "_assignment_state_from_adata",
        assignment_loader,
    )

    output_root = tmp_path / f"output-{confounding}-{enable_guidance}"
    argv = [
        "revise-benchmark",
        "--confounding",
        confounding,
        "--data-root",
        str(tmp_path / "data"),
        "--output-root",
        str(output_root),
        "--sample-name",
        "sample",
        "--config",
        "revise/revise.yaml",
        "--seed-scope",
        "run",
    ]
    if enable_guidance:
        argv.extend(
            [
                "--local-refinement-guidance",
                "prefer",
                "--local-refinement-compatibility-mode",
                "cost",
            ]
        )
    monkeypatch.setattr(sys, "argv", argv)

    cli.main()

    provenance_paths = list(output_root.rglob("provenance.json"))
    assert len(provenance_paths) == 1
    manifest = json.loads(provenance_paths[0].read_text())
    return manifest, captured, observed


@pytest.mark.parametrize(
    ("profile", "confounding"),
    [
        ("benchmark_impute_panel", "gene_panel"),
        ("benchmark_impute_dropout", "gene_dropout"),
    ],
)
def test_public_cli_route_default_off_reaches_manifest_without_assignment_reads(
    impute_module,
    monkeypatch,
    tmp_path,
    profile,
    confounding,
):
    manifest, captured, observed = _run_public_cli_route(
        impute_module,
        monkeypatch,
        tmp_path,
        confounding=confounding,
        enable_guidance=False,
    )

    assert manifest["profile"] == profile
    assert observed["registry_strategy"] == "ScSvcImputeBenchmarkStrategy"
    assert observed["prepare_routes"] == [f"sim2real:{confounding}"]
    assert observed["solve_routes"] == [f"sim2real:{confounding}"]
    assert observed["runner_configs"][0].posterior_conditioning_enabled is False
    np.testing.assert_allclose(captured["solve_calls"][0]["cost"], 0.0)
    assert captured["assignment_reads"] == 0
    guidance = manifest["assignment_guidance"]
    assert guidance["resolved"]["guidance"] == "off"
    assert guidance["events"][0]["outcome"] == "off"
    assert guidance["events"][0]["route"] == f"sim2real:{confounding}"


@pytest.mark.parametrize(
    ("profile", "confounding"),
    [
        ("benchmark_impute_panel", "gene_panel"),
        ("benchmark_impute_dropout", "gene_dropout"),
    ],
)
def test_public_cli_guidance_reaches_manifest_with_bilateral_lineage(
    impute_module,
    monkeypatch,
    tmp_path,
    profile,
    confounding,
):
    manifest, captured, observed = _run_public_cli_route(
        impute_module,
        monkeypatch,
        tmp_path,
        confounding=confounding,
        enable_guidance=True,
    )

    assert manifest["profile"] == profile
    assert observed["registry_strategy"] == "ScSvcImputeBenchmarkStrategy"
    assert observed["prepare_routes"] == [f"sim2real:{confounding}"]
    assert observed["solve_routes"] == [f"sim2real:{confounding}"]
    assert observed["runner_configs"][0].posterior_conditioning_enabled is True
    assert captured["assignment_reads"] == 2
    event = manifest["assignment_guidance"]["events"][0]
    assert event["route"] == f"sim2real:{confounding}"
    assert event["outcome"] == "applied"
    assert event["operator"] == "imputation_ot"
    assert event["left_assignment"]["lineage"][0] == {
        "operation": "load",
        "container": "obsm",
        "key": "Level1",
    }
    assert event["right_assignment"]["source"] == "aggregate(obsm[Level1])"
    assert event["right_assignment"]["lineage"][-1]["operation"] == "aggregate"
    assert np.count_nonzero(captured["solve_calls"][0]["cost"]) == 4


@pytest.mark.parametrize(
    ("mode", "solver", "message"),
    [
        ("application", "pot", "application local refinement"),
        ("benchmark", "tacco", "TACCO local refinement"),
    ],
)
def test_reference_mode_is_rejected_by_unsupported_preflight_surfaces(
    mode,
    solver,
    message,
):
    from revise.backend.policies import ModeValidationPolicy

    ctx = SimpleNamespace(
        merged_config={
            "local_refinement": {
                "guidance": "prefer",
                "compatibility": {"mode": "reference"},
            },
            "ot": {"lr": {"solver": solver}},
        },
        runtime={"mode": mode, "task": "sc_svc_impute"},
    )

    with pytest.raises(ValueError, match=message):
        ModeValidationPolicy._validate_solver_compatibility(ctx)
