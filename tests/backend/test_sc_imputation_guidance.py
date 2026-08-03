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


def _config():
    return SimpleNamespace(
        cell_type_col="Level1",
        rec_ot_method="pot",
        rec_impute_pot_reg=0.1,
        rec_impute_pot_reg_m=0.0,
        rec_impute_pot_reg_type="kl",
        rec_merge_subcluster_method="mean",
        rec_impute_prune_flag=False,
    )


def _inputs():
    spatial = AnnData(
        X=sparse.csr_matrix([[2.0, 1.0], [1.0, 2.0]]),
        obs=pd.DataFrame(
            {"Level1": ["A", "A"]},
            index=["spot-1", "spot-2"],
        ),
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
            index=["ref-1", "ref-2", "ref-3", "ref-4"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    spatial.obsm["Level1"] = np.full((spatial.n_obs, 1), np.nan)
    reference.obsm["Level1"] = np.full((reference.n_obs, 1), np.nan)
    return spatial, reference


def _runner(module, spatial, reference):
    runner = module.ScSVCImpute.__new__(module.ScSVCImpute)
    runner.st_adata = spatial
    runner.sc_ref_adata = reference
    runner.config = _config()
    runner.logger = logging.getLogger("test-sc-imputation-base")
    runner.svc = {}
    return runner


def _patch_problem(
    module,
    monkeypatch,
    runner,
    *,
    solver_error=None,
    assembly_error=None,
):
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

    def impute(spatial_block, _reference_block, **kwargs):
        captured["weights"].append(kwargs["neighbor_weights"].copy())
        if assembly_error is not None:
            raise assembly_error
        return spatial_block.copy()

    runner.gene_impute = SimpleNamespace(run=impute)
    return captured


def test_imputation_solves_the_base_problem(
    impute_module,
    monkeypatch,
):
    spatial, reference = _inputs()
    runner = _runner(impute_module, spatial, reference)
    captured = _patch_problem(impute_module, monkeypatch, runner)

    result = runner.local_impute(reference, "leiden_3")

    assert len(captured["solve_calls"]) == 1
    np.testing.assert_allclose(captured["solve_calls"][0]["cost"], 0.0)
    np.testing.assert_allclose(result.X.toarray(), spatial.X.toarray())


@pytest.mark.parametrize(
    ("spatial_values", "reference_values"),
    [
        ([[0.0, 0.0], [1.0, 1.0]], [[1.0, 1.0], [1.0, 1.0]]),
        ([[1.0, 1.0], [1.0, 1.0]], [[0.0, 0.0], [0.0, 0.0]]),
    ],
)
def test_zero_observed_marginals_return_empty_base_result(
    impute_module,
    monkeypatch,
    spatial_values,
    reference_values,
):
    spatial = AnnData(
        X=sparse.csr_matrix(spatial_values),
        obs=pd.DataFrame({"Level1": ["A", "A"]}, index=["s1", "s2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    reference = AnnData(
        X=sparse.csr_matrix(reference_values),
        obs=pd.DataFrame(
            {"Level1": ["A", "A"], "leiden_3": ["c1", "c1"]},
            index=["r1", "r2"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    runner = _runner(impute_module, spatial, reference)
    captured = _patch_problem(impute_module, monkeypatch, runner)

    result = runner.local_impute(reference, "leiden_3")

    assert captured["solve_calls"] == []
    assert result.n_obs == 0


def test_disjoint_categorical_subclusters_preserve_base_output(
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
    runner = _runner(impute_module, spatial, reference)
    captured = _patch_problem(impute_module, monkeypatch, runner)

    result = runner.local_impute(reference, "leiden_3")

    assert len(captured["solve_calls"]) == 2
    assert result.n_obs == 4


def test_reference_only_cell_type_is_skipped(
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
    runner = impute_module.ScSVCImpute.__new__(impute_module.ScSVCImpute)
    runner.st_adata = spatial
    runner.sc_ref_adata = reference
    runner.config = SimpleNamespace(
        cell_type_col="Level1",
        rec_ot_method="pot",
        rec_impute_pot_reg=0.1,
        rec_impute_pot_reg_m=0.0,
        rec_impute_pot_reg_type="kl",
        rec_merge_subcluster_method="mean",
        rec_impute_prune_flag=False,
    )
    runner.logger = logging.getLogger("test-reference-only-cell-type")
    captured = _patch_problem(impute_module, monkeypatch, runner)

    result = runner.local_impute(reference, "leiden_3")

    assert len(captured["solve_calls"]) == 1
    assert result.obs["Level1"].tolist() == ["A", "A"]
    assert result.obs_names.tolist() == ["spot-a1", "spot-a2"]


@pytest.mark.parametrize("failure_stage", ["solver", "assembly"])
def test_base_failures_still_propagate(
    impute_module,
    monkeypatch,
    failure_stage,
):
    spatial, reference = _inputs()
    runner = _runner(impute_module, spatial, reference)
    error = RuntimeError(f"{failure_stage} boom")
    _patch_problem(
        impute_module,
        monkeypatch,
        runner,
        solver_error=error if failure_stage == "solver" else None,
        assembly_error=error if failure_stage == "assembly" else None,
    )

    with pytest.raises(RuntimeError, match=f"{failure_stage} boom"):
        runner.local_impute(reference, "leiden_3")


def test_coupling_alignment_preserves_merged_reference_order(impute_module):
    coupling = pd.DataFrame(
        [[0.2, 0.8]],
        index=["spot-1"],
        columns=["s2", "s1"],
    )
    merged = AnnData(
        X=np.ones((2, 1)),
        obs=pd.DataFrame(
            {"leiden_3": ["s1", "s2"]},
            index=["merged-1", "merged-2"],
        ),
        var=pd.DataFrame(index=["g1"]),
    )

    aligned = impute_module._align_coupling_to_merged_subclusters(
        coupling,
        merged,
        "leiden_3",
    )

    assert aligned.columns.tolist() == ["merged-1", "merged-2"]
    np.testing.assert_allclose(aligned.to_numpy(), [[0.8, 0.2]])


def test_coupling_alignment_rejects_label_mismatch(impute_module):
    coupling = pd.DataFrame([[1.0]], columns=["s1"])
    merged = AnnData(
        X=np.ones((1, 1)),
        obs=pd.DataFrame({"leiden_3": ["s2"]}, index=["merged-1"]),
        var=pd.DataFrame(index=["g1"]),
    )

    with pytest.raises(ValueError, match="labels differ"):
        impute_module._align_coupling_to_merged_subclusters(
            coupling,
            merged,
            "leiden_3",
        )
