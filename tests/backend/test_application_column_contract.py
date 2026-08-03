from __future__ import annotations

import importlib
import logging
import sys
import types
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, concat as anndata_concat
from scipy import sparse


class _ReachedConfiguredColumn(RuntimeError):
    pass


ISOLATED_RUNNER_MODULE_PREFIXES = (
    "scanpy",
    "squidpy",
    "numba",
    "revise.backend.adapters",
    "revise.backend.kernels",
    "revise.backend.ops",
    "revise.backend.runners",
)
_MISSING = object()


def _module_names_for_prefixes(prefixes=ISOLATED_RUNNER_MODULE_PREFIXES):
    return tuple(
        name
        for name in sys.modules
        if any(
            name == prefix or name.startswith(f"{prefix}.")
            for prefix in prefixes
        )
    )


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


def _remove_modules(module_names) -> None:
    for module_name in module_names:
        sys.modules.pop(module_name, None)
    for module_name in module_names:
        parent_name, separator, attribute = module_name.rpartition(".")
        parent = sys.modules.get(parent_name)
        if separator and parent is not None and hasattr(parent, attribute):
            delattr(parent, attribute)


@pytest.fixture(autouse=True)
def lightweight_runner_dependencies():
    snapshot = _snapshot_modules(_module_names_for_prefixes())
    _remove_modules(_module_names_for_prefixes())
    if "scanpy" not in sys.modules:
        scanpy = types.ModuleType("scanpy")
        scanpy.AnnData = AnnData
        scanpy.concat = anndata_concat
        scanpy.pp = SimpleNamespace(
            filter_cells=lambda *_args, **_kwargs: None,
            filter_genes=lambda *_args, **_kwargs: None,
            normalize_total=lambda *_args, **_kwargs: None,
        )
        scanpy.pl = SimpleNamespace()
        scanpy.tl = SimpleNamespace()
        sys.modules["scanpy"] = scanpy
    if "squidpy" not in sys.modules:
        squidpy = types.ModuleType("squidpy")
        squidpy.gr = SimpleNamespace()
        sys.modules["squidpy"] = squidpy
    if "numba" not in sys.modules:
        numba = types.ModuleType("numba")

        def njit(*args, **_kwargs):
            if args and callable(args[0]):
                return args[0]
            return lambda function: function

        numba.njit = njit
        numba.prange = range
        sys.modules["numba"] = numba
    try:
        yield
    finally:
        _remove_modules(_module_names_for_prefixes())
        _restore_modules(snapshot)


@pytest.fixture
def load_runner():
    snapshot = _snapshot_modules(_module_names_for_prefixes())
    _remove_modules(_module_names_for_prefixes())

    def normalize_total(adata, target_sum=1e4, **_kwargs):
        values = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)
        row_sums = values.sum(axis=1, keepdims=True)
        normalized = values / np.maximum(row_sums, 1e-12) * float(target_sum)
        adata.X = sparse.csr_matrix(normalized) if sparse.issparse(adata.X) else normalized

    scanpy = types.ModuleType("scanpy")
    scanpy.AnnData = AnnData
    scanpy.concat = anndata_concat
    scanpy.pp = SimpleNamespace(normalize_total=normalize_total)
    scanpy.pl = SimpleNamespace()
    scanpy.tl = SimpleNamespace()
    squidpy = types.ModuleType("squidpy")
    squidpy.gr = SimpleNamespace()
    numba = types.ModuleType("numba")

    def njit(*args, **_kwargs):
        if args and callable(args[0]):
            return args[0]
        return lambda function: function

    numba.njit = njit
    numba.prange = range
    sys.modules["scanpy"] = scanpy
    sys.modules["squidpy"] = squidpy
    sys.modules["numba"] = numba

    try:
        yield lambda module_name: importlib.import_module(
            f"revise.backend.runners.{module_name}"
        )
    finally:
        _remove_modules(_module_names_for_prefixes())
        _restore_modules(snapshot)


def test_lightweight_fixture_restores_relevant_prefix_keys_and_identities():
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
    probe_name = "revise.backend.ops._column_contract_probe"
    probe = types.ModuleType(probe_name)
    fixture = lightweight_runner_dependencies.__wrapped__()
    next(fixture)
    try:
        importlib.import_module("revise.backend.ops.assignment")
        importlib.import_module("revise.backend.ops.local_ot")
        importlib.import_module("revise.backend.ops.posterior_conditioning")
        importlib.import_module("revise.backend.ops.sr_allocation")
        sys.modules[probe_name] = probe
        ops_parent = importlib.import_module("revise.backend.ops")
        setattr(ops_parent, "_column_contract_probe", probe)
        with pytest.raises(StopIteration):
            next(fixture)

        after = relevant_modules()
        assert after.keys() == before.keys()
        assert all(after[name] is module for name, module in before.items())
        for name in before:
            parent_name, separator, attribute = name.rpartition(".")
            if not separator or parent_name not in sys.modules:
                continue
            expected = before_parent_attributes[name]
            if expected is _MISSING:
                assert not hasattr(sys.modules[parent_name], attribute)
            else:
                assert getattr(sys.modules[parent_name], attribute) is expected
    finally:
        sys.modules.pop(probe_name, None)
        ops_parent = sys.modules.get("revise.backend.ops")
        if (
            ops_parent is not None
            and getattr(ops_parent, "_column_contract_probe", None) is probe
        ):
            delattr(ops_parent, "_column_contract_probe")
        try:
            fixture.close()
        except RuntimeError:
            pass


def _sp_argmax_runner(module):
    n_obs = 51
    obs_names = [f"cell-{idx}" for idx in range(n_obs)]
    st = AnnData(
        X=sparse.csr_matrix(np.ones((n_obs, 2))),
        obs=pd.DataFrame({"major_type": ["A"] * n_obs}, index=obs_names),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    reference = AnnData(
        X=sparse.csr_matrix(np.ones((2, 2))),
        obs=pd.DataFrame({"major_type": ["A", "A"]}, index=["ref-1", "ref-2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    runner = module.SpSVC.__new__(module.SpSVC)
    runner.st_adata = st
    runner.sc_ref_adata = reference
    runner.config = SimpleNamespace(
        plot_flag=False,
        cell_type_col="major_type",
        rec_graph_method="pca",
        rec_graph_alpha=0.0,
        rec_graph_exp_neighbor_num=1,
        rec_graph_spatial_neighbor_num=1,
        rec_graph_n_neighbors=1,
        rec_ot_method="pot",
        rec_pot_reg=0.1,
        rec_pot_reg_m=0.0,
        rec_pot_reg_type="kl",
        local_refinement_strength=1.0,
    )
    runner.logger = logging.getLogger("test-sp-argmax-fallback")
    runner.graph_aggregate = SimpleNamespace(run=lambda *, adata, **_kwargs: adata)
    runner.svc = {}
    return runner


def test_sp_adapter_normalizes_configured_reference_labels(monkeypatch, tmp_path):
    from revise.backend import adapters

    st = AnnData(
        X=sparse.csr_matrix(np.ones((2, 2))),
        obs=pd.DataFrame(index=["s1", "s2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    reference = AnnData(
        X=sparse.csr_matrix(np.ones((2, 2))),
        obs=pd.DataFrame(
            {
                "Patient": ["sample", "sample"],
                "major_type": ["T/NK", "T/NK"],
                "minor_type": ["T/1", "T/2"],
            },
            index=["c1", "c2"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )

    class InputService:
        def read_st_adata(self, path):
            return st.copy()

        def read_sc_ref_adata(self, path):
            return reference.copy()

    monkeypatch.setattr(adapters, "_input_service", lambda ctx: InputService())
    monkeypatch.setattr(
        adapters.sc.pp,
        "filter_cells",
        lambda *args, **kwargs: None,
        raising=False,
    )
    monkeypatch.setattr(
        adapters.sc.pp,
        "filter_genes",
        lambda *args, **kwargs: None,
        raising=False,
    )
    ctx = SimpleNamespace(
        merged_config={
            "ot": {
                "ga": {
                    "solver": "pot",
                    "pot": {"reg": 0.1, "reg_m": 0.0, "reg_type": "entropy"},
                },
                "lr": {
                    "solver": "pot",
                    "pot": {"reg": 0.1, "reg_m": 0.0, "reg_type": "kl"},
                },
                "impute": {"reg": 5.0, "reg_m": 0.0, "reg_type": "kl"},
            },
                "preprocess": {},
                "graph": {},
                "plot": {},
                "local_refinement": {"strength": 0.2},
        },
        io={
            "sample_name": "sample",
            "data_root": str(tmp_path),
            "st_file": "st.h5ad",
            "sc_ref_file": "sc.h5ad",
            "patient_key": "Patient",
        },
        columns={
            "cell_type_col": "major_type",
            "sub_cell_type_col": "minor_type",
            "confidence_col": "Confidence",
            "unknown_key": "Unknown",
        },
        run_dir=tmp_path,
        logger=logging.getLogger("test-sp-label-normalization"),
    )

    adapters.SpSvcApplicationStrategy().prepare_context(ctx)

    assert ctx.runner.sc_ref_adata.obs["major_type"].tolist() == ["T_NK", "T_NK"]
    assert ctx.runner.sc_ref_adata.obs["minor_type"].tolist() == ["T_1", "T_2"]


def test_reference_label_normalization_rejects_category_collisions():
    from revise.backend.adapters import _replace_slash_labels

    reference = AnnData(
        X=np.ones((2, 1)),
        obs=pd.DataFrame({"major_type": ["A/B", "A_B"]}, index=["c1", "c2"]),
        var=pd.DataFrame(index=["g1"]),
    )

    with pytest.raises(ValueError, match="collide after slash normalization"):
        _replace_slash_labels(reference, ["major_type"])


def test_sp_local_refinement_trims_by_the_configured_cell_type(monkeypatch):
    from revise.backend.runners import sp_svc_application as module

    st = AnnData(
        X=sparse.csr_matrix(np.ones((2, 2))),
        obs=pd.DataFrame({"major_type": ["A", "A"]}, index=["s1", "s2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    reference = AnnData(
        X=sparse.csr_matrix(np.ones((2, 2))),
        obs=pd.DataFrame({"major_type": ["A", "A"]}, index=["c1", "c2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    runner = module.SpSVC.__new__(module.SpSVC)
    runner.st_adata = st
    runner.sc_ref_adata = reference
    runner.config = SimpleNamespace(plot_flag=False, cell_type_col="major_type")
    runner.logger = logging.getLogger("test-sp-custom-column")

    def capture_column(_st, _reference, column):
        assert column == "major_type"
        raise _ReachedConfiguredColumn

    monkeypatch.setattr(module, "trim_sp_adata", capture_column)

    with pytest.raises(_ReachedConfiguredColumn):
        runner.local_refinement()


def test_sp_argmax_without_soft_posterior_fails_before_local_solver(
    monkeypatch,
    load_runner,
):
    module = load_runner("sp_svc_application")

    runner = _sp_argmax_runner(module)
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
    monkeypatch.setattr(
        module,
        "_compute_topk_expression",
        lambda **_kwargs: (
            np.ones((51, 1)),
            np.ones(1),
            np.zeros((51, 1), dtype=np.int32),
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
        "solve_local_ot",
        lambda *_args, **_kwargs: pytest.fail(
            "missing soft Q must fail before local solve"
        ),
    )

    with pytest.raises(ValueError, match=r"obsm\[major_type\]"):
        runner.local_refinement()


def test_sr_reference_profiles_ignore_an_unrelated_clusters_column(monkeypatch):
    from revise.backend.runners import sc_svc_sr_application as module

    st = AnnData(
        X=np.ones((1, 2)),
        obs=pd.DataFrame(index=["spot-1"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    st.obsm["major_type"] = pd.DataFrame(
        [[1.0]], index=st.obs_names, columns=["A"]
    )
    st.obs["major_type"] = "A"
    reference = AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame(
            {"major_type": ["A", "A"], "clusters": ["0", "1"]},
            index=["c1", "c2"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    runner = module.ScSVCSr.__new__(module.ScSVCSr)
    runner.st_adata = st
    runner.sc_ref_adata = reference
    runner.config = SimpleNamespace(cell_type_col="major_type")
    runner.spot_sr = SimpleNamespace(run=lambda _runner: None)
    runner.logger = logging.getLogger("test-sr-custom-column")

    def capture_column(_adata, *, key_type, type_list):
        assert key_type == "major_type"
        assert type_list == ["A"]
        raise _ReachedConfiguredColumn

    monkeypatch.setattr(module.sc.pp, "normalize_total", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(module, "construct_sc_ref", capture_column)

    with pytest.raises(_ReachedConfiguredColumn):
        runner.local_refinement()


def test_sr_benchmark_custom_broad_column_drives_mandatory_allocation(
    monkeypatch,
    load_runner,
):
    module = load_runner("sc_svc_sr_benchmark")

    st = AnnData(
        X=np.ones((1, 2)),
        obs=pd.DataFrame(index=["spot-1"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    st.obsm["major_type"] = pd.DataFrame(
        [[1.0]],
        index=st.obs_names,
        columns=["A"],
    )
    st.obs["major_type"] = "A"
    reference = AnnData(
        X=np.ones((1, 2)),
        obs=pd.DataFrame({"major_type": ["A"]}, index=["cell-1"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    runner = module.ScSVCSr.__new__(module.ScSVCSr)
    runner.st_adata = st
    runner.sc_ref_adata = reference
    runner.config = SimpleNamespace(
        cell_type_col="major_type",
        rec_graph_agg_enabled=False,
        local_refinement_strength=0.0,
    )
    runner.svc_obs = pd.DataFrame(
        {
            "spot_name": ["spot-1"],
            "cell_id": ["virtual-1"],
            "x": [0.0],
            "y": [0.0],
            "true_cell_type": ["Unknown"],
        }
    )
    runner.spot_sr = SimpleNamespace(
        run=lambda target: target.svc_obs.__setitem__("cell_type", "A")
    )
    runner.logger = logging.getLogger("test-sr-benchmark-literal-level1")
    runner.svc = {}

    monkeypatch.setattr(module.sc.pp, "normalize_total", lambda *_args, **_kwargs: None)

    runner.local_refinement()

    assert runner.svc["sc_svc_dec"].uns["sr_allocation"] == {
        "broad_key": "major_type",
        "posterior_key": "major_type",
        "operator": "closed_form_reference_allocation",
        "beta": 1.0,
    }
    assert runner.svc_obs["cell_type"].tolist() == ["A"]


def test_sr_benchmark_reference_uses_configured_broad_column_without_clusters(
    monkeypatch,
    load_runner,
):
    module = load_runner("sc_svc_sr_benchmark")
    st = AnnData(
        X=np.ones((1, 1)),
        obs=pd.DataFrame(index=["spot-1"]),
        var=pd.DataFrame(index=["g1"]),
    )
    st.obsm["major_type"] = pd.DataFrame(
        [[1.0]],
        index=st.obs_names,
        columns=["A"],
    )
    st.obs["major_type"] = "A"
    reference = AnnData(
        X=np.ones((1, 1)),
        obs=pd.DataFrame({"major_type": ["A"]}, index=["cell-1"]),
        var=pd.DataFrame(index=["g1"]),
    )
    runner = module.ScSVCSr.__new__(module.ScSVCSr)
    runner.st_adata = st
    runner.sc_ref_adata = reference
    runner.config = SimpleNamespace(
        cell_type_col="major_type",
        rec_graph_agg_enabled=False,
        local_refinement_strength=0.0,
    )
    runner.svc_obs = pd.DataFrame(
        {
            "spot_name": ["spot-1"],
            "cell_id": ["virtual-1"],
            "x": [0.0],
            "y": [0.0],
            "true_cell_type": ["Unknown"],
        }
    )
    reached_assignment = {"value": False}
    def assign_rows(target):
        reached_assignment["value"] = True
        target.svc_obs["cell_type"] = "A"

    runner.spot_sr = SimpleNamespace(run=assign_rows)
    runner.logger = logging.getLogger("test-sr-benchmark-literal-clusters")
    runner.svc = {}
    monkeypatch.setattr(module.sc.pp, "normalize_total", lambda *_args, **_kwargs: None)

    runner.local_refinement()

    assert reached_assignment["value"] is True
    assert runner.svc["sc_svc_dec"].obs["cell_type"].tolist() == ["A"]


def test_sr_zero_strength_preserves_quota_row_and_expression_allocation(
    monkeypatch,
    load_runner,
):
    module = load_runner("sc_svc_sr_application")
    from revise.backend.kernels.spot_sr import SpotSrKernel
    from revise.backend.ops.sr_allocation import (
        mandatory_reference_allocation as real_reference_allocation,
    )

    st = AnnData(
        X=np.array([[2.0, 8.0]]),
        obs=pd.DataFrame(index=["spot-1"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    st.obsm["major_type"] = pd.DataFrame(
        [[0.75, 0.25]],
        index=st.obs_names,
        columns=["A", "B"],
    )
    st.obs["major_type"] = "A"
    reference = AnnData(
        X=np.array([[9.0, 1.0], [1.0, 9.0]]),
        obs=pd.DataFrame(
            {"major_type": ["A", "B"]},
            index=["ref-a", "ref-b"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    config = SimpleNamespace(
        pm_on_cell_file="/path/that/does/not/exist.csv",
        svc_completeness=True,
        sr_assignment_seed=17,
        cell_type_col="major_type",
        local_refinement_strength=0.0,
        rec_match_spot_sum=True,
    )
    runner = module.ScSVCSr.__new__(module.ScSVCSr)
    runner.st_adata = st
    runner.sc_ref_adata = reference
    runner.config = config
    runner.svc_obs = pd.DataFrame(
        {
            "spot_name": ["spot-1"] * 3,
            "cell_id": ["virtual-1", "virtual-2", "virtual-3"],
            "x": [0.0] * 3,
            "y": [0.0] * 3,
            "true_cell_type": ["Unknown"] * 3,
        }
    )
    runner.spot_sr = SpotSrKernel(config, logging.getLogger("test-sr-assignment"))
    runner.graph_aggregate = SimpleNamespace()
    runner.logger = logging.getLogger("test-sr-allocation")
    runner.svc = {}
    captured = {}

    def capture_allocation(spot_expression, posterior, reference_profiles, **kwargs):
        captured["spot_expression"] = np.asarray(spot_expression).copy()
        captured["posterior"] = posterior.copy()
        allocated = real_reference_allocation(
            spot_expression,
            posterior,
            reference_profiles,
            **kwargs,
        )
        captured["allocated"] = allocated.copy()
        return allocated

    monkeypatch.setattr(module, "mandatory_reference_allocation", capture_allocation)

    runner.local_refinement()

    output = runner.svc["sc_svc_dec"]
    assert output.obs["cell_type"].value_counts().to_dict() == {"A": 2, "B": 1}
    np.testing.assert_allclose(captured["posterior"].to_numpy(), [[0.75, 0.25]])
    np.testing.assert_allclose(captured["spot_expression"], [[2000.0, 8000.0]])
    np.testing.assert_allclose(
        captured["allocated"],
        [[[1928.57142857, 71.42857143], [2000.0, 6000.0]]],
    )

    type_index = {"A": 0, "B": 1}
    replicated = np.stack(
        [
            captured["allocated"][0, :, type_index[cell_type]]
            for cell_type in output.obs["cell_type"]
        ]
    )
    ratio = captured["spot_expression"][0] / replicated.sum(axis=0)
    expected_rescaled = replicated * ratio

    # CHARACTERIZATION: allocation first creates type-level pools, row
    # assignment replicates the A pool twice, and only then does rescale make
    # the three virtual rows sum back to the normalized spot expression.
    np.testing.assert_allclose(output.X, expected_rescaled)
    np.testing.assert_allclose(np.asarray(output.X).sum(axis=0), [2000.0, 8000.0])
    assert not np.allclose(np.asarray(output.X).sum(axis=1), 10000.0)
