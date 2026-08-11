from __future__ import annotations

import importlib
import itertools
import logging
import sys
import types
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise.backend.kernels.spot_sr import SpotSrKernel
from revise.config import merge_unified_config, resolve_semantic_route
from revise.config.authority import _authority_document


def _merged(profile, *, runtime_overrides=None, io_overrides=None):
    raw = _authority_document()
    namespace, selector = next(
        (namespace, selector)
        for namespace, routes in raw["router"].items()
        for selector, spec in routes.items()
        if spec["profile"] == profile
    )
    runtime = resolve_semantic_route(
        raw,
        **({"svc_type": selector} if namespace == "application" else {"cf": selector}),
    )
    selected_profile = runtime.pop("profile")
    runtime.pop("warning")
    runtime.update(runtime_overrides or {})
    return merge_unified_config(
        raw_config=raw,
        profile=selected_profile,
        runtime_overrides=runtime,
        io_overrides=io_overrides or {},
        algorithm_overrides={},
    )


def _kernel(*, pm=None, seed=42, cell_type_col="Level1"):
    config = SimpleNamespace(
        pm_on_cell=pm,
        svc_completeness=True,
        sr_assignment_seed=seed,
        cell_type_col=cell_type_col,
    )
    return SpotSrKernel(config, logging.getLogger("test-spot-sr-assignment"))


def _svc(spot_cells):
    rows = [
        {"spot_name": spot, "cell_id": cell}
        for spot, cells in spot_cells.items()
        for cell in cells
    ]
    return pd.DataFrame(rows)


def _quota(index=("spot-a",), columns=("A", "B"), values=((1, 1),)):
    return pd.DataFrame(values, index=index, columns=columns, dtype=np.int64)


def test_pm_assignment_matches_brute_force_global_optimum_and_exact_quota():
    svc_obs = _svc({"spot-a": ["c1", "c2", "c3"]})
    quota = _quota(values=((2, 1),))
    pm = pd.DataFrame(
        [[9.0, 8.0], [7.0, 1.0], [6.0, 5.0]],
        index=["c1", "c2", "c3"],
        columns=["A", "B"],
    )

    result = _kernel(pm=pm).assign_cell_types(svc_obs, quota)

    assigned = tuple(result["cell_type"])
    assigned_score = sum(pm.loc[cell, cell_type] for cell, cell_type in zip(result["cell_id"], assigned))
    candidates = set(itertools.permutations(("A", "A", "B")))
    oracle = max(
        sum(pm.loc[cell, cell_type] for cell, cell_type in zip(svc_obs["cell_id"], candidate))
        for candidate in candidates
    )
    assert assigned_score == oracle
    assert result["cell_type"].value_counts().to_dict() == {"A": 2, "B": 1}


def test_pm_assignment_optimizes_raw_scores_not_row_normalized_scores():
    svc_obs = _svc({"spot-a": ["c1", "c2"]})
    quota = _quota()
    pm = pd.DataFrame(
        [[505.0, 495.0], [2.5, 0.5]],
        index=["c1", "c2"],
        columns=["A", "B"],
    )

    result = _kernel(pm=pm).assign_cell_types(svc_obs, quota)

    assert result.set_index("cell_id")["cell_type"].to_dict() == {"c1": "A", "c2": "B"}


def test_pm_is_aligned_by_string_cell_id_and_normalized_class_without_changing_values():
    svc_obs = _svc({"spot-a": [2, 1]})
    quota = _quota(columns=("Mono/Macro", "T"), values=((1, 1),))
    pm = pd.DataFrame(
        [[0.25, -3.0], [4.5, 0.75]],
        index=["1", "2"],
        columns=["T", "Mono_Macro"],
    )
    kernel = _kernel(pm=pm)

    result = kernel.assign_cell_types(svc_obs, quota)

    assert result["cell_id"].tolist() == ["2", "1"]
    expected = pd.DataFrame(
        [[0.75, 4.5], [-3.0, 0.25]],
        index=pd.Index(["2", "1"]),
        columns=["Mono_Macro", "T"],
    )
    pd.testing.assert_frame_equal(kernel.pm_on_cell, expected)
    assert result.set_index("cell_id")["cell_type"].to_dict() == {
        "2": "T",
        "1": "Mono_Macro",
    }


def test_run_with_pm_writes_raw_score_optimal_assignment_back_to_sc_svc():
    svc_obs = _svc({"spot-a": ["c1", "c2"]})
    svc_obs["true_cell_type"] = ["A", "B"]
    sc_svc = SimpleNamespace(
        svc_obs=svc_obs,
        st_adata=SimpleNamespace(
            obsm={"Level1": pd.DataFrame([[0.5, 0.5]], index=["spot-a"], columns=["A", "B"])}
        ),
    )
    pm = pd.DataFrame(
        [[505.0, 495.0], [2.5, 0.5]],
        index=["c1", "c2"],
        columns=["A", "B"],
    )

    _kernel(pm=pm).run(sc_svc)

    assert sc_svc.svc_obs.set_index("cell_id")["cell_type"].to_dict() == {"c1": "A", "c2": "B"}
    assert sc_svc.svc_obs["match"].tolist() == [True, True]


def test_run_without_pm_writes_seeded_exact_quota_assignment_back_to_sc_svc():
    svc_obs = _svc({"spot-a": ["c1", "c2", "c3", "c4"]})
    svc_obs["true_cell_type"] = "Unknown"
    sc_svc = SimpleNamespace(
        svc_obs=svc_obs,
        st_adata=SimpleNamespace(
            obsm={"Level1": pd.DataFrame([[0.25, 0.75]], index=["spot-a"], columns=["A", "B"])}
        ),
    )

    _kernel(seed=19).run(sc_svc)

    expected = np.random.default_rng(19).permutation(["A", "B", "B", "B"]).tolist()
    assert sc_svc.svc_obs["cell_type"].tolist() == expected
    assert sc_svc.svc_obs["cell_type"].value_counts().to_dict() == {"B": 3, "A": 1}


def test_run_reads_cell_contributions_from_the_configured_column():
    svc_obs = _svc({"spot-a": ["c1", "c2"]})
    svc_obs["true_cell_type"] = "Unknown"
    sc_svc = SimpleNamespace(
        svc_obs=svc_obs,
        st_adata=SimpleNamespace(
            obsm={
                "major_type": pd.DataFrame(
                    [[0.5, 0.5]], index=["spot-a"], columns=["A", "B"]
                )
            }
        ),
    )

    _kernel(seed=7, cell_type_col="major_type").run(sc_svc)

    assert sc_svc.svc_obs["cell_type"].value_counts().to_dict() == {"A": 1, "B": 1}


@pytest.mark.parametrize(
    ("pm", "message"),
    [
        (pd.DataFrame([[1.0, 0.0]], index=["c1"], columns=["A", "B"]), "cell ID"),
        (pd.DataFrame([[1.0], [0.0]], index=["c1", "c2"], columns=["A"]), "classes"),
    ],
)
def test_pm_allocation_rejects_missing_active_axes(pm, message):
    svc_obs = _svc({"spot-a": ["c1", "c2"]})

    with pytest.raises(ValueError, match=message):
        _kernel(pm=pm).assign_cell_types(svc_obs, _quota())

    assert "cell_type" not in svc_obs


def test_pm_allocation_rejects_extra_patient_level_rows():
    svc_obs = _svc({"spot-a": ["c1", "c2"]})
    pm = pd.DataFrame(
        [[1.0, 0.0], [0.0, 1.0], [0.2, 0.8]],
        index=["c1", "c2", "unrelated-case-cell"],
        columns=["A", "B"],
    )
    with pytest.raises(ValueError, match="extra"):
        _kernel(pm=pm).assign_cell_types(svc_obs, _quota())


def test_pm_allocation_rejects_extra_global_classes():
    svc_obs = _svc({"spot-a": ["c1", "c2"]})
    pm = pd.DataFrame(
        [[1.0, 0.0, 0.3], [0.0, 1.0, 0.7]],
        index=["c1", "c2"],
        columns=["A", "B", "unrelated-reference-class"],
    )
    with pytest.raises(ValueError, match="extra"):
        _kernel(pm=pm).assign_cell_types(svc_obs, _quota())


@pytest.mark.parametrize(
    "cell_ids",
    [
        ["c1", "c1"],
        [1, "1"],
    ],
)
def test_pm_contract_rejects_nonunique_svc_cell_ids_after_string_conversion(cell_ids):
    svc_obs = _svc({"spot-a": cell_ids})
    pm = pd.DataFrame([[1.0, 0.0], [0.0, 1.0]], index=["c1", "c2"], columns=["A", "B"])

    with pytest.raises(ValueError, match="SVC.*cell_id"):
        _kernel(pm=pm).assign_cell_types(svc_obs, _quota())


@pytest.mark.parametrize(
    "invalid_quota",
    [
        np.array([[1, 1]]),
        pd.DataFrame([[1, 1], [1, 1]], index=["spot-a", "spot-a"], columns=["A", "B"]),
        pd.DataFrame([[1, 1]], index=["spot-a"], columns=["A", "A"]),
        pd.DataFrame([[1, 1]], index=["spot-b"], columns=["A", "B"]),
        pd.DataFrame([[1, 1]], index=["spot-a"], columns=["A/B", "A_B"]),
        pd.DataFrame([["not-a-number", 1]], index=["spot-a"], columns=["A", "B"]),
        pd.DataFrame([[np.inf, 1]], index=["spot-a"], columns=["A", "B"]),
        pd.DataFrame([[-1, 3]], index=["spot-a"], columns=["A", "B"]),
        pd.DataFrame([[0.5, 1.5]], index=["spot-a"], columns=["A", "B"]),
        pd.DataFrame([[1, 0]], index=["spot-a"], columns=["A", "B"]),
    ],
)
def test_assignment_rejects_invalid_quota_without_mutating_svc(invalid_quota):
    svc_obs = _svc({"spot-a": ["c1", "c2"]})
    original = svc_obs.copy(deep=True)

    with pytest.raises((TypeError, ValueError)):
        _kernel().assign_cell_types_random(svc_obs, invalid_quota)

    pd.testing.assert_frame_equal(svc_obs, original)


def test_invalid_contributions_fail_before_run_writes_assignment():
    svc_obs = _svc({"spot-a": ["c1", "c2"]})
    sc_svc = SimpleNamespace(
        svc_obs=svc_obs,
        st_adata=SimpleNamespace(
            obsm={"Level1": pd.DataFrame([[0.4, 0.4]], index=["spot-a"], columns=["A", "B"])}
        ),
    )

    with pytest.raises(ValueError, match="sum"):
        _kernel().run(sc_svc)

    assert sc_svc.svc_obs is svc_obs
    assert "cell_type" not in sc_svc.svc_obs


def test_random_assignment_is_seeded_once_across_spots_and_preserves_exact_composition():
    svc_obs = _svc(
        {
            "spot-a": [f"a{i}" for i in range(6)],
            "spot-b": [f"b{i}" for i in range(6)],
        }
    )
    quota = _quota(
        index=("spot-b", "spot-a"),
        columns=("A", "B", "C"),
        values=((2, 2, 2), (1, 2, 3)),
    )

    first = _kernel(seed=17).assign_cell_types_random(svc_obs, quota)
    repeated = _kernel(seed=17).assign_cell_types_random(svc_obs, quota)
    different = _kernel(seed=18).assign_cell_types_random(svc_obs, quota)

    assert first["cell_type"].tolist() == repeated["cell_type"].tolist()
    assert first["cell_type"].tolist() != different["cell_type"].tolist()
    observed = first.groupby("spot_name", sort=False)["cell_type"].value_counts().unstack(fill_value=0)
    pd.testing.assert_frame_equal(
        observed.loc[quota.index, quota.columns].astype(np.int64),
        quota,
        check_names=False,
    )


def test_single_cell_random_assignment_obeys_quota_instead_of_top_type_shortcut():
    svc_obs = _svc({"spot-a": ["c1"]})
    quota = _quota(values=((0, 1),))

    result = _kernel(seed=3).assign_cell_types_random(svc_obs, quota)

    assert result.loc[0, "cell_type"] == "B"


def test_random_assignment_writes_by_row_position_when_dataframe_index_is_not_unique():
    svc_obs = _svc({"spot-a": ["c1", "c2"]})
    svc_obs.index = [7, 7]

    result = _kernel(seed=3).assign_cell_types_random(svc_obs, _quota())

    assert result["cell_type"].value_counts().to_dict() == {"A": 1, "B": 1}


@pytest.fixture
def adapters(monkeypatch):
    scanpy = types.ModuleType("scanpy")
    scanpy.pp = SimpleNamespace()
    scanpy.pl = SimpleNamespace()
    scanpy.tl = SimpleNamespace()
    monkeypatch.setitem(sys.modules, "scanpy", scanpy)
    monkeypatch.delitem(sys.modules, "revise.backend.adapters", raising=False)
    return importlib.import_module("revise.backend.adapters")


def test_true_cell_type_keeps_randomly_allocated_cells_as_unknown(adapters):
    meta = importlib.import_module("revise.backend.ops.meta")
    svc_obs = pd.DataFrame(
        {
            "cell_id": ["known", "synthetic"],
            "x": [1.0, 9.0],
            "y": [2.0, 8.0],
        }
    )
    gt = SimpleNamespace(
        obs=pd.DataFrame(
            {
                "cell_id": ["known"],
                "clusters": pd.Categorical(["T"]),
                "x": [3.0],
                "y": [4.0],
            }
        )
    )

    result = meta.get_true_cell_type(svc_obs, gt)

    assert result["true_cell_type"].tolist() == ["T", "Unknown"]
    assert result.loc[0, ["x", "y"]].tolist() == [3.0, 4.0]
    assert result.loc[1, ["x", "y"]].tolist() == [9.0, 8.0]


def test_true_cell_type_default_falls_back_to_historical_level1(adapters):
    meta = importlib.import_module("revise.backend.ops.meta")
    svc_obs = pd.DataFrame(
        {"cell_id": ["known"], "x": [1.0], "y": [2.0]}
    )
    gt = SimpleNamespace(
        obs=pd.DataFrame(
            {
                "cell_id": ["known"],
                "Level1": ["T"],
                "x": [3.0],
                "y": [4.0],
            }
        )
    )

    result = meta.get_true_cell_type(svc_obs, gt)

    assert result["true_cell_type"].tolist() == ["T"]


def test_true_cell_type_explicit_label_key_does_not_fallback(adapters):
    meta = importlib.import_module("revise.backend.ops.meta")
    svc_obs = pd.DataFrame(
        {"cell_id": ["known"], "x": [1.0], "y": [2.0]}
    )
    gt = SimpleNamespace(
        obs=pd.DataFrame(
            {
                "cell_id": ["known"],
                "Level1": ["T"],
                "x": [3.0],
                "y": [4.0],
            }
        )
    )

    with pytest.raises(KeyError, match="clusters"):
        meta.get_true_cell_type(svc_obs, gt, label_key="clusters")


@pytest.mark.parametrize(
    ("strategy_name", "profile", "runner_module", "runner_class"),
    [
        ("ScSvcSrApplicationStrategy", "application_sc_sr", "sc_svc_sr_application", "ScSVCSr"),
        ("ScSvcSrBenchmarkStrategy", "benchmark_sr_batch", "sc_svc_sr_benchmark", "ScSVCSr"),
    ],
)
def test_sr_adapters_propagate_runtime_seed_to_assignment_config(
    adapters,
    monkeypatch,
    tmp_path,
    strategy_name,
    profile,
    runner_module,
    runner_class,
):
    merged = _merged(profile)
    merged["io"]["data_root"] = str(tmp_path)
    merged["io"]["output_root"] = str(tmp_path)
    runtime = dict(merged["runtime"])
    runtime["seed"] = 731

    runner_stub = types.ModuleType(f"revise.backend.runners.{runner_module}")
    setattr(runner_stub, runner_class, object)
    monkeypatch.setitem(sys.modules, f"revise.backend.runners.{runner_module}", runner_stub)

    captured = {}

    conf_type = (
        adapters.ApplicationScSrConf
        if profile == "application_sc_sr"
        else adapters.BenchmarkSrConf
    )

    def capture_conf(**kwargs):
        conf = conf_type(**kwargs)
        captured["conf"] = conf
        return conf

    class StopAfterConfig(Exception):
        pass

    def stop_before_io(_ctx):
        raise StopAfterConfig

    monkeypatch.setattr(adapters, conf_type.__name__, capture_conf)
    monkeypatch.setattr(adapters, "_input_service", stop_before_io)
    ctx = SimpleNamespace(
        merged_config=merged,
        io=merged["io"],
        columns=merged["columns"],
        runtime=runtime,
        route_key=(
            f"{runtime['mode']}:"
            f"{runtime.get('application_route') or runtime.get('confounding')}"
        ),
        run_dir=tmp_path,
        logger=logging.getLogger(f"test-{strategy_name}"),
        compatibility_mode=False,
        pm_on_cell=pd.DataFrame([[1.0]], index=["c1"], columns=["A"]),
    )

    with pytest.raises(StopAfterConfig):
        getattr(adapters, strategy_name)().prepare_context(ctx)

    assert captured["conf"].sr_assignment_seed == 731
    assert captured["conf"].pm_on_cell is ctx.pm_on_cell


def test_sr_benchmark_subsample_restricts_pm_to_active_cells(
    adapters,
    monkeypatch,
    tmp_path,
):
    merged = _merged(
        "benchmark_sr_batch",
        runtime_overrides={"seed": 17},
        io_overrides={"sample_size": 1},
    )
    merged["io"]["data_root"] = str(tmp_path)
    merged["io"]["output_root"] = str(tmp_path / "output")

    st = AnnData(
        X=np.ones((2, 1)),
        obs=pd.DataFrame(index=["spot-1", "spot-2"]),
        var=pd.DataFrame(index=["g1"]),
    )
    st.uns["all_cells_in_spot"] = {
        "spot-1": ["c1", "c2"],
        "spot-2": ["c3"],
    }
    real = AnnData(
        X=np.ones((3, 1)),
        obs=pd.DataFrame(
            {"cell_id": ["c1", "c2", "c3"]},
            index=["row-1", "row-2", "row-3"],
        ),
        var=pd.DataFrame(index=["g1"]),
    )
    reference = AnnData(
        X=np.ones((2, 1)),
        obs=pd.DataFrame({"Level1": ["A", "B"]}, index=["r1", "r2"]),
        var=pd.DataFrame(index=["g1"]),
    )

    input_service = SimpleNamespace(
        read_st_adata=lambda _path: st.copy(),
        read_real_adata=lambda _path: real.copy(),
        read_sc_ref_adata=lambda _path: reference.copy(),
    )
    monkeypatch.setattr(adapters, "_input_service", lambda _ctx: input_service)
    monkeypatch.setattr(adapters, "ensure_all_cells_in_spot", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        adapters,
        "_subsample_obs",
        lambda adata, _size, _seed: (adata[["spot-1"], :].copy(), ["spot-1"]),
    )

    runner_module = types.ModuleType("revise.backend.runners.sc_svc_sr_benchmark")

    class Runner:
        def __init__(self, _st, _sc, conf, _real, _logger):
            self.conf = conf

    runner_module.ScSVCSr = Runner
    monkeypatch.setitem(
        sys.modules,
        "revise.backend.runners.sc_svc_sr_benchmark",
        runner_module,
    )
    pm = pd.DataFrame(
        [[1.0, 0.0], [0.5, 0.5], [0.0, 1.0]],
        index=["c1", "c2", "c3"],
        columns=["A", "B"],
    )
    ctx = SimpleNamespace(
        merged_config=merged,
        io=merged["io"],
        columns=merged["columns"],
        runtime=merged["runtime"],
        route_key="sim2real:spot_size",
        run_dir=tmp_path,
        logger=logging.getLogger("test-sr-subsampled-pm"),
        compatibility_mode=False,
        pm_on_cell=pm,
        input_specs=(
            SimpleNamespace(role="st", path="/resolved/st.h5ad"),
            SimpleNamespace(role="sc_ref", path="/resolved/sc_ref.h5ad"),
            SimpleNamespace(role="gt", path="/resolved/gt.h5ad"),
        ),
    )

    adapters.ScSvcSrBenchmarkStrategy().prepare_context(ctx)

    assert ctx.runner.conf.pm_on_cell.index.tolist() == ["c1", "c2"]
    assert ctx.real_st_adata.obs["cell_id"].tolist() == ["c1", "c2"]
    assert ctx.pm_on_cell.index.tolist() == ["c1", "c2", "c3"]


def test_sr_benchmark_derives_assignment_seed_from_process_rng_when_runtime_seed_is_none(
    adapters,
    monkeypatch,
    tmp_path,
):
    merged = _merged(
        "benchmark_sr_batch",
        runtime_overrides={"seed": None},
    )
    merged["io"]["data_root"] = str(tmp_path)
    merged["io"]["output_root"] = str(tmp_path)

    runner_stub = types.ModuleType("revise.backend.runners.sc_svc_sr_benchmark")
    runner_stub.ScSVCSr = object
    monkeypatch.setitem(
        sys.modules,
        "revise.backend.runners.sc_svc_sr_benchmark",
        runner_stub,
    )

    captured = {}

    conf_type = adapters.BenchmarkSrConf

    def capture_conf(**kwargs):
        conf = conf_type(**kwargs)
        captured["conf"] = conf
        return conf

    class StopAfterConfig(Exception):
        pass

    def stop_before_io(_ctx):
        raise StopAfterConfig

    monkeypatch.setattr(adapters, "BenchmarkSrConf", capture_conf)
    monkeypatch.setattr(adapters, "_input_service", stop_before_io)
    ctx = SimpleNamespace(
        merged_config=merged,
        io=merged["io"],
        columns=merged["columns"],
        runtime=merged["runtime"],
        route_key="sim2real:batch_effect",
        run_dir=tmp_path,
        logger=logging.getLogger("test-process-scope-sr-seed"),
        compatibility_mode=True,
    )

    np.random.seed(731)
    expected_seed = int(
        np.random.RandomState(731).randint(0, np.iinfo(np.int32).max)
    )
    with pytest.raises(StopAfterConfig):
        adapters.ScSvcSrBenchmarkStrategy().prepare_context(ctx)

    assert captured["conf"].sr_assignment_seed == expected_seed
