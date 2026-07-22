from __future__ import annotations

import importlib
import itertools
import logging
import sys
import types
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from revise.backend.kernels.spot_sr import SpotSrKernel
from revise.config import load_raw_config, merge_unified_config


CONFIG_PATH = Path(__file__).parents[2] / "revise" / "revise.yaml"


def _kernel(*, pm=None, seed=42, cell_type_col="Level1"):
    config = SimpleNamespace(
        pm_on_cell_file="/path/that/does/not/exist.csv",
        svc_completeness=True,
        sr_assignment_seed=seed,
        cell_type_col=cell_type_col,
    )
    kernel = SpotSrKernel(config, logging.getLogger("test-spot-sr-assignment"))
    kernel.pm_on_cell = pm
    return kernel


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
        index=[1, 2],
        columns=["T", "Mono/Macro"],
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
        (pd.DataFrame([[1.0, 0.0]], index=["c1"], columns=["A", "B"]), "cell IDs"),
        (
            pd.DataFrame([[1.0, 0.0], [0.0, 1.0], [0.2, 0.8]], index=["c1", "c2", "c3"], columns=["A", "B"]),
            "cell IDs",
        ),
        (pd.DataFrame([[1.0, 0.0], [0.0, 1.0]], index=["c1", "c1"], columns=["A", "B"]), "duplicate row"),
        (pd.DataFrame([[1.0, 0.0], [0.0, 1.0]], index=[1, "1"], columns=["A", "B"]), "string conversion"),
        (pd.DataFrame([[1.0], [0.0]], index=["c1", "c2"], columns=["A"]), "classes"),
        (pd.DataFrame([[1.0, 0.0, 2.0], [0.0, 1.0, 2.0]], index=["c1", "c2"], columns=["A", "B", "C"]), "classes"),
        (pd.DataFrame([[1.0, 0.0], [0.0, 1.0]], index=["c1", "c2"], columns=["A", "A"]), "duplicate column"),
        (
            pd.DataFrame([[1.0, 0.0], [0.0, 1.0]], index=["c1", "c2"], columns=["A/B", "A_B"]),
            "normalization",
        ),
        (pd.DataFrame([[np.nan, 0.0], [0.0, 1.0]], index=["c1", "c2"], columns=["A", "B"]), "finite"),
        (pd.DataFrame([[np.inf, 0.0], [0.0, 1.0]], index=["c1", "c2"], columns=["A", "B"]), "finite"),
    ],
)
def test_pm_contract_rejects_misalignment_collisions_and_nonfinite_values(pm, message):
    svc_obs = _svc({"spot-a": ["c1", "c2"]})

    with pytest.raises(ValueError, match=message):
        _kernel(pm=pm).assign_cell_types(svc_obs, _quota())

    assert "cell_type" not in svc_obs


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
                "clusters": ["T"],
                "x": [3.0],
                "y": [4.0],
            }
        )
    )

    result = meta.get_true_cell_type(svc_obs, gt)

    assert result["true_cell_type"].tolist() == ["T", "Unknown"]
    assert result.loc[0, ["x", "y"]].tolist() == [3.0, 4.0]
    assert result.loc[1, ["x", "y"]].tolist() == [9.0, 8.0]


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
    raw = load_raw_config(CONFIG_PATH)
    merged = merge_unified_config(
        raw_config=raw,
        profile=profile,
        runtime_overrides={},
        io_overrides={},
        set_overrides=(),
    )
    merged["io"]["data_root"] = str(tmp_path)
    merged["io"]["output_root"] = str(tmp_path)
    runtime = dict(merged["runtime"])
    runtime["seed"] = 731

    runner_stub = types.ModuleType(f"revise.backend.runners.{runner_module}")
    setattr(runner_stub, runner_class, object)
    monkeypatch.setitem(sys.modules, f"revise.backend.runners.{runner_module}", runner_stub)

    captured = {}

    def capture_conf(conf, _cfg):
        captured["conf"] = conf

    class StopAfterConfig(Exception):
        pass

    def stop_before_io(_ctx):
        raise StopAfterConfig

    monkeypatch.setattr(adapters, "_attach_posterior_conditioning_conf", capture_conf)
    monkeypatch.setattr(adapters, "_input_service", stop_before_io)
    ctx = SimpleNamespace(
        merged_config=merged,
        io=merged["io"],
        columns=merged["columns"],
        runtime=runtime,
        route_key=f"{runtime['platform']}:{runtime['confounding']}",
        run_dir=tmp_path,
        logger=logging.getLogger(f"test-{strategy_name}"),
        compatibility_mode=False,
    )

    with pytest.raises(StopAfterConfig):
        getattr(adapters, strategy_name)().prepare_context(ctx)

    assert captured["conf"].sr_assignment_seed == 731
