from __future__ import annotations

import ast
import importlib
import importlib.metadata
import logging
import sys
import types
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise.backend.ops.local_ot import solve_local_ot
from revise.recon.context import PipelineContext
from revise.svc import SVC


ROOT = Path(__file__).resolve().parents[1]
LOCAL_OT_CALLERS = {
    "revise/backend/kernels/local_anchoring.py": 1,
    "revise/backend/runners/sp_svc_application.py": 1,
    "revise/backend/runners/sp_svc_benchmark.py": 1,
    "revise/backend/runners/sc_svc_sr_application.py": 1,
    "revise/backend/runners/sc_svc_sr_benchmark.py": 1,
    "revise/backend/runners/sc_svc_impute_benchmark.py": 1,
}


def _context(tmp_path: Path, *, ga: str = "pot", lr: str = "tacco"):
    return PipelineContext(
        merged_config={
            "ot": {
                "ga": {"solver": ga},
                "lr": {"solver": lr},
            }
        },
        raw_config={},
        config_path="revise/revise.yaml",
        profile="application_sc",
        runtime={},
        route_key="sc_svc:segmentation",
        run_dir=tmp_path,
        logger=logging.getLogger("test-ot-events"),
    )


def _fake_ot(monkeypatch, coupling=None, error: Exception | None = None):
    module = types.ModuleType("ot")
    if coupling is None:
        coupling = np.array([[0.5, 0.0], [0.0, 0.5]])

    def solve(*args, **kwargs):
        if error is not None:
            raise error
        return np.asarray(coupling)

    module.unbalanced = SimpleNamespace(sinkhorn_unbalanced=solve)
    monkeypatch.setitem(sys.modules, "ot", module)
    return module


def _fake_tacco(monkeypatch, coupling=None, error: Exception | None = None):
    module = types.ModuleType("tacco")
    if coupling is None:
        coupling = np.array([[0.5, 0.0], [0.0, 0.5]])

    def solve(*args, **kwargs):
        if error is not None:
            raise error
        return np.asarray(coupling)

    module.__version__ = "0.5.0"
    module.utils = SimpleNamespace(solve_OT=solve)
    monkeypatch.setitem(sys.modules, "tacco", module)
    monkeypatch.setattr(
        importlib.metadata,
        "version",
        lambda package: "0.5.0" if package == "tacco" else "unknown",
    )
    return module


def test_pipeline_context_accepts_ordered_ga_and_repeated_lr_events(tmp_path):
    ctx = _context(tmp_path)

    assert ctx.ot_events == [
        {"phase": "ga", "solver": "pot", "status": "requested", "call": 0},
        {"phase": "lr", "solver": "tacco", "status": "requested", "call": 0},
    ]

    for event in (
        ("ga", "pot", "attempted"),
        ("ga", "pot", "completed"),
        ("lr", "tacco", "attempted"),
        ("lr", "tacco", "completed"),
        ("lr", "tacco", "attempted"),
        ("lr", "tacco", "completed"),
    ):
        ctx.record_ot_event(*event)

    assert ctx.ot_events == [
        {"phase": phase, "solver": solver, "status": status, "call": call}
        for phase, solver, status, call in (
            ("ga", "pot", "requested", 0),
            ("lr", "tacco", "requested", 0),
            ("ga", "pot", "attempted", 1),
            ("ga", "pot", "completed", 1),
            ("lr", "tacco", "attempted", 1),
            ("lr", "tacco", "completed", 1),
            ("lr", "tacco", "attempted", 2),
            ("lr", "tacco", "completed", 2),
        )
    ]


@pytest.mark.parametrize(
    "events, message",
    [
        ([("ga", "tacco", "attempted")], "configured solver"),
        ([("ga", "pot", "completed")], "attempted"),
        ([("ga", "pot", "requested")], "exactly once"),
        (
            [("lr", "tacco", "attempted"), ("lr", "tacco", "attempted")],
            "cannot repeat",
        ),
        (
            [
                ("ga", "pot", "attempted"),
                ("ga", "pot", "completed"),
                ("ga", "pot", "attempted"),
            ],
            "only once",
        ),
    ],
)
def test_pipeline_context_rejects_solver_mismatch_and_invalid_transition(
    tmp_path, events, message
):
    ctx = _context(tmp_path)

    with pytest.raises(ValueError, match=message):
        for event in events:
            ctx.record_ot_event(*event)


def test_local_pot_records_attempted_and_completed_without_importing_tacco(
    monkeypatch,
):
    _fake_ot(monkeypatch)
    real_import_module = importlib.import_module

    def guarded_import(name, package=None):
        if name == "tacco":
            raise AssertionError("POT path imported TACCO")
        return real_import_module(name, package)

    monkeypatch.setattr(importlib, "import_module", guarded_import)
    events = []

    coupling = solve_local_ot(
        [0.5, 0.5],
        [0.5, 0.5],
        [[0.0, 1.0], [1.0, 0.0]],
        method="pot",
        pot_reg=0.1,
        pot_reg_m=0.0,
        event_callback=lambda *event: events.append(event),
    )

    np.testing.assert_allclose(coupling, [[0.5, 0.0], [0.0, 0.5]])
    assert events == [("lr", "pot", "attempted"), ("lr", "pot", "completed")]


def test_local_tacco_records_attempted_and_completed(monkeypatch):
    _fake_tacco(monkeypatch)
    events = []

    solve_local_ot(
        [0.5, 0.5],
        [0.5, 0.5],
        [[0.0, 1.0], [1.0, 0.0]],
        method="tacco",
        event_callback=lambda *event: events.append(event),
    )

    assert events == [
        ("lr", "tacco", "attempted"),
        ("lr", "tacco", "completed"),
    ]


def test_local_failure_records_attempted_without_completed(monkeypatch):
    _fake_tacco(monkeypatch, error=RuntimeError("solver exploded"))
    events = []

    with pytest.raises(RuntimeError, match="solver exploded"):
        solve_local_ot(
            [0.5, 0.5],
            [0.5, 0.5],
            [[0.0, 1.0], [1.0, 0.0]],
            method="tacco",
            event_callback=lambda *event: events.append(event),
        )

    assert events == [("lr", "tacco", "attempted")]


def test_successful_events_are_persisted_to_svc_and_provenance_json(tmp_path):
    import json

    from revise.framework import REVISEPipeline

    ctx = _context(tmp_path, ga="pot", lr="pot")
    for event in (
        ("ga", "pot", "attempted"),
        ("ga", "pot", "completed"),
        ("lr", "pot", "attempted"),
        ("lr", "pot", "completed"),
    ):
        ctx.record_ot_event(*event)
    ctx.runner_config = SimpleNamespace(
        st_file_path=None,
        sc_ref_file_path=None,
        gt_svc_file_path=None,
    )
    ctx.svc = SVC(expr=None, spatial=None, svc_kind="sc")
    pipeline = REVISEPipeline.__new__(REVISEPipeline)

    pipeline._write_final_metadata(ctx)

    persisted = json.loads((tmp_path / "provenance.json").read_text())
    assert ctx.svc.provenance["ot_events"] == ctx.ot_events
    assert persisted["ot_events"] == ctx.ot_events
    assert "tacco" in persisted["packages"]


def test_local_empty_support_does_not_fabricate_branch_events(monkeypatch):
    events = []

    coupling = solve_local_ot(
        [0.5, 0.5],
        [0.5, 0.5],
        [[0.0, 1.0], [1.0, 0.0]],
        method="tacco",
        valid_support_mask=np.zeros((2, 2), dtype=bool),
        event_callback=lambda *event: events.append(event),
    )

    np.testing.assert_array_equal(coupling, np.zeros((2, 2)))
    assert events == []


def test_local_multiple_invocations_record_each_actual_call(monkeypatch):
    _fake_ot(monkeypatch)
    events = []
    kwargs = {
        "method": "pot",
        "pot_reg": 0.1,
        "pot_reg_m": 0.0,
        "event_callback": lambda *event: events.append(event),
    }

    for _ in range(2):
        solve_local_ot(
            [0.5, 0.5],
            [0.5, 0.5],
            [[0.0, 1.0], [1.0, 0.0]],
            **kwargs,
        )

    assert events == [
        ("lr", "pot", "attempted"),
        ("lr", "pot", "completed"),
        ("lr", "pot", "attempted"),
        ("lr", "pot", "completed"),
    ]


def test_lr_second_failed_call_keeps_distinct_attempted_call_id(tmp_path):
    ctx = _context(tmp_path)
    ctx.record_ot_event("lr", "tacco", "attempted")
    ctx.record_ot_event("lr", "tacco", "completed")
    ctx.record_ot_event("lr", "tacco", "attempted")

    assert ctx.ot_events[-2:] == [
        {"phase": "lr", "solver": "tacco", "status": "completed", "call": 1},
        {"phase": "lr", "solver": "tacco", "status": "attempted", "call": 2},
    ]
    assert not any(
        event["phase"] == "lr"
        and event["status"] == "completed"
        and event["call"] == 2
        for event in ctx.ot_events
    )


@pytest.mark.parametrize(
    "coupling, message",
    [
        (np.zeros((2, 2)), "positive total mass"),
        (np.array([[0.5, 0.0], [0.0, 0.0]]), "positive transported row mass"),
        (np.array([[0.5, 0.0], [0.5, 0.0]]), "positive transported column mass"),
    ],
)
def test_local_pot_rejects_unusable_finite_coupling_before_completed(
    monkeypatch, coupling, message
):
    _fake_ot(monkeypatch, coupling=coupling)
    events = []

    with pytest.raises(ValueError, match=message):
        solve_local_ot(
            [0.5, 0.5],
            [0.5, 0.5],
            [[0.0, 1.0], [1.0, 0.0]],
            method="pot",
            pot_reg=0.1,
            pot_reg_m=0.0,
            event_callback=lambda *event: events.append(event),
        )

    assert events == [("lr", "pot", "attempted")]


def test_missing_tacco_is_actionable_and_does_not_fallback(monkeypatch):
    monkeypatch.delitem(sys.modules, "tacco", raising=False)
    real_import_module = importlib.import_module

    def missing(name, package=None):
        if name == "tacco":
            raise ModuleNotFoundError("No module named 'tacco'", name="tacco")
        return real_import_module(name, package)

    monkeypatch.setattr(importlib, "import_module", missing)
    events = []

    with pytest.raises(
        ModuleNotFoundError,
        match=r'python -m pip install "tacco==0\.5\.0"',
    ):
        solve_local_ot(
            [0.5, 0.5],
            [0.5, 0.5],
            [[0.0, 1.0], [1.0, 0.0]],
            method="tacco",
            event_callback=lambda *event: events.append(event),
        )

    assert events == [("lr", "tacco", "attempted")]


def test_missing_tacco_transitive_dependency_is_not_reported_as_missing_tacco(
    monkeypatch,
):
    monkeypatch.delitem(sys.modules, "tacco", raising=False)

    def missing_dependency(name, package=None):
        raise ModuleNotFoundError(
            "No module named 'transitive_dep'", name="transitive_dep"
        )

    monkeypatch.setattr(importlib, "import_module", missing_dependency)

    with pytest.raises(
        ImportError,
        match=r'transitive_dep.*python -m pip install "tacco==0\.5\.0"',
    ) as caught:
        solve_local_ot(
            [0.5, 0.5],
            [0.5, 0.5],
            [[0.0, 1.0], [1.0, 0.0]],
            method="tacco",
        )

    assert "No module named 'tacco'" not in str(caught.value)


def test_unsupported_tacco_version_is_actionable_and_does_not_fallback(monkeypatch):
    module = types.ModuleType("tacco")
    module.utils = SimpleNamespace(solve_OT=lambda *args: None)
    monkeypatch.setitem(sys.modules, "tacco", module)
    monkeypatch.setattr(importlib.metadata, "version", lambda package: "0.5.1")
    events = []

    with pytest.raises(
        RuntimeError,
        match=r'requires tacco==0\.5\.0.*0\.5\.1.*python -m pip install',
    ):
        solve_local_ot(
            [0.5, 0.5],
            [0.5, 0.5],
            [[0.0, 1.0], [1.0, 0.0]],
            method="tacco",
            event_callback=lambda *event: events.append(event),
        )

    assert events == [("lr", "tacco", "attempted")]


def test_every_physical_local_ot_caller_passes_the_explicit_event_callback():
    found = {}
    for relative, expected_count in LOCAL_OT_CALLERS.items():
        tree = ast.parse((ROOT / relative).read_text())
        calls = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "solve_local_ot"
        ]
        found[relative] = len(calls)
        assert len(calls) == expected_count
        for call in calls:
            callback = next(
                (kw.value for kw in call.keywords if kw.arg == "event_callback"),
                None,
            )
            assert callback is not None, f"{relative} does not wire event_callback"
            assert ast.unparse(callback) == (
                "getattr(self.config, 'ot_event_callback', None)"
            )

    assert found == LOCAL_OT_CALLERS


def test_runner_strategy_attaches_production_ot_recorder_and_persists_calls(
    monkeypatch, tmp_path
):
    import json

    scanpy = types.ModuleType("scanpy")
    scanpy.pp = SimpleNamespace()
    scanpy.pl = SimpleNamespace()
    scanpy.tl = SimpleNamespace()
    monkeypatch.setitem(sys.modules, "scanpy", scanpy)
    import revise.backend.kernels as kernels
    from revise.backend.adapters import RunnerBackedStrategy
    from revise.framework import REVISEPipeline

    ctx = _context(tmp_path, ga="pot", lr="pot")
    ctx.runner_config = SimpleNamespace(
        annotate_mode="pot",
        rec_ot_method="pot",
        st_file_path=None,
        sc_ref_file_path=None,
        gt_svc_file_path=None,
    )
    ctx.runner = SimpleNamespace(st_adata=object(), sc_ref_adata=object())

    class FakeGlobalKernel:
        def run(self, target, reference, **kwargs):
            callback = ctx.runner_config.ot_event_callback
            callback("ga", "pot", "attempted")
            callback("ga", "pot", "completed")
            return target

    class ConcreteStrategy(RunnerBackedStrategy):
        def prepare_context(self, ctx):
            raise NotImplementedError

        def finalize_svc(self, ctx):
            raise NotImplementedError

    monkeypatch.setattr(kernels, "build_kernel", lambda *args, **kwargs: FakeGlobalKernel())
    ConcreteStrategy().global_anchoring(ctx)
    _fake_ot(monkeypatch)
    solve_local_ot(
        [0.5, 0.5],
        [0.5, 0.5],
        [[0.0, 1.0], [1.0, 0.0]],
        method="pot",
        pot_reg=0.1,
        pot_reg_m=0.0,
        event_callback=ctx.runner_config.ot_event_callback,
    )
    ctx.svc = SVC(expr=None, spatial=None, svc_kind="sc")

    REVISEPipeline.__new__(REVISEPipeline)._write_final_metadata(ctx)

    persisted = json.loads((tmp_path / "provenance.json").read_text())
    assert persisted["ot_events"] == [
        {"phase": "ga", "solver": "pot", "status": "requested", "call": 0},
        {"phase": "lr", "solver": "pot", "status": "requested", "call": 0},
        {"phase": "ga", "solver": "pot", "status": "attempted", "call": 1},
        {"phase": "ga", "solver": "pot", "status": "completed", "call": 1},
        {"phase": "lr", "solver": "pot", "status": "attempted", "call": 1},
        {"phase": "lr", "solver": "pot", "status": "completed", "call": 1},
    ]


def test_ci_has_mandatory_exact_tacco_smoke_job():
    ci = (ROOT / ".github/workflows/ci.yml").read_text()

    assert "tacco-smoke:" in ci
    assert "needs: test" in ci
    assert "tacco==0.5.0" in ci
    assert "POT==0.9.5" in ci
    assert "tests/integration/test_tacco_solver_smoke.py" in ci
    smoke = (ROOT / "tests/integration/test_tacco_solver_smoke.py").read_text()
    assert "importorskip" not in smoke
    assert "skipif" not in smoke
    assert "pytestmark" not in smoke
    assert "REVISE_TACCO_SMOKE" not in ci
    assert smoke.count("def test_real_tacco_050_completes_global_and_local_smoke") == 1
