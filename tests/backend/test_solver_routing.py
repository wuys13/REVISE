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
from anndata import AnnData, concat as anndata_concat
from scipy import sparse

from revise.backend.kernels.ot import OTKernel
from revise.recon.context import PipelineContext
from revise.svc import SVC


ROOT = Path(__file__).resolve().parents[2]
LOCAL_OT_CALLERS = {
    "revise/backend/runners/sp_svc_application.py": 1,
    "revise/backend/runners/sp_svc_benchmark.py": 1,
    "revise/backend/runners/sc_svc_sr_application.py": 1,
    "revise/backend/runners/sc_svc_sr_benchmark.py": 1,
    "revise/backend/runners/sc_svc_impute_benchmark.py": 1,
}
GUIDED_LOCAL_OT_CALLERS = {
    "revise/backend/runners/sp_svc_application.py",
    "revise/backend/runners/sp_svc_benchmark.py",
    "revise/backend/runners/sc_svc_sr_application.py",
    "revise/backend/runners/sc_svc_sr_benchmark.py",
}


def _load_runner_method(relative_path, class_name, method_name, namespace):
    tree = ast.parse((ROOT / relative_path).read_text())
    class_node = next(
        node
        for node in tree.body
        if isinstance(node, ast.ClassDef) and node.name == class_name
    )
    method_node = next(
        node
        for node in class_node.body
        if isinstance(node, ast.FunctionDef) and node.name == method_name
    )
    exec(
        compile(
            ast.fix_missing_locations(ast.Module(body=[method_node], type_ignores=[])),
            ROOT / relative_path,
            "exec",
        ),
        namespace,
    )
    return namespace[method_name]


def _context(
    tmp_path: Path,
    *,
    ga: str = "pot",
    lr: str = "tacco",
    task: str | None = None,
    strength: float | None = None,
):
    merged_config = {
        "ot": {
            "ga": {"solver": ga},
            "lr": {"solver": lr},
        }
    }
    if strength is not None:
        merged_config["local_refinement"] = {"strength": strength}
    return PipelineContext(
        merged_config=merged_config,
        profile="application_sc",
        runtime={"task": task} if task is not None else {},
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


def test_local_pot_routes_without_importing_tacco(monkeypatch):
    _fake_ot(monkeypatch)
    real_import_module = importlib.import_module

    def guarded_import(name, package=None):
        if name == "tacco":
            raise AssertionError("POT path imported TACCO")
        return real_import_module(name, package)

    monkeypatch.setattr(importlib, "import_module", guarded_import)
    coupling = OTKernel.couple(
        [0.5, 0.5],
        [0.5, 0.5],
        [[0.0, 1.0], [1.0, 0.0]],
        method="pot",
        pot_reg=0.1,
        pot_reg_m=0.0,
    )

    np.testing.assert_allclose(coupling, [[0.5, 0.0], [0.0, 0.5]])


def test_manifest_has_ot_config_but_no_ot_events(tmp_path):
    import json

    from revise.framework import REVISEPipeline

    ctx = _context(tmp_path, ga="pot", lr="pot")
    ctx.runner_config = SimpleNamespace(
        st_file_path=None,
        sc_ref_file_path=None,
        gt_svc_file_path=None,
    )
    ctx.svc = SVC(expr=None, spatial=None, svc_kind="sc")
    pipeline = REVISEPipeline.__new__(REVISEPipeline)

    pipeline._write_final_metadata(ctx)

    persisted = json.loads((tmp_path / "provenance.json").read_text())
    assert persisted["ot_config"] == ctx.merged_config["ot"]
    assert "ot_events" not in persisted
    assert "ot_events" not in ctx.svc.provenance


def test_local_empty_support_returns_zero_coupling():
    coupling = OTKernel.couple(
        [0.5, 0.5],
        [0.5, 0.5],
        [[0.0, 1.0], [1.0, 0.0]],
        method="tacco",
        valid_support_mask=np.zeros((2, 2), dtype=bool),
    )

    np.testing.assert_array_equal(coupling, np.zeros((2, 2)))


def test_sr_benchmark_singleton_short_circuits_before_assignment_validation():
    apply_graph_aggregation = _load_runner_method(
        "revise/backend/runners/sc_svc_sr_benchmark.py",
        "ScSVCSr",
        "_apply_graph_aggregation",
        {"np": np},
    )
    runner = SimpleNamespace(
        config=SimpleNamespace(),
        logger=SimpleNamespace(info=lambda *_args, **_kwargs: None),
        _get_graphagg_posterior_matrix=lambda: pytest.fail(
            "singleton has no local problem to validate"
        ),
    )

    result = types.MethodType(apply_graph_aggregation, runner)(
        np.array([[3.0, 7.0]])
    )

    np.testing.assert_array_equal(result, [[3.0, 7.0]])


def test_sp_benchmark_validates_ga_before_insufficient_units_short_circuit():
    assignment_loaded = {"value": False}

    def load_assignment(_adata, *, key, expected_categories):
        assert key == "major_type"
        assert expected_categories.equals(pd.Index(["A"]))
        assignment_loaded["value"] = True

    local_refinement = _load_runner_method(
        "revise/backend/runners/sp_svc_benchmark.py",
        "SpSVC",
        "local_refinement",
        {
            "np": np,
            "pd": pd,
            "scipy": SimpleNamespace(sparse=sparse),
            "sc": SimpleNamespace(concat=anndata_concat),
            "tqdm": lambda values, **_kwargs: values,
            "global_assignment_from_adata": load_assignment,
        },
    )
    st = AnnData(
        X=sparse.csr_matrix(np.ones((2, 2))),
        obs=pd.DataFrame(
            {"major_type": ["A", "A"], "no_effect": [False, False]},
            index=["cell-1", "cell-2"],
        ),
    )
    runner = SimpleNamespace(
        st_adata=st,
        sc_ref_adata=AnnData(
            X=sparse.csr_matrix(np.ones((1, 2))),
            obs=pd.DataFrame({"major_type": ["A"]}, index=["ref-1"]),
        ),
        config=SimpleNamespace(
            cell_type_col="major_type",
            local_refinement_strength=0.2,
        ),
        logger=SimpleNamespace(
            info=lambda *_args, **_kwargs: None,
            warning=lambda *_args, **_kwargs: None,
        ),
        svc={},
    )

    applied = types.MethodType(local_refinement, runner)()

    assert assignment_loaded["value"] is True
    assert applied is False
    assert runner.svc["sp_svc"].n_obs == 2


@pytest.mark.parametrize(
    "coupling, message",
    [
        (np.zeros((2, 2)), "positive total mass"),
        (np.array([[0.5, 0.0], [0.0, 0.0]]), "positive transported row mass"),
        (np.array([[0.5, 0.0], [0.5, 0.0]]), "positive transported column mass"),
    ],
)
def test_local_pot_rejects_unusable_finite_coupling(
    monkeypatch, coupling, message
):
    _fake_ot(monkeypatch, coupling=coupling)
    with pytest.raises(ValueError, match=message):
        OTKernel.couple(
            [0.5, 0.5],
            [0.5, 0.5],
            [[0.0, 1.0], [1.0, 0.0]],
            method="pot",
            pot_reg=0.1,
            pot_reg_m=0.0,
        )


def test_missing_tacco_is_actionable_and_does_not_fallback(monkeypatch):
    monkeypatch.delitem(sys.modules, "tacco", raising=False)
    real_import_module = importlib.import_module

    def missing(name, package=None):
        if name == "tacco":
            raise ModuleNotFoundError("No module named 'tacco'", name="tacco")
        return real_import_module(name, package)

    monkeypatch.setattr(importlib, "import_module", missing)
    with pytest.raises(
        ModuleNotFoundError,
        match=r'python -m pip install "tacco==0\.5\.0"',
    ) as caught:
        OTKernel.couple(
            [0.5, 0.5],
            [0.5, 0.5],
            [[0.0, 1.0], [1.0, 0.0]],
            method="tacco",
        )

    assert "algorithm.ot_method: pot" in str(caught.value)
    assert "does not fall back automatically" in str(caught.value)


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
        OTKernel.couple(
            [0.5, 0.5],
            [0.5, 0.5],
            [[0.0, 1.0], [1.0, 0.0]],
            method="tacco",
        )

    assert "No module named 'tacco'" not in str(caught.value)
    assert "algorithm.ot_method: pot" in str(caught.value)


def test_unsupported_tacco_version_is_actionable_and_does_not_fallback(monkeypatch):
    module = types.ModuleType("tacco")
    module.utils = SimpleNamespace(solve_OT=lambda *args: None)
    monkeypatch.setitem(sys.modules, "tacco", module)
    monkeypatch.setattr(importlib.metadata, "version", lambda package: "0.5.1")
    with pytest.raises(
        RuntimeError,
        match=r'requires tacco==0\.5\.0.*0\.5\.1.*python -m pip install',
    ) as caught:
        OTKernel.couple(
            [0.5, 0.5],
            [0.5, 0.5],
            [[0.0, 1.0], [1.0, 0.0]],
            method="tacco",
        )

    assert "algorithm.ot_method: pot" in str(caught.value)


def test_physical_local_ot_callers_use_ot_kernel_without_event_callback_wiring():
    found = {}
    for relative, expected_count in LOCAL_OT_CALLERS.items():
        source = (ROOT / relative).read_text()
        tree = ast.parse(source)
        calls = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id == "OTKernel"
            and node.func.attr == "couple"
        ]
        found[relative] = len(calls)
        assert len(calls) == expected_count
        assert "solve_local_ot" not in source
        assert all(
            keyword.arg != "event_callback"
            for call in calls
            for keyword in call.keywords
        )
        assert source.count('"local_refinement_applied_callback"') == (
            1 if relative in GUIDED_LOCAL_OT_CALLERS else 0
        )

    assert found == LOCAL_OT_CALLERS


def test_annotation_solver_and_ga_lr_boundaries_are_static():
    solver_patterns = (
        "tacco.tl.annotate",
        "tc.tl.annotate",
        "tacco.utils.solve_OT",
        "ot.unbalanced.sinkhorn_unbalanced",
    )
    offenders = []
    for path in (ROOT / "revise").rglob("*.py"):
        if path == ROOT / "revise/backend/kernels/ot.py":
            continue
        source = path.read_text()
        if any(pattern in source for pattern in solver_patterns):
            offenders.append(str(path.relative_to(ROOT)))
    assert offenders == []

    local_source = (ROOT / "revise/backend/kernels/local_anchoring.py").read_text()
    assert "GlobalAnchoringKernel" not in local_source
    for relative in (
        "revise/backend/runners/base_svc.py",
        "revise/backend/runners/application_svc.py",
        "revise/backend/runners/benchmark_svc.py",
    ):
        source = (ROOT / relative).read_text()
        assert "GlobalAnchoringKernel" not in source
        assert "annotate_method" not in source
        assert "def global_anchoring" not in source
    assert not (ROOT / "revise/backend/runners/base_svc_anchor.py").exists()


def test_runner_strategy_records_completed_conditioning_before_later_failure(
    monkeypatch, tmp_path
):
    scanpy = types.ModuleType("scanpy")
    scanpy.pp = SimpleNamespace()
    scanpy.pl = SimpleNamespace()
    scanpy.tl = SimpleNamespace()
    monkeypatch.setitem(sys.modules, "scanpy", scanpy)
    import revise.backend.kernels as kernels
    from revise.backend.adapters import RunnerBackedStrategy

    ctx = _context(
        tmp_path,
        ga="pot",
        lr="pot",
        task="sp_svc",
        strength=0.2,
    )
    ctx.runner_config = SimpleNamespace()

    class FakeGlobalKernel:
        def run(self, target, reference, **kwargs):
            return target

    class ConcreteStrategy(RunnerBackedStrategy):
        def prepare_context(self, ctx):
            raise NotImplementedError

        def finalize_svc(self, ctx):
            raise NotImplementedError

    ctx.runner = SimpleNamespace(st_adata=object(), sc_ref_adata=object())
    monkeypatch.setattr(
        kernels,
        "build_kernel",
        lambda *args, **kwargs: FakeGlobalKernel(),
    )
    strategy = ConcreteStrategy()
    strategy.global_anchoring(ctx)

    def local_refinement():
        ctx.runner_config.local_refinement_applied_callback()
        raise RuntimeError("later cell type failed")

    ctx.runner.local_refinement = local_refinement

    with pytest.raises(RuntimeError, match="later cell type failed"):
        strategy.solve_ot(ctx)

    assert ctx.local_refinement_record["applied"] is True


@pytest.mark.parametrize(
    ("task", "strength", "callback_expected"),
    [("sp_svc", 0.0, True), ("sc_svc", None, False), ("impute", None, False)],
)
def test_refinement_callback_is_independent_of_strength(
    monkeypatch, tmp_path, task, strength, callback_expected
):
    import revise.backend.kernels as kernels
    from revise.backend.adapters import RunnerBackedStrategy

    ctx = _context(tmp_path, task=task, strength=strength)
    ctx.runner_config = SimpleNamespace()
    ctx.runner = SimpleNamespace(st_adata=object(), sc_ref_adata=object())

    class FakeGlobalKernel:
        def run(self, target, reference, **kwargs):
            return target

    class ConcreteStrategy(RunnerBackedStrategy):
        def prepare_context(self, ctx):
            raise NotImplementedError

        def finalize_svc(self, ctx):
            raise NotImplementedError

    monkeypatch.setattr(
        kernels,
        "build_kernel",
        lambda *args, **kwargs: FakeGlobalKernel(),
    )

    ConcreteStrategy().global_anchoring(ctx)

    assert hasattr(ctx.runner_config, "local_refinement_applied_callback") is callback_expected
    assert ctx.local_refinement_record["applied"] is False


def test_ci_has_mandatory_exact_tacco_smoke_job():
    ci = (ROOT / ".github/workflows/ci.yml").read_text()

    assert "tacco-smoke:" in ci
    assert "needs: test" in ci
    assert "tacco==0.5.0" in ci
    assert "POT==0.9.5" in ci
    assert "tests/integration/solvers/test_tacco_solver_smoke.py" in ci
    smoke = (ROOT / "tests/integration/solvers/test_tacco_solver_smoke.py").read_text()
    assert "importorskip" not in smoke
    assert "skipif" not in smoke
    assert "pytestmark" not in smoke
    assert "REVISE_TACCO_SMOKE" not in ci
    assert smoke.count("def test_real_tacco_050_completes_global_and_local_smoke") == 1
