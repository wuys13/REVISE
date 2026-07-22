from __future__ import annotations

import ast
import importlib
import inspect
import os
from pathlib import Path
import subprocess
import sys
import textwrap
import types
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse


RUNNER_MODULE_NAMES = (
    "revise.backend.runners.application_svc",
    "revise.backend.runners.benchmark_svc",
    "revise.backend.runners.sc_svc_sr_application",
    "revise.backend.runners.sc_svc_sr_benchmark",
    "revise.backend.runners.sp_svc_application",
)
ISOLATED_MODULE_NAMES = (
    "scanpy",
    "matplotlib",
    "matplotlib.pyplot",
    "revise.analysis",
    "revise.analysis.metrics",
    "revise.backend.adapters",
    "revise.backend.kernels",
    "revise.backend.kernels.global_anchoring",
    "revise.backend.kernels.graph_aggregate",
    "revise.backend.kernels.spot_sr",
    "revise.backend.ops.distance",
    "revise.backend.ops.meta",
    "revise.backend.ops.shaver",
    "revise.backend.ops.topology",
    "revise.backend.runners.base_svc_anchor",
    *RUNNER_MODULE_NAMES,
)
_MISSING = object()


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


def _install_import_stubs() -> None:
    scanpy = types.ModuleType("scanpy")
    scanpy.pp = SimpleNamespace()
    scanpy.pl = SimpleNamespace()
    scanpy.tl = SimpleNamespace()
    scanpy._revise_test_stub = True
    sys.modules["scanpy"] = scanpy

    distance = types.ModuleType("revise.backend.ops.distance")
    distance.bhattacharyya_distance = lambda *args, **kwargs: None
    distance.similarity_to_distance = lambda *args, **kwargs: None
    sys.modules["revise.backend.ops.distance"] = distance

    matplotlib = types.ModuleType("matplotlib")
    matplotlib.__path__ = []
    pyplot = types.ModuleType("matplotlib.pyplot")
    pyplot.close = lambda *args, **kwargs: None
    pyplot.savefig = lambda *args, **kwargs: None
    sys.modules["matplotlib"] = matplotlib
    sys.modules["matplotlib.pyplot"] = pyplot

    analysis = types.ModuleType("revise.analysis")
    analysis.__path__ = []
    metrics = types.ModuleType("revise.analysis.metrics")
    metrics.compute_clustering_metrics = lambda *args, **kwargs: None
    sys.modules["revise.analysis"] = analysis
    sys.modules["revise.analysis.metrics"] = metrics


@pytest.fixture
def runner_modules():
    snapshot = _snapshot_modules(ISOLATED_MODULE_NAMES)
    for module_name in ISOLATED_MODULE_NAMES:
        sys.modules.pop(module_name, None)
    _install_import_stubs()
    try:
        yield SimpleNamespace(
            application=importlib.import_module(RUNNER_MODULE_NAMES[0]),
            benchmark=importlib.import_module(RUNNER_MODULE_NAMES[1]),
            sr_application=importlib.import_module(RUNNER_MODULE_NAMES[2]),
            sr_benchmark=importlib.import_module(RUNNER_MODULE_NAMES[3]),
            sp_application=importlib.import_module(RUNNER_MODULE_NAMES[4]),
        )
    finally:
        _restore_modules(snapshot)


def _adata(gene: str) -> AnnData:
    return AnnData(
        X=np.ones((1, 1), dtype=np.float64),
        obs=pd.DataFrame(index=["obs-1"]),
        var=pd.DataFrame(index=[gene]),
    )


def _check_message(exc_info, *fragments: str) -> None:
    message = str(exc_info.value)
    missing = [fragment for fragment in fragments if fragment not in message]
    if missing:
        raise AssertionError(f"missing message fragments {missing}: {message}")


def test_application_overlap_validation_is_contextual_value_error(runner_modules):
    target = SimpleNamespace(
        st_adata=_adata("st-gene"),
        sc_ref_adata=_adata("sc-gene"),
        config=SimpleNamespace(
            st_file_path="/data/sample_st.h5ad",
            sc_ref_file_path="/data/reference_sc.h5ad",
        ),
    )

    with pytest.raises(ValueError) as exc_info:
        runner_modules.application.ApplicationSVC._adata_validate(target)

    _check_message(
        exc_info,
        "st_file_path=/data/sample_st.h5ad",
        "sc_ref_file_path=/data/reference_sc.h5ad",
        "field=var_names_overlap",
        "expected=>=1",
        "actual=0",
    )


def test_benchmark_requires_contextual_ground_truth_value_error(runner_modules):
    target = SimpleNamespace(
        st_adata=_adata("gene"),
        sc_ref_adata=_adata("gene"),
        real_st_adata=None,
        config=SimpleNamespace(gt_svc_file_path="/data/ground_truth.h5ad"),
    )

    with pytest.raises(ValueError) as exc_info:
        runner_modules.benchmark.BenchmarkSVC._adata_validate(target)

    _check_message(
        exc_info,
        "gt_svc_file_path=/data/ground_truth.h5ad",
        "field=real_st_adata",
        "expected=loaded",
        "actual=None",
    )


def test_benchmark_overlap_validation_is_contextual_value_error(runner_modules):
    target = SimpleNamespace(
        st_adata=_adata("st-gene"),
        sc_ref_adata=_adata("sc-gene"),
        real_st_adata=_adata("ground-truth-gene"),
        config=SimpleNamespace(
            st_file_path="/data/benchmark_st.h5ad",
            sc_ref_file_path="/data/benchmark_sc.h5ad",
        ),
    )

    with pytest.raises(ValueError) as exc_info:
        runner_modules.benchmark.BenchmarkSVC._adata_validate(target)

    _check_message(
        exc_info,
        "st_file_path=/data/benchmark_st.h5ad",
        "sc_ref_file_path=/data/benchmark_sc.h5ad",
        "field=var_names_overlap",
        "expected=>=1",
        "actual=0",
    )


@pytest.mark.parametrize(
    ("module_attr", "st_path"),
    [
        ("sr_application", "/data/application_st.h5ad"),
        ("sr_benchmark", "/data/benchmark_st.h5ad"),
    ],
)
def test_sr_mapping_validation_is_contextual_key_error(
    runner_modules, module_attr, st_path
):
    target = SimpleNamespace(
        st_adata=_adata("gene"),
        config=SimpleNamespace(st_file_path=st_path),
    )

    with pytest.raises(KeyError) as exc_info:
        getattr(runner_modules, module_attr).ScSVCSr._adata_validate_dec(target)

    _check_message(
        exc_info,
        f"st_file_path={st_path}",
        "field=uns['all_cells_in_spot']",
        "expected=present",
        "actual=missing",
    )


@pytest.mark.parametrize(
    ("current", "original", "actual"),
    [
        (
            np.array([[1.0, 9.0], [3.0, 4.0]]),
            np.array([[1.0, 2.0], [3.0, 4.0]]),
            "actual=mismatches=1",
        ),
        (
            sparse.csr_matrix([[1.0, 9.0], [3.0, 4.0]]),
            sparse.csr_matrix([[1.0, 2.0], [3.0, 4.0]]),
            "actual=mismatches=1",
        ),
        (
            np.ones((1, 2), dtype=np.float64),
            np.ones((2, 1), dtype=np.float64),
            "actual=shape=(1, 2), original_shape=(2, 1)",
        ),
        (
            np.array([[1.0, 2.0]], dtype=np.float64),
            sparse.csr_matrix([[1.0, 2.0]]),
            "actual=current_representation=dense, original_representation=sparse",
        ),
    ],
    ids=("dense", "sparse", "shape", "mixed-representation"),
)
def test_sp_expression_invariant_raises_contextual_runtime_error(
    runner_modules, current, original, actual
):
    with pytest.raises(RuntimeError) as exc_info:
        runner_modules.sp_application._validate_expression_unchanged(
            current,
            original,
            cell_type="T-cell",
        )

    _check_message(
        exc_info,
        "cell_type=T-cell",
        "field=X",
        "expected=exactly unchanged",
        actual,
    )


@pytest.mark.parametrize(
    ("current", "original"),
    [
        (
            np.array([[1.0, 2.0], [3.0, 4.0]]),
            np.array([[1.0, 2.0], [3.0, 4.0]]),
        ),
        (
            sparse.csr_matrix([[1.0, 2.0], [3.0, 4.0]]),
            sparse.csr_matrix([[1.0, 2.0], [3.0, 4.0]]),
        ),
    ],
    ids=("dense", "sparse"),
)
def test_sp_expression_invariant_accepts_exactly_equal_representations(
    runner_modules, current, original
):
    runner_modules.sp_application._validate_expression_unchanged(
        current,
        original,
        cell_type="T-cell",
    )


def test_sp_local_refinement_calls_expression_invariant(runner_modules):
    source = textwrap.dedent(
        inspect.getsource(runner_modules.sp_application.SpSVC.local_refinement)
    )
    tree = ast.parse(source)
    calls = [
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    ]
    if "_validate_expression_unchanged" not in calls:
        raise AssertionError(
            "SpSVC.local_refinement must call _validate_expression_unchanged"
        )


def test_runner_import_stubs_do_not_leak_to_other_tests():
    consumers = (
        "revise.backend.adapters",
        "revise.backend.ops.meta",
        "revise.backend.ops.topology",
        "revise.backend.ops.shaver",
        *RUNNER_MODULE_NAMES,
    )
    leaked = [
        module_name
        for module_name in consumers
        if getattr(
            getattr(sys.modules.get(module_name), "sc", None),
            "_revise_test_stub",
            False,
        )
    ]
    if leaked:
        raise AssertionError(f"scanpy import stubs leaked through modules: {leaked}")


def test_runtime_validations_survive_python_optimized_mode(tmp_path):
    repo_root = Path(__file__).parents[2]
    env = os.environ.copy()
    env["MPLCONFIGDIR"] = str(tmp_path / "matplotlib")
    env["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] = "1"
    env["PYTHONPATH"] = os.pathsep.join(
        [str(repo_root), env.get("PYTHONPATH", "")]
    ).rstrip(os.pathsep)

    result = subprocess.run(
        [
            sys.executable,
            "-O",
            "-m",
            "pytest",
            "-p",
            "no:capture",
            "-q",
            "tests/recon/test_runtime_validation_errors.py",
            "-k",
            "not runtime_validations_survive_python_optimized_mode",
        ],
        cwd=repo_root,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    if result.returncode != 0:
        raise AssertionError(result.stdout + result.stderr)
