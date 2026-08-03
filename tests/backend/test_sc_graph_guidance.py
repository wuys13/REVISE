from __future__ import annotations

import importlib
import inspect
import logging
import sys
import types
from contextlib import contextmanager
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData


_MISSING = object()


def _snapshot_modules(names):
    modules = {name: sys.modules.get(name, _MISSING) for name in names}
    parent_attributes = {}
    for name in names:
        parent_name, separator, attribute = name.rpartition(".")
        if separator and parent_name in sys.modules:
            parent_attributes[(parent_name, attribute)] = getattr(
                sys.modules[parent_name], attribute, _MISSING
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
def _isolated_runner_module():
    names = (
        "scanpy",
        "revise.backend.runners.sc_svc_application",
        "revise.backend.runners.application_svc",
        "revise.backend.kernels",
        "revise.analysis",
        "revise.analysis.bio",
    )
    snapshot = _snapshot_modules(names)
    _remove_modules(names)

    scanpy = types.ModuleType("scanpy")
    scanpy.pp = SimpleNamespace()
    scanpy.pl = SimpleNamespace()
    sys.modules["scanpy"] = scanpy

    runners_base = types.ModuleType("revise.backend.runners.application_svc")
    runners_base.ApplicationSVC = object
    sys.modules[runners_base.__name__] = runners_base

    kernels = types.ModuleType("revise.backend.kernels")
    kernels.GraphClusterKernel = object
    kernels.LocalAnchoringKernel = object
    sys.modules[kernels.__name__] = kernels

    analysis = types.ModuleType("revise.analysis")
    analysis.__path__ = []
    sys.modules[analysis.__name__] = analysis
    bio = types.ModuleType("revise.analysis.bio")
    bio.get_degs = lambda *_args, **_kwargs: None
    bio.conclusions_write = lambda *_args, **_kwargs: None
    bio.plot_volcano = lambda *_args, **_kwargs: None
    sys.modules[bio.__name__] = bio

    try:
        yield importlib.import_module("revise.backend.runners.sc_svc_application")
    finally:
        _remove_modules(names)
        _restore_modules(snapshot)


def _spatial(level1_posterior):
    names = ["sp-a1", "sp-a2", "sp-b1", "sp-b2"]
    adata = AnnData(
        X=np.ones((4, 2), dtype=np.float64),
        obs=pd.DataFrame(
            {
                "Level1": ["A", "A", "B", "B"],
                "Confidence": np.max(level1_posterior, axis=1),
            },
            index=names,
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    adata.obsm["Level1"] = pd.DataFrame(
        level1_posterior,
        index=names,
        columns=["A", "B"],
    )
    return adata


def _reference():
    return AnnData(
        X=np.ones((4, 2), dtype=np.float64),
        obs=pd.DataFrame(
            {
                "Level1": ["A", "A", "B", "B"],
                "Level2": ["a1", "a2", "b1", "b2"],
            },
            index=["ref-a1", "ref-a2", "ref-b1", "ref-b2"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )


def _run_a_cohort(module, level1_posterior):
    runner = module.ScSVC.__new__(module.ScSVC)
    runner.st_adata = _spatial(level1_posterior)
    runner.sc_ref_adata = _reference()
    runner.config = SimpleNamespace(
        cell_type_col="Level1",
        confidence_col="Confidence",
        unknown_key="Unknown",
    )
    runner.logger = logging.getLogger("test-sc-svc-store-only")
    runner.cluster_col = "SVC_cluster"
    captured = {"local": [], "graph": []}

    class LocalProducer:
        def run(self, target, reference, **kwargs):
            captured["local"].append(
                {
                    "target": tuple(target.obs_names),
                    "reference": tuple(reference.obs_names),
                    "column": kwargs["cell_type_col"],
                }
            )
            result = target.copy()
            if kwargs["cell_type_col"] == "Level2":
                result.obsm["Level2"] = pd.DataFrame(
                    [[0.8, 0.2], [0.1, 0.9]],
                    index=result.obs_names,
                    columns=["a1", "a2"],
                )
                result.obs["Level2"] = ["a1", "a2"]
            return result

    runner.local_annotate_method = LocalProducer()

    def graph_run(adata, resolutions, evaluation_label):
        captured["graph"].append(
            {
                "obs_names": tuple(adata.obs_names),
                "level2": tuple(adata.obs["Level2"]),
                "level2_q": adata.obsm["Level2"].copy(),
                "resolutions": tuple(resolutions),
                "evaluation_label": evaluation_label,
            }
        )
        result = adata.copy()
        result.obs["leiden_0.5"] = pd.Categorical(["0", "1"])
        return result, pd.DataFrame(), 0.5

    runner.graph_cluster = SimpleNamespace(run=graph_run)
    spatial_result, expression_result = runner.local_refinement(
        "A", "Level2", [0.5], select_res=0.5
    )
    return spatial_result, expression_result, captured


def test_graph_cluster_public_api_is_minimal():
    names = (
        "scanpy",
        "squidpy",
        "revise.backend.kernels.graph_cluster",
    )
    snapshot = _snapshot_modules(names)
    _remove_modules(names)
    scanpy = types.ModuleType("scanpy")
    scanpy.pp = SimpleNamespace()
    scanpy.tl = SimpleNamespace()
    scanpy.pl = SimpleNamespace()
    squidpy = types.ModuleType("squidpy")
    squidpy.gr = SimpleNamespace()
    sys.modules["scanpy"] = scanpy
    sys.modules["squidpy"] = squidpy
    try:
        graph_module = importlib.import_module(
            "revise.backend.kernels.graph_cluster"
        )
        parameters = inspect.signature(
            graph_module.GraphClusterKernel.run
        ).parameters
    finally:
        _remove_modules(names)
        _restore_modules(snapshot)

    assert tuple(parameters) == ("self", "adata", "resolution", "label")


def test_same_argmax_different_ga_q_preserves_cohort_and_level2_consumption():
    first_q = np.array(
        [[0.9, 0.1], [0.6, 0.4], [0.2, 0.8], [0.1, 0.9]]
    )
    second_q = np.array(
        [[0.51, 0.49], [0.99, 0.01], [0.49, 0.51], [0.3, 0.7]]
    )

    with _isolated_runner_module() as module:
        first_spatial, _first_expr, first = _run_a_cohort(module, first_q)
        second_spatial, _second_expr, second = _run_a_cohort(module, second_q)

    assert first["local"] == second["local"] == [
        {
            "target": ("sp-a1", "sp-a2"),
            "reference": ("ref-a1", "ref-a2"),
            "column": "Level2",
        },
        {
            "target": ("ref-a1", "ref-a2"),
            "reference": ("sp-a1", "sp-a2"),
            "column": "SVC_cluster",
        },
    ]
    assert first["graph"][0]["obs_names"] == second["graph"][0]["obs_names"]
    assert first["graph"][0]["level2"] == second["graph"][0]["level2"]
    pd.testing.assert_frame_equal(
        first["graph"][0]["level2_q"], second["graph"][0]["level2_q"]
    )
    assert first_spatial.obs["SVC_cluster"].tolist() == ["0", "1"]
    assert second_spatial.obs["SVC_cluster"].tolist() == ["0", "1"]


def test_sc_svc_retains_ga_q_and_graph_cluster_consumes_level2_results():
    level1_q = np.array(
        [[0.9, 0.1], [0.6, 0.4], [0.2, 0.8], [0.1, 0.9]]
    )

    with _isolated_runner_module() as module:
        spatial_result, _expression_result, captured = _run_a_cohort(
            module, level1_q
        )

    expected = pd.DataFrame(
        level1_q[:2],
        index=["sp-a1", "sp-a2"],
        columns=["A", "B"],
    )
    pd.testing.assert_frame_equal(spatial_result.obsm["Level1"], expected)
    np.testing.assert_allclose(expected.max(axis=1), [0.9, 0.6])
    assert set(captured["graph"][0]) == {
        "obs_names",
        "level2",
        "level2_q",
        "resolutions",
        "evaluation_label",
    }


def test_missing_reference_and_singleton_fail_before_local_anchoring():
    with _isolated_runner_module() as module:
        runner = module.ScSVC.__new__(module.ScSVC)
        runner.config = SimpleNamespace(cell_type_col="Level1")
        runner.logger = logging.getLogger("test-sc-svc-invalid-cohort")
        runner.cluster_col = "SVC_cluster"
        runner.local_annotate_method = SimpleNamespace(
            run=lambda *_args, **_kwargs: pytest.fail(
                "invalid cohorts must bypass local anchoring"
            )
        )

        runner.st_adata = _spatial(
            np.array([[0.9, 0.1], [0.6, 0.4], [0.2, 0.8], [0.1, 0.9]])
        )
        runner.sc_ref_adata = _reference()[
            _reference().obs["Level1"] == "B"
        ].copy()
        with pytest.raises(ValueError, match="no reference cells"):
            runner.local_refinement("A", "Level2", [0.5])

        runner.st_adata = runner.st_adata[["sp-a1"]].copy()
        runner.sc_ref_adata = _reference()
        with pytest.raises(ValueError, match="at least two spatial cells"):
            runner.local_refinement("A", "Level2", [0.5])
