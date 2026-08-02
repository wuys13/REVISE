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

from revise.backend.ops.assignment import (
    AssignmentState,
    AssignmentStateError,
    one_hot_assignment,
)
from revise.backend.ops.assignment_guidance import (
    AssignmentGuidanceCollector,
    FallbackReason,
    NotApplicableReason,
)


_MISSING = object()
_TEST_ISOLATION_PREFIXES = (
    "scanpy",
    "squidpy",
    "revise.backend.adapters",
    "revise.backend.kernels.graph_cluster",
    "revise.backend.runners.application_svc",
    "revise.backend.runners.sc_svc_application",
    "revise.analysis.bio",
)
_GRAPH_MODULES = (
    "scanpy",
    "squidpy",
    "revise.backend.kernels.graph_cluster",
)


def _module_names_for_prefixes(prefixes):
    return tuple(
        name
        for name in sys.modules
        if any(
            name == prefix or name.startswith(f"{prefix}.")
            for prefix in prefixes
        )
    )


@pytest.fixture(autouse=True)
def _restore_sc_graph_test_modules():
    names = _module_names_for_prefixes(_TEST_ISOLATION_PREFIXES)
    snapshot = _snapshot_modules(names)
    yield
    current = _module_names_for_prefixes(_TEST_ISOLATION_PREFIXES)
    _remove_modules(current)
    _restore_modules(snapshot)


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


@pytest.fixture
def graph_module():
    snapshot = _snapshot_modules(_GRAPH_MODULES)
    _remove_modules(_GRAPH_MODULES)
    scanpy = types.ModuleType("scanpy")
    scanpy.pp = SimpleNamespace()
    scanpy.tl = SimpleNamespace()
    scanpy.pl = SimpleNamespace()
    squidpy = types.ModuleType("squidpy")
    squidpy.gr = SimpleNamespace()
    sys.modules["scanpy"] = scanpy
    sys.modules["squidpy"] = squidpy
    try:
        yield importlib.import_module("revise.backend.kernels.graph_cluster")
    finally:
        _remove_modules(_GRAPH_MODULES)
        _restore_modules(snapshot)


def _config(
    collector,
    *,
    guidance="prefer",
    compatibility_mode="cost",
    strength=2.0,
):
    return SimpleNamespace(
        rec_random_state=0,
        rec_graph_alpha=0.0,
        rec_graph_method="pca",
        rec_ot_method="pot",
        plot_flag=False,
        assignment_guidance_policy=guidance,
        posterior_conditioning_enabled=guidance != "off",
        posterior_conditioning_mode=compatibility_mode,
        posterior_conditioning_strict=guidance == "require",
        posterior_conditioning_beta=1.0,
        posterior_conditioning_min_affinity=0.1,
        posterior_conditioning_cost_strength=strength,
        assignment_guidance_callback=collector.callback,
        assignment_guidance_route="real2real:cellular",
    )


def _adata():
    names = ["cell-1", "cell-2", "cell-3", "cell-4"]
    adata = AnnData(
        X=np.ones((4, 3), dtype=np.float64),
        obs=pd.DataFrame(index=names),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
    )
    adata.obsm["Level1"] = pd.DataFrame(
        [
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ],
        index=names,
        columns=["broad-a", "broad-b"],
    )
    return adata


def _soft_level2(adata):
    return AssignmentState(
        values=np.array(
            [
                [0.9, 0.1],
                [0.8, 0.2],
                [0.1, 0.9],
                [0.2, 0.8],
            ]
        ),
        observation_labels=adata.obs_names,
        category_labels=("sub-a", "sub-b"),
        source="local_anchoring:obsm[Level2]",
        level="Level2",
        value_semantics="soft",
        lineage=[
            {
                "operation": "local_anchoring",
                "container": "obsm",
                "key": "Level2",
            }
        ],
    )


def _patch_graph_runtime(
    module,
    monkeypatch,
    *,
    fail_base=False,
    fail_guided=False,
    guided_labels_by_resolution=None,
):
    base_graph = sparse.csr_matrix(
        [
            [0.0, 1.0, 0.5, 0.0],
            [1.0, 0.0, 0.0, 0.5],
            [0.5, 0.0, 0.0, 1.0],
            [0.0, 0.5, 1.0, 0.0],
        ]
    )
    spatial_graph = sparse.eye(4, format="csr")
    captured = {"leiden": []}

    def highly_variable_genes(adata, **_kwargs):
        adata.var["highly_variable"] = True

    def neighbors(adata, **_kwargs):
        adata.obsp["connectivities"] = base_graph.copy()

    def spatial_neighbors(adata):
        adata.obsp["spatial_connectivities"] = spatial_graph.copy()

    def leiden(adata, *, adjacency, resolution, key_added, **_kwargs):
        graph = adjacency.copy()
        is_guided = not np.allclose(graph.toarray(), base_graph.toarray())
        captured["leiden"].append(
            {
                "resolution": resolution,
                "key_added": key_added,
                "guided": is_guided,
                "adjacency": graph,
            }
        )
        if fail_base and not is_guided:
            raise RuntimeError("base resolution selection boom")
        if fail_guided and is_guided:
            raise RuntimeError("guided clustering boom")
        guided_labels = (guided_labels_by_resolution or {}).get(
            float(resolution)
        )
        if is_guided and guided_labels is not None:
            adata.obs[key_added] = pd.Categorical(guided_labels)
        elif float(resolution) == 1.0:
            adata.obs[key_added] = pd.Categorical(["0", "0", "1", "1"])
        else:
            adata.obs[key_added] = pd.Categorical(["0", "1", "0", "1"])

    module.sc.pp = SimpleNamespace(
        filter_genes=lambda *_args, **_kwargs: None,
        normalize_total=lambda *_args, **_kwargs: None,
        log1p=lambda *_args, **_kwargs: None,
        highly_variable_genes=highly_variable_genes,
        pca=lambda *_args, **_kwargs: None,
        neighbors=neighbors,
    )
    module.sc.tl = SimpleNamespace(leiden=leiden)
    module.sc.pl = SimpleNamespace(scatter=lambda *_args, **_kwargs: None)
    module.sq.gr = SimpleNamespace(spatial_neighbors=spatial_neighbors)

    from revise.backend.ops import coefficients

    monkeypatch.setattr(
        coefficients,
        "get_weighted_align_score",
        lambda _adata, *, res, label: 0.9 if float(res) == 1.0 else 0.2,
    )
    return base_graph, captured


def _run_kernel(
    module,
    monkeypatch,
    *,
    guidance="prefer",
    state=True,
    fail_guided=False,
):
    collector = AssignmentGuidanceCollector()
    adata = _adata()
    base_graph, captured = _patch_graph_runtime(
        module,
        monkeypatch,
        fail_guided=fail_guided,
    )
    guidance_state = _soft_level2(adata) if state is True else state
    kernel = module.GraphClusterKernel(
        _config(collector, guidance=guidance),
        logging.getLogger("test-standard-sc-graph-guidance"),
    )
    result = kernel.run(
        adata,
        resolution=[0.5, 1.0],
        label="Level2",
        guidance_state=guidance_state,
        problem_key="standard-sc:A",
    )
    return result, collector, captured, base_graph


@pytest.mark.parametrize(
    ("canonical", "legacy", "expected"),
    [("require", "off", "require"), ("off", "require", "off")],
)
def test_graph_route_prefers_canonical_guidance_over_conflicting_legacy_fields(
    graph_module,
    canonical,
    legacy,
    expected,
):
    collector = AssignmentGuidanceCollector()
    config = _config(collector, guidance=legacy)
    config.assignment_guidance_policy = canonical
    kernel = graph_module.GraphClusterKernel(
        config,
        logging.getLogger("test-canonical-graph-guidance"),
    )

    assert kernel._guidance_mode() == expected


def test_default_off_does_not_read_or_validate_level2_state(
    graph_module,
    monkeypatch,
):
    monkeypatch.setattr(
        graph_module,
        "resolve_assignment_guidance",
        lambda *_args, **_kwargs: pytest.fail(
            "off guidance must not read Assignment State"
        ),
        raising=False,
    )

    result, collector, captured, base_graph = _run_kernel(
        graph_module,
        monkeypatch,
        guidance="off",
        state=AssignmentState(
            values=np.full((4, 2), np.nan),
            observation_labels=("bad",) * 4,
            category_labels=("bad", "bad"),
            source="invalid",
            level="Level2",
            value_semantics="soft",
            lineage=[],
        ),
    )

    _adata_out, _metrics, best_res = result
    assert best_res == 1.0
    assert [call["guided"] for call in captured["leiden"]] == [False, False]
    for call in captured["leiden"]:
        np.testing.assert_allclose(call["adjacency"].toarray(), base_graph.toarray())
    [event] = collector.events
    assert event["outcome"] == "off"
    assert event["availability"] == "not_checked"
    assert event["left_assignment"] is None


def test_prefer_uses_explicit_level2_not_conflicting_global_level1(
    graph_module,
    monkeypatch,
):
    observed = {}
    real_compatibility = graph_module.assignment_compatibility

    def capture(left, right, **kwargs):
        observed["left"] = left
        observed["right"] = right
        observed["support"] = kwargs["support"]
        return real_compatibility(left, right, **kwargs)

    monkeypatch.setattr(graph_module, "assignment_compatibility", capture)

    (_adata_out, metrics, best_res), collector, captured, _base = _run_kernel(
        graph_module,
        monkeypatch,
    )

    assert best_res == 1.0
    assert metrics.set_index("resolution").loc[1.0, "align_score"] == 0.9
    assert [call["guided"] for call in captured["leiden"]] == [
        False,
        False,
        True,
    ]
    assert [call["resolution"] for call in captured["leiden"]] == [0.5, 1.0, 1.0]
    assert observed["left"].level == "Level2"
    assert observed["left"].category_labels == ("sub-a", "sub-b")
    assert observed["right"].source == "local_anchoring:obsm[Level2]"
    rows, columns = observed["support"]
    assert rows.ndim == columns.ndim == 1
    assert len(rows) < 4 * 4
    [event] = collector.events
    assert event["operator"] == "graph_edge"
    assert event["outcome"] == "applied"
    assert event["left_assignment"]["level"] == "Level2"
    assert event["right_assignment"]["level"] == "Level2"


def test_resolution_search_is_unguided_then_fixed_best_resolution_is_guided(
    graph_module,
    monkeypatch,
):
    (_adata_out, _metrics, best_res), _collector, captured, _base = _run_kernel(
        graph_module,
        monkeypatch,
    )

    assert best_res == 1.0
    assert [
        (call["resolution"], call["guided"])
        for call in captured["leiden"]
    ] == [(0.5, False), (1.0, False), (1.0, True)]


@pytest.mark.parametrize("semantics", ["soft", "one_hot"])
def test_soft_and_one_hot_level2_use_the_same_graph_edge_interface(
    graph_module,
    monkeypatch,
    semantics,
):
    adata = _adata()
    state = _soft_level2(adata)
    if semantics == "one_hot":
        state = one_hot_assignment(
            ("sub-a", "sub-a", "sub-b", "sub-b"),
            observation_labels=adata.obs_names,
            category_labels=("sub-a", "sub-b"),
            source="local_anchoring:argmax[Level2]",
            level="Level2",
            lineage=state.lineage,
        )

    (_out, _metrics, _best), collector, captured, _base = _run_kernel(
        graph_module,
        monkeypatch,
        state=state,
    )

    assert [call["guided"] for call in captured["leiden"]] == [
        False,
        False,
        True,
    ]
    [event] = collector.events
    assert event["operator"] == "graph_edge"
    assert event["outcome"] == "applied"
    assert event["left_assignment"]["value_semantics"] == semantics


def test_missing_level2_falls_back_under_prefer(
    graph_module,
    monkeypatch,
):
    (_out, _metrics, _best), collector, captured, _base = _run_kernel(
        graph_module,
        monkeypatch,
        state=None,
    )

    assert [call["guided"] for call in captured["leiden"]] == [False, False]
    [event] = collector.events
    assert event["outcome"] == "fallback"
    assert event["availability"] == "unavailable"
    assert event["reason"] == FallbackReason.ASSIGNMENT_MISSING.value


def test_missing_level2_fails_require_after_unguided_resolution_selection(
    graph_module,
    monkeypatch,
):
    collector = AssignmentGuidanceCollector()
    adata = _adata()
    _base, captured = _patch_graph_runtime(graph_module, monkeypatch)
    kernel = graph_module.GraphClusterKernel(
        _config(collector, guidance="require"),
        logging.getLogger("test-standard-sc-require"),
    )

    with pytest.raises(AssignmentStateError, match="assignment_state_unavailable"):
        kernel.run(
            adata,
            resolution=[0.5, 1.0],
            label="Level2",
            guidance_state=None,
            problem_key="standard-sc:A",
        )

    assert [call["guided"] for call in captured["leiden"]] == [False, False]
    [event] = collector.events
    assert event["outcome"] == "failed"
    assert event["attempted"] is False
    assert event["reason"] is None


@pytest.mark.parametrize(
    ("guidance", "expected_outcome", "raises"),
    [
        ("prefer", "fallback", False),
        ("require", "failed", True),
    ],
)
def test_invalid_level2_obeys_prefer_require_policy(
    graph_module,
    monkeypatch,
    guidance,
    expected_outcome,
    raises,
):
    collector = AssignmentGuidanceCollector()
    adata = _adata()
    _base, captured = _patch_graph_runtime(graph_module, monkeypatch)
    invalid = _soft_level2(adata)
    invalid.values[0, 0] = np.nan
    kernel = graph_module.GraphClusterKernel(
        _config(collector, guidance=guidance),
        logging.getLogger("test-standard-sc-invalid-level2"),
    )

    if raises:
        with pytest.raises(AssignmentStateError, match="values_nan"):
            kernel.run(
                adata,
                resolution=[0.5, 1.0],
                label="Level2",
                guidance_state=invalid,
                problem_key="standard-sc:A",
            )
        assert [call["guided"] for call in captured["leiden"]] == [False, False]
    else:
        kernel.run(
            adata,
            resolution=[0.5, 1.0],
            label="Level2",
            guidance_state=invalid,
            problem_key="standard-sc:A",
        )
        assert [call["guided"] for call in captured["leiden"]] == [False, False]
    [event] = collector.events
    assert event["outcome"] == expected_outcome
    assert event["reason"] == (
        None
        if raises
        else FallbackReason.ASSIGNMENT_INVALID.value
    )


def test_base_resolution_failure_leaves_guidance_provenance_unstarted(
    graph_module,
    monkeypatch,
):
    collector = AssignmentGuidanceCollector()
    adata = _adata()
    _base, captured = _patch_graph_runtime(
        graph_module,
        monkeypatch,
        fail_base=True,
    )
    kernel = graph_module.GraphClusterKernel(
        _config(collector),
        logging.getLogger("test-standard-sc-base-selection-failure"),
    )

    with pytest.raises(RuntimeError, match="base resolution selection boom"):
        kernel.run(
            adata,
            resolution=[0.5, 1.0],
            label="Level2",
            guidance_state=_soft_level2(adata),
            problem_key="standard-sc:A",
        )

    assert [call["guided"] for call in captured["leiden"]] == [False]
    assert collector.manifest()["events"] == []


def test_guided_clustering_failure_is_terminal_and_rethrown(
    graph_module,
    monkeypatch,
):
    collector = AssignmentGuidanceCollector()
    adata = _adata()
    _base, captured = _patch_graph_runtime(
        graph_module,
        monkeypatch,
        fail_guided=True,
    )
    kernel = graph_module.GraphClusterKernel(
        _config(collector),
        logging.getLogger("test-standard-sc-guided-failure"),
    )

    with pytest.raises(RuntimeError, match="guided clustering boom"):
        kernel.run(
            adata,
            resolution=[0.5, 1.0],
            label="Level2",
            guidance_state=_soft_level2(adata),
            problem_key="standard-sc:A",
        )

    assert [call["guided"] for call in captured["leiden"]] == [
        False,
        False,
        True,
    ]
    [event] = collector.events
    assert event["attempted"] is True
    assert event["outcome"] == "failed"
    assert event["reason"] is None


def test_graph_compatibility_failure_is_terminal_and_rethrown(
    graph_module,
    monkeypatch,
):
    collector = AssignmentGuidanceCollector()
    adata = _adata()
    _base, _captured = _patch_graph_runtime(graph_module, monkeypatch)
    monkeypatch.setattr(
        graph_module,
        "assignment_compatibility",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            RuntimeError("compatibility boom")
        ),
    )
    kernel = graph_module.GraphClusterKernel(
        _config(collector),
        logging.getLogger("test-standard-sc-compatibility-failure"),
    )

    with pytest.raises(RuntimeError, match="compatibility boom"):
        kernel.run(
            adata,
            resolution=[0.5, 1.0],
            label="Level2",
            guidance_state=_soft_level2(adata),
            problem_key="standard-sc:A",
        )

    [event] = collector.events
    assert event["attempted"] is True
    assert event["outcome"] == "failed"
    assert event["reason"] is None


def test_guided_clustering_interrupt_is_terminal_and_rethrown(
    graph_module,
    monkeypatch,
):
    collector = AssignmentGuidanceCollector()
    adata = _adata()
    _base, _captured = _patch_graph_runtime(graph_module, monkeypatch)
    original_leiden = graph_module.sc.tl.leiden

    def interrupt_guided(adata, *, adjacency, **kwargs):
        if not np.allclose(
            adjacency.toarray(),
            _base.toarray(),
        ):
            raise KeyboardInterrupt("guided clustering interrupted")
        return original_leiden(adata, adjacency=adjacency, **kwargs)

    graph_module.sc.tl.leiden = interrupt_guided
    kernel = graph_module.GraphClusterKernel(
        _config(collector),
        logging.getLogger("test-standard-sc-guided-interrupt"),
    )

    with pytest.raises(KeyboardInterrupt, match="guided clustering interrupted"):
        kernel.run(
            adata,
            resolution=[0.5, 1.0],
            label="Level2",
            guidance_state=_soft_level2(adata),
            problem_key="standard-sc:A",
        )

    [event] = collector.events
    assert event["attempted"] is True
    assert event["outcome"] == "interrupted"
    assert event["reason"] is None


def test_singleton_records_not_applicable_without_assignment_fallback(
    graph_module,
):
    collector = AssignmentGuidanceCollector()
    adata = AnnData(
        X=np.ones((1, 1)),
        obs=pd.DataFrame(index=["cell-1"]),
        var=pd.DataFrame(index=["g1"]),
    )
    kernel = graph_module.GraphClusterKernel(
        _config(collector),
        logging.getLogger("test-standard-sc-singleton"),
    )

    with pytest.raises(ValueError, match="at least two spatial cells"):
        kernel.run(
            adata,
            resolution=[0.5],
            label="Level2",
            guidance_state=None,
            problem_key="standard-sc:A",
        )

    [event] = collector.events
    assert event["applicability"] == "not_applicable"
    assert event["outcome"] == "not_applicable"
    assert event["reason"] == NotApplicableReason.INSUFFICIENT_UNITS.value
    assert event["reason_details"] == {
        "unit": "spatial_cell",
        "observed": 1,
        "required": 2,
    }


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

    base = types.ModuleType("revise.backend.runners.application_svc")

    class ApplicationSVC:
        pass

    base.ApplicationSVC = ApplicationSVC
    sys.modules[base.__name__] = base

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


def test_runner_routes_level1_cohort_but_passes_first_level2_state_explicitly():
    with _isolated_runner_module() as module:
        collector = AssignmentGuidanceCollector()
        config = _config(collector, guidance="off")
        config.assignment_guidance_policy = "require"
        config.cell_type_col = "Level1"
        config.confidence_col = "Confidence"
        config.unknown_key = "Unknown"
        spatial = AnnData(
            X=np.ones((4, 2)),
            obs=pd.DataFrame(
                {"Level1": ["A", "A", "B", "B"]},
                index=["sp-a1", "sp-a2", "sp-b1", "sp-b2"],
            ),
            var=pd.DataFrame(index=["g1", "g2"]),
        )
        spatial.obsm["Level1"] = pd.DataFrame(
            [[0.0, 1.0], [0.0, 1.0], [1.0, 0.0], [1.0, 0.0]],
            index=spatial.obs_names,
            columns=["A", "B"],
        )
        reference = AnnData(
            X=np.ones((4, 2)),
            obs=pd.DataFrame(
                {
                    "Level1": ["A", "A", "B", "B"],
                    "Level2": ["a1", "a2", "b1", "b2"],
                },
                index=["ref-a1", "ref-a2", "ref-b1", "ref-b2"],
            ),
            var=pd.DataFrame(index=["g1", "g2"]),
        )
        runner = module.ScSVC.__new__(module.ScSVC)
        runner.st_adata = spatial
        runner.sc_ref_adata = reference
        runner.config = config
        runner.logger = logging.getLogger("test-standard-sc-runner")
        runner.cluster_col = "SVC_cluster"
        local_calls = []

        def local_run(target, ref, **kwargs):
            local_calls.append(
                {
                    "target": tuple(target.obs_names),
                    "reference": tuple(ref.obs_names),
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

        def assignment_state(annotated, key):
            posterior = annotated.obsm[key]
            return AssignmentState(
                values=posterior.to_numpy(),
                observation_labels=posterior.index,
                category_labels=posterior.columns,
                source=f"local_anchoring:obsm[{key}]",
                level=key,
                value_semantics="soft",
                lineage=[
                    {
                        "operation": "local_anchoring",
                        "container": "obsm",
                        "key": key,
                    }
                ],
            )

        runner.local_annotate_method = SimpleNamespace(
            run=local_run,
            assignment_state=assignment_state,
        )
        captured = {}

        def graph_run(
            adata,
            resolutions,
            evaluation_label,
            *,
            guidance_state,
            problem_key,
            selected_resolution,
        ):
            captured.update(
                adata=adata,
                resolutions=resolutions,
                evaluation_label=evaluation_label,
                guidance_state=guidance_state,
                problem_key=problem_key,
                selected_resolution=selected_resolution,
            )
            result = adata.copy()
            result.obs["leiden_0.5"] = pd.Categorical(["0", "1"])
            metrics = pd.DataFrame(
                {"resolution": [0.5], "cluster_num": [2], "align_score": [1.0]}
            )
            return result, metrics, 0.5

        runner.graph_cluster = SimpleNamespace(run=graph_run)

        runner.local_refinement("A", "Level2", [0.5], select_res=0.5)

        assert local_calls == [
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
        assert tuple(captured["adata"].obs_names) == ("sp-a1", "sp-a2")
        assert captured["evaluation_label"] == "Level2"
        assert captured["guidance_state"].level == "Level2"
        assert captured["guidance_state"].source == "local_anchoring:obsm[Level2]"
        assert tuple(captured["guidance_state"].category_labels) == ("a1", "a2")
        assert captured["problem_key"] == "standard-sc:A"
        assert captured["selected_resolution"] == 0.5


def test_runner_canonical_off_does_not_read_level2_state_despite_legacy_cost():
    with _isolated_runner_module() as module:
        collector = AssignmentGuidanceCollector()
        config = _config(collector, guidance="require")
        config.assignment_guidance_policy = "off"
        config.cell_type_col = "Level1"
        config.confidence_col = "Confidence"
        config.unknown_key = "Unknown"
        spatial = AnnData(
            X=np.ones((2, 2)),
            obs=pd.DataFrame(
                {"Level1": ["A", "A"]},
                index=["sp-a1", "sp-a2"],
            ),
            var=pd.DataFrame(index=["g1", "g2"]),
        )
        reference = AnnData(
            X=np.ones((2, 2)),
            obs=pd.DataFrame(
                {
                    "Level1": ["A", "A"],
                    "Level2": ["a1", "a2"],
                },
                index=["ref-a1", "ref-a2"],
            ),
            var=pd.DataFrame(index=["g1", "g2"]),
        )
        runner = module.ScSVC.__new__(module.ScSVC)
        runner.st_adata = spatial
        runner.sc_ref_adata = reference
        runner.config = config
        runner.logger = logging.getLogger("test-standard-sc-canonical-off")
        runner.cluster_col = "SVC_cluster"
        runner.local_annotate_method = SimpleNamespace(
            run=lambda target, _reference, **_kwargs: target.copy(),
            assignment_state=lambda *_args, **_kwargs: pytest.fail(
                "canonical off must not read optional Level2 assignment state"
            ),
        )
        captured = {}

        def graph_run(
            adata,
            _resolutions,
            _evaluation_label,
            *,
            guidance_state,
            **_kwargs,
        ):
            captured["guidance_state"] = guidance_state
            result = adata.copy()
            result.obs["leiden_0.5"] = pd.Categorical(["0", "1"])
            return result, pd.DataFrame(), 0.5

        runner.graph_cluster = SimpleNamespace(run=graph_run)

        runner.local_refinement("A", "Level2", [0.5], select_res=0.5)

        assert captured["guidance_state"] is None


def test_runner_missing_reference_records_not_applicable_before_local_solver():
    with _isolated_runner_module() as module:
        collector = AssignmentGuidanceCollector()
        config = _config(collector)
        config.cell_type_col = "Level1"
        spatial = AnnData(
            X=np.ones((2, 1)),
            obs=pd.DataFrame(
                {"Level1": ["A", "A"]},
                index=["sp-1", "sp-2"],
            ),
            var=pd.DataFrame(index=["g1"]),
        )
        reference = AnnData(
            X=np.ones((1, 1)),
            obs=pd.DataFrame(
                {"Level1": ["B"], "Level2": ["b1"]},
                index=["ref-1"],
            ),
            var=pd.DataFrame(index=["g1"]),
        )
        runner = module.ScSVC.__new__(module.ScSVC)
        runner.st_adata = spatial
        runner.sc_ref_adata = reference
        runner.config = config
        runner.logger = logging.getLogger("test-standard-sc-missing-reference")
        runner.cluster_col = "SVC_cluster"
        runner.local_annotate_method = SimpleNamespace(
            run=lambda *_args, **_kwargs: pytest.fail(
                "missing reference must bypass local solver"
            )
        )
        recorded = {}
        runner.graph_cluster = SimpleNamespace(
            record_not_applicable=lambda **kwargs: recorded.update(kwargs)
        )

        with pytest.raises(ValueError, match="no reference cells"):
            runner.local_refinement("A", "Level2", [0.5])

        assert recorded == {
            "problem_key": "standard-sc:A",
            "reason": NotApplicableReason.REFERENCE_UNAVAILABLE,
            "reason_details": {"role": "reference_cell"},
        }


def test_reference_compatibility_is_rejected_for_standard_sc_preflight():
    from revise.backend.policies import ModeValidationPolicy
    from revise.config import load_raw_config, merge_unified_config

    merged = merge_unified_config(
        raw_config=load_raw_config("revise/revise.yaml"),
        profile="application_sc",
        runtime_overrides={},
        io_overrides={},
        algorithm_overrides={
            "local_refinement": {
                "guidance": "prefer",
                "compatibility": {"mode": "reference"},
            }
        },
    )

    with pytest.raises(ValueError, match="graph_edge"):
        ModeValidationPolicy._validate_solver_compatibility(
            SimpleNamespace(
                merged_config=merged,
                runtime=merged["runtime"],
            )
        )


def test_graph_kernel_rejects_reference_before_graph_or_guidance_attempt(
    graph_module,
    monkeypatch,
):
    collector = AssignmentGuidanceCollector()
    config = _config(
        collector,
        compatibility_mode="reference",
    )
    monkeypatch.setattr(
        graph_module.sc,
        "pp",
        SimpleNamespace(
            normalize_total=lambda *_args, **_kwargs: pytest.fail(
                "reference must fail before graph preprocessing"
            )
        ),
    )
    kernel = graph_module.GraphClusterKernel(
        config,
        logging.getLogger("test-standard-sc-reference-defense"),
    )

    with pytest.raises(
        ValueError,
        match="graph_edge local refinement does not support reference",
    ):
        kernel.run(
            _adata(),
            resolution=[0.5],
            label="Level2",
            guidance_state=_soft_level2(_adata()),
            problem_key="standard-sc:A",
        )

    assert collector.events == []


def test_evaluation_label_and_assignment_level_are_independent(
    graph_module,
    monkeypatch,
):
    collector = AssignmentGuidanceCollector()
    adata = _adata()
    adata.obsm["evaluation_soft"] = pd.DataFrame(
        np.full((4, 2), 0.5),
        index=adata.obs_names,
        columns=["eval-a", "eval-b"],
    )
    _base, captured = _patch_graph_runtime(graph_module, monkeypatch)
    kernel = graph_module.GraphClusterKernel(
        _config(collector),
        logging.getLogger("test-standard-sc-independent-evaluation-label"),
    )

    kernel.run(
        adata,
        resolution=[0.5, 1.0],
        label="evaluation_soft",
        guidance_state=_soft_level2(adata),
        problem_key="standard-sc:A",
    )

    assert [call["guided"] for call in captured["leiden"]] == [
        False,
        False,
        True,
    ]
    [event] = collector.events
    assert event["outcome"] == "applied"
    assert event["left_assignment"]["level"] == "Level2"


def test_guided_graph_preserves_duplicate_self_loop_and_directed_edge_semantics(
    graph_module,
):
    collector = AssignmentGuidanceCollector()
    config = _config(collector, strength=2.0)
    kernel = graph_module.GraphClusterKernel(
        config,
        logging.getLogger("test-standard-sc-sparse-edge-semantics"),
    )
    graph = sparse.coo_matrix(
        (
            np.array([1.0, 2.0, 4.0, 5.0]),
            (
                np.array([0, 0, 1, 2]),
                np.array([1, 1, 1, 0]),
            ),
        ),
        shape=(3, 3),
    )
    state = one_hot_assignment(
        ("A", "A", "B"),
        observation_labels=("cell-0", "cell-1", "cell-2"),
        category_labels=("A", "B"),
        source="local_anchoring:argmax[Level2]",
        level="Level2",
    )

    guided = kernel._guided_graph(graph, state)

    assert sparse.isspmatrix_csr(guided)
    assert guided.nnz == 3
    np.testing.assert_allclose(
        guided.toarray(),
        [
            [0.0, 3.0, 0.0],
            [0.0, 4.0, 0.0],
            [0.05, 0.0, 0.0],
        ],
    )


def test_runner_logs_final_guided_cluster_count_at_user_selected_resolution(
    graph_module,
    monkeypatch,
):
    with _isolated_runner_module() as module:
        collector = AssignmentGuidanceCollector()
        config = _config(collector)
        config.cell_type_col = "Level1"
        config.confidence_col = "Confidence"
        config.unknown_key = "Unknown"
        spatial = AnnData(
            X=np.ones((4, 3)),
            obs=pd.DataFrame(
                {"Level1": ["A"] * 4},
                index=["sp-1", "sp-2", "sp-3", "sp-4"],
            ),
            var=pd.DataFrame(index=["g1", "g2", "g3"]),
        )
        reference = AnnData(
            X=np.ones((4, 3)),
            obs=pd.DataFrame(
                {
                    "Level1": ["A"] * 4,
                    "Level2": ["a1", "a1", "a2", "a2"],
                },
                index=["ref-1", "ref-2", "ref-3", "ref-4"],
            ),
            var=pd.DataFrame(index=["g1", "g2", "g3"]),
        )
        runner = module.ScSVC.__new__(module.ScSVC)
        runner.st_adata = spatial
        runner.sc_ref_adata = reference
        runner.config = config
        runner.cluster_col = "SVC_cluster"
        messages = []
        runner.logger = SimpleNamespace(
            info=lambda message, *args: messages.append(
                message % args if args else message
            ),
            warning=lambda *_args, **_kwargs: None,
        )

        class LocalProducer:
            def run(self, target, _reference, **kwargs):
                result = target.copy()
                if kwargs["cell_type_col"] == "Level2":
                    result.obsm["Level2"] = pd.DataFrame(
                        [
                            [0.9, 0.1],
                            [0.8, 0.2],
                            [0.1, 0.9],
                            [0.2, 0.8],
                        ],
                        index=result.obs_names,
                        columns=["a1", "a2"],
                    )
                    result.obs["Level2"] = ["a1", "a1", "a2", "a2"]
                return result

            @staticmethod
            def assignment_state(annotated, key):
                posterior = annotated.obsm[key]
                return AssignmentState(
                    values=posterior.to_numpy(),
                    observation_labels=posterior.index,
                    category_labels=posterior.columns,
                    source=f"local_anchoring:obsm[{key}]",
                    level=key,
                    value_semantics="soft",
                    lineage=[],
                )

        runner.local_annotate_method = LocalProducer()
        _base, captured = _patch_graph_runtime(
            graph_module,
            monkeypatch,
            guided_labels_by_resolution={
                0.5: ["0", "1", "2", "2"],
            },
        )
        runner.graph_cluster = graph_module.GraphClusterKernel(
            config,
            runner.logger,
        )

        spatial_result, _expression_result = runner.local_refinement(
            "A",
            "Level2",
            [0.5, 1.0],
            select_res=0.5,
        )

        assert [
            (call["resolution"], call["guided"])
            for call in captured["leiden"]
        ] == [(0.5, False), (1.0, False), (0.5, True)]
        assert spatial_result.obs["SVC_cluster"].nunique() == 3
        assert "resolution 0.5 got cluster number 3" in messages


@pytest.mark.parametrize(
    ("enable_guidance", "expected_outcome", "expected_reads"),
    [
        (False, "off", 0),
        (True, "applied", 1),
    ],
)
def test_public_standard_sc_route_reaches_framework_manifest_with_level2_outcome(
    graph_module,
    monkeypatch,
    tmp_path,
    enable_guidance,
    expected_outcome,
    expected_reads,
):
    from revise.backend import adapters
    from revise.backend.policies import (
        ModeEvaluationPolicy,
        ModeValidationPolicy,
    )
    from revise.backend.registry import StrategyRegistry
    from revise.framework import REVISEPipeline
    from revise.recon.pipeline import UnifiedReconstructionPipeline

    application_svc = types.ModuleType(
        "revise.backend.runners.application_svc"
    )

    class ApplicationSVC:
        pass

    application_svc.ApplicationSVC = ApplicationSVC
    monkeypatch.setitem(
        sys.modules,
        application_svc.__name__,
        application_svc,
    )
    bio = types.ModuleType("revise.analysis.bio")
    bio.get_degs = lambda *_args, **_kwargs: None
    bio.conclusions_write = lambda *_args, **_kwargs: None
    bio.plot_volcano = lambda *_args, **_kwargs: None
    monkeypatch.setitem(sys.modules, bio.__name__, bio)
    kernels_package = importlib.import_module("revise.backend.kernels")
    monkeypatch.setitem(
        kernels_package.__dict__,
        "LocalAnchoringKernel",
        object,
    )
    from revise.backend.runners import sc_svc_application as runner_module

    spatial = AnnData(
        X=np.ones((4, 3)),
        obs=pd.DataFrame(
            {
                "transcript_counts": [100, 100, 100, 100],
                "Level1": ["A", "A", "A", "A"],
            },
            index=["sp-1", "sp-2", "sp-3", "sp-4"],
        ),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
    )
    spatial.obsm["Level1"] = pd.DataFrame(
        [[0.0, 1.0]] * 4,
        index=spatial.obs_names,
        columns=["A", "B"],
    )
    reference = AnnData(
        X=np.ones((4, 3)),
        obs=pd.DataFrame(
            {
                "Patient": ["sample"] * 4,
                "Level1": ["A"] * 4,
                "Level2": ["a1", "a1", "a2", "a2"],
            },
            index=["ref-1", "ref-2", "ref-3", "ref-4"],
        ),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
    )
    captured = {
        "registry_strategy": None,
        "routes": [],
        "assignment_reads": 0,
    }

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
    monkeypatch.setattr(adapters, "_install_safe_topology_patch", lambda: None)
    monkeypatch.setattr(
        adapters,
        "_input_service",
        lambda _ctx: SimpleNamespace(
            read_st_adata=lambda _path: spatial.copy(),
            read_sc_ref_adata=lambda _path: reference.copy(),
        ),
    )
    monkeypatch.setattr(
        adapters.sc.pp,
        "filter_genes",
        lambda *_args, **_kwargs: None,
        raising=False,
    )

    real_registry_get = StrategyRegistry.get

    def registry_get(registry, strategy_id):
        captured["registry_strategy"] = strategy_id
        return real_registry_get(registry, strategy_id)

    monkeypatch.setattr(StrategyRegistry, "get", registry_get)

    def global_anchoring(strategy, ctx):
        strategy._attach_assignment_guidance_callback(ctx)
        captured["routes"].append(ctx.route_key)

    monkeypatch.setattr(
        adapters.ScSvcApplicationStrategy,
        "global_anchoring",
        global_anchoring,
    )

    _base_graph, graph_calls = _patch_graph_runtime(graph_module, monkeypatch)

    class LocalProducer:
        def run(self, target, _reference, **kwargs):
            result = target.copy()
            if kwargs["cell_type_col"] == "Level2":
                result.obsm["Level2"] = pd.DataFrame(
                    [
                        [0.9, 0.1],
                        [0.8, 0.2],
                        [0.1, 0.9],
                        [0.2, 0.8],
                    ],
                    index=result.obs_names,
                    columns=["a1", "a2"],
                )
                result.obs["Level2"] = ["a1", "a1", "a2", "a2"]
            return result

        def assignment_state(self, annotated, key):
            captured["assignment_reads"] += 1
            posterior = annotated.obsm[key]
            return AssignmentState(
                values=posterior.to_numpy(),
                observation_labels=posterior.index,
                category_labels=posterior.columns,
                source=f"local_anchoring:obsm[{key}]",
                level=key,
                value_semantics="soft",
                lineage=[{"operation": "local_anchoring"}],
            )

    def runner_init(self, st_adata, sc_ref_adata, config, logger):
        self.st_adata = st_adata.copy()
        self.sc_ref_adata = sc_ref_adata.copy()
        self.config = config
        self.logger = logger
        self.cluster_col = "SVC_cluster"
        self.local_annotate_method = LocalProducer()
        self.graph_cluster = graph_module.GraphClusterKernel(config, logger)

    monkeypatch.setattr(runner_module.ScSVC, "__init__", runner_init)

    output_root = tmp_path / f"output-{enable_guidance}"
    algorithm_overrides = {
        "preprocess": {"st_min_transcripts": 0, "st_min_cells": 0, "sc_min_cells": 0},
        "graph": {"method": "pca", "alpha": 0.0},
        "sc": {
            "select_ct": "A",
            "resolutions": [0.5, 1.0],
            "select_resolution": None,
        },
    }
    if enable_guidance:
        algorithm_overrides["local_refinement"] = {
            "guidance": "prefer",
            "compatibility": {"mode": "cost"},
        }

    pipeline = REVISEPipeline("revise/revise.yaml")
    pipeline._run_with_algorithm_overrides(
        profile="application_sc",
        runtime_overrides={},
        io_overrides={
            "data_root": str(tmp_path / "data"),
            "output_root": str(output_root),
            "sample_name": "sample",
        },
        algorithm_overrides=algorithm_overrides,
        dry_run=False,
    )

    paths = list(output_root.rglob("provenance.json"))
    assert len(paths) == 1
    manifest = json.loads(paths[0].read_text())
    assert captured["registry_strategy"] == "ScSvcApplicationStrategy"
    assert captured["routes"] == ["sc_svc:segmentation"]
    assert captured["assignment_reads"] == expected_reads
    [event] = manifest["assignment_guidance"]["events"]
    assert event["route"] == "sc_svc:segmentation"
    assert event["operator"] == "graph_edge"
    assert event["outcome"] == expected_outcome
    assert manifest["assignment_guidance"]["resolved"]["guidance"] == (
        "prefer" if enable_guidance else "off"
    )
    if enable_guidance:
        assert event["left_assignment"]["level"] == "Level2"
        assert [call["guided"] for call in graph_calls["leiden"]] == [
            False,
            False,
            True,
        ]
    else:
        assert event["left_assignment"] is None
        assert [call["guided"] for call in graph_calls["leiden"]] == [
            False,
            False,
        ]


def test_all_cell_adapter_records_not_applicable_for_bypassed_cohorts(
    graph_module,
    monkeypatch,
):
    from revise.backend import adapters

    collector = AssignmentGuidanceCollector()
    config = _config(collector)
    config.cell_type_col = "Level1"
    spatial = AnnData(
        X=np.ones((3, 1)),
        obs=pd.DataFrame(
            {"Level1": ["A", "B", "B"]},
            index=["sp-a", "sp-b1", "sp-b2"],
        ),
        var=pd.DataFrame(index=["g1"]),
    )
    reference = AnnData(
        X=np.ones((1, 1)),
        obs=pd.DataFrame(
            {"Level1": ["A"], "Level2": ["a1"]},
            index=["ref-a"],
        ),
        var=pd.DataFrame(index=["g1"]),
    )
    runner = SimpleNamespace(
        st_adata=spatial,
        sc_ref_adata=reference,
        config=config,
        graph_cluster=graph_module.GraphClusterKernel(
            config,
            logging.getLogger("test-standard-sc-all-cell-bypass"),
        ),
        local_refinement=lambda *_args, **_kwargs: pytest.fail(
            "singleton and missing-reference cohorts must bypass local refinement"
        ),
    )

    def singleton_outputs(_runner, candidate, cell_type_col, sub_cell_type_col):
        part = spatial[spatial.obs[cell_type_col] == candidate].copy()
        part.obs[sub_cell_type_col] = [f"{candidate}-singleton"]
        part.obs["SVC_cluster"] = pd.Categorical(["0"])
        expr_part = part.copy()
        expr_part.obsm["SVC_cluster"] = pd.DataFrame(
            [[1.0]],
            index=expr_part.obs_names,
            columns=["0"],
        )
        return part, expr_part

    monkeypatch.setattr(
        adapters,
        "_singleton_sc_svc_outputs",
        singleton_outputs,
    )
    monkeypatch.setattr(
        adapters.sc,
        "concat",
        anndata_concat,
        raising=False,
    )
    ctx = SimpleNamespace(
        merged_config={
            "sc": {
                "select_ct": "all",
                "resolutions": [0.5],
                "select_resolution": None,
            }
        },
        columns={
            "cell_type_col": "Level1",
            "sub_cell_type_col": "Level2",
        },
        runner=runner,
        logger=logging.getLogger("test-standard-sc-all-cell-strategy"),
        artifacts={},
    )

    adapters.ScSvcApplicationStrategy().solve_ot(ctx)

    assert {
        event["problem_key"]: (event["outcome"], event["reason"])
        for event in collector.events
    } == {
        "standard-sc:A": (
            "not_applicable",
            NotApplicableReason.INSUFFICIENT_UNITS.value,
        ),
        "standard-sc:B": (
            "not_applicable",
            NotApplicableReason.REFERENCE_UNAVAILABLE.value,
        ),
    }
