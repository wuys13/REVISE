from __future__ import annotations

import logging
import sys
import types
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData


_ISOLATED_PREFIXES = (
    "scanpy",
    "revise.backend.adapters",
    "revise.backend.runners.sc_svc_application",
)
_MISSING = object()


def _isolated_module_names():
    return tuple(
        name
        for name in sys.modules
        if any(
            name == prefix or name.startswith(f"{prefix}.")
            for prefix in _ISOLATED_PREFIXES
        )
    )


@pytest.fixture(autouse=True)
def _restore_sc_test_modules():
    names = _isolated_module_names()
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
    yield
    current = _isolated_module_names()
    for name in current:
        sys.modules.pop(name, None)
    for name in current:
        parent_name, separator, attribute = name.rpartition(".")
        parent = sys.modules.get(parent_name)
        if separator and parent is not None and hasattr(parent, attribute):
            delattr(parent, attribute)
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


def _import_sc_svc(monkeypatch):
    scanpy = types.ModuleType("scanpy")
    scanpy.pp = SimpleNamespace()
    scanpy.pl = SimpleNamespace()
    scanpy.tl = SimpleNamespace()
    monkeypatch.setitem(sys.modules, "scanpy", scanpy)
    if "revise.backend.adapters" in sys.modules:
        monkeypatch.setattr(sys.modules["revise.backend.adapters"], "sc", scanpy)

    ot = types.ModuleType("ot")
    monkeypatch.setitem(sys.modules, "ot", ot)

    distance = types.ModuleType("revise.backend.ops.distance")
    distance.bhattacharyya_distance = lambda *args, **kwargs: None
    monkeypatch.setitem(sys.modules, "revise.backend.ops.distance", distance)

    graph_cluster = types.ModuleType("revise.backend.kernels.graph_cluster")

    class GraphClusterKernel:
        def __init__(self, config, logger):
            pass

    graph_cluster.GraphClusterKernel = GraphClusterKernel
    monkeypatch.setitem(
        sys.modules, "revise.backend.kernels.graph_cluster", graph_cluster
    )

    analysis = types.ModuleType("revise.analysis")
    analysis.__path__ = []
    monkeypatch.setitem(sys.modules, "revise.analysis", analysis)
    bio = types.ModuleType("revise.analysis.bio")
    bio.get_degs = lambda *args, **kwargs: None
    bio.conclusions_write = lambda *args, **kwargs: None
    bio.plot_volcano = lambda *args, **kwargs: None
    monkeypatch.setitem(sys.modules, "revise.analysis.bio", bio)

    from revise.backend.runners.sc_svc_application import ScSVC

    return ScSVC


def _adata(obs_names, level1, cell_type_col="Level1"):
    return AnnData(
        X=np.ones((len(obs_names), 2), dtype=np.float64),
        obs=pd.DataFrame({cell_type_col: level1}, index=obs_names),
        var=pd.DataFrame(index=["g1", "g2"]),
    )


@pytest.mark.parametrize("cell_type_col", ["Level1", "custom_level1"])
def test_ist_local_refinement_uses_configured_columns_and_local_ot(
    monkeypatch,
    cell_type_col,
):
    ScSVC = _import_sc_svc(monkeypatch)
    config = SimpleNamespace(
        annotate_mode="pot",
        rec_ot_method="tacco",
        cell_type_col=cell_type_col,
        confidence_col="Confidence",
        unknown_key="Unknown",
    )
    runner = ScSVC(
        _adata(["sp1", "sp2"], ["A", "A"], cell_type_col),
        _adata(["sc1", "sc2"], ["A", "A"], cell_type_col),
        config,
        logger=None,
    )

    def fail_global(*args, **kwargs):
        raise AssertionError("sc-SVC local refinement reused the global annotation kernel")

    runner.annotate_method.run = fail_global
    local_calls = []

    def run_local(target, reference, **kwargs):
        local_calls.append(kwargs["cell_type_col"])
        return target.copy()

    runner.local_annotate_method.run = run_local

    clustered = _adata(["sp1", "sp2"], ["A", "A"], cell_type_col)
    clustered.obs["leiden_0.5"] = pd.Categorical(["0", "1"])
    merge_df = pd.DataFrame({"resolution": [0.5], "cluster_num": [2]})
    runner.graph_cluster.run = lambda *args, **kwargs: (clustered, merge_df, 0.5)

    runner.local_refinement("A", "Level2", [0.5], select_res=0.5)

    assert local_calls == ["Level2", "SVC_cluster"]
    assert runner.local_annotate_method.method == "tacco"


def test_application_sc_config_carries_local_ot_method():
    from revise.config.runner_conf import ApplicationScConf

    config = ApplicationScConf(
        sample_name="sample",
        raw_data_path="data",
        result_root_path="output",
        cell_type_col="Level1",
        confidence_col="Confidence",
        unknown_key="Unknown",
        st_file="sp.h5ad",
        sc_ref_file="sc.h5ad",
        annotate_mode="pot",
        rec_ot_method="tacco",
        tacco_annotate_multi_center=1,
        tacco_annotate_lamb=0.001,
    )

    assert config.annotate_mode == "pot"
    assert config.rec_ot_method == "tacco"


def test_application_sc_passes_configured_tacco_parameters_to_all_three_calls(
    monkeypatch,
):
    from revise.backend.kernels import global_anchoring, local_anchoring
    from revise.config.runner_conf import ApplicationScConf

    config = ApplicationScConf(
        sample_name="sample",
        raw_data_path="data",
        result_root_path="output",
        cell_type_col="Level1",
        confidence_col="Confidence",
        unknown_key="Unknown",
        st_file="sp.h5ad",
        sc_ref_file="sc.h5ad",
        annotate_mode="tacco",
        rec_ot_method="tacco",
        tacco_annotate_multi_center=1,
        tacco_annotate_lamb=0.001,
    )
    calls = []

    def annotate(
        adata,
        reference,
        annotation_key,
        *,
        result_key,
        return_reference,
        multi_center,
        lamb,
    ):
        calls.append(
            {
                "annotation_key": annotation_key,
                "return_reference": return_reference,
                "multi_center": multi_center,
                "lamb": lamb,
            }
        )
        categories = pd.Index(
            pd.unique(reference.obs[annotation_key].astype(str))
        )
        adata.obsm[result_key] = pd.DataFrame(
            np.full((adata.n_obs, len(categories)), 1.0 / len(categories)),
            index=adata.obs_names,
            columns=categories,
        )
        return adata, reference

    monkeypatch.setattr(
        global_anchoring,
        "require_tacco",
        lambda: SimpleNamespace(tl=SimpleNamespace(annotate=annotate)),
    )

    target = AnnData(
        X=np.array([[2.0, 1.0], [1.0, 2.0]]),
        obs=pd.DataFrame(index=["sp1", "sp2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    level1_reference = AnnData(
        X=np.array([[2.0, 0.0], [0.0, 2.0]]),
        obs=pd.DataFrame({"Level1": ["A", "B"]}, index=["sc1", "sc2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    global_anchoring.GlobalAnchoringKernel(
        config,
        logging.getLogger("test"),
    ).run(target, level1_reference, cell_type_col="Level1")

    level2_reference = level1_reference.copy()
    level2_reference.obs["Level2"] = ["A1", "A2"]
    local = local_anchoring.LocalAnchoringKernel(
        config,
        logging.getLogger("test"),
    )
    spatial = local.run(target, level2_reference, cell_type_col="Level2")
    spatial.obs["SVC_cluster"] = pd.Categorical(["0", "1"])
    local.run(level1_reference, spatial, cell_type_col="SVC_cluster")

    assert [call["annotation_key"] for call in calls] == [
        "Level1",
        "Level2",
        "SVC_cluster",
    ]
    assert all(call["return_reference"] is True for call in calls)
    assert all(call["multi_center"] == 1 for call in calls)
    assert all(call["lamb"] == 0.001 for call in calls)


@pytest.mark.parametrize("method", ["pot", "tacco"])
def test_local_anchoring_routes_normalized_problem_to_shared_solver(
    monkeypatch, method
):
    _import_sc_svc(monkeypatch)
    from revise.backend.kernels import local_anchoring

    config = SimpleNamespace(
        rec_ot_method=method,
        rec_pot_reg=0.2,
        rec_pot_reg_m=0.3,
        rec_pot_reg_type="kl",
        cell_type_col="Level1",
        confidence_col="Confidence",
        unknown_key="Unknown",
    )
    kernel = local_anchoring.LocalAnchoringKernel(config)
    target = AnnData(
        X=np.array([[2.0, 0.0], [0.0, 2.0]]),
        obs=pd.DataFrame(index=["sp1", "sp2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    reference = AnnData(
        X=np.array([[1.0, 0.0], [0.0, 1.0]]),
        obs=pd.DataFrame({"Level2": ["A", "B"]}, index=["sc1", "sc2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    monkeypatch.setattr(
        local_anchoring,
        "bhattacharyya_distance",
        lambda profiles, expression: np.array([[0.0, 2.0], [1.0, 3.0]]),
    )
    captured = {}

    def solve(source, target_mass, cost, **kwargs):
        captured.update(
            source=np.asarray(source),
            target=np.asarray(target_mass),
            cost=np.asarray(cost),
            kwargs=kwargs,
        )
        return np.array([[0.5, 0.0], [0.0, 0.5]])

    monkeypatch.setattr(local_anchoring, "solve_local_ot", solve)

    result = kernel.run(target, reference, cell_type_col="Level2")

    np.testing.assert_allclose(captured["source"], [0.5, 0.5])
    np.testing.assert_allclose(captured["target"], [0.5, 0.5])
    np.testing.assert_allclose(
        captured["cost"],
        [[0.0, 1.0 / 3.0], [2.0 / 3.0, 1.0]],
    )
    assert captured["kwargs"] == {
        "method": method,
        "pot_reg": 0.2,
        "pot_reg_m": 0.3,
        "pot_reg_type": "kl",
        "pot_verbose": False,
        "pot_num_iter_max": 5000,
    }
    assert result.obs["Level2"].tolist() == ["A", "B"]
    assert tuple(result.obsm["Level2"].index) == ("sp1", "sp2")
    assert tuple(result.obsm["Level2"].columns) == ("A", "B")
    np.testing.assert_allclose(
        result.obsm["Level2"].to_numpy(), [[1.0, 0.0], [0.0, 1.0]]
    )
    np.testing.assert_allclose(result.obs["Confidence"], [1.0, 1.0])


@pytest.mark.parametrize(
    ("cell_type_col", "sub_cell_type_col"),
    [("Level1", "Level2"), ("custom_level1", "custom_level2")],
)
def test_ist_adapter_propagates_configured_columns_and_local_ot(
    monkeypatch,
    tmp_path,
    cell_type_col,
    sub_cell_type_col,
):
    _import_sc_svc(monkeypatch)
    import scanpy
    from revise.backend import adapters

    scanpy.pp.filter_genes = lambda *args, **kwargs: None
    st_adata = AnnData(
        X=np.ones((2, 2), dtype=np.float64),
        obs=pd.DataFrame(index=["sp1", "sp2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    sc_adata = AnnData(
        X=np.ones((2, 2), dtype=np.float64),
        obs=pd.DataFrame(
            {
                "Patient": ["sample", "sample"],
                cell_type_col: ["A/B", "A/B"],
                sub_cell_type_col: ["A/1", "A/2"],
            },
            index=["sc1", "sc2"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    if cell_type_col != "Level1":
        sc_adata.obs["Level1"] = ["legacy", "legacy"]
        sc_adata.obs["Level2"] = ["legacy1", "legacy2"]

    class InputService:
        def read_st_adata(self, path):
            return st_adata.copy()

        def read_sc_ref_adata(self, path):
            return sc_adata.copy()

    monkeypatch.setattr(adapters, "_input_service", lambda ctx: InputService())
    ctx = SimpleNamespace(
        merged_config={
            "ot": {
                "ga": {
                    "solver": "pot",
                    "pot": {"reg": 0.1, "reg_m": 0.0, "reg_type": "entropy"},
                },
                "lr": {
                    "solver": "tacco",
                    "pot": {"reg": 0.2, "reg_m": 0.3, "reg_type": "kl"},
                },
                "impute": {"reg": 5.0, "reg_m": 0.0, "reg_type": "kl"},
            },
            "preprocess": {},
            "graph": {},
            "sc": {
                "tacco_annotate": {
                    "multi_center": 1,
                    "lamb": 0.001,
                }
            },
        },
        io={
            "sample_name": "sample",
            "data_root": str(tmp_path),
            "st_file": "sp.h5ad",
            "sc_ref_file": "sc.h5ad",
            "patient_key": "Patient",
        },
        columns={
            "cell_type_col": cell_type_col,
            "sub_cell_type_col": sub_cell_type_col,
            "confidence_col": "Confidence",
            "unknown_key": "Unknown",
        },
        run_dir=tmp_path,
        logger=None,
    )

    adapters.ScSvcApplicationStrategy().prepare_context(ctx)

    assert ctx.runner_config.annotate_mode == "pot"
    assert ctx.runner_config.rec_ot_method == "tacco"
    assert ctx.runner_config.tacco_annotate_multi_center == 1
    assert ctx.runner_config.tacco_annotate_lamb == 0.001
    assert ctx.runner.local_annotate_method.method == "tacco"
    assert list(ctx.runner.sc_ref_adata.obs.columns) == [
        cell_type_col,
        sub_cell_type_col,
    ]
    assert ctx.runner.sc_ref_adata.obs[cell_type_col].tolist() == ["A_B", "A_B"]
    assert ctx.runner.sc_ref_adata.obs[sub_cell_type_col].tolist() == ["A_1", "A_2"]


@pytest.mark.parametrize(
    "select_ct",
    [None, "", " ", "all", "*", "__all__", "all_cell_types"],
)
def test_ist_adapter_requires_one_concrete_cell_type(
    monkeypatch,
    select_ct,
):
    _import_sc_svc(monkeypatch)
    from revise.backend import adapters

    ctx = SimpleNamespace(
        merged_config={
            "sc": {
                "select_ct": select_ct,
                "resolutions": [0.5],
            }
        },
        columns={"cell_type_col": "Level1", "sub_cell_type_col": "Level2"},
        runner=SimpleNamespace(
            st_adata=_adata(["sp1", "sp2"], ["A", "A"]),
            local_refinement=lambda *_args, **_kwargs: pytest.fail(
                "invalid selection reached local refinement"
            ),
        ),
        logger=logging.getLogger("test-concrete-sc-selection"),
    )

    strategy = adapters.ScSvcApplicationStrategy()
    with pytest.raises(
        ValueError,
        match="route.select_cell_type must name one concrete broad cell type",
    ):
        strategy.solve_ot(ctx)


def test_ist_adapter_refines_only_the_selected_cell_type(
    monkeypatch,
):
    _import_sc_svc(monkeypatch)
    from revise.backend import adapters

    spatial = _adata(["sp1", "sp2"], ["T", "T"])
    expression = _adata(["sc1", "sc2"], ["T", "T"])
    calls = []

    def local_refinement(select_ct, sub_cell_type_col, resolutions, select_res=None):
        calls.append((select_ct, sub_cell_type_col, resolutions, select_res))
        return spatial, expression

    applied = []
    ctx = SimpleNamespace(
        merged_config={
            "sc": {
                "select_ct": "T",
                "resolutions": [0.5],
                "select_resolution": 0.5,
            }
        },
        columns={"cell_type_col": "Level1", "sub_cell_type_col": "Level2"},
        runner=SimpleNamespace(local_refinement=local_refinement),
        logger=logging.getLogger("test-single-sc-selection"),
        artifacts={},
        record_local_refinement=applied.append,
    )

    adapters.ScSvcApplicationStrategy().solve_ot(ctx)

    assert calls == [("T", "Level2", [0.5], 0.5)]
    assert applied == [True]
    assert ctx.artifacts["outputs"] == {
        "sc_svc_spatial": spatial,
        "sc_svc_expr": expression,
    }
    assert ctx.artifacts["selected_cell_type"] == "T"


@pytest.mark.parametrize("method", ["pot", "tacco"])
def test_application_ot_method_switches_global_and_local_together(method):
    from reconstruct import _engine_overrides

    request = SimpleNamespace(
        svc_type="sc-SVC",
        ot_method=method,
        local_refinement_strength=None,
        select_cell_type="T",
        broad_column="Level1",
        subtype_column="Level2",
    )

    overrides = _engine_overrides(
        SimpleNamespace(
            ot_method=request.ot_method,
            broad_column=request.broad_column,
            subtype_column=request.subtype_column,
            select_cell_type=request.select_cell_type,
            local_refinement_strength=request.local_refinement_strength,
            seed=None,
            st_path=Path("st"),
            reference_path=Path("ref"),
            pm_on_cell_path=None,
            output_dir=Path("out"),
            output_name="sample",
            st_format="h5ad",
            spatialdata_table=None,
            spatialdata_element=None,
        )
    )[2]

    assert overrides["ot"]["ga"]["solver"] == method
    assert overrides["ot"]["lr"]["solver"] == method
