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
    "ot",
    "scanpy",
    "revise.analysis",
    "revise.application.preprocess",
    "revise.backend.adapters",
    "revise.backend.kernels",
    "revise.backend.ops",
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
    kernels_package = sys.modules.get("revise.backend.kernels")
    kernel_exports = {}
    if kernels_package is not None:
        declared_exports = vars(kernels_package).get("_KERNEL_EXPORTS", {})
        kernel_exports = {
            name: vars(kernels_package).get(name, _MISSING)
            for name in declared_exports
        }
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
    if kernels_package is not None:
        for name, value in kernel_exports.items():
            if value is _MISSING:
                vars(kernels_package).pop(name, None)
            else:
                vars(kernels_package)[name] = value


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
    reference = _adata(["sc1", "sc2"], ["A", "A"], cell_type_col)
    reference.obs["Level2"] = ["A1", "A2"]
    runner = ScSVC(
        _adata(["sp1", "sp2"], ["A", "A"], cell_type_col),
        reference,
        config,
        logger=None,
    )

    assert not hasattr(runner, "annotate_method")
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
        annotate_pot_reg=0.1,
        annotate_pot_reg_m=0.0,
        annotate_pot_reg_type="entropy",
        tacco_annotate_multi_center=1,
        tacco_annotate_lamb=0.001,
        rec_graph_n_neighbors=10,
        rec_graph_exp_neighbor_num=15,
        rec_graph_spatial_neighbor_num=6,
        rec_graph_method="joint",
        rec_graph_alpha=0.2,
        rec_random_state=42,
        rec_pot_reg=0.1,
        rec_pot_reg_m=0.0,
        rec_pot_reg_type="entropy",
        rec_alpha=0.5,
    )

    assert config.annotate_mode == "pot"
    assert config.rec_ot_method == "tacco"


def test_application_sc_passes_configured_tacco_parameters_to_all_three_calls(
    monkeypatch,
):
    from revise.backend.kernels import global_anchoring, local_anchoring, ot
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
        annotate_pot_reg=0.1,
        annotate_pot_reg_m=0.0,
        annotate_pot_reg_type="entropy",
        tacco_annotate_multi_center=1,
        tacco_annotate_lamb=0.001,
        rec_graph_n_neighbors=10,
        rec_graph_exp_neighbor_num=15,
        rec_graph_spatial_neighbor_num=6,
        rec_graph_method="joint",
        rec_graph_alpha=0.2,
        rec_random_state=42,
        rec_pot_reg=0.1,
        rec_pot_reg_m=0.0,
        rec_pot_reg_type="entropy",
        rec_alpha=0.5,
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
        ot,
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


def test_application_sc_pot_local_annotation_uses_annotation_contract(monkeypatch):
    from revise.backend.kernels import local_anchoring
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
        rec_ot_method="pot",
        annotate_pot_reg=0.1,
        annotate_pot_reg_m=0.0,
        annotate_pot_reg_type="entropy",
        tacco_annotate_multi_center=1,
        tacco_annotate_lamb=0.001,
        rec_graph_n_neighbors=10,
        rec_graph_exp_neighbor_num=15,
        rec_graph_spatial_neighbor_num=6,
        rec_graph_method="joint",
        rec_graph_alpha=0.2,
        rec_random_state=0,
        rec_pot_reg=0.2,
        rec_pot_reg_m=0.3,
        rec_pot_reg_type="kl",
        rec_alpha=0.5,
    )
    captured = {}

    def annotate(target, reference, **kwargs):
        captured.update(kwargs)
        result = target.copy()
        result.obs[kwargs["annotation_key"]] = "A"
        result.obs[kwargs["confidence_key"]] = 1.0
        return result

    monkeypatch.setattr(
        local_anchoring,
        "OTKernel",
        SimpleNamespace(annotate=annotate),
    )
    target = _adata(["sp1", "sp2"], ["A", "A"])
    reference = _adata(["sc1", "sc2"], ["A", "B"])
    reference.obs["Level2"] = ["A1", "B1"]

    result = local_anchoring.LocalAnchoringKernel(config).run(
        target,
        reference,
        cell_type_col="Level2",
    )

    assert captured == {
        "method": "pot",
        "annotation_key": "Level2",
        "confidence_key": "Confidence",
        "pot_reg": 0.2,
        "pot_reg_m": 0.3,
        "pot_reg_type": "kl",
        "pot_verbose": False,
        "pot_num_iter_max": 5000,
        "multi_center": None,
        "lamb": None,
        "unknown_key": "Unknown",
    }
    assert result.obs["Level2"].tolist() == ["A", "A"]


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
        X=np.array([[1.0, 1.0], [0.0, 0.0]], dtype=np.float64),
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

    from revise.application.preprocess import prepare_sc_svc_pair

    prepared_st, prepared_sc = prepare_sc_svc_pair(
        st_adata,
        sc_adata,
        broad_column=cell_type_col,
        subtype_column=sub_cell_type_col,
    )
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
                "graph": {
                    "method": "joint",
                    "alpha": 0.2,
                    "n_neighbors": 10,
                    "exp_neighbors": 15,
                    "spatial_neighbors": 6,
                    "random_state": 0,
                },
            "reconstruct": {"alpha": 0.5},
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
        runtime={"seed": 42},
        input_specs=(
            SimpleNamespace(role="st", path="sp.h5ad"),
            SimpleNamespace(role="sc_ref", path="sc.h5ad"),
        ),
        logger=None,
        st_adata=prepared_st,
        sc_ref_adata=prepared_sc,
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
    ["", " ", "all", "*", "__all__", "all_cell_types"],
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
    spatial.obs["Level2"] = ["T1", "T2"]
    spatial.obs["SVC_cluster"] = pd.Categorical(["0", "1"])
    expression = _adata(["sc1", "sc2"], ["T", "T"])
    expression.obs["SVC_cluster"] = pd.Categorical(["0", "1"])
    reference = _adata(["ref1", "ref2"], ["T", "T"])
    reference.obs["Level2"] = ["T1", "T2"]
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
        runner=SimpleNamespace(
            st_adata=_adata(["sp1", "sp2"], ["T", "T"]),
            sc_ref_adata=reference,
            local_refinement=local_refinement,
        ),
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
    assert ctx.artifacts["processed_cell_types"] == ["T"]
    assert ctx.artifacts["skipped_cell_types"] == []
    assert ctx.artifacts["label_assignments"]["Level2"].tolist() == ["T1", "T2"]


def test_explicit_single_type_keeps_legacy_one_subtype_refinement(monkeypatch):
    _import_sc_svc(monkeypatch)
    from revise.backend import adapters

    st = _adata(["sp1", "sp2"], ["T", "T"])
    reference = _adata(["ref1", "ref2"], ["T", "T"])
    reference.obs["Level2"] = ["T1", "T1"]
    spatial = st.copy()
    spatial.obs["Level2"] = ["T1", "T1"]
    spatial.obs["SVC_cluster"] = pd.Categorical(["0", "0"])
    expression = reference.copy()
    expression.obs["SVC_cluster"] = pd.Categorical(["0", "0"])
    calls = []
    ctx = SimpleNamespace(
        merged_config={"sc": {"select_ct": "T", "resolutions": [0.5]}},
        columns={"cell_type_col": "Level1", "sub_cell_type_col": "Level2"},
        runner=SimpleNamespace(
            st_adata=st,
            sc_ref_adata=reference,
            local_refinement=lambda *_args, **_kwargs: (
                calls.append("T") or (spatial, expression)
            ),
        ),
        logger=logging.getLogger("test-legacy-one-subtype"),
        artifacts={},
        record_local_refinement=lambda _value: None,
    )

    adapters.ScSvcApplicationStrategy().solve_ot(ctx)

    assert calls == ["T"]
    assert ctx.artifacts["processed_cell_types"] == ["T"]


def test_ist_adapter_runs_all_eligible_actual_types_once_and_namespaces_clusters(
    monkeypatch,
):
    _import_sc_svc(monkeypatch)
    from revise.backend import adapters

    st = _adata(["b1", "a1", "c1", "a2", "b2"], ["B", "A", "C", "A", "B"])
    reference = _adata(
        ["ar1", "ar2", "br1", "br2", "cr1", "missing"],
        ["A", "A", "B", "B", "C", "A"],
    )
    reference.obs["Level2"] = ["A1", "A2", "B1", "B2", "C1", None]
    from revise.application.preprocess import prepare_sc_svc_pair

    st, reference = prepare_sc_svc_pair(
        st,
        reference,
        broad_column="Level1",
        subtype_column="Level2",
    )
    assert pd.isna(reference.obs.loc["missing", "Level2"])
    calls = []

    def local_refinement(cell_type, *_args, **_kwargs):
        calls.append(cell_type)
        spatial_ids = st.obs_names[st.obs["Level1"] == cell_type]
        spatial = st[spatial_ids, :].copy()
        spatial.obs["Level2"] = [f"{cell_type}1", f"{cell_type}2"]
        spatial.obs["SVC_cluster"] = pd.Categorical(["0", "1"])
        expression = reference[reference.obs["Level1"] == cell_type, :].copy()
        expression = expression[expression.obs["Level2"].notna(), :].copy()
        expression.obs["SVC_cluster"] = pd.Categorical(
            ["0", "1"]
        )
        return spatial, expression

    applied = []
    ctx = SimpleNamespace(
        merged_config={
            "sc": {
                "select_ct": None,
                "resolutions": [0.5],
                "select_resolution": None,
            }
        },
        columns={"cell_type_col": "Level1", "sub_cell_type_col": "Level2"},
        runner=SimpleNamespace(
            st_adata=st,
            sc_ref_adata=reference,
            local_refinement=local_refinement,
        ),
        logger=logging.getLogger("test-all-sc-types"),
        artifacts={},
        record_local_refinement=applied.append,
    )

    adapters.ScSvcApplicationStrategy().solve_ot(ctx)

    assert calls == ["A", "B"]
    assert applied == [True]
    assert ctx.artifacts["processed_cell_types"] == ["A", "B"]
    assert ctx.artifacts["skipped_cell_types"] == [
        {
            "cell_type": "C",
            "reason": "insufficient_valid_reference_subtypes",
            "valid_subtype_count": 1,
        }
    ]
    spatial = ctx.artifacts["outputs"]["sc_svc_spatial"]
    expression = ctx.artifacts["outputs"]["sc_svc_expr"]
    assert spatial.obs_names.tolist() == ["a1", "a2", "b1", "b2"]
    assert expression.obs_names.tolist() == ["ar1", "ar2", "br1", "br2"]
    assert spatial.obs["SVC_cluster"].tolist() == [
        '["A","0"]', '["A","1"]', '["B","0"]', '["B","1"]'
    ]
    assert expression.obs["SVC_cluster"].tolist() == [
        '["A","0"]', '["A","1"]', '["B","0"]', '["B","1"]'
    ]
    assignments = ctx.artifacts["label_assignments"]
    assert assignments.loc[["a1", "a2", "b1", "b2"], "Level2"].tolist() == [
        "A1", "A2", "B1", "B2"
    ]
    assert pd.isna(assignments.loc["c1", "Level2"])


def test_cluster_namespace_is_injective_when_labels_contain_separator(monkeypatch):
    _import_sc_svc(monkeypatch)
    from revise.backend import adapters

    left = _adata(["left"], ["A::B"])
    left.obs["SVC_cluster"] = pd.Categorical(["C"])
    right = _adata(["right"], ["A"])
    right.obs["SVC_cluster"] = pd.Categorical(["B::C"])

    left_key = adapters._namespace_clusters(left, "A::B").obs["SVC_cluster"].iloc[0]
    right_key = adapters._namespace_clusters(right, "A").obs["SVC_cluster"].iloc[0]

    assert left_key == '["A::B","C"]'
    assert right_key == '["A","B::C"]'
    assert left_key != right_key


def test_ist_adapter_propagates_eligible_type_failure(monkeypatch):
    _import_sc_svc(monkeypatch)
    from revise.backend import adapters

    st = _adata(["a1", "a2"], ["A", "A"])
    reference = _adata(["r1", "r2"], ["A", "A"])
    reference.obs["Level2"] = ["A1", "A2"]
    ctx = SimpleNamespace(
        merged_config={"sc": {"select_ct": None, "resolutions": [0.5]}},
        columns={"cell_type_col": "Level1", "sub_cell_type_col": "Level2"},
        runner=SimpleNamespace(
            st_adata=st,
            sc_ref_adata=reference,
            local_refinement=lambda *_args, **_kwargs: (_ for _ in ()).throw(
                RuntimeError("LR exploded")
            ),
        ),
        logger=logging.getLogger("test-failing-sc-type"),
        artifacts={},
        record_local_refinement=lambda _value: None,
    )

    with pytest.raises(RuntimeError, match="LR exploded"):
        adapters.ScSvcApplicationStrategy().solve_ot(ctx)


def test_full_sample_fails_when_every_actual_type_is_ineligible(monkeypatch):
    _import_sc_svc(monkeypatch)
    from revise.backend import adapters

    st = _adata(["a1", "a2", "b1", "b2"], ["A", "A", "B", "B"])
    reference = _adata(["ar1", "ar2", "br1", "br2"], ["A", "A", "B", "B"])
    reference.obs["Level2"] = [None, "  ", "B1", "B1"]
    ctx = SimpleNamespace(
        merged_config={"sc": {"select_ct": None, "resolutions": [0.5]}},
        columns={"cell_type_col": "Level1", "sub_cell_type_col": "Level2"},
        runner=SimpleNamespace(
            st_adata=st,
            sc_ref_adata=reference,
            local_refinement=lambda *_args, **_kwargs: pytest.fail(
                "ineligible type reached LR"
            ),
        ),
        logger=logging.getLogger("test-all-ineligible-types"),
        artifacts={},
        record_local_refinement=lambda _value: None,
    )

    with pytest.raises(ValueError, match="no eligible broad cell types"):
        adapters.ScSvcApplicationStrategy().solve_ot(ctx)

    assert ctx.artifacts["processed_cell_types"] == []
    assert ctx.artifacts["skipped_cell_types"] == [
        {
            "cell_type": "A",
            "reason": "insufficient_valid_reference_subtypes",
            "valid_subtype_count": 0,
        },
        {
            "cell_type": "B",
            "reason": "insufficient_valid_reference_subtypes",
            "valid_subtype_count": 1,
        },
    ]


def test_full_sample_pipeline_runs_ga_once_before_all_local_types(monkeypatch, tmp_path):
    _import_sc_svc(monkeypatch)
    from revise.backend import adapters
    from revise.recon.context import PipelineContext
    from revise.recon.pipeline import UnifiedReconstructionPipeline

    st = _adata(["a1", "a2", "b1", "b2"], ["A", "A", "B", "B"])
    reference = _adata(["ar1", "ar2", "br1", "br2"], ["A", "A", "B", "B"])
    reference.obs["Level2"] = ["A1", "A2", "B1", "B2"]
    local_calls = []

    def local_refinement(cell_type, *_args, **_kwargs):
        local_calls.append(cell_type)
        spatial = st[st.obs["Level1"] == cell_type, :].copy()
        spatial.obs["Level2"] = [f"{cell_type}1", f"{cell_type}2"]
        spatial.obs["SVC_cluster"] = pd.Categorical(["0", "1"])
        expression = reference[reference.obs["Level1"] == cell_type, :].copy()
        expression.obs["SVC_cluster"] = pd.Categorical(["0", "1"])
        return spatial, expression

    runner = SimpleNamespace(
        st_adata=st,
        sc_ref_adata=reference,
        local_refinement=local_refinement,
    )

    class PipelineStrategy(adapters.ScSvcApplicationStrategy):
        def __init__(self):
            self.ga_calls = 0

        def prepare_context(self, ctx):
            ctx.runner = runner

        def global_anchoring(self, ctx):
            self.ga_calls += 1
            ctx.st_adata = runner.st_adata
            ctx.artifacts["ga_spatial"] = runner.st_adata

    strategy = PipelineStrategy()
    ctx = PipelineContext(
        merged_config={
            "io": {"save_outputs": False},
            "columns": {
                "cell_type_col": "Level1",
                "sub_cell_type_col": "Level2",
                "confidence_col": "Confidence",
            },
            "sc": {
                "select_ct": None,
                "resolutions": [0.5],
                "select_resolution": None,
            },
        },
        profile="test",
        runtime={
            "mode": "application",
            "task": "sc_svc",
            "svc_kind": "sc",
            "strategy": strategy.strategy_id,
            "application_route": "sc-SVC",
            "application_mode": "cluster",
            "compatibility_mode": False,
        },
        route_key="application:sc-SVC:cluster",
        run_dir=tmp_path,
        logger=logging.getLogger("test-full-sample-pipeline"),
    )
    validation = SimpleNamespace(validate=lambda _ctx: None)
    evaluation = SimpleNamespace(should_evaluate=lambda _ctx: False)

    UnifiedReconstructionPipeline(strategy, validation, evaluation).run(ctx)

    assert strategy.ga_calls == 1
    assert local_calls == ["A", "B"]
    assert ctx.svc.provenance["processed_cell_types"] == ["A", "B"]


@pytest.mark.parametrize("method", ["pot", "tacco"])
def test_application_ot_method_switches_global_and_local_together(method):
    from revise.application.config import _compile_engine_config

    request = SimpleNamespace(
        svc_type="sc-SVC",
        ot_method=method,
        local_refinement_strength=None,
        select_cell_type="T",
        broad_column="Level1",
        subtype_column="Level2",
    )

    overrides = _compile_engine_config(
        SimpleNamespace(
            svc_type=request.svc_type,
            ot_method=request.ot_method,
            broad_column=request.broad_column,
            subtype_column=request.subtype_column,
            select_cell_type=request.select_cell_type,
            local_refinement_strength=request.local_refinement_strength,
            local_refinement_alpha=0.2,
            local_refinement_resolutions=(0.6, 0.7, 0.8),
            local_refinement_graph_method=None,
            local_refinement_graph_alpha=None,
            local_refinement_graph_n_neighbors=None,
            local_refinement_graph_exp_neighbors=None,
            local_refinement_graph_spatial_neighbors=None,
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
