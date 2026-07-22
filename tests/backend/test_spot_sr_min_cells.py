from __future__ import annotations

import importlib
import logging
from pathlib import Path
import sys
import types
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise.config import load_raw_config, merge_unified_config


CONFIG_PATH = Path(__file__).parents[2] / "revise" / "revise.yaml"
ISOLATED_MODULE_NAMES = (
    "scanpy",
    "revise.backend.adapters",
    "revise.backend.kernels",
    "revise.backend.kernels.global_anchoring",
    "revise.backend.kernels.graph_aggregate",
    "revise.backend.kernels.spot_sr",
    "revise.backend.ops.distance",
    "revise.backend.ops.meta",
    "revise.backend.ops.shaver",
    "revise.backend.ops.topology",
    "revise.backend.runners.application_svc",
    "revise.backend.runners.base_svc_anchor",
    "revise.backend.runners.sc_svc_sr_application",
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


def _install_scanpy_preprocessing_stub():
    def filter_cells(adata, *, min_counts=None, **_kwargs):
        if min_counts is None:
            return
        keep = np.asarray(adata.X.sum(axis=1)).ravel() >= min_counts
        adata._inplace_subset_obs(keep)

    def filter_genes(adata, *, min_cells=None, min_counts=None, **_kwargs):
        if min_cells is not None:
            keep = np.asarray((adata.X != 0).sum(axis=0)).ravel() >= min_cells
        elif min_counts is not None:
            keep = np.asarray(adata.X.sum(axis=0)).ravel() >= min_counts
        else:
            return
        adata._inplace_subset_var(keep)

    scanpy = types.ModuleType("scanpy")
    scanpy.pp = SimpleNamespace(
        filter_cells=filter_cells,
        filter_genes=filter_genes,
    )
    scanpy.pl = SimpleNamespace()
    scanpy.tl = SimpleNamespace()
    scanpy.AnnData = AnnData
    scanpy._revise_test_stub = True
    sys.modules["scanpy"] = scanpy

    distance = types.ModuleType("revise.backend.ops.distance")
    distance.bhattacharyya_distance = lambda *args, **kwargs: None
    distance.similarity_to_distance = lambda *args, **kwargs: None
    sys.modules["revise.backend.ops.distance"] = distance


@pytest.fixture
def adapters():
    snapshot = _snapshot_modules(ISOLATED_MODULE_NAMES)
    for module_name in ISOLATED_MODULE_NAMES:
        sys.modules.pop(module_name, None)
    _install_scanpy_preprocessing_stub()
    try:
        yield importlib.import_module("revise.backend.adapters")
    finally:
        _restore_modules(snapshot)


@pytest.fixture
def sst_adatas():
    expression = np.array(
        [
            [1.0, 2.0],
            [1.0, 0.0],
        ],
        dtype=np.float64,
    )
    genes = ["widespread", "single_high"]
    annotations = {
        "Patient": ["sample", "sample"],
        "Level1": ["A", "A"],
        "Level2": ["A1", "A2"],
    }
    st_adata = AnnData(
        X=expression.copy(),
        obs=pd.DataFrame(annotations, index=["spot-1", "spot-2"]),
        var=pd.DataFrame(index=genes),
    )
    st_adata.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 1.0]])
    st_adata.uns["all_cells_in_spot"] = {
        "spot-1": ["cell-1"],
        "spot-2": ["cell-2"],
    }
    sc_ref_adata = AnnData(
        X=expression.copy(),
        obs=pd.DataFrame(annotations, index=["cell-1", "cell-2"]),
        var=pd.DataFrame(index=genes),
    )
    return st_adata, sc_ref_adata


def test_application_spot_sr_prepare_context_filters_genes_by_observation_count(
    adapters, monkeypatch, tmp_path, sst_adatas
):
    merged = merge_unified_config(
        raw_config=load_raw_config(CONFIG_PATH),
        profile="application_sc_sr",
        runtime_overrides={},
        io_overrides={},
        set_overrides=(),
    )
    merged["io"].update(
        sample_name="sample",
        data_root=str(tmp_path),
        output_root=str(tmp_path),
        st_file="st.h5ad",
        sc_ref_file="sc.h5ad",
        patient_key="Patient",
    )
    merged["preprocess"].update(
        st_min_transcripts=0,
        st_min_cells=2,
        sc_min_cells=2,
    )
    st_adata, sc_ref_adata = sst_adatas

    class InputService:
        def read_st_adata(self, _path):
            return st_adata.copy()

        def read_sc_ref_adata(self, _path):
            return sc_ref_adata.copy()

    monkeypatch.setattr(adapters, "_input_service", lambda _ctx: InputService())
    ctx = SimpleNamespace(
        merged_config=merged,
        io=merged["io"],
        columns=merged["columns"],
        runtime=merged["runtime"],
        route_key="sc_svc_sr:spot_size",
        run_dir=tmp_path,
        logger=logging.getLogger("test-spot-sr-min-cells"),
        compatibility_mode=False,
    )

    adapters.ScSvcSrApplicationStrategy().prepare_context(ctx)

    assert type(ctx.runner).__name__ == "ScSVCSr"
    assert ctx.st_adata.var_names.tolist() == ["widespread"]
    assert ctx.sc_ref_adata.var_names.tolist() == ["widespread"]


def test_application_spot_sr_prepare_context_honors_configured_annotation_columns(
    adapters, monkeypatch, tmp_path
):
    merged = merge_unified_config(
        raw_config=load_raw_config(CONFIG_PATH),
        profile="application_sc_sr",
        runtime_overrides={},
        io_overrides={},
        set_overrides=(
            "columns.cell_type_col=major_type",
            "columns.sub_cell_type_col=minor_type",
        ),
    )
    merged["io"].update(
        sample_name="sample",
        data_root=str(tmp_path),
        output_root=str(tmp_path),
        st_file="st.h5ad",
        sc_ref_file="sc.h5ad",
        patient_key="Patient",
    )
    merged["preprocess"].update(
        st_min_transcripts=0,
        st_min_cells=1,
        sc_min_cells=1,
    )
    st = AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame(index=["spot-1", "spot-2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    st.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 1.0]])
    st.uns["all_cells_in_spot"] = {
        "spot-1": ["cell-1"],
        "spot-2": ["cell-2"],
    }
    reference = AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame(
            {
                "Patient": ["sample", "sample"],
                "major_type": ["A/B", "A/B"],
                "minor_type": ["A1", "A2"],
            },
            index=["cell-1", "cell-2"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )

    class InputService:
        def read_st_adata(self, _path):
            return st.copy()

        def read_sc_ref_adata(self, _path):
            return reference.copy()

    monkeypatch.setattr(adapters, "_input_service", lambda _ctx: InputService())
    ctx = SimpleNamespace(
        merged_config=merged,
        io=merged["io"],
        columns=merged["columns"],
        runtime=merged["runtime"],
        route_key="sc_svc_sr:spot_size",
        run_dir=tmp_path,
        logger=logging.getLogger("test-spot-sr-custom-columns"),
        compatibility_mode=False,
    )

    adapters.ScSvcSrApplicationStrategy().prepare_context(ctx)

    assert ctx.runner.config.cell_type_col == "major_type"
    assert ctx.runner.sc_ref_adata.obs["major_type"].tolist() == ["A_B", "A_B"]


def test_application_spot_sr_validates_overlap_after_gene_filtering(
    adapters, monkeypatch, tmp_path
):
    merged = merge_unified_config(
        raw_config=load_raw_config(CONFIG_PATH),
        profile="application_sc_sr",
        runtime_overrides={},
        io_overrides={},
        set_overrides=(),
    )
    merged["io"].update(
        sample_name="sample",
        data_root=str(tmp_path),
        output_root=str(tmp_path),
        st_file="st.h5ad",
        sc_ref_file="sc.h5ad",
        patient_key="Patient",
    )
    merged["preprocess"].update(
        st_min_transcripts=0,
        st_min_cells=2,
        sc_min_cells=2,
    )
    annotations = {
        "Patient": ["sample", "sample"],
        "Level1": ["A", "A"],
        "Level2": ["A1", "A2"],
    }
    st_adata = AnnData(
        X=np.array([[1.0, 1.0], [0.0, 1.0]]),
        obs=pd.DataFrame(annotations, index=["spot-1", "spot-2"]),
        var=pd.DataFrame(index=["shared", "st-only"]),
    )
    st_adata.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 1.0]])
    st_adata.uns["all_cells_in_spot"] = {
        "spot-1": ["cell-1"],
        "spot-2": ["cell-2"],
    }
    sc_ref_adata = AnnData(
        X=np.array([[1.0, 1.0], [0.0, 1.0]]),
        obs=pd.DataFrame(annotations, index=["cell-1", "cell-2"]),
        var=pd.DataFrame(index=["shared", "sc-only"]),
    )

    class InputService:
        def read_st_adata(self, _path):
            return st_adata.copy()

        def read_sc_ref_adata(self, _path):
            return sc_ref_adata.copy()

    monkeypatch.setattr(adapters, "_input_service", lambda _ctx: InputService())
    ctx = SimpleNamespace(
        merged_config=merged,
        io=merged["io"],
        columns=merged["columns"],
        runtime=merged["runtime"],
        route_key="sc_svc_sr:spot_size",
        run_dir=tmp_path,
        logger=logging.getLogger("test-spot-sr-post-filter-overlap"),
        compatibility_mode=False,
    )

    with pytest.raises(ValueError) as exc_info:
        adapters.ScSvcSrApplicationStrategy().prepare_context(ctx)

    message = str(exc_info.value)
    expected = (
        f"st_file_path={tmp_path / 'sample_st.h5ad'}",
        f"sc_ref_file_path={tmp_path / 'sc.h5ad'}",
        "field=var_names_overlap",
        "expected=>=1",
        "actual=0",
    )
    missing = [fragment for fragment in expected if fragment not in message]
    if missing:
        raise AssertionError(f"missing message fragments {missing}: {message}")
