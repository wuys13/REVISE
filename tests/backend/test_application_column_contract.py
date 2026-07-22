from __future__ import annotations

import logging
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse


class _ReachedConfiguredColumn(RuntimeError):
    pass


def test_sp_adapter_normalizes_configured_reference_labels(monkeypatch, tmp_path):
    from revise.backend import adapters

    st = AnnData(
        X=sparse.csr_matrix(np.ones((2, 2))),
        obs=pd.DataFrame(index=["s1", "s2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    reference = AnnData(
        X=sparse.csr_matrix(np.ones((2, 2))),
        obs=pd.DataFrame(
            {
                "Patient": ["sample", "sample"],
                "major_type": ["T/NK", "T/NK"],
                "minor_type": ["T/1", "T/2"],
            },
            index=["c1", "c2"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )

    class InputService:
        def read_st_adata(self, path):
            return st.copy()

        def read_sc_ref_adata(self, path):
            return reference.copy()

    monkeypatch.setattr(adapters, "_input_service", lambda ctx: InputService())
    monkeypatch.setattr(adapters.sc.pp, "filter_cells", lambda *args, **kwargs: None)
    monkeypatch.setattr(adapters.sc.pp, "filter_genes", lambda *args, **kwargs: None)
    ctx = SimpleNamespace(
        merged_config={
            "ot": {
                "ga": {
                    "solver": "pot",
                    "pot": {"reg": 0.1, "reg_m": 0.0, "reg_type": "entropy"},
                },
                "lr": {
                    "solver": "pot",
                    "pot": {"reg": 0.1, "reg_m": 0.0, "reg_type": "kl"},
                },
                "impute": {"reg": 5.0, "reg_m": 0.0, "reg_type": "kl"},
            },
            "preprocess": {},
            "graph": {},
            "plot": {},
            "posterior_conditioning": {},
        },
        io={
            "sample_name": "sample",
            "data_root": str(tmp_path),
            "st_file": "st.h5ad",
            "sc_ref_file": "sc.h5ad",
            "patient_key": "Patient",
        },
        columns={
            "cell_type_col": "major_type",
            "sub_cell_type_col": "minor_type",
            "confidence_col": "Confidence",
            "unknown_key": "Unknown",
        },
        run_dir=tmp_path,
        logger=logging.getLogger("test-sp-label-normalization"),
    )

    adapters.SpSvcApplicationStrategy().prepare_context(ctx)

    assert ctx.runner.sc_ref_adata.obs["major_type"].tolist() == ["T_NK", "T_NK"]
    assert ctx.runner.sc_ref_adata.obs["minor_type"].tolist() == ["T_1", "T_2"]


def test_reference_label_normalization_rejects_category_collisions():
    from revise.backend.adapters import _replace_slash_labels

    reference = AnnData(
        X=np.ones((2, 1)),
        obs=pd.DataFrame({"major_type": ["A/B", "A_B"]}, index=["c1", "c2"]),
        var=pd.DataFrame(index=["g1"]),
    )

    with pytest.raises(ValueError, match="collide after slash normalization"):
        _replace_slash_labels(reference, ["major_type"])


def test_sp_local_refinement_trims_by_the_configured_cell_type(monkeypatch):
    from revise.backend.runners import sp_svc_application as module

    st = AnnData(
        X=sparse.csr_matrix(np.ones((2, 2))),
        obs=pd.DataFrame({"major_type": ["A", "A"]}, index=["s1", "s2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    reference = AnnData(
        X=sparse.csr_matrix(np.ones((2, 2))),
        obs=pd.DataFrame({"major_type": ["A", "A"]}, index=["c1", "c2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    runner = module.SpSVC.__new__(module.SpSVC)
    runner.st_adata = st
    runner.sc_ref_adata = reference
    runner.config = SimpleNamespace(plot_flag=False, cell_type_col="major_type")
    runner.logger = logging.getLogger("test-sp-custom-column")

    def capture_column(_st, _reference, column):
        assert column == "major_type"
        raise _ReachedConfiguredColumn

    monkeypatch.setattr(module, "trim_sp_adata", capture_column)

    with pytest.raises(_ReachedConfiguredColumn):
        runner.local_refinement()


def test_sr_reference_profiles_ignore_an_unrelated_clusters_column(monkeypatch):
    from revise.backend.runners import sc_svc_sr_application as module

    st = AnnData(
        X=np.ones((1, 2)),
        obs=pd.DataFrame(index=["spot-1"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    st.obsm["major_type"] = pd.DataFrame(
        [[1.0]], index=st.obs_names, columns=["A"]
    )
    reference = AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame(
            {"major_type": ["A", "A"], "clusters": ["0", "1"]},
            index=["c1", "c2"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    runner = module.ScSVCSr.__new__(module.ScSVCSr)
    runner.st_adata = st
    runner.sc_ref_adata = reference
    runner.config = SimpleNamespace(cell_type_col="major_type")
    runner.spot_sr = SimpleNamespace(run=lambda _runner: None)
    runner.logger = logging.getLogger("test-sr-custom-column")

    def capture_column(_adata, *, key_type, type_list):
        assert key_type == "major_type"
        assert type_list == ["A"]
        raise _ReachedConfiguredColumn

    monkeypatch.setattr(module.sc.pp, "normalize_total", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(module, "construct_sc_ref", capture_column)

    with pytest.raises(_ReachedConfiguredColumn):
        runner.local_refinement()
