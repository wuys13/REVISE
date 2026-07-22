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
