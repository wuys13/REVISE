from __future__ import annotations

import logging
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise.backend import adapters


def _adata(obs_names, *, n_vars=2, obs=None):
    frame = pd.DataFrame(index=pd.Index(obs_names)) if obs is None else obs
    return AnnData(
        X=np.ones((len(obs_names), n_vars)),
        obs=frame,
        var=pd.DataFrame(index=[f"g{i}" for i in range(n_vars)]),
    )


@pytest.mark.parametrize("route", ["segmentation", "bin2cell"])
def test_segmentation_preprocess_aligns_sampled_st_and_gt(route, monkeypatch):
    st = _adata(["c1", "c2", "c3"])
    gt = _adata(["c3", "c1", "other"])
    reference = _adata(["r1", "r2"])
    monkeypatch.setattr(
        adapters,
        "_subsample_obs",
        lambda adata, size, seed: (adata[["c2", "c3", "c1"], :].copy(), ["c2", "c3", "c1"]),
    )

    prepared_st, prepared_ref, prepared_gt = adapters.preprocess_data(
        route,
        st,
        reference,
        gt,
        SimpleNamespace(),
        sample_size=3,
        sample_seed=17,
        logger=logging.getLogger("benchmark-preprocess-seg"),
    )

    assert prepared_st.obs_names.tolist() == ["c3", "c1"]
    assert prepared_gt.obs_names.tolist() == ["c3", "c1"]
    assert prepared_ref is reference


@pytest.mark.parametrize("route", ["batch_effect", "spot_size"])
def test_sr_preprocess_ensures_mapping_before_sampling_and_filters_gt_pm(
    route,
    monkeypatch,
):
    calls = []
    st = _adata(["spot-1", "spot-2"])
    st.uns["all_cells_in_spot"] = {
        "spot-1": ["c1", "c2"],
        "spot-2": ["c3"],
    }
    gt_obs = pd.DataFrame(
        {"cell_id": ["c1", "c2", "c3"]},
        index=["row-1", "row-2", "row-3"],
    )
    gt = _adata(gt_obs.index, obs=gt_obs)
    reference = _adata(["r1", "r2"])
    conf = SimpleNamespace(
        pm_on_cell=pd.DataFrame(
            np.ones((3, 1)),
            index=["c1", "c2", "c3"],
            columns=["A"],
        )
    )
    monkeypatch.setattr(
        adapters,
        "ensure_all_cells_in_spot",
        lambda *_args, **_kwargs: calls.append("ensure"),
    )

    def sample(adata, size, seed):
        calls.append("sample")
        return adata[["spot-1"], :].copy(), ["spot-1"]

    monkeypatch.setattr(adapters, "_subsample_obs", sample)

    prepared_st, prepared_ref, prepared_gt = adapters.preprocess_data(
        route,
        st,
        reference,
        gt,
        conf,
        sample_size=1,
        sample_seed=17,
        logger=logging.getLogger("benchmark-preprocess-sr"),
    )

    assert calls == ["ensure", "sample"]
    assert prepared_st.obs_names.tolist() == ["spot-1"]
    assert prepared_gt.obs["cell_id"].tolist() == ["c1", "c2"]
    assert conf.pm_on_cell.index.tolist() == ["c1", "c2"]
    assert prepared_ref is reference


@pytest.mark.parametrize("route", ["gene_panel", "gene_dropout"])
def test_impute_preprocess_aligns_inputs_filters_singletons_and_caps_pcs(
    route,
    monkeypatch,
):
    st = _adata(["c1", "c2", "c3"], n_vars=4)
    gt = _adata(["c1", "c2", "c3"], n_vars=4)
    ref_obs = pd.DataFrame(
        {"Level1": ["A", "A", "B", "C"]},
        index=["r1", "r2", "r3", "r4"],
    )
    reference = _adata(ref_obs.index, n_vars=4, obs=ref_obs)
    conf = SimpleNamespace(cell_type_col="Level1", rec_graph_n_pcs=50)
    calls = []

    def sample(adata, size, seed):
        calls.append(adata.obs_names.tolist())
        keep = adata.obs_names[:size]
        return adata[keep, :].copy(), keep.to_numpy()

    monkeypatch.setattr(adapters, "_subsample_obs", sample)

    prepared_st, prepared_ref, prepared_gt = adapters.preprocess_data(
        route,
        st,
        reference,
        gt,
        conf,
        sample_size=3,
        sample_seed=17,
        logger=logging.getLogger("benchmark-preprocess-impute"),
    )

    assert calls == [["c1", "c2", "c3"], ["r1", "r2", "r3", "r4"]]
    assert prepared_st.obs_names.tolist() == prepared_gt.obs_names.tolist() == ["c1", "c2", "c3"]
    assert prepared_ref.obs["Level1"].tolist() == ["A", "A"]
    assert conf.rec_graph_n_pcs == 1


def test_benchmark_preprocess_rejects_unknown_route():
    data = _adata(["c1"])
    with pytest.raises(ValueError, match="Unsupported benchmark route"):
        adapters.preprocess_data(
            "noise",
            data,
            data,
            data,
            SimpleNamespace(),
            sample_size=None,
            sample_seed=42,
            logger=logging.getLogger("benchmark-preprocess-unknown"),
        )
