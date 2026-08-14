from __future__ import annotations

import builtins
import sys
from types import ModuleType, SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from revise.analysis.advanced.trajectory import infer_palantir


def _adata() -> ad.AnnData:
    return ad.AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame({"Level1": ["SMC", "Fibroblast"]}, index=["c1", "c2"]),
    )


def test_infer_palantir_runs_the_current_omicverse_sequence_without_x_umap(
    monkeypatch,
):
    calls = []

    class FakeTrajectory:
        def __init__(self, adata):
            self.adata = adata

        def set_origin_cells(self, origin):
            calls.append(("origin", origin))

        def set_terminal_cells(self, terminal):
            calls.append(("terminal", terminal))

        def inference(self, **kwargs):
            calls.append(("inference", kwargs))
            self.adata.obs["palantir_pseudotime"] = [0.0, 1.0]

    def fake_traj_infer(adata, **kwargs):
        calls.append(("TrajInfer", kwargs))
        return FakeTrajectory(adata)

    omicverse = ModuleType("omicverse")
    omicverse.single = SimpleNamespace(TrajInfer=fake_traj_infer)
    monkeypatch.setitem(sys.modules, "omicverse", omicverse)
    adata = _adata()
    assert "X_umap" not in adata.obsm

    pseudotime = infer_palantir(
        adata,
        groupby="Level1",
        basis="X_umap",
        use_rep="scaled|original|X_pca",
        n_comps=50,
        origin_cells="SMC",
        terminal_cells=["Fibroblast"],
        num_waypoints=500,
    )

    assert calls == [
        (
            "TrajInfer",
            {
                "basis": "X_umap",
                "groupby": "Level1",
                "use_rep": "scaled|original|X_pca",
                "n_comps": 50,
            },
        ),
        ("origin", "SMC"),
        ("terminal", ["Fibroblast"]),
        ("inference", {"method": "palantir", "num_waypoints": 500}),
    ]
    assert pseudotime.index.tolist() == ["c1", "c2"]
    assert pseudotime.tolist() == [0.0, 1.0]


@pytest.mark.parametrize(
    ("groupby", "origin_cells", "terminal_cells", "message"),
    [
        ("missing", "SMC", ["Fibroblast"], "groupby"),
        ("Level1", "Missing", ["Fibroblast"], "origin_cells"),
        ("Level1", "SMC", ["Missing"], "terminal_cells"),
    ],
)
def test_infer_palantir_validates_group_labels_before_importing_omicverse(
    groupby, origin_cells, terminal_cells, message
):
    with pytest.raises(ValueError, match=message):
        infer_palantir(
            _adata(),
            groupby=groupby,
            basis="X_umap",
            use_rep="scaled|original|X_pca",
            n_comps=50,
            origin_cells=origin_cells,
            terminal_cells=terminal_cells,
            num_waypoints=500,
        )


def test_infer_palantir_missing_dependency_names_trajectory_extra(
    monkeypatch,
):
    real_import = builtins.__import__

    def blocked_import(name, *args, **kwargs):
        if name == "omicverse":
            raise ModuleNotFoundError("No module named 'omicverse'", name=name)
        return real_import(name, *args, **kwargs)

    monkeypatch.delitem(sys.modules, "omicverse", raising=False)
    monkeypatch.setattr(builtins, "__import__", blocked_import)

    with pytest.raises(ImportError, match=r"revise-svc\[trajectory\]"):
        infer_palantir(
            _adata(),
            groupby="Level1",
            basis="X_umap",
            use_rep="scaled|original|X_pca",
            n_comps=50,
            origin_cells="SMC",
            terminal_cells=["Fibroblast"],
            num_waypoints=500,
        )


def test_infer_palantir_rejects_a_missing_pseudotime_result(monkeypatch):
    class FakeTrajectory:
        def set_origin_cells(self, origin):
            pass

        def set_terminal_cells(self, terminal):
            pass

        def inference(self, **kwargs):
            pass

    omicverse = ModuleType("omicverse")
    omicverse.single = SimpleNamespace(TrajInfer=lambda adata, **kwargs: FakeTrajectory())
    monkeypatch.setitem(sys.modules, "omicverse", omicverse)

    with pytest.raises(RuntimeError, match="palantir_pseudotime"):
        infer_palantir(
            _adata(),
            groupby="Level1",
            basis="X_umap",
            use_rep="scaled|original|X_pca",
            n_comps=50,
            origin_cells="SMC",
            terminal_cells=["Fibroblast"],
            num_waypoints=500,
        )
