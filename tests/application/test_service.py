"""Application result assembly contract.

Covers: sc-SVC pair publication and canonical publication metadata.
Proof limit: uses synthetic AnnData and does not run scientific kernels.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, read_h5ad

from revise.svc import SVC


class PublicationContext(SimpleNamespace):
    def set_pending_publication(self, *, commit, rollback):
        self._publication_commit = commit
        self._publication_rollback = rollback

    def rollback_pending_publication(self):
        rollback = getattr(self, "_publication_rollback", None)
        self._publication_commit = None
        self._publication_rollback = None
        if rollback is not None:
            rollback()


def _sc_outputs():
    spatial = AnnData(
        X=np.zeros((3, 1)),
        obs=pd.DataFrame(
            {"SVC_cluster": ["A", "B", "A"], "position_owner": [1, 2, 3]},
            index=["position-1", "position-2", "position-3"],
        ),
        var=pd.DataFrame(index=["spatial-placeholder"]),
    )
    spatial.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]])
    expression = AnnData(
        X=np.array([[2.0, 4.0], [6.0, 8.0], [10.0, 12.0]]),
        obs=pd.DataFrame(
            {"SVC_cluster": ["A", "A", "B"]},
            index=["cell-1", "cell-2", "cell-3"],
        ),
        var=pd.DataFrame({"gene_owner": ["x", "y"]}, index=["g1", "g2"]),
    )
    return spatial, expression


def test_sc_svc_publication_preserves_the_notebook_spatial_and_expression_pair(
    tmp_path,
):
    from revise.application.service import _build_sc_public_results

    spatial, expression = _sc_outputs()
    records = []
    ctx = PublicationContext(
        svc=SVC(
            expr=expression,
            spatial=spatial,
            svc_kind="sc",
            artifacts={
                "outputs": {
                    "sc_svc_spatial": spatial,
                    "sc_svc_expr": expression,
                }
            },
        ),
        run_dir=tmp_path / "run",
        merged_config={"ot": {"ga": {"solver": "pot"}, "lr": {"solver": "pot"}}},
        provenance={},
        artifact_records=records,
        record_artifact=records.append,
    )
    args = SimpleNamespace(
        svc_type="sc-SVC",
        output_root=str(tmp_path / "out"),
        sample_name="sample",
        seed=17,
        ot_method="pot",
        select_ct="T",
    )

    results, paths = _build_sc_public_results(args, "application_sc", ctx)

    output_dir = tmp_path / "out" / "sample" / "sc-SVC" / "T"
    assert paths == {
        "spatial": output_dir / "sc_SVC_spatial.h5ad",
        "expression": output_dir / "sc_SVC_expr.h5ad",
    }
    assert set(results) == {"spatial", "expression"}
    assert not (output_dir / "SVC.h5ad").exists()
    assert set(ctx.provenance["results"]) == {"spatial", "expression"}

    published_spatial = read_h5ad(paths["spatial"])
    published_expression = read_h5ad(paths["expression"])
    assert published_spatial.obs_names.tolist() == spatial.obs_names.tolist()
    assert published_spatial.var_names.tolist() == spatial.var_names.tolist()
    assert np.array_equal(published_spatial.X, spatial.X)
    assert np.array_equal(published_spatial.obsm["spatial"], spatial.obsm["spatial"])
    assert published_expression.obs_names.tolist() == expression.obs_names.tolist()
    assert published_expression.var_names.tolist() == expression.var_names.tolist()
    assert np.array_equal(published_expression.X, expression.X)
    assert published_spatial.uns["revise_reconstruction"]["output_role"] == "spatial"
    assert (
        published_expression.uns["revise_reconstruction"]["output_role"]
        == "expression"
    )
    assert "ot_events" not in published_spatial.uns["revise_reconstruction"]


def test_single_file_publication_reloads_staged_h5ad_before_replace(
    monkeypatch,
    tmp_path,
):
    from revise.application import service

    output = AnnData(
        X=np.ones((1, 1)),
        obs=pd.DataFrame(index=["cell-1"]),
        var=pd.DataFrame(index=["g1"]),
    )
    output_path = tmp_path / "out" / "sample" / "SVC.h5ad"
    output_path.parent.mkdir(parents=True)
    output_path.write_bytes(b"previous-valid-result")
    ctx = PublicationContext(
        svc=SVC(
            expr=output,
            spatial=output,
            svc_kind="sp",
            artifacts={"outputs": {"sp_svc": output}},
        ),
        run_dir=tmp_path / "run",
        merged_config={"ot": {"ga": {"solver": "pot"}, "lr": {"solver": "pot"}}},
        provenance={},
        record_artifact=lambda artifact: None,
    )
    args = SimpleNamespace(
        svc_type="sp-SVC",
        output_root=str(tmp_path / "out"),
        sample_name="sample",
        seed=17,
        ot_method="pot",
    )

    def corrupt_write(self, path, *args, **kwargs):
        Path(path).write_bytes(b"not-an-h5ad")

    monkeypatch.setattr(AnnData, "write_h5ad", corrupt_write)

    with pytest.raises(OSError):
        service._build_public_result(args, "application_sp", "sp_svc", ctx)

    assert output_path.read_bytes() == b"previous-valid-result"
    assert list(output_path.parent.iterdir()) == [output_path]


def test_sc_pair_publication_rolls_back_both_existing_files_on_replace_failure(
    monkeypatch,
    tmp_path,
):
    from revise.application import service

    spatial, expression = _sc_outputs()
    output_dir = tmp_path / "out" / "sample" / "sc-SVC" / "T"
    output_dir.mkdir(parents=True)
    spatial_path = output_dir / "sc_SVC_spatial.h5ad"
    expression_path = output_dir / "sc_SVC_expr.h5ad"
    legacy_path = output_dir / "SVC.h5ad"
    spatial_path.write_bytes(b"old-spatial")
    expression_path.write_bytes(b"old-expression")
    legacy_path.write_bytes(b"old-legacy")
    records = []
    ctx = PublicationContext(
        svc=SVC(
            expr=expression,
            spatial=spatial,
            svc_kind="sc",
            artifacts={
                "outputs": {
                    "sc_svc_spatial": spatial,
                    "sc_svc_expr": expression,
                }
            },
        ),
        run_dir=tmp_path / "run",
        runtime={"seed": 17},
        merged_config={
            "runtime": {"seed": 17},
            "sc": {"tacco_annotate": {"multi_center": 1, "lamb": 0.001}},
            "ot": {"ga": {"solver": "tacco"}, "lr": {"solver": "tacco"}},
        },
        provenance={},
        artifact_records=records,
        record_artifact=records.append,
    )
    args = SimpleNamespace(
        svc_type="sc-SVC",
        output_root=str(tmp_path / "out"),
        sample_name="sample",
        seed=17,
        ot_method=None,
        select_ct="T",
    )
    original_replace = service.os.replace
    failed = False

    def replace(source, destination):
        nonlocal failed
        destination = Path(destination)
        if destination == expression_path and not failed:
            failed = True
            raise OSError("injected second-result replace failure")
        return original_replace(source, destination)

    monkeypatch.setattr(service.os, "replace", replace)

    with pytest.raises(OSError, match="second-result"):
        service._build_sc_public_results(args, "application_sc", ctx)

    assert spatial_path.read_bytes() == b"old-spatial"
    assert expression_path.read_bytes() == b"old-expression"
    assert legacy_path.read_bytes() == b"old-legacy"
    assert records == []
    assert "results" not in ctx.provenance


def test_sc_pair_delayed_rollback_restores_pair_without_touching_unrelated_file(
    tmp_path,
):
    from revise.application import service

    spatial, expression = _sc_outputs()
    output_dir = tmp_path / "out" / "sample" / "sc-SVC" / "T"
    output_dir.mkdir(parents=True)
    spatial_path = output_dir / "sc_SVC_spatial.h5ad"
    expression_path = output_dir / "sc_SVC_expr.h5ad"
    legacy_path = output_dir / "SVC.h5ad"
    spatial_path.write_bytes(b"old-spatial")
    expression_path.write_bytes(b"old-expression")
    legacy_path.write_bytes(b"old-legacy")
    records = []
    ctx = PublicationContext(
        svc=SVC(
            expr=expression,
            spatial=spatial,
            svc_kind="sc",
            artifacts={
                "outputs": {
                    "sc_svc_spatial": spatial,
                    "sc_svc_expr": expression,
                }
            },
        ),
        run_dir=tmp_path / "run",
        runtime={"seed": 17},
        merged_config={
            "runtime": {"seed": 17},
            "ot": {"ga": {"solver": "pot"}, "lr": {"solver": "pot"}},
        },
        provenance={},
        artifact_records=records,
        record_artifact=records.append,
    )
    args = SimpleNamespace(
        svc_type="sc-SVC",
        output_root=str(tmp_path / "out"),
        sample_name="sample",
        seed=17,
        ot_method=None,
        select_ct="T",
    )

    service._build_sc_public_results(args, "application_sc", ctx)
    assert legacy_path.read_bytes() == b"old-legacy"
    assert set(ctx.provenance["results"]) == {"spatial", "expression"}

    ctx.rollback_pending_publication()

    assert spatial_path.read_bytes() == b"old-spatial"
    assert expression_path.read_bytes() == b"old-expression"
    assert legacy_path.read_bytes() == b"old-legacy"
    assert records == []
    assert "results" not in ctx.provenance
