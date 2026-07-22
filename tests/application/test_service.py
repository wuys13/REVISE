"""Application result assembly contract.

Covers: sc-SVC spatial/expression merge and canonical publication metadata.
Proof limit: uses synthetic AnnData and does not run scientific kernels.
"""

from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
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


def test_merge_sc_svc_uses_spatial_rows_and_expression_genes():
    from revise.application.service import merge_sc_svc

    spatial, expression = _sc_outputs()

    merged = merge_sc_svc(spatial, expression, mode="mean", seed=17)

    assert merged.obs_names.tolist() == spatial.obs_names.tolist()
    assert merged.var_names.tolist() == expression.var_names.tolist()
    assert merged.obs["position_owner"].tolist() == [1, 2, 3]
    assert merged.var["gene_owner"].tolist() == ["x", "y"]
    assert np.array_equal(merged.obsm["spatial"], spatial.obsm["spatial"])
    assert np.array_equal(merged.X, np.array([[4.0, 6.0], [10.0, 12.0], [4.0, 6.0]]))
    assert merged.uns["revise_reconstruction"]["mapping_mode"] == "mean"
    assert "platform" not in merged.uns["revise_reconstruction"]


def test_sc_svc_publication_records_only_the_confirmed_result_type(tmp_path):
    from revise.application.service import _build_public_result

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
        ot_events=[],
        provenance={},
        record_artifact=records.append,
    )
    args = SimpleNamespace(
        svc_type="sc-SVC",
        output_root=str(tmp_path / "out"),
        sample_name="sample",
        seed=17,
        ot_method="pot",
        sc_mapping="mean",
    )

    _, path = _build_public_result(args, "application_sc", None, ctx)

    assert path == tmp_path / "out" / "sample" / "SVC.h5ad"
    assert ctx.provenance["result"] == {
        "filename": "SVC.h5ad",
        "type": "sc-SVC",
    }
    published = read_h5ad(path)
    metadata = published.uns["revise_reconstruction"]
    assert metadata["svc_type"] == "sc-SVC"
    assert "platform" not in metadata
    assert json.loads(metadata["ot_events"]) == []


def test_sc_svc_publication_uses_the_resolved_yaml_seed_when_cli_omits_it(tmp_path):
    from revise.application.service import _build_public_result

    spatial, expression = _sc_outputs()
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
        ot_events=[],
        provenance={},
        record_artifact=lambda artifact: None,
    )
    args = SimpleNamespace(
        svc_type="sc-SVC",
        output_root=str(tmp_path / "out"),
        sample_name="sample",
        seed=None,
        ot_method=None,
        sc_mapping="mean",
    )

    result, _ = _build_public_result(args, "application_sc", None, ctx)

    assert result.uns["revise_reconstruction"]["mapping_seed"] == 17
