from types import SimpleNamespace

import numpy as np
import pandas as pd
from anndata import AnnData, read_h5ad

from revise.svc import SVC


class PublicationContext(SimpleNamespace):
    def set_pending_publication(self, *, commit, rollback):
        self._commit = commit
        self._rollback = rollback

    def rollback_pending_publication(self):
        self._rollback()

    def record_artifact(self, artifact):
        self.artifact_records.append(artifact)


def _ot_config():
    return {"ot": {"ga": {"solver": "pot"}, "lr": {"solver": "pot"}}}


def _sc_pair():
    spatial = AnnData(
        X=np.ones((2, 1)),
        obs=pd.DataFrame({"SVC_cluster": ["A", "B"]}, index=["p1", "p2"]),
        var=pd.DataFrame(index=["spatial-placeholder"]),
    )
    spatial.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 1.0]])
    expression = AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame({"SVC_cluster": ["A", "B"]}, index=["c1", "c2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    return spatial, expression


def test_publish_single_writes_one_stable_h5ad(tmp_path):
    from revise.application.publication import publish_single

    output = AnnData(
        X=np.ones((1, 1)),
        obs=pd.DataFrame(index=["cell-1"]),
        var=pd.DataFrame(index=["g1"]),
    )
    ctx = PublicationContext(
        svc=SVC(
            expr=output,
            spatial=output,
            svc_kind="sp",
            artifacts={"outputs": {"sp_svc": output}},
            provenance={"primary_output_key": "sp_svc"},
        ),
        runtime={"application_route": "sp-SVC", "svc_kind": "sp"},
        profile="application_sp",
        run_dir=tmp_path / "run",
        merged_config=_ot_config(),
        provenance={},
        artifact_records=[],
    )
    request = SimpleNamespace(
        svc_type="sp-SVC",
        output_root=tmp_path / "out",
        sample_name="sample",
        seed=17,
        ot_method="pot",
    )

    result, path = publish_single(request, ctx)

    assert path == tmp_path / "out" / "sample" / "SVC.h5ad"
    assert read_h5ad(path).shape == result.shape
    assert result.uns["revise_reconstruction"]["profile"] == "application_sp"


def test_publish_pair_writes_the_sc_svc_spatial_and_expression_outputs(tmp_path):
    from revise.application.publication import publish_pair

    spatial, expression = _sc_pair()
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
            provenance={"selected_cell_type": "T"},
        ),
        runtime={"application_route": "sc-SVC", "svc_kind": "sc"},
        profile="application_sc",
        run_dir=tmp_path / "run",
        merged_config={**_ot_config(), "sc": {}},
        provenance={},
        artifact_records=[],
    )
    request = SimpleNamespace(
        svc_type="sc-SVC",
        output_root=tmp_path / "out",
        sample_name="sample",
        select_cell_type="T",
        ot_method="pot",
    )

    results, paths = publish_pair(request, ctx)

    assert set(results) == {"spatial", "expression"}
    assert all(path.exists() for path in paths.values())
    assert paths["spatial"].name == "sc_SVC_spatial.h5ad"
    assert paths["expression"].name == "sc_SVC_expr.h5ad"


def test_pair_rollback_restores_both_previous_files(tmp_path):
    from revise.application.publication import publish_pair

    spatial, expression = _sc_pair()
    output_dir = tmp_path / "out" / "sample" / "sc-SVC" / "T"
    output_dir.mkdir(parents=True)
    previous = {
        "spatial": output_dir / "sc_SVC_spatial.h5ad",
        "expression": output_dir / "sc_SVC_expr.h5ad",
    }
    previous["spatial"].write_bytes(b"previous-spatial")
    previous["expression"].write_bytes(b"previous-expression")
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
            provenance={"selected_cell_type": "T"},
        ),
        runtime={"application_route": "sc-SVC", "svc_kind": "sc"},
        profile="application_sc",
        run_dir=tmp_path / "run",
        merged_config={**_ot_config(), "sc": {}},
        provenance={},
        artifact_records=[],
    )
    request = SimpleNamespace(
        svc_type="sc-SVC",
        output_root=tmp_path / "out",
        sample_name="sample",
        select_cell_type="T",
        ot_method="pot",
    )

    publish_pair(request, ctx)
    ctx.rollback_pending_publication()

    assert previous["spatial"].read_bytes() == b"previous-spatial"
    assert previous["expression"].read_bytes() == b"previous-expression"
