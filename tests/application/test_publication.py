from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from anndata import AnnData, read_h5ad


def _adata():
    return AnnData(
        X=np.ones((1, 1)),
        obs=pd.DataFrame(index=["cell-1"]),
        var=pd.DataFrame(index=["g1"]),
    )


def _config(tmp_path, svc_type, *, output_name="sample_result"):
    return SimpleNamespace(
        svc_type=svc_type,
        output_dir=tmp_path / "out",
        output_name=output_name,
        select_cell_type="T" if svc_type == "sc-SVC" else None,
    )


def _ctx(tmp_path, outputs):
    return SimpleNamespace(
        svc=SimpleNamespace(
            artifacts={"outputs": outputs},
        ),
        profile="application_test",
        run_dir=tmp_path / "run",
        merged_config={"ot": {"ga": {"solver": "pot"}, "lr": {"solver": "pot"}}},
    )


def test_single_output_uses_logical_output_name(tmp_path):
    from revise.application.publication import output_paths, publish_outputs

    config = _config(tmp_path, "sp-SVC")
    paths = output_paths(config)
    adata = _adata()
    before_x = adata.X.copy()
    before_obs = adata.obs_names.copy()
    before_var = adata.var_names.copy()
    ctx = _ctx(tmp_path, {"sp_svc": adata})
    ctx.svc.provenance = {"primary_output_key": "sp_svc"}

    result = publish_outputs(config, paths, ctx)

    assert result is adata
    assert paths["svc"].name == "sample_result.h5ad"
    assert read_h5ad(paths["svc"]).uns["revise_reconstruction"]["output_role"] == "svc"
    np.testing.assert_array_equal(adata.X, before_x)
    assert adata.obs_names.equals(before_obs)
    assert adata.var_names.equals(before_var)


def test_single_output_without_name_uses_svc_filename(tmp_path):
    from revise.application.publication import output_paths

    assert output_paths(_config(tmp_path, "sp-SVC", output_name=None))["svc"].name == "svc.h5ad"


def test_sc_pair_writes_both_fixed_artifact_roles(tmp_path):
    from revise.application.publication import output_paths, publish_outputs

    config = _config(tmp_path, "sc-SVC", output_name=None)
    paths = output_paths(config)
    spatial = _adata()
    expression = _adata()
    ctx = _ctx(tmp_path, {"sc_svc_spatial": spatial, "sc_svc_expr": expression})
    ctx.svc.provenance = {"primary_output_key": "sc_svc_expr"}

    result = publish_outputs(config, paths, ctx)

    assert result[0] is spatial
    assert result[1] is expression
    assert paths["spatial"].exists()
    assert paths["expression"].exists()
    assert paths["spatial"].name == "spatial.h5ad"
    assert paths["expression"].name == "expr.h5ad"


def test_named_sc_pair_uses_name_as_filename_prefix(tmp_path):
    from revise.application.publication import output_paths

    paths = output_paths(_config(tmp_path, "sc-SVC", output_name="foo"))

    assert paths["spatial"].name == "foo_spatial.h5ad"
    assert paths["expression"].name == "foo_expr.h5ad"


def test_sr_returns_and_writes_graphagg_when_it_is_primary(tmp_path):
    from revise.application.publication import output_paths, publish_outputs

    config = _config(tmp_path, "sc-SVC-sr", output_name="foo")
    raw = _adata()
    graphagg = _adata()
    ctx = _ctx(
        tmp_path,
        {"sc_svc_dec": raw, "sc_svc_dec_graphagg": graphagg},
    )
    ctx.svc.provenance = {"primary_output_key": "sc_svc_dec_graphagg"}
    paths = output_paths(config)

    result = publish_outputs(config, paths, ctx)

    assert result is graphagg
    assert paths["svc"].name == "foo.h5ad"
    assert read_h5ad(paths["svc"]).shape == graphagg.shape
