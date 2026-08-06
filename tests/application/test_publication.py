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


def _config(tmp_path, svc_type):
    return SimpleNamespace(
        svc_type=svc_type,
        output_dir=tmp_path / "out",
        output_name="sample_result",
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
    import reconstruct

    config = _config(tmp_path, "sp-SVC")
    paths = reconstruct._output_paths(config)
    adata = _adata()
    written = reconstruct._write_outputs(config, paths, _ctx(tmp_path, {"sp_svc": adata}))

    assert paths["svc"].name == "sample_result.h5ad"
    assert read_h5ad(paths["svc"]).uns["revise_reconstruction"]["output_role"] == "svc"
    assert written["svc"].shape == adata.shape


def test_sc_pair_writes_both_fixed_artifact_roles(tmp_path):
    import reconstruct

    config = _config(tmp_path, "sc-SVC")
    paths = reconstruct._output_paths(config)
    spatial = _adata()
    expression = _adata()
    reconstruct._write_outputs(
        config,
        paths,
        _ctx(tmp_path, {"sc_svc_spatial": spatial, "sc_svc_expr": expression}),
    )

    assert paths["spatial"].exists()
    assert paths["expression"].exists()
    assert paths["spatial"].name.endswith("_spatial.h5ad")
    assert paths["expression"].name.endswith("_expr.h5ad")
