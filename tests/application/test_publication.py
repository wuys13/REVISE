from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, read_h5ad


def _adata():
    return AnnData(
        X=np.ones((1, 1)),
        obs=pd.DataFrame(index=["cell-1"]),
        var=pd.DataFrame(index=["g1"]),
    )


def _config(tmp_path, svc_type, *, mode=None, output_name="sample_result"):
    if svc_type == "sc-SVC" and mode is None:
        mode = "cluster"
    output_root = tmp_path / "out-root"
    output_dir = output_root / "T" if mode == "cluster" else output_root
    return SimpleNamespace(
        svc_type=svc_type,
        mode=mode,
        output_root=output_root,
        output_dir=output_dir,
        output_name=output_name,
        select_cell_type="T" if mode == "cluster" else None,
    )


def _ctx(tmp_path, outputs, *, application_config_metadata=None):
    context = SimpleNamespace(
        svc=SimpleNamespace(
            artifacts={"outputs": outputs},
        ),
        profile="application_test",
        run_dir=tmp_path / "run",
        merged_config={"ot": {"ga": {"solver": "pot"}, "lr": {"solver": "pot"}}},
        application_config_metadata=application_config_metadata or {},
        artifact_records=[],
    )
    context.set_pending_publication = lambda *, commit, rollback: setattr(
        context, "pending_publication", (commit, rollback)
    )
    context.record_artifact = context.artifact_records.append
    return context


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


def test_sc_pair_rolls_back_both_public_files_when_second_replace_fails(
    monkeypatch,
    tmp_path,
):
    import revise.application.publication as publication

    config = _config(tmp_path, "sc-SVC", output_name=None)
    paths = publication.output_paths(config)
    config.output_dir.mkdir(parents=True)
    paths["spatial"].write_bytes(b"previous-spatial")
    paths["expression"].write_bytes(b"previous-expression")
    ctx = _ctx(
        tmp_path,
        {"sc_svc_spatial": _adata(), "sc_svc_expr": _adata()},
    )
    ctx.svc.provenance = {"primary_output_key": "sc_svc_expr"}
    original_replace = publication.os.replace

    def fail_expression_replace(source, target):
        if Path(target) == paths["expression"] and str(source).endswith(".tmp.h5ad"):
            raise OSError("expression publication failed")
        return original_replace(source, target)

    monkeypatch.setattr(publication.os, "replace", fail_expression_replace)

    with pytest.raises(OSError, match="expression publication failed"):
        publication.publish_outputs(config, paths, ctx)

    assert paths["spatial"].read_bytes() == b"previous-spatial"
    assert paths["expression"].read_bytes() == b"previous-expression"
    assert ctx.artifact_records == []


def test_named_sc_pair_uses_name_as_filename_prefix(tmp_path):
    from revise.application.publication import output_paths

    paths = output_paths(_config(tmp_path, "sc-SVC", output_name="foo"))

    assert paths["spatial"].name == "foo_spatial.h5ad"
    assert paths["expression"].name == "foo_expr.h5ad"


def test_sr_returns_and_writes_graphagg_when_it_is_primary(tmp_path):
    from revise.application.publication import output_paths, publish_outputs

    config = _config(tmp_path, "sc-SVC", mode="sr", output_name="foo")
    raw = AnnData(
        X=np.array([[1.0, 2.0], [3.0, 4.0]]),
        obs=pd.DataFrame(index=["raw-1", "raw-2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    graphagg = AnnData(
        X=np.array([[10.0, 20.0], [30.0, 40.0], [50.0, 60.0]]),
        obs=pd.DataFrame(index=["graphagg-1", "graphagg-2", "graphagg-3"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    ctx = _ctx(
        tmp_path,
        {"sc_svc_dec": raw, "sc_svc_dec_graphagg": graphagg},
    )
    ctx.svc.provenance = {"primary_output_key": "sc_svc_dec_graphagg"}
    paths = output_paths(config)

    result = publish_outputs(config, paths, ctx)

    assert result is graphagg
    assert paths["svc"].name == "foo.h5ad"
    published = read_h5ad(paths["svc"])
    assert published.shape == graphagg.shape
    assert published.obs_names.tolist() == ["graphagg-1", "graphagg-2", "graphagg-3"]
    np.testing.assert_array_equal(published.X, graphagg.X)
    assert published.shape != raw.shape


def test_h5ad_reuses_pipeline_application_metadata_and_adds_publication_fields(tmp_path):
    from revise.application.publication import output_paths, publish_outputs
    from revise.utils.provenance import hash_jsonable

    config = _config(tmp_path, "sc-SVC", mode="cluster", output_name="foo")
    request = {
        "application_route": "sc-SVC",
        "application_mode": "cluster",
        "selected_cell_type": "T",
        "output": {
            "root": str(config.output_root),
            "dir": str(config.output_dir),
            "name": "foo",
        },
    }
    application_config_metadata = {
        "source_path": "application.yaml",
        "output_root": str(config.output_root),
        "output_dir": str(config.output_dir),
        "effective_request": request,
        "effective_request_hash": hash_jsonable(request),
    }
    spatial = _adata()
    expression = _adata()
    ctx = _ctx(
        tmp_path,
        {"sc_svc_spatial": spatial, "sc_svc_expr": expression},
        application_config_metadata=application_config_metadata,
    )
    ctx.svc.provenance = {"primary_output_key": "sc_svc_expr"}

    publish_outputs(config, output_paths(config), ctx)

    metadata = read_h5ad(config.output_dir / "foo_spatial.h5ad").uns[
        "revise_reconstruction"
    ]
    for key, value in application_config_metadata.items():
        assert metadata[key] == value
    assert metadata["profile"] == "application_test"
    assert metadata["run_manifest"] == str(tmp_path / "run" / "provenance.json")
    assert metadata["output_role"] == "spatial"
