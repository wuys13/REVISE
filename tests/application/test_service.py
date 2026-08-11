from pathlib import Path


def test_cli_overrides_are_recorded_and_applied_after_yaml():
    from revise.application.config import apply_application_overrides

    document = {
        "application": {"svc_type": "sp-SVC"},
        "algorithm": {"ot_method": "pot"},
        "output": {"dir": "output", "name": "old"},
    }
    overrides = apply_application_overrides(
        document,
        {"svc_type": "sc-SVC", "ot_method": "tacco", "output_name": "new"},
    )

    assert document["application"]["svc_type"] == "sc-SVC"
    assert document["algorithm"]["ot_method"] == "tacco"
    assert document["output"]["name"] == "new"
    assert overrides == {"svc_type": "sc-SVC", "ot_method": "tacco", "output_name": "new"}


def test_engine_overrides_carry_sc_alpha_and_resolutions():
    import reconstruct
    from types import SimpleNamespace

    config = SimpleNamespace(
        svc_type="sc-SVC",
        st_path=Path("st.h5ad"),
        reference_path=Path("ref.h5ad"),
        reference_filter_column=None,
        reference_filter_value=None,
        pm_on_cell_path=None,
        output_dir=Path("output"),
        output_name="run",
        st_format="h5ad",
        spatialdata_table=None,
        spatialdata_element=None,
        broad_column="Level1",
        subtype_column="Level2",
        select_cell_type="T",
        local_refinement_strength=None,
        local_refinement_alpha=0.2,
        local_refinement_resolutions=(0.6, 0.7, 0.8),
        local_refinement_graph_method=None,
        local_refinement_graph_alpha=None,
        local_refinement_graph_n_neighbors=None,
        local_refinement_graph_exp_neighbors=None,
        local_refinement_graph_spatial_neighbors=None,
        local_refinement_match_spot_sum=None,
        ot_method="tacco",
        seed=42,
    )

    _runtime, _io, algorithm = reconstruct._engine_overrides(config)

    assert algorithm["graph"]["alpha"] == 0.2
    assert algorithm["sc"]["resolutions"] == [0.6, 0.7, 0.8]


def test_engine_sample_name_uses_route_identity_for_named_and_unnamed_outputs():
    import reconstruct
    from types import SimpleNamespace

    config = SimpleNamespace(
        svc_type="sp-SVC",
        st_path=Path("st.h5ad"),
        reference_path=Path("ref.h5ad"),
        reference_filter_column="Donor",
        reference_filter_value="D1",
        pm_on_cell_path=None,
        output_dir=Path("output"),
        output_name=None,
        st_format="h5ad",
        spatialdata_table=None,
        spatialdata_element=None,
        broad_column="Level1",
        subtype_column=None,
        select_cell_type=None,
        local_refinement_strength=0.2,
        local_refinement_alpha=None,
        local_refinement_resolutions=None,
        local_refinement_graph_method=None,
        local_refinement_graph_alpha=None,
        local_refinement_graph_n_neighbors=None,
        local_refinement_graph_exp_neighbors=None,
        local_refinement_graph_spatial_neighbors=None,
        local_refinement_match_spot_sum=None,
        ot_method="pot",
        seed=42,
    )

    for output_name in (None, "foo"):
        config.output_name = output_name
        _runtime, io, _algorithm = reconstruct._engine_overrides(config)

        assert io["sample_name"] == "sp-SVC"


def test_public_run_returns_pipeline_artifact_and_dry_run_returns_none(
    monkeypatch,
    tmp_path,
):
    import reconstruct
    from types import SimpleNamespace

    document = {
        "schema_version": 1,
        "application": {"svc_type": "sp-SVC"},
        "paths": {"root_dir": str(tmp_path)},
        "algorithm": {"ot_method": "pot"},
        "inputs": {
            "st": {"path": "st.h5ad", "format": "h5ad"},
            "reference": {"path": "ref.h5ad", "format": "h5ad"},
        },
        "preprocessing": {
            "spatial": {"min_transcript_counts": None, "min_cell_counts": 0},
            "reference": {"min_transcript_counts": None, "min_cell_counts": 0},
        },
        "global_anchoring": {"broad_column": "Level1"},
        "local_refinement": {"strength": 0.2},
        "output": {"dir": "output", "name": None},
        "execution": {"seed": 42},
    }
    config_path = tmp_path / "run.yaml"
    import yaml

    config_path.write_text(yaml.safe_dump(document), encoding="utf-8")
    artifact = _fake_adata()

    class FakePipeline:
        def run(self, **kwargs):
            if kwargs["dry_run"]:
                return SimpleNamespace(
                    provenance={"dry_run": True, "run_dir": str(tmp_path / "run")},
                    artifacts={},
                    summary=lambda: {},
                )
            svc = SimpleNamespace(
                provenance={"primary_output_key": "sp_svc", "run_dir": str(tmp_path / "run")},
                artifacts={"outputs": {"sp_svc": artifact}},
                summary=lambda: {},
            )
            kwargs["finalize_callback"](
                SimpleNamespace(
                    svc=svc,
                    profile="application_test",
                    run_dir=tmp_path / "run",
                    merged_config={"ot": {}},
                )
            )
            return svc

    monkeypatch.setattr(reconstruct, "REVISEPipeline", FakePipeline)

    result = reconstruct.run_application(config_path)
    dry_result = reconstruct.run_application(config_path, dry_run=True)

    assert result is artifact
    assert dry_result is None
    assert not hasattr(reconstruct, "ApplicationExecution")


def _fake_adata():
    import numpy as np
    import pandas as pd
    from anndata import AnnData

    return AnnData(
        X=np.ones((1, 1)),
        obs=pd.DataFrame(index=["cell"]),
        var=pd.DataFrame(index=["gene"]),
    )
