from pathlib import Path


def test_cli_overrides_are_recorded_and_applied_after_yaml():
    import reconstruct

    document = {
        "application": {"svc_type": "sp-SVC"},
        "algorithm": {"ot_method": "pot"},
        "output": {"dir": "output", "name": "old"},
    }
    overrides = reconstruct._apply_overrides(
        document,
        {"svc_type": "sc-SVC", "ot_method": "tacco", "output_name": "new"},
    )

    assert document["application"]["svc_type"] == "sc-SVC"
    assert document["algorithm"]["ot_method"] == "tacco"
    assert document["output"]["name"] == "new"
    assert overrides == {"svc_type": "sc-SVC", "ot_method": "tacco", "output_name": "new"}


def test_output_expansion_is_owned_by_reconstruct():
    import reconstruct
    from types import SimpleNamespace

    config = SimpleNamespace(
        svc_type="sc-SVC",
        output_dir=Path("output"),
        output_name="sample_sc-SVC",
    )
    paths = reconstruct._output_paths(config)
    assert paths["spatial"].name == "sample_sc-SVC_spatial.h5ad"
    assert paths["expression"].name == "sample_sc-SVC_expr.h5ad"
