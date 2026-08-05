from __future__ import annotations

from pathlib import Path

import pytest
import yaml


def _document(svc_type: str = "sp-SVC") -> dict:
    local_refinement = (
        {"subtype_column": "Level2", "select_cell_type": "T cell"}
        if svc_type == "sc-SVC"
        else {"strength": 0.2}
    )
    return {
        "schema_version": 1,
        "application": {"svc_type": svc_type, "sample_name": "sample"},
        "paths": {"root_dir": "."},
        "algorithm": {"ot_method": "pot"},
        "inputs": {
            "mode": "direct",
            "st": {"path": "data/sample_st.h5ad", "format": "h5ad"},
            "reference": {
                "path": "data/sc_ref.h5ad",
                "format": "h5ad",
                "patient_key": "Patient",
            },
        },
        "global_anchoring": {"broad_column": "Level1"},
        "local_refinement": local_refinement,
        "output": {"path": "output"},
        "execution": {"action": "run", "seed": 42},
    }


def _write_config(tmp_path: Path, document: dict, name: str = "run.yaml") -> Path:
    path = tmp_path / "configs" / "application" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(document, sort_keys=False), encoding="utf-8")
    return path


def test_dot_root_resolves_against_runtime_cwd(monkeypatch, tmp_path):
    from revise.application.request import load_application_request

    runtime_root = tmp_path / "runtime"
    runtime_root.mkdir()
    path = _write_config(tmp_path, _document())
    monkeypatch.chdir(runtime_root)

    request = load_application_request(path)

    assert request.config_path == path.resolve()
    assert request.declared_root == "."
    assert request.resolved_root == runtime_root.resolve()
    assert request.cwd == runtime_root.resolve()
    assert request.st_path == runtime_root / "data/sample_st.h5ad"
    assert request.reference_path == runtime_root / "data/sc_ref.h5ad"
    assert request.output_root == runtime_root / "output"
    assert not hasattr(request, "base_config")


def test_existing_absolute_root_is_accepted(tmp_path):
    from revise.application.request import load_application_request

    root = tmp_path / "workspace"
    root.mkdir()
    document = _document()
    document["paths"]["root_dir"] = str(root)
    request = load_application_request(_write_config(tmp_path, document))

    assert request.resolved_root == root.resolve()
    assert request.st_path == root / "data/sample_st.h5ad"


def test_root_must_be_dot_or_an_existing_absolute_directory(tmp_path):
    from revise.application.request import ApplicationConfigError, load_application_request

    document = _document()
    document["paths"]["root_dir"] = "relative"

    with pytest.raises(ApplicationConfigError, match=r"paths\.root_dir"):
        load_application_request(_write_config(tmp_path, document))


def test_nonexistent_absolute_root_is_rejected(tmp_path):
    from revise.application.request import ApplicationConfigError, load_application_request

    document = _document()
    document["paths"]["root_dir"] = str(tmp_path / "missing")

    with pytest.raises(ApplicationConfigError, match="existing absolute directory"):
        load_application_request(_write_config(tmp_path, document))


def test_direct_child_paths_cannot_escape_the_root(tmp_path):
    from revise.application.request import ApplicationConfigError, load_application_request

    document = _document()
    document["inputs"]["st"]["path"] = "../outside.h5ad"

    with pytest.raises(ApplicationConfigError, match=r"inputs\."):
        load_application_request(_write_config(tmp_path, document))


def test_input_paths_may_be_missing_and_output_is_not_created_by_compilation(tmp_path):
    from revise.application.request import load_application_request

    root = tmp_path / "root"
    root.mkdir()
    document = _document()
    document["paths"]["root_dir"] = str(root)

    request = load_application_request(_write_config(tmp_path, document))

    assert not request.st_path.exists()
    assert not request.reference_path.exists()
    assert not request.output_root.exists()


def test_dry_run_truth_table_only_marks_run_override(tmp_path):
    from revise.application.request import load_application_request

    run_path = _write_config(tmp_path, _document(), "run.yaml")
    preflight = _document()
    preflight["execution"]["action"] = "preflight"
    preflight_path = _write_config(tmp_path, preflight, "preflight.yaml")

    cases = [
        (run_path, False, "run", False),
        (run_path, True, "preflight", True),
        (preflight_path, False, "preflight", False),
        (preflight_path, True, "preflight", False),
    ]
    for path, flag, effective, override in cases:
        request = load_application_request(path, dry_run=flag)
        assert request.effective_action == effective
        assert request.dry_run_override is override


def test_legacy_layout_preserves_sample_prefix_and_data_root(tmp_path):
    from revise.application.request import load_application_request

    root = tmp_path / "root"
    root.mkdir()
    document = _document("sc-SVC-sr")
    document["paths"]["root_dir"] = str(root)
    document["application"]["sample_name"] = "donor1"
    document["inputs"] = {
        "mode": "legacy_layout",
        "data_root": "data",
        "st": {"file": "st.h5ad", "format": "h5ad"},
        "reference": {
            "file": "sc_ref.h5ad",
            "format": "h5ad",
            "patient_key": "Patient",
        },
    }

    request = load_application_request(_write_config(tmp_path, document))

    assert request.data_root == root / "data"
    assert request.st_path == root / "data/donor1_st.h5ad"
    assert request.reference_path == root / "data/sc_ref.h5ad"


def test_direct_and_legacy_layout_fields_are_strictly_exclusive(tmp_path):
    from revise.application.request import ApplicationConfigError, load_application_request

    document = _document()
    document["inputs"]["st"]["file"] = "st.h5ad"

    with pytest.raises(ApplicationConfigError, match="unknown field"):
        load_application_request(_write_config(tmp_path, document))


@pytest.mark.parametrize(
    ("svc_type", "field", "value"),
    [
        ("sc-SVC", "strength", 0.2),
        ("sp-SVC", "select_cell_type", "T cell"),
        ("sc-SVC-sr", "select_cell_type", "T cell"),
    ],
)
def test_local_refinement_route_fields_fail_closed(tmp_path, svc_type, field, value):
    from revise.application.request import ApplicationConfigError, load_application_request

    document = _document(svc_type)
    document["local_refinement"][field] = value

    with pytest.raises(ApplicationConfigError, match="not allowed|unknown field"):
        load_application_request(_write_config(tmp_path, document))


def test_sc_svc_requires_one_concrete_safe_cell_type(tmp_path):
    from revise.application.request import ApplicationConfigError, load_application_request

    document = _document("sc-SVC")
    document["local_refinement"]["select_cell_type"] = None

    with pytest.raises(ApplicationConfigError, match="concrete broad cell type"):
        load_application_request(_write_config(tmp_path, document))


def test_output_components_cannot_contain_path_separators(tmp_path):
    from revise.application.request import ApplicationConfigError, load_application_request

    document = _document()
    document["application"]["sample_name"] = "sample/name"

    with pytest.raises(ApplicationConfigError, match="sample_name"):
        load_application_request(_write_config(tmp_path, document))


def test_strength_must_be_non_negative(tmp_path):
    from revise.application.request import ApplicationConfigError, load_application_request

    document = _document()
    document["local_refinement"]["strength"] = -0.1

    with pytest.raises(ApplicationConfigError, match="non-negative finite"):
        load_application_request(_write_config(tmp_path, document))


def test_seed_must_fit_numpy_random_state_range(tmp_path):
    from revise.application.request import ApplicationConfigError, load_application_request

    document = _document()
    document["execution"]["seed"] = -1

    with pytest.raises(ApplicationConfigError, match="integer|between"):
        load_application_request(_write_config(tmp_path, document))


def test_unknown_fields_are_not_silently_ignored(tmp_path):
    from revise.application.request import ApplicationConfigError, load_application_request

    document = _document()
    document["extra"] = {}

    with pytest.raises(ApplicationConfigError, match="unknown field.*extra"):
        load_application_request(_write_config(tmp_path, document))


def test_full_engine_config_explains_that_engine_yaml_is_internal(tmp_path):
    from revise.application.request import ApplicationConfigError, load_application_request

    path = _write_config(
        tmp_path,
        {"version": 2, "defaults": {}, "router": {}, "profiles": {}},
    )

    with pytest.raises(ApplicationConfigError, match="full engine config.*internally managed"):
        load_application_request(path)
