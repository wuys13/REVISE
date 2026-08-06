from __future__ import annotations

from pathlib import Path

import pytest
import yaml


def _document(svc_type: str = "sp-SVC") -> dict:
    local_refinement = (
        {"subtype_column": "Level2", "select_cell_type": "T"}
        if svc_type == "sc-SVC"
        else {"strength": 0.2}
    )
    return {
        "schema_version": 1,
        "application": {"svc_type": svc_type},
        "paths": {"root_dir": "."},
        "algorithm": {"ot_method": "pot"},
        "inputs": {
            "st": {"path": "data/sample_st.h5ad", "format": "h5ad"},
            "reference": {"path": "data/sc_ref.h5ad", "format": "h5ad"},
        },
        "global_anchoring": {"broad_column": "Level1"},
        "local_refinement": local_refinement,
        "output": {"dir": "output", "name": "sample_sp-SVC"},
        "execution": {"seed": 42},
    }


def _write_config(tmp_path: Path, document: dict, name: str = "run.yaml") -> Path:
    path = tmp_path / name
    path.write_text(yaml.safe_dump(document, sort_keys=False), encoding="utf-8")
    return path


def test_dot_root_resolves_against_runtime_cwd(monkeypatch, tmp_path):
    from revise.application.config import compile_application_config, load_application_yaml

    runtime_root = tmp_path / "runtime"
    runtime_root.mkdir()
    path = _write_config(tmp_path, _document())
    monkeypatch.chdir(runtime_root)
    source, document = load_application_yaml(path)
    config = compile_application_config(document, source=source)

    assert config.resolved_root == runtime_root.resolve()
    assert config.st_path == runtime_root / "data/sample_st.h5ad"
    assert config.reference_path == runtime_root / "data/sc_ref.h5ad"
    assert config.output_dir == runtime_root / "output"


def test_input_paths_may_be_missing_and_compilation_does_not_create_output(tmp_path):
    from revise.application.config import compile_application_config, load_application_yaml

    root = tmp_path / "root"
    root.mkdir()
    document = _document()
    document["paths"]["root_dir"] = str(root)
    source, effective = load_application_yaml(_write_config(tmp_path, document))
    config = compile_application_config(effective, source=source)

    assert not config.st_path.exists()
    assert not config.reference_path.exists()
    assert not config.output_dir.exists()


def test_route_specific_fields_fail_closed(tmp_path):
    from revise.application.config import ApplicationConfigError, compile_application_config, load_application_yaml

    document = _document("sp-SVC")
    document["local_refinement"]["select_cell_type"] = "T"
    source, effective = load_application_yaml(_write_config(tmp_path, document))
    with pytest.raises(ApplicationConfigError, match="unknown field"):
        compile_application_config(effective, source=source)


def test_sc_svc_requires_one_concrete_cell_type(tmp_path):
    from revise.application.config import ApplicationConfigError, compile_application_config, load_application_yaml

    document = _document("sc-SVC")
    document["local_refinement"]["select_cell_type"] = None
    source, effective = load_application_yaml(_write_config(tmp_path, document))
    with pytest.raises(ApplicationConfigError, match="concrete broad cell type"):
        compile_application_config(effective, source=source)


def test_removed_action_is_rejected(tmp_path):
    from revise.application.config import ApplicationConfigError, compile_application_config, load_application_yaml

    document = _document()
    document["execution"]["action"] = "preflight"
    source, effective = load_application_yaml(_write_config(tmp_path, document))
    with pytest.raises(ApplicationConfigError, match="execution.action was removed"):
        compile_application_config(effective, source=source)
