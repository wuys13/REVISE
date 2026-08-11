from __future__ import annotations

from pathlib import Path

import pytest
import yaml


def _document(svc_type: str = "sp-SVC") -> dict:
    local_refinement = (
        {
            "subtype_column": "Level2",
            "select_cell_type": "T",
            "alpha": 0.2,
            "resolutions": [0.6, 0.7, 0.8],
        }
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
            "reference": {
                "path": "data/sc_ref.h5ad",
                "format": "h5ad",
                "filter_column": None,
                "filter_value": None,
            },
        },
        "preprocessing": {
            "spatial": {"min_transcript_counts": 60, "min_cell_counts": 100},
            "reference": {"min_transcript_counts": None, "min_cell_counts": 100},
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


def test_output_name_may_be_missing_or_null(tmp_path):
    from revise.application.config import compile_application_config, load_application_yaml

    for output in ({"dir": "output"}, {"dir": "output", "name": None}):
        document = _document()
        document["output"] = output
        source, effective = load_application_yaml(_write_config(tmp_path, document))

        config = compile_application_config(effective, source=source)

        assert config.output_name is None


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


def test_reference_filter_fields_must_be_supplied_together(tmp_path):
    from revise.application.config import ApplicationConfigError, compile_application_config, load_application_yaml

    document = _document()
    document["inputs"]["reference"]["filter_column"] = "Patient"
    source, effective = load_application_yaml(_write_config(tmp_path, document))

    with pytest.raises(ApplicationConfigError, match="filter_column.*filter_value"):
        compile_application_config(effective, source=source)


def test_sc_preprocessing_and_local_refinement_parameters_are_compiled(tmp_path):
    from revise.application.config import compile_application_config, load_application_yaml

    document = _document("sc-SVC")
    document["inputs"]["reference"].update(
        filter_column="Patient",
        filter_value="P2CRC",
    )
    document["local_refinement"].update(
        alpha=0.2,
        resolutions=[0.6, 0.7, 0.8],
    )
    source, effective = load_application_yaml(_write_config(tmp_path, document))
    config = compile_application_config(effective, source=source)

    assert config.reference_filter_column == "Patient"
    assert config.reference_filter_value == "P2CRC"
    assert config.spatial_min_transcript_counts == 60
    assert config.spatial_min_cell_counts == 100
    assert config.reference_min_transcript_counts is None
    assert config.reference_min_cell_counts == 100
    assert config.local_refinement_alpha == 0.2
    assert config.local_refinement_resolutions == (0.6, 0.7, 0.8)


def test_sc_sr_local_refinement_graph_and_match_spot_sum_are_compiled(tmp_path):
    from revise.application.config import compile_application_config, load_application_yaml

    document = _document("sc-SVC-sr")
    document["local_refinement"] = {
        "strength": 0.0,
        "graph": {
            "method": "pca",
            "alpha": 0.2,
            "n_neighbors": 10,
            "exp_neighbors": 10,
            "spatial_neighbors": 10,
        },
        "match_spot_sum": True,
    }
    source, effective = load_application_yaml(_write_config(tmp_path, document))

    config = compile_application_config(effective, source=source)

    assert config.local_refinement_graph_method == "pca"
    assert config.local_refinement_graph_alpha == 0.2
    assert config.local_refinement_graph_n_neighbors == 10
    assert config.local_refinement_graph_exp_neighbors == 10
    assert config.local_refinement_graph_spatial_neighbors == 10
    assert config.local_refinement_match_spot_sum is True


def test_sp_preprocessing_count_and_gene_thresholds_are_compiled(tmp_path):
    from revise.application.config import compile_application_config, load_application_yaml

    document = _document("sp-SVC")
    document["preprocessing"] = {
        "spatial": {
            "min_transcript_counts": None,
            "min_counts": 20,
            "min_cell_counts": 30,
        },
        "reference": {
            "min_transcript_counts": None,
            "min_genes": 20,
            "min_cell_counts": 50,
        },
    }
    source, effective = load_application_yaml(_write_config(tmp_path, document))

    config = compile_application_config(effective, source=source)

    assert config.spatial_min_counts == 20
    assert config.reference_min_genes == 20
    assert config.spatial_min_transcript_counts is None
    assert config.reference_min_transcript_counts is None


def test_effective_request_metadata_hash_covers_application_algorithm_inputs(tmp_path):
    from revise.application.config import compile_application_config, load_application_yaml
    from revise.application.publication import application_metadata, output_paths
    from revise.utils.provenance import hash_jsonable

    document = _document("sc-SVC")
    document["inputs"]["reference"].update(
        filter_column="Patient",
        filter_value="P2CRC",
    )
    source, effective_document = load_application_yaml(_write_config(tmp_path, document))
    config = compile_application_config(effective_document, source=source)

    metadata = application_metadata(
        config,
        cli_overrides={},
        paths=output_paths(config),
        dry_run=True,
    )

    request = metadata["effective_request"]
    assert request["preprocessing"] == {
        "spatial": {
            "min_transcript_counts": 60,
            "min_counts": None,
            "min_cell_counts": 100,
        },
        "reference": {
            "min_transcript_counts": None,
            "min_genes": None,
            "min_cell_counts": 100,
        },
    }
    assert request["inputs"]["reference_filter"] == {
        "column": "Patient",
        "value": "P2CRC",
    }
    assert request["local_refinement"] == document["local_refinement"]
    assert metadata["effective_request_hash"] == hash_jsonable(request)


def test_application_null_seed_resolves_to_authority_default(tmp_path):
    from revise.application.config import compile_application_config, load_application_yaml
    from revise.application.publication import application_metadata, output_paths

    document = _document()
    document["execution"] = {"seed": None}
    source, effective = load_application_yaml(_write_config(tmp_path, document))
    config = compile_application_config(effective, source=source)

    metadata = application_metadata(
        config,
        cli_overrides={},
        paths=output_paths(config),
        dry_run=False,
    )

    assert config.seed == 42
    assert metadata["effective_request"]["execution"] == {"seed": 42}


def test_removed_action_is_rejected(tmp_path):
    from revise.application.config import ApplicationConfigError, compile_application_config, load_application_yaml

    document = _document()
    document["execution"]["action"] = "preflight"
    source, effective = load_application_yaml(_write_config(tmp_path, document))
    with pytest.raises(ApplicationConfigError, match="execution.action was removed"):
        compile_application_config(effective, source=source)
