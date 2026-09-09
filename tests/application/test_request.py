from __future__ import annotations

from pathlib import Path

import pytest
import yaml


def _document(svc_type: str = "sp-SVC", mode: str | None = None) -> dict:
    local_refinement = (
        {
            "subtype_column": "Level2",
            "select_cell_type": "T",
            "alpha": 0.2,
            "resolutions": [0.6, 0.7, 0.8],
        }
        if svc_type == "sc-SVC" and mode != "sr"
        else {"strength": 0.2}
    )
    application = {"svc_type": svc_type}
    if mode is not None:
        application["mode"] = mode
    return {
        "schema_version": 1,
        "application": application,
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

    document = _document("sc-SVC", "cluster")
    document["local_refinement"]["select_cell_type"] = None
    source, effective = load_application_yaml(_write_config(tmp_path, document))
    with pytest.raises(ApplicationConfigError, match="concrete broad cell type"):
        compile_application_config(effective, source=source)


@pytest.mark.parametrize(
    ("selected", "normalized"),
    [
        (" T ", "T"),
        ("Fibroblast", "Fibroblast"),
        ("Mono/Macro", "Mono_Macro"),
        (" Mono/Macro ", "Mono_Macro"),
    ],
)
def test_sc_cluster_output_dir_is_derived_from_output_root_and_normalized_selection(
    tmp_path,
    selected,
    normalized,
):
    from revise.application.config import compile_application_config, load_application_yaml

    document = _document("sc-SVC", "cluster")
    document["paths"]["root_dir"] = str(tmp_path)
    document["local_refinement"]["select_cell_type"] = selected
    document["output"]["dir"] = "output/P2CRC_Xenium"
    source, effective = load_application_yaml(_write_config(tmp_path, document))

    config = compile_application_config(effective, source=source)

    assert config.select_cell_type == normalized
    assert config.output_root == tmp_path / "output/P2CRC_Xenium"
    assert config.output_dir == config.output_root / normalized


@pytest.mark.parametrize(
    "selected",
    ["", "all", "*", ".", "..", "../T", r"T\\x", "T\x00", "T\x7f", "T\x85"],
)
def test_sc_cluster_rejects_unsafe_cell_type_for_output_directory(tmp_path, selected):
    from revise.application.config import ApplicationConfigError, compile_application_config, load_application_yaml

    document = _document("sc-SVC", "cluster")
    document["local_refinement"]["select_cell_type"] = selected
    source, effective = load_application_yaml(_write_config(tmp_path, document))

    with pytest.raises(ApplicationConfigError, match="concrete broad cell type"):
        compile_application_config(effective, source=source)


def test_cluster_override_rederives_output_dir_from_original_output_root(tmp_path):
    from revise.application.config import (
        compile_application_config,
        load_application_yaml,
        override_select_cell_type,
    )

    document = _document("sc-SVC", "cluster")
    document["output"]["dir"] = "output/P2CRC_Xenium"
    source, effective = load_application_yaml(_write_config(tmp_path, document))
    config = compile_application_config(effective, source=source)

    fib = override_select_cell_type(config, "Fibroblast")
    mono = override_select_cell_type(fib, "Mono/Macro")

    assert fib.output_dir == config.output_root / "Fibroblast"
    assert mono.select_cell_type == "Mono_Macro"
    assert mono.output_dir == config.output_root / "Mono_Macro"
    assert mono.output_dir != config.output_dir / "T"


def test_reference_filter_fields_must_be_supplied_together(tmp_path):
    from revise.application.config import ApplicationConfigError, compile_application_config, load_application_yaml

    document = _document()
    document["inputs"]["reference"]["filter_column"] = "Patient"
    source, effective = load_application_yaml(_write_config(tmp_path, document))

    with pytest.raises(ApplicationConfigError, match="filter_column.*filter_value"):
        compile_application_config(effective, source=source)


def test_sc_preprocessing_and_local_refinement_parameters_are_compiled(tmp_path):
    from revise.application.config import compile_application_config, load_application_yaml

    document = _document("sc-SVC", "cluster")
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


def test_sc_sr_mode_local_refinement_graph_and_match_spot_sum_are_compiled(tmp_path):
    from revise.application.config import compile_application_config, load_application_yaml

    document = _document("sc-SVC", "sr")
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

    assert config.mode == "sr"
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

    document = _document("sc-SVC", "cluster")
    document["inputs"]["reference"].update(
        filter_column="Patient",
        filter_value="P2CRC",
    )
    source, effective_document = load_application_yaml(_write_config(tmp_path, document))
    config = compile_application_config(effective_document, source=source)

    metadata = application_metadata(
        config,
        paths=output_paths(config),
    )

    request = metadata["effective_request"]
    assert request["mode"] == "cluster"
    assert request["application_route"] == "sc-SVC"
    assert request["application_mode"] == "cluster"
    assert request["selected_cell_type"] == "T"
    assert request["output"] == {
        "root": str(config.output_root),
        "dir": str(config.output_dir),
        "name": "sample_sp-SVC",
    }
    assert metadata["output_root"] == str(config.output_root)
    assert metadata["output_dir"] == str(config.output_dir)
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
        paths=output_paths(config),
    )

    assert config.seed == 42
    assert metadata["effective_request"]["execution"] == {"seed": 42}


def test_removed_action_is_rejected(tmp_path):
    from revise.application.config import ApplicationConfigError, compile_application_config, load_application_yaml

    document = _document()
    document["execution"]["action"] = "preflight"
    source, effective = load_application_yaml(_write_config(tmp_path, document))
    with pytest.raises(ApplicationConfigError, match="execution.action is not supported"):
        compile_application_config(effective, source=source)


def test_sc_svc_requires_an_explicit_mode_with_migration_guidance(tmp_path):
    from revise.application.config import ApplicationConfigError, compile_application_config, load_application_yaml

    source, effective = load_application_yaml(_write_config(tmp_path, _document("sc-SVC")))

    with pytest.raises(ApplicationConfigError, match="application.mode is required"):
        compile_application_config(effective, source=source)


def test_legacy_sc_svc_sr_type_is_rejected_with_hard_cut_guidance(tmp_path):
    from revise.application.config import ApplicationConfigError, compile_application_config, load_application_yaml

    source, effective = load_application_yaml(_write_config(tmp_path, _document("sc-SVC-sr")))

    with pytest.raises(ApplicationConfigError, match="sc-SVC-sr.*mode: sr"):
        compile_application_config(effective, source=source)


def test_sp_svc_rejects_application_mode(tmp_path):
    from revise.application.config import ApplicationConfigError, compile_application_config, load_application_yaml

    source, effective = load_application_yaml(_write_config(tmp_path, _document("sp-SVC", "cluster")))

    with pytest.raises(ApplicationConfigError, match="application.mode is only valid for sc-SVC"):
        compile_application_config(effective, source=source)
