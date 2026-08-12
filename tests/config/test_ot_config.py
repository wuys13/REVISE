from __future__ import annotations

import copy
import importlib
import inspect
import logging
import math
import sys
import types
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from revise.config import (
    ConfigError,
    merge_unified_config,
    resolve_semantic_route,
)
from revise.config.authority import _authority_document
from revise.config.loader import _validate_raw_config
from revise.framework import REVISEPipeline
from revise.recon.context import PipelineContext
from revise.svc import SVC


ROOT = Path(__file__).parents[2]


def test_semantic_router_resolves_each_mode_without_cross_domain_selector():
    raw = _authority_document()

    application = resolve_semantic_route(raw, svc_type="sc-SVC")
    benchmark = resolve_semantic_route(raw, cf="segmentation")

    assert application == {
        "mode": "application",
        "application_route": "sc-SVC",
        "profile": "application_sc",
        "task": "sc_svc",
        "svc_kind": "sc",
        "strategy": "ScSvcApplicationStrategy",
        "warning": None,
    }
    assert "confounding" not in application
    assert benchmark == {
        "mode": "benchmark",
        "confounding": "segmentation",
        "profile": "benchmark_seg",
        "task": "sp_svc",
        "svc_kind": "sp",
        "strategy": "SpSvcBenchmarkSegStrategy",
        "warning": None,
    }
    assert "application_route" not in benchmark


def test_semantic_router_application_wins_and_invalid_winner_never_falls_back():
    raw = _authority_document()

    selected = resolve_semantic_route(
        raw,
        svc_type="sp-SVC",
        cf="segmentation",
    )

    assert selected["application_route"] == "sp-SVC"
    assert selected["profile"] == "application_sp"
    assert "segmentation" in selected["warning"]
    assert "sp-SVC" in selected["warning"]
    assert "ignored" in selected["warning"]
    with pytest.raises(ConfigError, match="invalid.*svc_type"):
        resolve_semantic_route(raw, svc_type="invalid", cf="segmentation")
    with pytest.raises(ConfigError, match="svc_type.*cf"):
        resolve_semantic_route(raw)


def test_public_run_uses_selectors_and_rejects_route_identity_overrides():
    pipeline = REVISEPipeline()

    with pytest.raises(ConfigError, match="route identity"):
        pipeline.run(
            svc_type="sp-SVC",
            runtime_overrides={"confounding": "segmentation"},
            dry_run=True,
        )


def _raw_config():
    return _authority_document()


def _runtime_for_profile(raw, profile):
    selected_profile = profile or "application_sp"
    for namespace, routes in raw["router"].items():
        for selector, route in routes.items():
            if route["profile"] != selected_profile:
                continue
            key = "application_route" if namespace == "application" else "confounding"
            return {
                "mode": namespace,
                key: selector,
                "task": route["task"],
                "svc_kind": route["svc_kind"],
                "strategy": route["strategy"],
            }
    raise AssertionError(f"No test route uses profile {selected_profile!r}")


def _merge(raw, profile=None, algorithm_overrides=None):
    selected_profile = profile or "application_sp"
    return merge_unified_config(
        raw_config=raw,
        profile=selected_profile,
        runtime_overrides=_runtime_for_profile(raw, selected_profile),
        io_overrides={},
        algorithm_overrides=algorithm_overrides or {},
    )


def test_application_sc_profile_uses_configured_notebook_tacco_defaults():
    merged = _merge(_raw_config(), "application_sc")

    assert merged["ot"]["ga"]["solver"] == "tacco"
    assert merged["ot"]["lr"]["solver"] == "tacco"
    assert merged["sc"]["resolutions"] == [0.6, 0.7, 0.8]
    assert merged["sc"]["tacco_annotate"] == {
        "multi_center": 1,
        "lamb": 0.001,
    }


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("multi_center", True),
        ("multi_center", 0),
        ("multi_center", 1.5),
        ("lamb", True),
        ("lamb", 0),
        ("lamb", math.nan),
    ],
)
def test_application_sc_tacco_annotation_parameters_are_strict(field, value):
    raw = copy.deepcopy(_raw_config())
    raw["profiles"]["application_sc"]["sc"]["tacco_annotate"][field] = value

    with pytest.raises(ConfigError, match=field):
        _merge(raw, "application_sc")


def test_application_sc_tacco_annotation_parameters_are_locked():
    with pytest.raises(ConfigError, match="locked parameter"):
        _merge(
            _raw_config(),
            "application_sc",
            {"sc": {"tacco_annotate": {"lamb": 0.1}}},
        )


def test_structured_algorithm_overrides_merge_only_algorithm_sections():
    merged = merge_unified_config(
        raw_config=_raw_config(),
        profile="application_sp",
        runtime_overrides=_runtime_for_profile(_raw_config(), "application_sp"),
        io_overrides={},
        algorithm_overrides={
            "graph": {"method": "pca", "n_neighbors": 7},
            "ot": {"ga": {"solver": "tacco"}},
        },
    )

    assert merged["graph"]["method"] == "pca"
    assert merged["graph"]["n_neighbors"] == 7
    assert merged["ot"]["ga"]["solver"] == "tacco"
    assert merged["ot"]["lr"]["solver"] == "pot"


@pytest.mark.parametrize("section", ["runtime", "io"])
def test_algorithm_overrides_reject_run_identity_sections(section):
    with pytest.raises(ConfigError, match="algorithm_overrides cannot modify"):
        merge_unified_config(
            raw_config=_raw_config(),
            profile="application_sp",
            runtime_overrides={},
            io_overrides={},
            algorithm_overrides={section: {"seed": 1}},
        )


def test_runtime_none_keeps_the_authority_seed():
    merged = merge_unified_config(
        raw_config=_raw_config(),
        profile="benchmark_seg",
        runtime_overrides={
            **_runtime_for_profile(_raw_config(), "benchmark_seg"),
            "seed": None,
        },
        io_overrides={},
        algorithm_overrides={},
    )

    assert merged["runtime"]["seed"] == 42


def test_non_seed_runtime_none_remains_omitted():
    merged = merge_unified_config(
        raw_config=_raw_config(),
        profile="application_sp",
        runtime_overrides={
            **_runtime_for_profile(_raw_config(), "application_sp"),
            "compatibility_mode": None,
        },
        io_overrides={},
        algorithm_overrides={},
    )

    assert merged["runtime"]["compatibility_mode"] is False


def test_pipeline_accepts_preloaded_application_data():
    from revise.framework import REVISEPipeline
    from reconstruct import run_application

    public_parameters = inspect.signature(REVISEPipeline.run).parameters

    assert "algorithm_overrides" in public_parameters
    assert "set_overrides" not in public_parameters
    assert "st_adata" in public_parameters
    assert list(inspect.signature(run_application).parameters) == ["config_path"]


@pytest.mark.parametrize("phase", ["ga", "lr"])
@pytest.mark.parametrize(
    ("parameter", "value"),
    [("reg", 0.2), ("reg_m", 0.2), ("reg_type", "kl")],
)
def test_algorithm_overrides_cannot_modify_locked_parameters(
    phase,
    parameter,
    value,
):
    path = f"ot.{phase}.pot.{parameter}"
    with pytest.raises(ConfigError, match=f"locked parameter '{path}'"):
        merge_unified_config(
            raw_config=_raw_config(),
            profile="application_sp",
            runtime_overrides={},
            io_overrides={},
            algorithm_overrides={"ot": {phase: {"pot": {parameter: value}}}},
        )


def test_legacy_expose_in_cli_is_rejected_and_cannot_unlock_algorithm_parameters():
    raw = _raw_config()
    raw["locked_params"]["expose_in_cli"] = True

    with pytest.raises(ConfigError, match="Unknown keys in locked_params"):
        _validate_raw_config(raw)

    with pytest.raises(ConfigError, match="locked parameter 'ot.ga.pot.reg'"):
        _merge(
            raw,
            algorithm_overrides={"ot": {"ga": {"pot": {"reg": 0.2}}}},
        )


def test_explicit_ot_method_overrides_conflicting_profile_solvers():
    from revise.application.config import _compile_engine_config

    raw = _raw_config()
    raw["profiles"]["application_sp"]["ot"] = {
        "ga": {"solver": "tacco"},
        "lr": {"solver": "tacco"},
    }
    config = SimpleNamespace(
        svc_type="sp-SVC", ot_method="pot", broad_column="Level1", subtype_column=None,
        select_cell_type=None, local_refinement_strength=None,
        local_refinement_alpha=None, local_refinement_resolutions=None,
        local_refinement_graph_method=None, local_refinement_graph_alpha=None,
        local_refinement_graph_n_neighbors=None, local_refinement_graph_exp_neighbors=None,
        local_refinement_graph_spatial_neighbors=None, local_refinement_match_spot_sum=None,
        seed=None, st_path=Path("st"), reference_path=Path("ref"),
        pm_on_cell_path=None, output_dir=Path("out"), output_name="sample",
        st_format="h5ad", spatialdata_table=None, spatialdata_element=None,
    )
    merged = _merge(raw, "application_sp", _compile_engine_config(config)[2])

    assert merged["ot"]["ga"]["solver"] == "pot"
    assert merged["ot"]["lr"]["solver"] == "pot"


@pytest.fixture
def adapters(monkeypatch):
    """Import production adapters while isolating the unavailable scanpy package."""
    scanpy = types.ModuleType("scanpy")
    scanpy.pp = SimpleNamespace()
    scanpy.pl = SimpleNamespace()
    scanpy.tl = SimpleNamespace()
    monkeypatch.setitem(sys.modules, "scanpy", scanpy)
    sys.modules.pop("revise.backend.adapters", None)
    return importlib.import_module("revise.backend.adapters")


def test_default_ot_schema_is_single_public_surface():
    merged = _merge(_raw_config())

    assert merged["ot"] == {
        "ga": {
            "solver": "pot",
            "pot": {"reg": 0.1, "reg_m": 0.0, "reg_type": "entropy"},
        },
        "lr": {
            "solver": "pot",
            "pot": {"reg": 0.05, "reg_m": 1.0, "reg_type": "kl"},
        },
        "impute": {"reg": 5.0, "reg_m": 0.0, "reg_type": "kl"},
    }
    assert "annotate" not in merged
    assert "local_ot" not in merged
    assert "ot_solver" not in merged["runtime"]


@pytest.mark.parametrize("profile", [
    "application_sp",
    "application_sc",
    "application_sc_sr",
    "benchmark_seg",
    "benchmark_bin2cell",
    "benchmark_sr_batch",
    "benchmark_sr_spot_size",
    "benchmark_impute_panel",
    "benchmark_impute_dropout",
])
@pytest.mark.parametrize(
    ("ga_solver", "lr_solver"),
    [("pot", "pot"), ("pot", "tacco"), ("tacco", "pot"), ("tacco", "tacco")],
)
def test_every_profile_and_ga_lr_combination_reaches_production_mapping(
    adapters, profile, ga_solver, lr_solver
):
    merged = _merge(
        _raw_config(),
        profile,
        {"ot": {"ga": {"solver": ga_solver}, "lr": {"solver": lr_solver}}},
    )

    kwargs = adapters._ot_runner_kwargs(merged)

    assert kwargs == {
        "annotate_mode": ga_solver,
        "annotate_pot_reg": float(merged["ot"]["ga"]["pot"]["reg"]),
        "annotate_pot_reg_m": float(merged["ot"]["ga"]["pot"]["reg_m"]),
        "annotate_pot_reg_type": str(merged["ot"]["ga"]["pot"]["reg_type"]),
        "rec_ot_method": lr_solver,
        "rec_pot_reg": float(merged["ot"]["lr"]["pot"]["reg"]),
        "rec_pot_reg_m": float(merged["ot"]["lr"]["pot"]["reg_m"]),
        "rec_pot_reg_type": str(merged["ot"]["lr"]["pot"]["reg_type"]),
    }


def test_router_profiles_exist():
    raw = _raw_config()
    routed_profiles = {
        route["profile"]
        for routes in raw["router"].values()
        for route in routes.values()
    }

    assert routed_profiles <= raw["profiles"].keys()


def test_impute_reuses_lr_solver_with_independent_numerics(adapters):
    merged = _merge(
        _raw_config(),
        "benchmark_impute_panel",
        {"ot": {"lr": {"solver": "tacco"}, "impute": {"reg": 7.0}}},
    )

    kwargs = adapters._ot_runner_kwargs(merged, impute=True)

    assert kwargs["rec_ot_method"] == "tacco"
    assert kwargs["rec_impute_pot_reg"] == 7.0
    assert kwargs["rec_impute_pot_reg_m"] == 0.0
    assert kwargs["rec_impute_pot_reg_type"] == "kl"
    assert "rec_pot_reg" not in kwargs


@pytest.mark.parametrize(
    ("location", "replacement"),
    [
        ("defaults.annotate", "annotate.mode -> ot.ga.solver"),
        ("profiles.application_sp.local_ot", "local_ot.method -> ot.lr.solver"),
        ("router.benchmark.bin2cell.ot_solver", "ot_solver -> ot.ga.solver + ot.lr.solver"),
        ("defaults.runtime.ot_solver", "ot_solver -> ot.ga.solver + ot.lr.solver"),
        ("defaults.ot.global", "ot.global -> ot.ga.pot"),
        ("profiles.application_sp.ot.local", "ot.local -> ot.lr.pot"),
    ],
)
def test_legacy_raw_profile_and_router_keys_report_replacements(location, replacement):
    raw = _raw_config()
    target = raw
    parts = location.split(".")
    for part in parts[:-1]:
        target = target.setdefault(part, {})
    leaf = parts[-1]
    if leaf == "annotate":
        target[leaf] = {"mode": "pot"}
    elif leaf == "local_ot":
        target[leaf] = {"method": "pot"}
    elif leaf == "ot_solver":
        target[leaf] = "pot"
    elif leaf in {"global", "local"}:
        target[leaf] = {"reg": 0.1, "reg_m": 0.0, "reg_type": "entropy"}

    with pytest.raises(ConfigError, match=replacement.replace("+", r"\+")):
        _validate_raw_config(raw)


def test_runtime_override_legacy_solver_reports_replacement():
    with pytest.raises(ConfigError, match=r"runtime/router ot_solver -> ot.ga.solver \+ ot.lr.solver"):
        merge_unified_config(
            raw_config=_raw_config(),
            profile=None,
            runtime_overrides={"ot_solver": "pot"},
            io_overrides={},
            algorithm_overrides={},
        )


@pytest.mark.parametrize(
    "algorithm_overrides",
    [
        {"ot": {"ga": {"solver": "invalid"}}},
        {"ot": {"lr": {"extra": True}}},
        {"unexpected": {"enabled": True}},
    ],
)
def test_algorithm_overrides_cannot_bypass_resolved_strict_validation(
    algorithm_overrides,
):
    raw = _raw_config()
    raw["locked_params"]["keys"] = []

    with pytest.raises(ConfigError):
        _merge(raw, algorithm_overrides=algorithm_overrides)


def test_locked_leaf_cannot_be_bypassed_by_overriding_its_parent():
    with pytest.raises(ConfigError, match="locked parameter 'ot.ga.pot.reg'"):
        _merge(
            _raw_config(),
            algorithm_overrides={
                "ot": {
                    "ga": {
                        "pot": {"reg": 0.2, "reg_m": 0.0, "reg_type": "entropy"}
                    }
                }
            },
        )


@pytest.mark.parametrize("value", [False, 0, 1, "true", None])
def test_svc_completeness_is_strictly_true(value):
    raw = _raw_config()
    raw["defaults"]["sc"]["svc_completeness"] = value

    with pytest.raises(ConfigError, match="sc.svc_completeness must be exactly true"):
        _merge(raw)


def test_false_completeness_fails_before_input_or_output_path_processing(tmp_path):
    raw = _raw_config()
    raw["defaults"]["sc"]["svc_completeness"] = False
    raw["defaults"]["io"]["output_root"] = str(tmp_path / "must-not-exist")
    raw["defaults"]["io"]["data_root"] = str(tmp_path / "missing-inputs")
    with pytest.raises(ConfigError, match="sc.svc_completeness must be exactly true"):
        _merge(raw)

    assert not (tmp_path / "must-not-exist").exists()


def test_noop_plugin_layer_is_removed():
    import revise.backend as backend
    from revise.backend import registry

    assert not hasattr(backend, "PluginRegistry")
    assert not hasattr(backend, "build_default_plugin_registry")
    assert not hasattr(registry, "PluginRegistry")
    assert not hasattr(registry, "build_default_plugin_registry")
    assert not (ROOT / "revise" / "backend" / "plugins.py").exists()


def test_runtime_and_context_route_have_no_ot_marker(tmp_path):
    merged = _merge(_raw_config())
    ctx = PipelineContext(
        merged_config=merged,
        profile=None,
        runtime=merged["runtime"],
        route_key="application:sp-SVC",
        run_dir=tmp_path,
        logger=logging.getLogger("test_ot_config"),
    )

    assert "ot_solver" not in merged["runtime"]
    assert "ot_solver" not in ctx.route


def _application_request(**overrides):
    values = {
        "svc_type": "sc-SVC",
        "ot_method": None,
        "select_cell_type": "T",
        "broad_column": "Level1",
        "subtype_column": "Level2",
        "local_refinement_strength": None,
    }
    values.update(overrides)
    return SimpleNamespace(**values)


def test_application_yaml_keeps_declared_ot_method():
    from revise.application.config import compile_application_config, load_application_yaml

    source, document = load_application_yaml(
        ROOT / "configs/application/VisiumHD.yaml"
    )
    config = compile_application_config(document, source=source)

    assert config.ot_method == "pot"


def test_structured_config_supports_mixed_ot_solvers():
    merged = _merge(
        _raw_config(),
        "application_sc",
        {"ot": {"ga": {"solver": "tacco"}, "lr": {"solver": "pot"}}},
    )

    assert merged["ot"]["ga"]["solver"] == "tacco"
    assert merged["ot"]["lr"]["solver"] == "pot"


@pytest.mark.parametrize("method", ["pot", "tacco"])
def test_explicit_application_ot_method_overrides_both_phases(method):
    from revise.application.config import _compile_engine_config

    config = SimpleNamespace(
        svc_type="sp-SVC", ot_method=method, broad_column="Level1", subtype_column=None,
        select_cell_type=None, local_refinement_strength=None,
        local_refinement_alpha=None, local_refinement_resolutions=None,
        local_refinement_graph_method=None, local_refinement_graph_alpha=None,
        local_refinement_graph_n_neighbors=None, local_refinement_graph_exp_neighbors=None,
        local_refinement_graph_spatial_neighbors=None, local_refinement_match_spot_sum=None,
        seed=None, st_path=Path("st"), reference_path=Path("ref"),
        pm_on_cell_path=None, output_dir=Path("out"), output_name="sample",
        st_format="h5ad", spatialdata_table=None, spatialdata_element=None,
    )
    overrides = _compile_engine_config(config)[2]

    assert overrides["ot"]["ga"]["solver"] == method
    assert overrides["ot"]["lr"]["solver"] == method


def test_framework_provenance_records_resolved_ot_config_without_events(tmp_path):
    from revise.framework import REVISEPipeline

    merged = _merge(_raw_config())
    svc = SVC(expr=None, spatial=None, svc_kind="sc")
    ctx = SimpleNamespace(
        runner_config=None,
        merged_config=merged,
        profile=None,
        route={"mode": "application", "application_route": "sp-SVC"},
        route_key="application:sp-SVC",
        run_dir=tmp_path,
        stage_trace=[],
        quality_metrics={},
        svc=svc,
        software_versions={},
        engine_defaults_hash=None,
        authority_hash=None,
        algorithm_config_hash=None,
        effective_config_hash=None,
        local_refinement_record={
            "route": "application:sp-SVC",
            "applied": False,
            "strength": 0.2,
        },
    )

    REVISEPipeline()._write_final_metadata(ctx)

    assert svc.provenance["ot_config"] == merged["ot"]
    assert "ot_events" not in svc.provenance
    assert not any("actual" in key or "completed" in key for key in svc.provenance)


def test_application_metadata_cannot_override_canonical_provenance(tmp_path):
    from revise.framework import REVISEPipeline

    merged = _merge(_raw_config())
    svc = SVC(expr=None, spatial=None, svc_kind="sc")
    ctx = SimpleNamespace(
        merged_config=merged,
        profile=None,
        route={"mode": "application", "application_route": "sp-SVC"},
        route_key="application:sp-SVC",
        run_dir=tmp_path,
        run_record={"status": "failed"},
        stage_records=[],
        artifact_records=[],
        quality_metrics={},
        svc=svc,
        software_versions={},
        engine_defaults_hash=None,
        authority_hash=None,
        algorithm_config_hash=None,
        effective_config_hash=None,
        local_refinement_record={},
        application_config_metadata={"run": {"status": "succeeded"}},
    )

    REVISEPipeline()._write_final_metadata(ctx)

    assert svc.provenance["run"] == {"status": "failed"}
    assert svc.provenance["application_config"] == {}


def test_application_sc_sr_config_fields_accept_production_mapping(adapters):
    from dataclasses import fields
    from revise.config.runner_conf import ApplicationScSrConf

    merged = _merge(_raw_config(), "application_sc_sr")
    mapping = adapters._ot_runner_kwargs(merged)
    field_names = {field.name for field in fields(ApplicationScSrConf)}

    assert set(mapping) <= field_names
    conf = ApplicationScSrConf(
        sample_name="sample",
        raw_data_path="data",
        result_root_path="output",
        cell_type_col="Level1",
        confidence_col="Confidence",
        unknown_key="Unknown",
        st_file="st.h5ad",
        sc_ref_file="sc.h5ad",
        rec_graph_n_neighbors=int(merged["graph"]["n_neighbors"]),
        rec_graph_method=str(merged["graph"]["method"]),
        rec_graph_alpha=float(merged["graph"]["alpha"]),
        rec_graph_exp_neighbor_num=int(merged["graph"]["exp_neighbors"]),
        rec_graph_spatial_neighbor_num=int(merged["graph"]["spatial_neighbors"]),
        rec_alpha=float(merged["reconstruct"]["alpha"]),
        rec_match_spot_sum=bool(merged["sc"]["match_spot_sum"]),
        rec_graph_agg_enabled=bool(merged["sc"]["sr_graph_agg_enabled"]),
        svc_completeness=bool(merged["sc"]["svc_completeness"]),
        sr_assignment_seed=int(merged["runtime"]["seed"]),
        local_refinement_strength=float(merged["local_refinement"]["strength"]),
        **mapping,
    )
    assert conf.rec_ot_method == merged["ot"]["lr"]["solver"]


@pytest.mark.parametrize("phase", ["ga", "lr"])
@pytest.mark.parametrize(
    "bad_solver",
    [["pot"], {"name": "pot"}, 1, True, None, "unknown"],
)
def test_raw_solver_must_be_a_supported_string(phase, bad_solver):
    raw = _raw_config()
    raw["defaults"]["ot"][phase]["solver"] = bad_solver

    with pytest.raises(ConfigError, match=rf"defaults\.ot\.{phase}\.solver"):
        _validate_raw_config(raw)


@pytest.mark.parametrize("phase", ["ga", "lr"])
@pytest.mark.parametrize(
    "encoded_solver",
    ["[pot]", "{name: pot}", "1", "true", "null", "unknown"],
)
def test_algorithm_solver_must_be_a_supported_string(phase, encoded_solver):
    with pytest.raises(ConfigError, match=rf"resolved\.ot\.{phase}\.solver"):
        _merge(
            _raw_config(),
            algorithm_overrides={
                "ot": {phase: {"solver": yaml.safe_load(encoded_solver)}}
            },
        )


OT_NUMERIC_SECTIONS = ["ga.pot", "lr.pot", "impute"]


@pytest.mark.parametrize("section", OT_NUMERIC_SECTIONS)
@pytest.mark.parametrize("bad_value", [0, -0.1, True, "0.1", [0.1], {"v": 0.1}, math.nan, math.inf])
def test_ot_reg_must_be_finite_positive_real(section, bad_value):
    raw = _raw_config()
    target = raw["defaults"]["ot"]
    for part in section.split("."):
        target = target[part]
    target["reg"] = bad_value

    with pytest.raises(ConfigError, match=rf"resolved\.ot\.{section}\.reg"):
        _merge(raw)


@pytest.mark.parametrize("section", OT_NUMERIC_SECTIONS)
@pytest.mark.parametrize("bad_value", [-0.1, True, "0.0", [0.0], {"v": 0.0}, math.nan, math.inf])
def test_ot_reg_m_must_be_finite_non_negative_real(section, bad_value):
    raw = _raw_config()
    target = raw["defaults"]["ot"]
    for part in section.split("."):
        target = target[part]
    target["reg_m"] = bad_value

    with pytest.raises(ConfigError, match=rf"resolved\.ot\.{section}\.reg_m"):
        _merge(raw)


@pytest.mark.parametrize("section", OT_NUMERIC_SECTIONS)
@pytest.mark.parametrize("bad_value", [1, True, None, ["kl"], "unsupported"])
def test_ot_reg_type_must_be_supported_string(section, bad_value):
    raw = _raw_config()
    target = raw["defaults"]["ot"]
    for part in section.split("."):
        target = target[part]
    target["reg_type"] = bad_value

    with pytest.raises(ConfigError, match=rf"resolved\.ot\.{section}\.reg_type"):
        _merge(raw)


@pytest.mark.parametrize("profile", [
    "application_sp",
    "application_sc",
    "application_sc_sr",
    "benchmark_seg",
    "benchmark_bin2cell",
    "benchmark_sr_batch",
    "benchmark_sr_spot_size",
    "benchmark_impute_panel",
    "benchmark_impute_dropout",
])
def test_current_profile_ot_numerics_remain_valid(profile):
    _merge(_raw_config(), profile)


@pytest.mark.parametrize(
    ("strategy_name", "profile", "runner_module", "runner_class", "impute"),
    [
        ("SpSvcApplicationStrategy", "application_sp", "sp_svc_application", "SpSVC", False),
        ("ScSvcApplicationStrategy", "application_sc", "sc_svc_application", "ScSVC", False),
        ("ScSvcSrApplicationStrategy", "application_sc_sr", "sc_svc_sr_application", "ScSVCSr", False),
        ("SpSvcBenchmarkSegStrategy", "benchmark_seg", "sp_svc_benchmark", "SpSVC", False),
        ("ScSvcSrBenchmarkStrategy", "benchmark_sr_batch", "sc_svc_sr_benchmark", "ScSVCSr", False),
        ("ScSvcImputeBenchmarkStrategy", "benchmark_impute_panel", "sc_svc_impute_benchmark", "ScSVCImpute", True),
    ],
)
def test_six_strategies_put_ot_mapping_on_actual_runner_config(
    adapters,
    monkeypatch,
    tmp_path,
    strategy_name,
    profile,
    runner_module,
    runner_class,
    impute,
):
    import numpy as np
    import pandas as pd
    from anndata import AnnData

    adapters.sc.pp.filter_cells = lambda *args, **kwargs: None
    adapters.sc.pp.filter_genes = lambda *args, **kwargs: None
    raw = _raw_config()
    raw["locked_params"]["keys"] = []
    merged = _merge(
        raw,
        profile,
        {
            "ot": {
                "ga": {
                    "solver": "tacco",
                    "pot": {"reg": 0.2, "reg_m": 0.3, "reg_type": "kl"},
                },
                "lr": {
                    "solver": "pot",
                    "pot": {"reg": 0.4, "reg_m": 0.5, "reg_type": "entropy"},
                },
                "impute": {"reg": 0.6, "reg_m": 0.7, "reg_type": "entropy"},
            }
        },
    )
    merged["io"]["data_root"] = str(tmp_path)
    merged["io"]["output_root"] = str(tmp_path)
    if profile.startswith("application_"):
        merged["io"]["sample_name"] = "1"

    obs_names = ["cell1", "cell2", "cell3", "cell4"]
    genes = ["g1", "g2", "g3"]
    st = AnnData(
        X=np.ones((4, 3)),
        obs=pd.DataFrame(index=obs_names),
        var=pd.DataFrame(index=genes),
    )
    st.obsm["spatial"] = np.arange(8, dtype=float).reshape(4, 2)
    st.uns["all_cells_in_spot"] = {name: [name] for name in obs_names}
    sc_ref = AnnData(
        X=np.ones((4, 3)),
        obs=pd.DataFrame(
            {
                "Patient": (
                    [1] * 4
                    if profile.startswith("application_")
                    else [merged["io"]["sample_name"]] * 4
                ),
                "Level1": ["A", "A", "B", "B"],
                "Level2": ["A1", "A2", "B1", "B2"],
            },
            index=obs_names,
        ),
        var=pd.DataFrame(index=genes),
    )
    real = st.copy()
    loaded_paths = []

    class InputService:
        def read_st_adata(self, path):
            loaded_paths.append(("st", path))
            return st.copy()

        def read_sc_ref_adata(self, path):
            loaded_paths.append(("sc_ref", path))
            return sc_ref.copy()

        def read_real_adata(self, path):
            loaded_paths.append(("gt", path))
            return real.copy()

    captured = {}

    class DummyRunner:
        def __init__(self, *args):
            captured["runner_config"] = next(
                arg for arg in args if hasattr(arg, "rec_ot_method")
            )

    runner_stub = types.ModuleType(f"revise.backend.runners.{runner_module}")
    setattr(runner_stub, runner_class, DummyRunner)
    monkeypatch.setitem(
        sys.modules,
        f"revise.backend.runners.{runner_module}",
        runner_stub,
    )
    monkeypatch.setattr(adapters, "_input_service", lambda _ctx: InputService())

    ctx = SimpleNamespace(
        merged_config=merged,
        io=merged["io"],
        columns=merged["columns"],
        route_key=(
            f"{merged['runtime']['mode']}:"
            f"{merged['runtime'].get('application_route') or merged['runtime'].get('confounding')}"
        ),
        run_dir=tmp_path,
        logger=logging.getLogger(f"test-{strategy_name}"),
        compatibility_mode=bool(merged["runtime"].get("compatibility_mode", False)),
        runtime=merged["runtime"],
        input_specs=tuple(
            SimpleNamespace(role=role, path=f"/resolved/{role}.h5ad")
            for role in (
                ("st", "sc_ref", "gt")
                if profile.startswith("benchmark_")
                else ("st", "sc_ref")
            )
        ),
        st_adata=st.copy() if profile.startswith("application_") else None,
        sc_ref_adata=sc_ref.copy() if profile.startswith("application_") else None,
    )
    getattr(adapters, strategy_name)().prepare_context(ctx)

    expected = adapters._ot_runner_kwargs(merged, impute=impute)
    assert captured["runner_config"] is ctx.runner_config
    assert {
        key: getattr(ctx.runner_config, key)
        for key in expected
    } == expected
    if profile.startswith("application_"):
        assert loaded_paths == []
    else:
        assert dict(loaded_paths) == {
            spec.role: spec.path for spec in ctx.input_specs
        }
    if profile.startswith("application_"):
        assert ctx.sc_ref_adata.n_obs == 4
    if "local_refinement" in merged:
        assert ctx.runner_config.local_refinement_strength == merged["local_refinement"]["strength"]


def test_active_docs_do_not_advertise_ot_solver_plugins_or_route_markers():
    repo_root = Path(__file__).parents[2]
    active_sources = [
        repo_root / "docs" / "source" / "architecture.rst",
        repo_root / "docs" / "source" / "configuration.rst",
    ]
    text = "\n".join(path.read_text(encoding="utf-8") for path in active_sources)

    assert "OT solver marker" not in text
    assert "OT solver route marker" not in text
    assert "CF/solver" not in text
    assert "confounding strategy, and OT solver" not in text
