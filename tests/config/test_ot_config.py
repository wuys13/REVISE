from __future__ import annotations

import argparse
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

from revise.config import ConfigError, load_raw_config, merge_unified_config
from revise.recon.context import PipelineContext
from revise.svc import SVC


CONFIG_PATH = Path(__file__).parents[2] / "revise" / "revise.yaml"


def _raw_config():
    return load_raw_config(CONFIG_PATH)


def _merge(raw, profile=None, algorithm_overrides=None):
    return merge_unified_config(
        raw_config=raw,
        profile=profile,
        runtime_overrides={},
        io_overrides={},
        algorithm_overrides=algorithm_overrides or {},
    )


def _write_config(tmp_path, raw):
    path = tmp_path / "revise.yaml"
    path.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")
    return path


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
        runtime_overrides={},
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


def test_runtime_none_is_an_explicit_override():
    merged = merge_unified_config(
        raw_config=_raw_config(),
        profile="benchmark_seg",
        runtime_overrides={"seed": None},
        io_overrides={},
        algorithm_overrides={},
    )

    assert merged["runtime"]["seed"] is None


def test_non_seed_runtime_none_remains_omitted():
    merged = merge_unified_config(
        raw_config=_raw_config(),
        profile="application_sp",
        runtime_overrides={"platform": None},
        io_overrides={},
        algorithm_overrides={},
    )

    assert merged["runtime"]["platform"] == "sp_svc"


def test_pipeline_public_api_does_not_expose_algorithm_overrides():
    from revise.framework import REVISEPipeline
    from revise.recon.facade import sc_svc, sp_svc

    public_parameters = inspect.signature(REVISEPipeline.run).parameters
    internal_parameters = inspect.signature(
        REVISEPipeline._run_with_algorithm_overrides
    ).parameters

    assert "algorithm_overrides" not in public_parameters
    assert "set_overrides" not in public_parameters
    assert "algorithm_overrides" in internal_parameters
    assert "algorithm_overrides" not in inspect.signature(sp_svc).parameters
    assert "algorithm_overrides" not in inspect.signature(sc_svc).parameters


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


def test_legacy_expose_in_cli_is_rejected_and_cannot_unlock_algorithm_parameters(
    tmp_path,
):
    raw = _raw_config()
    raw["locked_params"]["expose_in_cli"] = True

    with pytest.raises(ConfigError, match="Unknown keys in locked_params"):
        load_raw_config(_write_config(tmp_path, raw))

    with pytest.raises(ConfigError, match="locked parameter 'ot.ga.pot.reg'"):
        _merge(
            raw,
            algorithm_overrides={"ot": {"ga": {"pot": {"reg": 0.2}}}},
        )


def test_algorithm_overrides_cannot_change_strategy_through_hyperresolution():
    with pytest.raises(ConfigError, match="run identity.*sc.hyperresolution"):
        _merge(
            _raw_config(),
            "application_sc",
            {
                "sc": {
                    "hyperresolution": {
                        "enabled": True,
                        "strategy": "InjectedStrategy",
                    }
                }
            },
        )


def test_explicit_ot_method_overrides_conflicting_profile_solvers():
    from revise.application import service

    raw = _raw_config()
    raw["profiles"]["application_sp"]["ot"] = {
        "ga": {"solver": "tacco"},
        "lr": {"solver": "tacco"},
    }
    merged = _merge(
        raw,
        "application_sp",
        service._build_algorithm_overrides(_cli_args(ot_method="pot")),
    )

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
    "application_sc_hyper",
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


def test_profiles_cover_the_flattened_router_route_set():
    raw = _raw_config()
    routed_profiles = {
        (_merge(raw, profile)["runtime"]["platform"], _merge(raw, profile)["runtime"]["confounding"])
        for profile in raw["profiles"]
    }
    router_routes = {
        (platform, confounding)
        for platform, confoundings in raw["router"].items()
        for confounding in confoundings
    }

    assert routed_profiles == router_routes


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


def test_adapter_projects_only_local_refinement_strength(adapters):
    merged = _merge(
        _raw_config(),
        "application_sc_sr",
        {
            "local_refinement": {"strength": 3.0},
        },
    )
    conf = SimpleNamespace()

    adapters._attach_local_refinement_strength(conf, merged)

    assert vars(conf) == {"local_refinement_strength": 3.0}


@pytest.mark.parametrize(
    ("location", "replacement"),
    [
        ("defaults.annotate", "annotate.mode -> ot.ga.solver"),
        ("profiles.application_sp.local_ot", "local_ot.method -> ot.lr.solver"),
        ("router.sp_svc.bin2cell.ot_solver", "ot_solver -> ot.ga.solver + ot.lr.solver"),
        ("defaults.runtime.ot_solver", "ot_solver -> ot.ga.solver + ot.lr.solver"),
        ("defaults.ot.global", "ot.global -> ot.ga.pot"),
        ("profiles.application_sp.ot.local", "ot.local -> ot.lr.pot"),
    ],
)
def test_legacy_raw_profile_and_router_keys_report_replacements(tmp_path, location, replacement):
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
        load_raw_config(_write_config(tmp_path, raw))


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
    from revise.framework import REVISEPipeline

    raw = _raw_config()
    raw["defaults"]["sc"]["svc_completeness"] = False
    raw["defaults"]["io"]["output_root"] = str(tmp_path / "must-not-exist")
    raw["defaults"]["io"]["data_root"] = str(tmp_path / "missing-inputs")
    pipeline = REVISEPipeline(str(_write_config(tmp_path, raw)))

    with pytest.raises(ConfigError, match="sc.svc_completeness must be exactly true"):
        pipeline.run()

    assert not (tmp_path / "must-not-exist").exists()


def test_noop_plugin_layer_is_removed_and_hyper_profile_selects_strategy_directly():
    import revise.backend as backend
    from revise.backend import registry

    merged = _merge(_raw_config(), "application_sc_hyper")

    assert merged["runtime"]["strategy"] == "ScSvcHyperApplicationStrategy"
    assert not hasattr(backend, "PluginRegistry")
    assert not hasattr(backend, "build_default_plugin_registry")
    assert not hasattr(registry, "PluginRegistry")
    assert not hasattr(registry, "build_default_plugin_registry")
    assert not (CONFIG_PATH.parent / "backend" / "plugins.py").exists()


def test_runtime_and_context_route_have_no_ot_marker(tmp_path):
    merged = _merge(_raw_config())
    ctx = PipelineContext(
        merged_config=merged,
        raw_config=_raw_config(),
        config_path=str(CONFIG_PATH),
        profile=None,
        runtime=merged["runtime"],
        route_key="sp_svc:bin2cell",
        run_dir=tmp_path,
        logger=logging.getLogger("test_ot_config"),
    )

    assert "ot_solver" not in merged["runtime"]
    assert "ot_solver" not in ctx.route


def _cli_args(**overrides):
    values = {
        "svc_type": "sc-SVC",
        "ot_method": None,
        "select_ct": "T",
        "cell_type_col": "Level1",
        "sub_cell_type_col": "Level2",
    }
    values.update(overrides)
    return argparse.Namespace(**values)


def test_cli_omitted_ot_method_parses_as_none(monkeypatch):
    from revise.application import cli

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "reconstruct.py",
            "--svc-type", "sc-SVC",
            "--sample-name", "sample",
            "--st-file", "st.h5ad",
            "--sc-ref-file", "sc.h5ad",
            "--data-root", "data",
            "--select-ct", "T",
        ],
    )

    assert cli.parse_args().ot_method is None


def test_structured_config_supports_mixed_ot_solvers():
    merged = _merge(
        _raw_config(),
        "application_sc",
        {"ot": {"ga": {"solver": "tacco"}, "lr": {"solver": "pot"}}},
    )

    assert merged["ot"]["ga"]["solver"] == "tacco"
    assert merged["ot"]["lr"]["solver"] == "pot"


@pytest.mark.parametrize("method", ["pot", "tacco"])
def test_explicit_cli_ot_flag_overrides_both_phases(method):
    from revise.application import service

    overrides = service._build_algorithm_overrides(_cli_args(ot_method=method))

    assert overrides["ot"]["ga"]["solver"] == method
    assert overrides["ot"]["lr"]["solver"] == method


def test_framework_provenance_records_resolved_ot_config_without_events(tmp_path):
    from revise.framework import REVISEPipeline

    merged = _merge(_raw_config())
    svc = SVC(expr=None, spatial=None, svc_kind="sc")
    ctx = SimpleNamespace(
        runner_config=None,
        merged_config=merged,
        config_path=str(CONFIG_PATH),
        profile=None,
        route={"platform": "sp_svc"},
        route_key="sp_svc:bin2cell",
        run_dir=tmp_path,
        stage_trace=[],
        quality_metrics={},
        svc=svc,
        software_versions={},
        local_refinement_record={
            "route": "sp_svc:bin2cell",
            "applied": False,
            "strength": 0.2,
        },
    )

    REVISEPipeline()._write_final_metadata(ctx)

    assert svc.provenance["ot_config"] == merged["ot"]
    assert "ot_events" not in svc.provenance
    assert not any("actual" in key or "completed" in key for key in svc.provenance)


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
        **mapping,
    )
    assert conf.rec_ot_method == merged["ot"]["lr"]["solver"]


@pytest.mark.parametrize("phase", ["ga", "lr"])
@pytest.mark.parametrize(
    "bad_solver",
    [["pot"], {"name": "pot"}, 1, True, None, "unknown"],
)
def test_raw_solver_must_be_a_supported_string(tmp_path, phase, bad_solver):
    raw = _raw_config()
    raw["defaults"]["ot"][phase]["solver"] = bad_solver

    with pytest.raises(ConfigError, match=rf"defaults\.ot\.{phase}\.solver"):
        load_raw_config(_write_config(tmp_path, raw))


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
    "application_sc_hyper",
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
        route_key=f"{merged['runtime']['platform']}:{merged['runtime']['confounding']}",
        run_dir=tmp_path,
        logger=logging.getLogger(f"test-{strategy_name}"),
        compatibility_mode=bool(merged["runtime"].get("compatibility_mode", False)),
        input_specs=tuple(
            SimpleNamespace(role=role, path=f"/resolved/{role}.h5ad")
            for role in (
                ("st", "sc_ref", "gt")
                if profile.startswith("benchmark_")
                else ("st", "sc_ref")
            )
        ),
    )
    getattr(adapters, strategy_name)().prepare_context(ctx)

    expected = adapters._ot_runner_kwargs(merged, impute=impute)
    assert captured["runner_config"] is ctx.runner_config
    assert {
        key: getattr(ctx.runner_config, key)
        for key in expected
    } == expected
    assert dict(loaded_paths) == {
        spec.role: spec.path for spec in ctx.input_specs
    }
    if profile.startswith("application_"):
        assert ctx.sc_ref_adata.n_obs == 4


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
