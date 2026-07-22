from __future__ import annotations

import argparse
import importlib
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


CONFIG_PATH = Path(__file__).parents[1] / "revise" / "revise.yaml"


def _raw_config():
    return load_raw_config(CONFIG_PATH)


def _merge(raw, profile=None, set_overrides=()):
    return merge_unified_config(
        raw_config=raw,
        profile=profile,
        runtime_overrides={},
        io_overrides={},
        set_overrides=set_overrides,
    )


def _write_config(tmp_path, raw):
    path = tmp_path / "revise.yaml"
    path.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")
    return path


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
        [f"ot.ga.solver={ga_solver}", f"ot.lr.solver={lr_solver}"],
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
        ["ot.lr.solver=tacco", "ot.impute.reg=7.0"],
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


@pytest.mark.parametrize(
    ("override", "replacement"),
    [
        ("annotate.mode=pot", "annotate.mode -> ot.ga.solver"),
        ("local_ot.method=pot", "local_ot.method -> ot.lr.solver"),
        ("runtime.ot_solver=pot", "runtime/router ot_solver -> ot.ga.solver + ot.lr.solver"),
        ("ot.global={reg: 0.1}", "ot.global -> ot.ga.pot"),
        ("ot.local={reg: 0.1}", "ot.local -> ot.lr.pot"),
    ],
)
def test_legacy_set_keys_report_replacements_instead_of_key_errors(override, replacement):
    with pytest.raises(ConfigError, match=replacement.replace("+", r"\+")):
        _merge(_raw_config(), set_overrides=[override])


def test_runtime_override_legacy_solver_reports_replacement():
    with pytest.raises(ConfigError, match=r"runtime/router ot_solver -> ot.ga.solver \+ ot.lr.solver"):
        merge_unified_config(
            raw_config=_raw_config(),
            profile=None,
            runtime_overrides={"ot_solver": "pot"},
            io_overrides={},
            set_overrides=[],
        )


@pytest.mark.parametrize(
    "override",
    [
        "ot={ga: {solver: pot, pot: {reg: 0.1, reg_m: 0.0, reg_type: entropy}}}",
        "ot.ga={solver: invalid, pot: {reg: 0.1, reg_m: 0.0, reg_type: entropy}}",
        "ot.lr={solver: pot, pot: {reg: 0.1, reg_m: 0.0, reg_type: kl}, extra: true}",
    ],
)
def test_whole_parent_overrides_cannot_bypass_resolved_strict_validation(override):
    raw = _raw_config()
    raw["locked_params"]["expose_in_cli"] = True

    with pytest.raises(ConfigError):
        _merge(raw, set_overrides=[override])


def test_locked_leaf_cannot_be_bypassed_by_overriding_its_parent():
    with pytest.raises(ConfigError, match="locked parameter 'ot.ga.pot'"):
        _merge(
            _raw_config(),
            set_overrides=["ot.ga.pot={reg: 0.2, reg_m: 0.0, reg_type: entropy}"],
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


def test_plugin_registry_keeps_platform_and_cf_but_has_no_ot_api():
    from revise.backend.registry import build_default_plugin_registry

    registry = build_default_plugin_registry()
    payload = {
        "runtime": {"platform": "sp_svc", "confounding": "bin2cell"},
        "merged_config": {},
    }

    assert registry.get_platform_adapter("default").adapt(payload) is payload
    assert registry.get_cf_strategy("bin2cell").apply(payload) is payload
    assert set(registry._platform_adapters) == {"default", "sim2real"}
    assert not hasattr(registry, "register_ot_solver")
    assert not hasattr(registry, "get_ot_solver")
    assert not hasattr(registry, "_ot_solvers")


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
        "select_ct": "all",
        "cell_type_col": "Level1",
        "sub_cell_type_col": "Level2",
        "set_overrides": [],
    }
    values.update(overrides)
    return argparse.Namespace(**values)


def test_cli_omitted_ot_method_parses_as_none(monkeypatch):
    from revise.application import cli

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "application_reconstruct.py",
            "--svc-type", "sc-SVC",
            "--sample-name", "sample",
            "--st-file", "st.h5ad",
            "--sc-ref-file", "sc.h5ad",
            "--data-root", "data",
        ],
    )

    assert cli.parse_args().ot_method is None


def test_cli_without_ot_flag_preserves_mixed_set_solvers():
    from revise.application import service

    overrides = service._build_set_overrides(
        _cli_args(set_overrides=["ot.ga.solver=tacco", "ot.lr.solver=pot"])
    )

    assert "ot.ga.solver=tacco" in overrides
    assert "ot.lr.solver=pot" in overrides
    assert overrides.count("ot.ga.solver=tacco") == 1
    assert overrides.count("ot.lr.solver=pot") == 1


@pytest.mark.parametrize("method", ["pot", "tacco"])
def test_explicit_cli_ot_flag_overrides_both_phases(method):
    from revise.application import service

    overrides = service._build_set_overrides(_cli_args(ot_method=method))

    assert f"ot.ga.solver={method}" in overrides
    assert f"ot.lr.solver={method}" in overrides


@pytest.mark.parametrize(
    "override",
    [
        "ot.ga.solver=pot",
        "ot.lr.solver=pot",
        "ot.ga={solver: pot}",
        "ot={ga: {solver: pot}, lr: {solver: pot}}",
    ],
)
def test_explicit_cli_ot_flag_conflicts_with_overlapping_set(override):
    from revise.application import service

    with pytest.raises(ValueError, match="Conflicting high-level CLI option"):
        service._build_set_overrides(
            _cli_args(ot_method="tacco", set_overrides=[override])
        )


def test_framework_provenance_records_resolved_ot_config_and_events(tmp_path):
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
        ot_events=[
            {"phase": "ga", "solver": "pot", "status": "requested", "call": 0},
            {"phase": "lr", "solver": "pot", "status": "requested", "call": 0},
        ],
    )

    REVISEPipeline()._write_final_metadata(ctx)

    assert svc.provenance["ot_config"] == merged["ot"]
    assert svc.provenance["ot_events"] == ctx.ot_events
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
def test_set_solver_must_be_a_supported_string(phase, encoded_solver):
    with pytest.raises(ConfigError, match=rf"resolved\.ot\.{phase}\.solver"):
        _merge(
            _raw_config(),
            set_overrides=[f"ot.{phase}.solver={encoded_solver}"],
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
    raw["locked_params"]["expose_in_cli"] = True
    merged = _merge(
        raw,
        profile,
        [
            "ot.ga.solver=tacco",
            "ot.ga.pot.reg=0.2",
            "ot.ga.pot.reg_m=0.3",
            "ot.ga.pot.reg_type=kl",
            "ot.lr.solver=pot",
            "ot.lr.pot.reg=0.4",
            "ot.lr.pot.reg_m=0.5",
            "ot.lr.pot.reg_type=entropy",
            "ot.impute.reg=0.6",
            "ot.impute.reg_m=0.7",
            "ot.impute.reg_type=entropy",
        ],
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
    repo_root = Path(__file__).parents[1]
    active_sources = [
        repo_root / "docs" / "source" / "architecture.rst",
        repo_root / "docs" / "source" / "configuration.rst",
        repo_root / "revise" / "backend" / "plugins.py",
    ]
    text = "\n".join(path.read_text(encoding="utf-8") for path in active_sources)

    assert "OT solver marker" not in text
    assert "OT solver route marker" not in text
    assert "CF/solver" not in text
    assert "confounding strategy, and OT solver" not in text
