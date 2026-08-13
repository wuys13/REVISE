from __future__ import annotations

import inspect
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]


def test_typed_authority_owns_defaults_routes_and_locked_keys():
    from revise.config.authority import ENGINE_DEFAULTS, LOCKED_KEYS, ROUTES, RouteSpec

    assert ENGINE_DEFAULTS["runtime"]["seed"] == 42
    assert set(ROUTES["application"]) == {"sp-SVC", "sc-SVC", "sc-SVC-sr"}
    assert set(ROUTES["benchmark"]) == {
        "segmentation",
        "bin2cell",
        "batch_effect",
        "spot_size",
        "gene_panel",
        "gene_dropout",
    }
    assert all(
        isinstance(route, RouteSpec)
        for namespace in ROUTES.values()
        for route in namespace.values()
    )
    assert ROUTES["application"]["sp-SVC"].overrides["reconstruct"]["alpha"] == 0.5
    assert ROUTES["application"]["sc-SVC"].overrides["reconstruct"]["alpha"] == 0.5
    assert ROUTES["benchmark"]["segmentation"].overrides["local_refinement"] == {"strength": 0.2}
    assert ROUTES["benchmark"]["batch_effect"].overrides["local_refinement"] == {"strength": 0.0}
    assert "ot.ga.pot.reg" in LOCKED_KEYS


def test_pipeline_has_no_external_engine_config_path():
    from revise.framework import REVISEPipeline

    assert list(inspect.signature(REVISEPipeline).parameters) == []
    with pytest.raises(TypeError):
        REVISEPipeline(config_path="revise/revise.yaml")


def test_segmentation_evaluation_policy_is_engine_authority():
    from revise.config.authority import ENGINE_DEFAULTS

    assert ENGINE_DEFAULTS["benchmark"] == {
        "evaluate": False,
        "dropout_total_counts": 60,
        "swapping_total_counts": 300,
        "lower_ts": 0.2,
        "upper_ts": 0.8,
    }


@pytest.mark.parametrize("enabled", [True, False])
def test_evaluation_policy_reads_the_resolved_benchmark_flag(enabled):
    from types import SimpleNamespace

    from revise.backend.policies import ModeEvaluationPolicy

    ctx = SimpleNamespace(
        runtime={"mode": "benchmark"},
        merged_config={"benchmark": {"evaluate": enabled}},
    )

    assert ModeEvaluationPolicy().should_evaluate(ctx) is enabled


def test_sc_cluster_random_state_preserves_the_legacy_algorithm_value():
    from revise.config.authority import ENGINE_DEFAULTS

    assert ENGINE_DEFAULTS["graph"]["random_state"] == 0


def test_legacy_engine_yaml_and_public_loader_are_absent():
    import revise.config as config

    assert not (ROOT / "revise" / "revise.yaml").exists()
    assert not hasattr(config, "load_raw_config")


def test_adapters_do_not_supply_hidden_config_defaults():
    source = (ROOT / "revise" / "backend" / "adapters.py").read_text(
        encoding="utf-8"
    )

    assert "def _cfg_get" not in source
    assert "runtime.get(\"seed\", 42)" not in source
    assert "columns.get(\"cell_type_col\", \"Level1\")" not in source
    assert "sc_cfg.get(\"resolutions\", [0.6, 0.7, 0.8])" not in source
