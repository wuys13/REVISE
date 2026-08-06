from __future__ import annotations

import math
from pathlib import Path

import pytest

from revise.config import ConfigError, load_raw_config, merge_unified_config


CONFIG_PATH = Path(__file__).parents[2] / "revise" / "revise.yaml"


def _merge(profile: str, algorithm_overrides=None):
    return merge_unified_config(
        raw_config=load_raw_config(CONFIG_PATH),
        profile=profile,
        runtime_overrides={},
        io_overrides={},
        algorithm_overrides=algorithm_overrides or {},
    )


@pytest.mark.parametrize(
    ("profile", "expected"),
    [
        ("application_sp", 0.2),
        ("benchmark_seg", 0.2),
        ("application_sc_sr", 0.0),
        ("benchmark_sr_batch", 0.0),
    ],
)
def test_local_refinement_has_one_route_default(profile, expected):
    assert _merge(profile)["local_refinement"] == {"strength": expected}


@pytest.mark.parametrize("profile", ["application_sc", "benchmark_impute_panel"])
def test_routes_without_local_posterior_conditioning_omit_identity_value(profile):
    assert "local_refinement" not in _merge(profile)


@pytest.mark.parametrize("profile", ["application_sc", "benchmark_impute_panel"])
def test_routes_without_local_posterior_conditioning_reject_explicit_strength(profile):
    with pytest.raises(ConfigError, match="does not accept local_refinement.strength"):
        _merge(profile, {"local_refinement": {"strength": 0.0}})


@pytest.mark.parametrize("strength", [-0.1, math.inf, math.nan, True, "0.2"])
def test_local_refinement_strength_must_be_finite_non_negative_real(strength):
    with pytest.raises(ConfigError, match="local_refinement.strength"):
        _merge("application_sp", {"local_refinement": {"strength": strength}})


@pytest.mark.parametrize(
    "removed",
    [
        {"local_refinement": {"guidance": "prefer"}},
        {"local_refinement": {"compatibility": {"mode": "cost"}}},
        {"posterior_conditioning": {"mode": "cost"}},
    ],
)
def test_removed_assignment_guidance_yaml_uses_one_migration_error(removed):
    with pytest.raises(
        ConfigError,
        match=r"Assignment guidance options were removed; use local_refinement\.strength",
    ):
        _merge("application_sp", removed)
