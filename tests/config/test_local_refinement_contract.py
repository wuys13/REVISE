from __future__ import annotations

import math

import pytest

from revise.config import ConfigError, merge_unified_config, resolve_semantic_route
from revise.config.authority import _authority_document


def _merge(profile: str, algorithm_overrides=None):
    raw = _authority_document()
    route = next(
        (namespace, selector)
        for namespace, routes in raw["router"].items()
        for selector, spec in routes.items()
        if spec["profile"] == profile
    )
    selector = (
        {"svc_type": route[1]}
        if route[0] == "application"
        else {"cf": route[1]}
    )
    runtime = resolve_semantic_route(raw, **selector)
    runtime.pop("profile")
    runtime.pop("warning")
    return merge_unified_config(
        raw_config=raw,
        profile=profile,
        runtime_overrides=runtime,
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
