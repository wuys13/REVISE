from __future__ import annotations

import copy
import math
from numbers import Real
from pathlib import Path
from typing import Any, Dict, List, TypedDict

try:
    import yaml
except ImportError as exc:  # pragma: no cover
    raise ImportError("PyYAML is required for revise config loading") from exc


APPLICATION_DEFAULT_CF = {
    "sp_svc": "bin2cell",
    "sc_svc": "segmentation",
    "sc_svc_sr": "spot_size",
}

TOP_LEVEL_KEYS = {"version", "defaults", "router", "profiles", "locked_params", "schemas"}
DEFAULT_SECTION_KEYS = {
    "runtime",
    "io",
    "columns",
    "preprocess",
    "graph",
    "ot",
    "plot",
    "reconstruct",
    "benchmark",
    "sc",
    "impute",
    "local_refinement",
    "posterior_conditioning",
}
RUNTIME_KEYS = {
    "seed",
    "deterministic",
    "compatibility_mode",
    "platform",
    "confounding",
    "mode",
    "task",
    "svc_kind",
    "strategy",
}
IO_KEYS = {
    "data_root",
    "output_root",
    "sample_name",
    "st_file",
    "sc_ref_file",
    "gt_svc_file",
    "seg_method",
    "spot_size",
    "patient_key",
    "sample_size",
    "save_outputs",
    "input_format",
    "spatialdata_path",
    "spatialdata_reader",
    "spatialdata_table",
    "spatialdata_spatial_element",
    "spatialdata_coordinate_system",
}
COLUMNS_KEYS = {"cell_type_col", "sub_cell_type_col", "confidence_col", "unknown_key"}
PREPROCESS_KEYS = {"st_min_counts", "st_min_cells", "sc_min_counts", "sc_min_cells", "st_min_transcripts"}
GRAPH_KEYS = {"method", "alpha", "n_neighbors", "exp_neighbors", "spatial_neighbors"}
OT_KEYS = {"ga", "lr", "impute"}
OT_PHASE_KEYS = {"solver", "pot"}
OT_LEAF_KEYS = {"reg", "reg_m", "reg_type"}
OT_SOLVERS = {"pot", "tacco"}
OT_REG_TYPES = {"entropy", "kl"}
PLOT_KEYS = {"enabled", "cluster_resolutions", "min_genes", "min_cells", "sample_size"}
RECONSTRUCT_KEYS = {"alpha"}
BENCHMARK_KEYS = {"evaluate"}
SC_KEYS = {
    "select_ct",
    "resolutions",
    "select_resolution",
    "hyperresolution",
    "match_spot_sum",
    "svc_completeness",
    "sr_graph_agg_enabled",
    "sr_graph_agg_low_conf_only",
    "sr_graph_agg_low_conf_quantile",
    "sr_graph_agg_anchor_only",
    "sr_graph_agg_anchor_high_conf_quantile",
    "sr_graph_agg_confidence_mode",
    "sr_graph_agg_conf_weighted_alpha",
    "sr_graph_agg_conf_alpha_min",
    "sr_graph_agg_conf_alpha_max",
    "sr_graph_agg_conf_alpha_power",
    "sr_noise_enabled",
    "sr_noise_lambda",
    "sr_noise_k",
    "sr_noise_weight",
    "sr_noise_preserve_total_counts",
    "sr_noise_seed",
    "tacco_annotate",
}
SC_HYPER_KEYS = {"enabled", "strategy", "resolutions", "select_resolution"}
SC_TACCO_ANNOTATE_KEYS = {"multi_center", "lamb"}
IMPUTE_KEYS = {
    "merge_subcluster_method",
    "subcluster_resolution",
    "in_panel_subcluster_resolution",
    "prune",
    "n_neighbors",
    "method",
}
POSTERIOR_CONDITIONING_KEYS = {
    "enabled",
    "mode",
    "posterior_key",
    "beta",
    "min_affinity",
    "cost_strength",
    "strict",
}
LOCAL_REFINEMENT_KEYS = {"guidance", "compatibility"}
LOCAL_REFINEMENT_COMPATIBILITY_KEYS = {
    "mode",
    "beta",
    "min_affinity",
    "strength",
}
GUIDANCE_POLICIES = {"off", "prefer", "require"}
GUIDANCE_COMPATIBILITY_MODES = {"cost", "reference"}
ROUTE_LEAF_KEYS = {"mode", "task", "svc_kind", "strategy"}
LOCKED_PARAMS_KEYS = {"keys"}
ALGORITHM_IDENTITY_PATHS = {"sc.hyperresolution"}


class ConfigError(ValueError):
    """Unified configuration error."""


class AssignmentGuidanceRequestRecord(TypedDict):
    configured_guidance: str | None
    configured_compatibility_mode: str | None
    resolution_source: str
    deprecations: list[str]


class AssignmentGuidanceRequestEvidence(TypedDict, total=False):
    assignment_guidance: AssignmentGuidanceRequestRecord


class ResolvedConfig(dict[str, Any]):
    """Algorithm config with separately copied, non-canonical request evidence.

    Coercing this object with ``dict(config)`` intentionally keeps only the
    algorithm mapping; use ``copy()``/``copy.copy``/``copy.deepcopy`` when the
    request evidence must travel with it.
    """

    request_evidence: AssignmentGuidanceRequestEvidence

    def __init__(
        self,
        *args: Any,
        request_evidence: AssignmentGuidanceRequestEvidence | None = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(*args, **kwargs)
        self.request_evidence = copy.deepcopy(request_evidence or {})

    def copy(self) -> "ResolvedConfig":
        return type(self)(self, request_evidence=self.request_evidence)

    def __copy__(self) -> "ResolvedConfig":
        return self.copy()

    def __deepcopy__(self, memo: dict[int, Any]) -> "ResolvedConfig":
        copied = type(self)()
        memo[id(self)] = copied
        copied.update(copy.deepcopy(dict(self), memo))
        copied.request_evidence = copy.deepcopy(self.request_evidence, memo)
        return copied


def _ensure_mapping(value: Any, ctx: str) -> Dict[str, Any]:
    if not isinstance(value, dict):
        raise ConfigError(f"{ctx} must be a mapping")
    return value


def _reject_unknown_keys(mapping: Dict[str, Any], allowed: set[str], ctx: str) -> None:
    unknown = sorted(set(mapping.keys()) - allowed)
    if unknown:
        raise ConfigError(f"Unknown keys in {ctx}: {unknown}")


def _validate_solver(value: Any, ctx: str) -> None:
    if not isinstance(value, str) or value not in OT_SOLVERS:
        raise ConfigError(
            f"{ctx} must be a string in {sorted(OT_SOLVERS)}; got {value!r}"
        )


def _validate_pot_values(pot: Dict[str, Any], ctx: str) -> None:
    reg = pot["reg"]
    if isinstance(reg, bool) or not isinstance(reg, Real) or not math.isfinite(reg) or reg <= 0:
        raise ConfigError(f"{ctx}.reg must be a finite real number greater than 0; got {reg!r}")

    reg_m = pot["reg_m"]
    if (
        isinstance(reg_m, bool)
        or not isinstance(reg_m, Real)
        or not math.isfinite(reg_m)
        or reg_m < 0
    ):
        raise ConfigError(
            f"{ctx}.reg_m must be a finite real number greater than or equal to 0; "
            f"got {reg_m!r}"
        )

    reg_type = pot["reg_type"]
    if not isinstance(reg_type, str) or reg_type not in OT_REG_TYPES:
        raise ConfigError(
            f"{ctx}.reg_type must be a string in {sorted(OT_REG_TYPES)}; got {reg_type!r}"
        )


def _validate_ot_section(ot_cfg: Dict[str, Any], ctx: str, *, resolved: bool = False) -> None:
    if "global" in ot_cfg:
        raise ConfigError(f"{ctx}.global is no longer supported; replace ot.global -> ot.ga.pot")
    if "local" in ot_cfg:
        raise ConfigError(f"{ctx}.local is no longer supported; replace ot.local -> ot.lr.pot")
    _reject_unknown_keys(ot_cfg, OT_KEYS, ctx)
    for name in ("ga", "lr"):
        if name not in ot_cfg:
            if resolved:
                raise ConfigError(f"Missing required OT section: {ctx}.{name}")
            continue
        phase = _ensure_mapping(ot_cfg[name], f"{ctx}.{name}")
        _reject_unknown_keys(phase, OT_PHASE_KEYS, f"{ctx}.{name}")
        if resolved and set(phase) != OT_PHASE_KEYS:
            missing = sorted(OT_PHASE_KEYS - set(phase))
            raise ConfigError(f"Missing required OT keys in {ctx}.{name}: {missing}")
        if "solver" in phase:
            _validate_solver(phase["solver"], f"{ctx}.{name}.solver")
        if "pot" in phase:
            pot = _ensure_mapping(phase["pot"], f"{ctx}.{name}.pot")
            _reject_unknown_keys(pot, OT_LEAF_KEYS, f"{ctx}.{name}.pot")
            if resolved and set(pot) != OT_LEAF_KEYS:
                missing = sorted(OT_LEAF_KEYS - set(pot))
                raise ConfigError(f"Missing required OT keys in {ctx}.{name}.pot: {missing}")
            if resolved:
                _validate_pot_values(pot, f"{ctx}.{name}.pot")

    if "impute" not in ot_cfg:
        if resolved:
            raise ConfigError(f"Missing required OT section: {ctx}.impute")
    else:
        impute = _ensure_mapping(ot_cfg["impute"], f"{ctx}.impute")
        _reject_unknown_keys(impute, OT_LEAF_KEYS, f"{ctx}.impute")
        if resolved and set(impute) != OT_LEAF_KEYS:
            missing = sorted(OT_LEAF_KEYS - set(impute))
            raise ConfigError(f"Missing required OT keys in {ctx}.impute: {missing}")
        if resolved:
            _validate_pot_values(impute, f"{ctx}.impute")


def _validate_sc_section(sc_cfg: Dict[str, Any], ctx: str) -> None:
    _reject_unknown_keys(sc_cfg, SC_KEYS, ctx)
    hyper = sc_cfg.get("hyperresolution")
    if hyper is not None:
        hyper_map = _ensure_mapping(hyper, f"{ctx}.hyperresolution")
        _reject_unknown_keys(hyper_map, SC_HYPER_KEYS, f"{ctx}.hyperresolution")

    tacco_annotate = sc_cfg.get("tacco_annotate")
    if tacco_annotate is not None:
        tacco_map = _ensure_mapping(tacco_annotate, f"{ctx}.tacco_annotate")
        _reject_unknown_keys(
            tacco_map,
            SC_TACCO_ANNOTATE_KEYS,
            f"{ctx}.tacco_annotate",
        )
        missing = sorted(SC_TACCO_ANNOTATE_KEYS - set(tacco_map))
        if missing:
            raise ConfigError(
                f"Missing required TACCO annotation keys in "
                f"{ctx}.tacco_annotate: {missing}"
            )
        multi_center = tacco_map["multi_center"]
        if (
            isinstance(multi_center, bool)
            or not isinstance(multi_center, int)
            or multi_center <= 0
        ):
            raise ConfigError(
                f"{ctx}.tacco_annotate.multi_center must be a positive integer; "
                f"got {multi_center!r}"
            )
        lamb = tacco_map["lamb"]
        if (
            isinstance(lamb, bool)
            or not isinstance(lamb, Real)
            or not math.isfinite(lamb)
            or lamb <= 0
        ):
            raise ConfigError(
                f"{ctx}.tacco_annotate.lamb must be a finite real number "
                f"greater than 0; got {lamb!r}"
            )


def _validate_guidance_number(
    value: Any,
    ctx: str,
    *,
    allow_zero: bool = False,
    maximum: float | None = None,
) -> None:
    valid = (
        not isinstance(value, bool)
        and isinstance(value, Real)
        and math.isfinite(value)
        and (value >= 0 if allow_zero else value > 0)
        and (maximum is None or value <= maximum)
    )
    if not valid:
        bounds = "greater than or equal to 0" if allow_zero else "greater than 0"
        if maximum is not None:
            bounds += f" and less than or equal to {maximum:g}"
        raise ConfigError(f"{ctx} must be a finite real number {bounds}; got {value!r}")


def _validate_local_refinement(section: Dict[str, Any], ctx: str) -> None:
    _reject_unknown_keys(section, LOCAL_REFINEMENT_KEYS, ctx)
    if "guidance" in section:
        guidance = section["guidance"]
        if not isinstance(guidance, str) or guidance not in GUIDANCE_POLICIES:
            raise ConfigError(
                f"{ctx}.guidance must be one of {sorted(GUIDANCE_POLICIES)}; "
                f"got {guidance!r}"
            )
    compatibility = section.get("compatibility")
    if compatibility is None:
        return
    compatibility = _ensure_mapping(compatibility, f"{ctx}.compatibility")
    _reject_unknown_keys(
        compatibility,
        LOCAL_REFINEMENT_COMPATIBILITY_KEYS,
        f"{ctx}.compatibility",
    )
    if (
        section.get("guidance") == "off"
        and set(compatibility) == LOCAL_REFINEMENT_COMPATIBILITY_KEYS
        and all(
            compatibility[key] is None
            for key in LOCAL_REFINEMENT_COMPATIBILITY_KEYS
        )
    ):
        return
    if "mode" in compatibility:
        mode = compatibility["mode"]
        if not isinstance(mode, str) or mode not in GUIDANCE_COMPATIBILITY_MODES:
            raise ConfigError(
                f"{ctx}.compatibility.mode must be one of "
                f"{sorted(GUIDANCE_COMPATIBILITY_MODES)}; got {mode!r}"
            )
    if "beta" in compatibility:
        _validate_guidance_number(compatibility["beta"], f"{ctx}.compatibility.beta")
    if "min_affinity" in compatibility:
        _validate_guidance_number(
            compatibility["min_affinity"],
            f"{ctx}.compatibility.min_affinity",
            maximum=1,
        )
    if "strength" in compatibility and compatibility["strength"] is not None:
        _validate_guidance_number(
            compatibility["strength"],
            f"{ctx}.compatibility.strength",
            allow_zero=True,
        )


def _validate_legacy_posterior_conditioning(
    section: Dict[str, Any],
    ctx: str,
) -> None:
    _reject_unknown_keys(section, POSTERIOR_CONDITIONING_KEYS, ctx)
    posterior_key = section.get("posterior_key")
    if posterior_key not in (None, ""):
        raise ConfigError(
            f"{ctx}.posterior_key is no longer supported; Assignment State is "
            "provided explicitly by each local-refinement route"
        )
    if "mode" in section:
        mode = section["mode"]
        if mode not in {"off", "cost", "reference"}:
            raise ConfigError(
                f"{ctx}.mode must be one of ['cost', 'off', 'reference']; got {mode!r}"
            )
    for key in ("enabled", "strict"):
        if key in section and not isinstance(section[key], bool):
            raise ConfigError(f"{ctx}.{key} must be a boolean; got {section[key]!r}")
    if "beta" in section:
        _validate_guidance_number(section["beta"], f"{ctx}.beta")
    if "min_affinity" in section:
        _validate_guidance_number(
            section["min_affinity"],
            f"{ctx}.min_affinity",
            maximum=1,
        )
    if "cost_strength" in section:
        _validate_guidance_number(
            section["cost_strength"],
            f"{ctx}.cost_strength",
            allow_zero=True,
        )


def _validate_sections(section_map: Dict[str, Any], ctx: str) -> None:
    if "annotate" in section_map:
        raise ConfigError(
            f"{ctx}.annotate is no longer supported; replace annotate.mode -> ot.ga.solver"
        )
    if "local_ot" in section_map:
        raise ConfigError(
            f"{ctx}.local_ot is no longer supported; replace local_ot.method -> ot.lr.solver"
        )
    _reject_unknown_keys(section_map, DEFAULT_SECTION_KEYS, ctx)
    for section, value in section_map.items():
        value_map = _ensure_mapping(value, f"{ctx}.{section}")
        if section == "runtime":
            if "ot_solver" in value_map:
                raise ConfigError(
                    f"{ctx}.runtime.ot_solver is no longer supported; replace "
                    "ot_solver -> ot.ga.solver + ot.lr.solver"
                )
            _reject_unknown_keys(value_map, RUNTIME_KEYS, f"{ctx}.runtime")
        elif section == "io":
            _reject_unknown_keys(value_map, IO_KEYS, f"{ctx}.io")
        elif section == "columns":
            _reject_unknown_keys(value_map, COLUMNS_KEYS, f"{ctx}.columns")
        elif section == "preprocess":
            _reject_unknown_keys(value_map, PREPROCESS_KEYS, f"{ctx}.preprocess")
        elif section == "graph":
            _reject_unknown_keys(value_map, GRAPH_KEYS, f"{ctx}.graph")
        elif section == "ot":
            _validate_ot_section(value_map, f"{ctx}.ot")
        elif section == "plot":
            _reject_unknown_keys(value_map, PLOT_KEYS, f"{ctx}.plot")
        elif section == "reconstruct":
            _reject_unknown_keys(value_map, RECONSTRUCT_KEYS, f"{ctx}.reconstruct")
        elif section == "benchmark":
            _reject_unknown_keys(value_map, BENCHMARK_KEYS, f"{ctx}.benchmark")
        elif section == "sc":
            _validate_sc_section(value_map, f"{ctx}.sc")
        elif section == "impute":
            _reject_unknown_keys(value_map, IMPUTE_KEYS, f"{ctx}.impute")
        elif section == "local_refinement":
            _validate_local_refinement(value_map, f"{ctx}.local_refinement")
        elif section == "posterior_conditioning":
            _validate_legacy_posterior_conditioning(
                value_map,
                f"{ctx}.posterior_conditioning",
            )


def _validate_router(router: Dict[str, Any]) -> None:
    for platform, conf_map in router.items():
        conf_map = _ensure_mapping(conf_map, f"router.{platform}")
        for confounding, route in conf_map.items():
            route_map = _ensure_mapping(route, f"router.{platform}.{confounding}")
            if "ot_solver" in route_map:
                raise ConfigError(
                    f"router.{platform}.{confounding}.ot_solver is no longer supported; "
                    "replace ot_solver -> ot.ga.solver + ot.lr.solver"
                )
            _reject_unknown_keys(route_map, ROUTE_LEAF_KEYS, f"router.{platform}.{confounding}")
            required = {
                "mode",
                "task",
                "svc_kind",
                "strategy",
            }
            missing = sorted(k for k in required if route_map.get(k) in (None, ""))
            if missing:
                raise ConfigError(
                    f"Missing required route keys in router.{platform}.{confounding}: {missing}"
                )


def _validate_locked_params(locked: Dict[str, Any]) -> None:
    _reject_unknown_keys(locked, LOCKED_PARAMS_KEYS, "locked_params")
    if "keys" in locked and not isinstance(locked["keys"], list):
        raise ConfigError("locked_params.keys must be a list")


def _validate_raw_config(raw: Dict[str, Any]) -> None:
    _reject_unknown_keys(raw, TOP_LEVEL_KEYS, "config root")

    defaults = _ensure_mapping(raw.get("defaults", {}), "defaults")
    _validate_sections(defaults, "defaults")

    profiles = _ensure_mapping(raw.get("profiles", {}), "profiles")
    for name, profile_cfg in profiles.items():
        _validate_sections(_ensure_mapping(profile_cfg, f"profiles.{name}"), f"profiles.{name}")

    router = _ensure_mapping(raw.get("router", {}), "router")
    _validate_router(router)

    locked = _ensure_mapping(raw.get("locked_params", {}), "locked_params")
    _validate_locked_params(locked)


def load_raw_config(path: str | Path) -> Dict[str, Any]:
    config_path = Path(path)
    if not config_path.exists():
        raise ConfigError(f"Config file not found: {config_path}")
    with config_path.open("r", encoding="utf-8") as f:
        raw = yaml.safe_load(f) or {}
    if not isinstance(raw, dict):
        raise ConfigError("Config root must be a mapping")
    _validate_raw_config(raw)
    return raw


def deep_merge(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    out = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = deep_merge(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


def set_by_dotted_path(config: Dict[str, Any], dotted_key: str, value: Any, create_missing: bool = False) -> None:
    parts = dotted_key.split(".")
    cur: Dict[str, Any] = config
    for part in parts[:-1]:
        if part not in cur:
            if not create_missing:
                raise KeyError(dotted_key)
            cur[part] = {}
        if not isinstance(cur[part], dict):
            raise KeyError(dotted_key)
        cur = cur[part]

    leaf = parts[-1]
    if leaf not in cur and not create_missing:
        raise KeyError(dotted_key)
    cur[leaf] = value


def _resolve_runtime_route(raw_config: Dict[str, Any], merged: Dict[str, Any]) -> Dict[str, Any]:
    runtime = merged.setdefault("runtime", {})
    platform = runtime.get("platform")
    if not platform:
        raise ConfigError("runtime.platform is required")

    confounding = runtime.get("confounding")
    if not confounding and platform in APPLICATION_DEFAULT_CF:
        confounding = APPLICATION_DEFAULT_CF[platform]
        runtime["confounding"] = confounding

    if platform == "sim2real" and not confounding:
        raise ConfigError("runtime.confounding is required when runtime.platform=sim2real")

    router = raw_config.get("router", {})
    route = router.get(platform, {}).get(confounding)
    if route is None:
        available = sorted(router.get(platform, {}).keys())
        raise ConfigError(
            f"No route found for platform={platform}, confounding={confounding}. "
            f"Available confounding values for {platform}: {available}"
        )

    runtime.update(route)
    hyperresolution = merged.get("sc", {}).get("hyperresolution", {})
    if runtime.get("task") == "sc_svc" and hyperresolution.get("enabled"):
        strategy = hyperresolution.get("strategy")
        if not strategy:
            raise ConfigError("sc.hyperresolution.strategy is required when enabled")
        runtime["strategy"] = strategy
    return route


def _validate_runtime(merged: Dict[str, Any]) -> None:
    runtime = merged.get("runtime", {})
    required = ["platform", "confounding", "mode", "task", "svc_kind", "strategy"]
    missing = [k for k in required if runtime.get(k) in (None, "")]
    if missing:
        raise ConfigError(f"Missing runtime keys after router resolution: {missing}")


def _validate_resolved_config(merged: Dict[str, Any]) -> None:
    """Strictly validate the complete config after every merge and route update."""
    _validate_sections(merged, "resolved")
    _validate_ot_section(
        _ensure_mapping(merged.get("ot", {}), "resolved.ot"),
        "resolved.ot",
        resolved=True,
    )
    _validate_runtime(merged)
    sc_cfg = _ensure_mapping(merged.get("sc", {}), "resolved.sc")
    if sc_cfg.get("svc_completeness") is not True:
        raise ConfigError("sc.svc_completeness must be exactly true")
    runtime = merged["runtime"]
    if runtime.get("mode") == "application" and runtime.get("task") == "sc_svc":
        solvers = {str(merged["ot"][phase]["solver"]) for phase in ("ga", "lr")}
        if "tacco" in solvers and sc_cfg.get("tacco_annotate") is None:
            raise ConfigError(
                "sc.tacco_annotate is required when application sc-SVC uses "
                "TACCO for Global Anchoring or Local Refinement"
            )


def _guidance_route_defaults(merged: Dict[str, Any]) -> Dict[str, Any]:
    runtime = merged["runtime"]
    guidance = "prefer" if runtime.get("task") == "sp_svc" else "off"
    strength = 1.0 if runtime.get("task") == "sc_svc" else 0.2
    return {
        "guidance": guidance,
        "compatibility": {
            "mode": "cost",
            "beta": 1.0,
            "min_affinity": 0.05,
            "strength": strength,
        },
    }


def _translate_legacy_guidance(
    legacy: Dict[str, Any],
    *,
    defaults: Dict[str, Any],
) -> Dict[str, Any]:
    enabled = legacy.get("enabled")
    mode = legacy.get("mode")
    strict = legacy.get("strict")
    disabled = enabled is False or mode == "off"

    if enabled is False and mode in {"cost", "reference"}:
        raise ConfigError(
            "posterior_conditioning has conflicting legacy enablement and mode"
        )
    if mode == "off" and enabled is True:
        raise ConfigError(
            "posterior_conditioning has conflicting legacy enablement and mode"
        )
    if disabled and strict is True:
        raise ConfigError(
            "posterior_conditioning strict=true conflicts with disabled guidance"
        )

    translated = copy.deepcopy(defaults)
    if disabled:
        translated["guidance"] = "off"
    elif strict is True:
        translated["guidance"] = "require"
    elif enabled is True or mode in {"cost", "reference"}:
        translated["guidance"] = "prefer"

    compatibility = translated["compatibility"]
    if mode in {"cost", "reference"}:
        compatibility["mode"] = mode
    for old, new in (
        ("beta", "beta"),
        ("min_affinity", "min_affinity"),
        ("cost_strength", "strength"),
    ):
        if old in legacy:
            compatibility[new] = legacy[old]
    return translated


def _normalize_guidance_state(state: Dict[str, Any]) -> Dict[str, Any]:
    normalized = copy.deepcopy(state)
    if normalized["guidance"] == "off":
        normalized["compatibility"] = {
            "mode": None,
            "beta": None,
            "min_affinity": None,
            "strength": None,
        }
    return normalized


def _finalize_guidance_state(
    state: Dict[str, Any],
    *,
    defaults: Dict[str, Any],
) -> Dict[str, Any]:
    finalized = copy.deepcopy(state)
    if (
        finalized["guidance"] != "off"
        and finalized["compatibility"].get("strength") is None
    ):
        finalized["compatibility"]["strength"] = defaults["compatibility"]["strength"]
    return _normalize_guidance_state(finalized)


def _resolve_assignment_guidance(
    merged: Dict[str, Any],
) -> tuple[Dict[str, Any], Dict[str, Any]]:
    defaults = _guidance_route_defaults(merged)
    new = merged.pop("local_refinement", None)
    legacy = merged.pop("posterior_conditioning", None)
    if new is not None:
        _validate_local_refinement(
            _ensure_mapping(new, "resolved.local_refinement"),
            "resolved.local_refinement",
        )
    if legacy is not None:
        _validate_legacy_posterior_conditioning(
            _ensure_mapping(legacy, "resolved.posterior_conditioning"),
            "resolved.posterior_conditioning",
        )

    new_state = deep_merge(defaults, new or {})
    new_effective = _finalize_guidance_state(new_state, defaults=defaults)
    legacy_state = (
        _translate_legacy_guidance(legacy, defaults=defaults)
        if legacy is not None
        else None
    )
    legacy_effective = (
        _finalize_guidance_state(legacy_state, defaults=defaults)
        if legacy_state is not None
        else None
    )
    if new is not None and legacy_effective is not None:
        new_dimensions = set()
        if "guidance" in new:
            new_dimensions.add("guidance")
        new_dimensions.update((new.get("compatibility") or {}).keys())

        legacy_dimensions = set()
        if set(legacy) & {"enabled", "mode", "strict"}:
            legacy_dimensions.add("guidance")
        if legacy.get("mode") in {"cost", "reference", "off"}:
            legacy_dimensions.add("mode")
        for legacy_key, dimension in (
            ("beta", "beta"),
            ("min_affinity", "min_affinity"),
            ("cost_strength", "strength"),
        ):
            if legacy_key in legacy:
                legacy_dimensions.add(dimension)

        overlap = new_dimensions & legacy_dimensions
        conflicting = (
            "guidance" in overlap
            and new_effective["guidance"] != legacy_effective["guidance"]
        )
        composed = deep_merge(legacy_state, new)
        resolved = _finalize_guidance_state(composed, defaults=defaults)
        if resolved["guidance"] != "off":
            for dimension in overlap - {"guidance"}:
                if (
                    new_effective["compatibility"][dimension]
                    != legacy_effective["compatibility"][dimension]
                ):
                    conflicting = True
        if conflicting:
            raise ConfigError(
                "local_refinement conflicts with legacy posterior_conditioning"
            )
        source = "new"
    elif new is not None:
        resolved = new_effective
        source = "new"
    elif legacy_effective is not None:
        resolved = legacy_effective
        source = "legacy"
    else:
        resolved = _finalize_guidance_state(defaults, defaults=defaults)
        source = "route_default"

    if new is not None:
        configured_guidance = new.get(
            "guidance",
            legacy_state["guidance"]
            if legacy is not None
            and set(legacy) & {"enabled", "mode", "strict"}
            else None,
        )
        configured_mode = (new.get("compatibility") or {}).get(
            "mode",
            legacy.get("mode")
            if legacy is not None
            and legacy.get("mode") in GUIDANCE_COMPATIBILITY_MODES
            else None,
        )
    elif legacy is not None:
        configured_guidance = (
            legacy_state["guidance"]
            if set(legacy) & {"enabled", "mode", "strict"}
            else None
        )
        configured_mode = (
            legacy.get("mode")
            if legacy.get("mode") in GUIDANCE_COMPATIBILITY_MODES
            else None
        )
    else:
        configured_guidance = None
        configured_mode = None
    evidence = {
        "assignment_guidance": {
            "configured_guidance": configured_guidance,
            "configured_compatibility_mode": configured_mode,
            "resolution_source": source,
            "deprecations": (
                [
                    f"posterior_conditioning.{key}"
                    for key in sorted(legacy)
                ]
                if legacy is not None
                else []
            ),
        }
    }
    return resolved, evidence


def _paths_overlap(left: str, right: str) -> bool:
    return (
        left == right
        or left.startswith(f"{right}.")
        or right.startswith(f"{left}.")
    )


def _leaf_paths(values: Dict[str, Any], prefix: str = "") -> List[str]:
    paths: List[str] = []
    for key, value in values.items():
        path = f"{prefix}.{key}" if prefix else key
        if isinstance(value, dict):
            paths.extend(_leaf_paths(value, path))
        else:
            paths.append(path)
    return paths


def merge_unified_config(
    raw_config: Dict[str, Any],
    profile: str | None,
    runtime_overrides: Dict[str, Any],
    io_overrides: Dict[str, Any],
    algorithm_overrides: Dict[str, Any],
) -> ResolvedConfig:
    # Merge order is intentionally strict and explicit:
    # 1) defaults
    # 2) selected profile overrides
    # 3) runtime/io CLI overrides
    # 4) trusted, structured algorithm overrides
    #
    # This mirrors the design spec and makes provenance deterministic.
    defaults = raw_config.get("defaults", {})
    if not isinstance(defaults, dict):
        raise ConfigError("defaults must be a mapping")

    merged = copy.deepcopy(defaults)

    if profile:
        profiles = raw_config.get("profiles", {})
        if profile not in profiles:
            raise ConfigError(f"Unknown profile: {profile}")
        merged = deep_merge(merged, profiles[profile])

    if "ot_solver" in runtime_overrides:
        raise ConfigError(
            "runtime/router ot_solver -> ot.ga.solver + ot.lr.solver"
        )

    for key, value in runtime_overrides.items():
        if value is None and key != "seed":
            continue
        set_by_dotted_path(merged, f"runtime.{key}", value, create_missing=False)

    for key, value in io_overrides.items():
        if value is None:
            continue
        set_by_dotted_path(merged, f"io.{key}", value, create_missing=False)

    locked_config = raw_config.get("locked_params", {})
    locked_keys = set(locked_config.get("keys", []))
    forbidden_sections = sorted(set(algorithm_overrides) & {"runtime", "io"})
    if forbidden_sections:
        raise ConfigError(
            "algorithm_overrides cannot modify run identity sections: "
            + ", ".join(forbidden_sections)
        )
    override_paths = _leaf_paths(algorithm_overrides)
    identity_paths = sorted(
        guarded
        for guarded in ALGORITHM_IDENTITY_PATHS
        if any(_paths_overlap(key, guarded) for key in override_paths)
    )
    if identity_paths:
        raise ConfigError(
            "algorithm_overrides cannot modify run identity through: "
            + ", ".join(identity_paths)
        )
    for key in override_paths:
        if any(_paths_overlap(key, locked) for locked in locked_keys):
            raise ConfigError(
                f"Override rejected for locked parameter '{key}'. "
                "These low-level parameters are governed by the selected profile."
            )
    merged = deep_merge(merged, algorithm_overrides)

    _resolve_runtime_route(raw_config, merged)
    resolved_guidance, request_evidence = _resolve_assignment_guidance(merged)
    merged["local_refinement"] = resolved_guidance
    _validate_resolved_config(merged)
    return ResolvedConfig(merged, request_evidence=request_evidence)


def infer_default_profile(raw_config: Dict[str, Any], runtime_overrides: Dict[str, Any]) -> str | None:
    """Pick a profile based on explicit runtime route, if a direct match exists."""
    runtime = runtime_overrides
    if not runtime.get("platform"):
        return None

    profiles = raw_config.get("profiles", {})
    for name, profile_cfg in profiles.items():
        pr = profile_cfg.get("runtime", {})
        if runtime.get("platform") and pr.get("platform") != runtime.get("platform"):
            continue
        if runtime.get("confounding") and pr.get("confounding") != runtime.get("confounding"):
            continue
        return name
    return None
