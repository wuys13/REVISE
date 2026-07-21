from __future__ import annotations

from typing import Any, Dict

from revise.backend.contracts import CFStrategy
from revise.backend.contracts import PlatformAdapter


class BasePlatformAdapter(PlatformAdapter):
    """Default platform adapter base class.

    Platform adapters are intentionally lightweight in the unified runtime:
    they normalize/validate platform-level payload shape before confounding
    factor plugins are applied. OT selection is handled by the resolved
    ``ot.ga`` and ``ot.lr`` configuration in the runner adapters.
    """

    platform_id: str = "default"

    def adapt(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        runtime = payload.setdefault("runtime", {})
        platform = runtime.get("platform")
        if not platform:
            raise ValueError("runtime.platform is required before platform adapter resolution")
        return payload


class Sim2RealPlatformAdapter(BasePlatformAdapter):
    platform_id = "sim2real"

    def adapt(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        payload = super().adapt(payload)
        runtime = payload["runtime"]
        if runtime.get("confounding") is None:
            raise ValueError("runtime.confounding is required for platform=sim2real")
        return payload


class HSTPlatformAdapter(BasePlatformAdapter):
    platform_id = "hST"


class ISTPlatformAdapter(BasePlatformAdapter):
    platform_id = "iST"


class SSTPlatformAdapter(BasePlatformAdapter):
    platform_id = "sST"


class DefaultCFStrategy(CFStrategy):
    """No-op CF strategy used for routes without special plugin behavior."""

    cf_id: str = "default"

    def apply(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        return payload


class SegmentationCFStrategy(DefaultCFStrategy):
    """Segmentation CF strategy with nested sc-hyperresolution switching.

    When application sc-SVC enables `sc.hyperresolution.enabled`, we route to
    a dedicated hyperresolution strategy implementation instead of the standard
    sc strategy.
    """

    cf_id = "segmentation"

    def apply(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        payload = super().apply(payload)
        runtime = payload.setdefault("runtime", {})
        merged_config = payload.get("merged_config", {})
        sc_cfg = merged_config.get("sc", {})
        hyper_cfg = sc_cfg.get("hyperresolution", {})

        enabled = bool(hyper_cfg.get("enabled", False))
        is_sc_application = runtime.get("mode") == "application" and runtime.get("task") == "sc_svc"
        if enabled and is_sc_application:
            strategy = hyper_cfg.get("strategy", "ScSvcHyperApplicationStrategy")
            runtime["strategy"] = strategy
            payload["resolution_trace"] = {
                "hyperresolution_enabled": True,
                "strategy": strategy,
            }
        return payload
