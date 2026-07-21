from __future__ import annotations

from typing import Dict

from revise.backend.contracts import CFStrategy
from revise.backend.contracts import LocalRefinementStrategy
from revise.backend.contracts import PlatformAdapter


class StrategyRegistry:
    def __init__(self) -> None:
        self._strategies: Dict[str, LocalRefinementStrategy] = {}

    def register(self, strategy: LocalRefinementStrategy) -> None:
        sid = strategy.strategy_id
        if sid in self._strategies:
            raise KeyError(f"Strategy already registered: {sid}")
        self._strategies[sid] = strategy

    def get(self, strategy_id: str) -> LocalRefinementStrategy:
        if strategy_id not in self._strategies:
            raise KeyError(f"Strategy not registered: {strategy_id}")
        return self._strategies[strategy_id]

    def available(self):
        return sorted(self._strategies.keys())


class PluginRegistry:
    def __init__(self) -> None:
        self._platform_adapters: Dict[str, PlatformAdapter] = {}
        self._cf_strategies: Dict[str, CFStrategy] = {}

    def register_platform_adapter(self, adapter_id: str, adapter: PlatformAdapter) -> None:
        if adapter_id in self._platform_adapters:
            raise KeyError(f"Platform adapter already registered: {adapter_id}")
        self._platform_adapters[adapter_id] = adapter

    def get_platform_adapter(self, adapter_id: str) -> PlatformAdapter:
        if adapter_id in self._platform_adapters:
            return self._platform_adapters[adapter_id]
        if "default" in self._platform_adapters:
            return self._platform_adapters["default"]
        raise KeyError(f"Platform adapter not registered: {adapter_id}")

    def register_cf_strategy(self, cf_id: str, strategy: CFStrategy) -> None:
        if cf_id in self._cf_strategies:
            raise KeyError(f"CF strategy already registered: {cf_id}")
        self._cf_strategies[cf_id] = strategy

    def get_cf_strategy(self, cf_id: str) -> CFStrategy:
        if cf_id in self._cf_strategies:
            return self._cf_strategies[cf_id]
        if "default" in self._cf_strategies:
            return self._cf_strategies["default"]
        raise KeyError(f"CF strategy not registered: {cf_id}")

def build_default_registry() -> StrategyRegistry:
    from revise.backend.adapters import (
        ScSvcApplicationStrategy,
        ScSvcHyperApplicationStrategy,
        ScSvcImputeBenchmarkStrategy,
        ScSvcSrApplicationStrategy,
        ScSvcSrBenchmarkStrategy,
        SpSvcApplicationStrategy,
        SpSvcBenchmarkSegStrategy,
    )

    reg = StrategyRegistry()
    reg.register(SpSvcApplicationStrategy())
    reg.register(ScSvcApplicationStrategy())
    reg.register(ScSvcHyperApplicationStrategy())
    reg.register(ScSvcSrApplicationStrategy())
    reg.register(SpSvcBenchmarkSegStrategy())
    reg.register(ScSvcSrBenchmarkStrategy())
    reg.register(ScSvcImputeBenchmarkStrategy())
    return reg


def build_default_plugin_registry() -> PluginRegistry:
    from revise.backend.plugins import (
        BasePlatformAdapter,
        DefaultCFStrategy,
        HSTPlatformAdapter,
        ISTPlatformAdapter,
        SSTPlatformAdapter,
        SegmentationCFStrategy,
        Sim2RealPlatformAdapter,
    )

    reg = PluginRegistry()
    reg.register_platform_adapter("default", BasePlatformAdapter())
    reg.register_platform_adapter("sim2real", Sim2RealPlatformAdapter())
    reg.register_platform_adapter("hST", HSTPlatformAdapter())
    reg.register_platform_adapter("iST", ISTPlatformAdapter())
    reg.register_platform_adapter("sST", SSTPlatformAdapter())

    reg.register_cf_strategy("default", DefaultCFStrategy())
    reg.register_cf_strategy("segmentation", SegmentationCFStrategy())
    for cf in ["bin2cell", "batch_effect", "spot_size", "gene_panel", "gene_dropout"]:
        reg.register_cf_strategy(cf, DefaultCFStrategy())

    return reg
