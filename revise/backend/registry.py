from __future__ import annotations

from typing import Dict

from revise.backend.contracts import LocalRefinementStrategy


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


def build_default_registry() -> StrategyRegistry:
    from revise.backend.adapters import (
        ScSvcApplicationStrategy,
        ScSvcImputeBenchmarkStrategy,
        ScSvcSuperResolutionApplicationStrategy,
        ScSvcSrBenchmarkStrategy,
        SpSvcApplicationStrategy,
        SpSvcBenchmarkSegStrategy,
    )

    reg = StrategyRegistry()
    reg.register(SpSvcApplicationStrategy())
    reg.register(ScSvcApplicationStrategy())
    reg.register(ScSvcSuperResolutionApplicationStrategy())
    reg.register(SpSvcBenchmarkSegStrategy())
    reg.register(ScSvcSrBenchmarkStrategy())
    reg.register(ScSvcImputeBenchmarkStrategy())
    return reg
