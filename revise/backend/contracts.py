from __future__ import annotations

from abc import ABC
from abc import abstractmethod
if False:  # pragma: no cover
    from revise.recon.context import PipelineContext
    from revise.svc import SVC


class LocalRefinementStrategy(ABC):
    """Task-level reconstruction strategy contract."""

    strategy_id: str

    @abstractmethod
    def prepare_context(self, ctx: "PipelineContext") -> None:
        raise NotImplementedError

    @abstractmethod
    def global_anchoring(self, ctx: "PipelineContext") -> None:
        raise NotImplementedError

    def prepare_local_units(self, ctx: "PipelineContext") -> None:
        return None

    def build_graph(self, ctx: "PipelineContext") -> None:
        return None

    def build_ot_problem(self, ctx: "PipelineContext") -> None:
        return None

    @abstractmethod
    def solve_ot(self, ctx: "PipelineContext") -> None:
        raise NotImplementedError

    def update_expression(self, ctx: "PipelineContext") -> None:
        return None

    @abstractmethod
    def finalize_svc(self, ctx: "PipelineContext") -> "SVC":
        raise NotImplementedError


class InputValidationPolicy(ABC):
    @abstractmethod
    def validate(self, ctx: "PipelineContext") -> None:
        raise NotImplementedError


class EvaluationPolicy(ABC):
    @abstractmethod
    def should_evaluate(self, ctx: "PipelineContext") -> bool:
        raise NotImplementedError
