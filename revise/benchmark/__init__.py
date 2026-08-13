"""Sim2Real-ST benchmark entrypoints."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from revise.benchmark.cli import run_benchmark as run_benchmark

__all__ = ["run_benchmark"]


def __getattr__(name: str):
    if name == "run_benchmark":
        from revise.benchmark.cli import run_benchmark

        return run_benchmark
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
