from __future__ import annotations

from typing import Any, Dict

from revise.framework import REVISEPipeline


def sp_svc(
    *,
    pipeline: REVISEPipeline,
    profile: str = "application_sp",
    runtime_overrides: Dict[str, Any] | None = None,
    io_overrides: Dict[str, Any] | None = None,
):
    return pipeline.run(
        profile=profile,
        runtime_overrides=runtime_overrides or {},
        io_overrides=io_overrides or {},
    )


def sc_svc(
    *,
    pipeline: REVISEPipeline,
    profile: str = "application_sc",
    runtime_overrides: Dict[str, Any] | None = None,
    io_overrides: Dict[str, Any] | None = None,
):
    return pipeline.run(
        profile=profile,
        runtime_overrides=runtime_overrides or {},
        io_overrides=io_overrides or {},
    )
