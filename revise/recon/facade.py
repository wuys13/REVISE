from __future__ import annotations

from typing import Any, Dict, Iterable

from revise.framework import REVISEPipeline


def sp_svc(
    *,
    pipeline: REVISEPipeline,
    profile: str = "application_sp",
    runtime_overrides: Dict[str, Any] | None = None,
    io_overrides: Dict[str, Any] | None = None,
    set_overrides: Iterable[str] = (),
):
    return pipeline.run(
        profile=profile,
        runtime_overrides=runtime_overrides or {},
        io_overrides=io_overrides or {},
        set_overrides=list(set_overrides),
    )


def sc_svc(
    *,
    pipeline: REVISEPipeline,
    profile: str = "application_sc",
    runtime_overrides: Dict[str, Any] | None = None,
    io_overrides: Dict[str, Any] | None = None,
    set_overrides: Iterable[str] = (),
):
    return pipeline.run(
        profile=profile,
        runtime_overrides=runtime_overrides or {},
        io_overrides=io_overrides or {},
        set_overrides=list(set_overrides),
    )
