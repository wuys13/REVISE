from __future__ import annotations

from typing import Any, Dict

from revise.framework import REVISEPipeline


def sp_svc(
    *,
    pipeline: REVISEPipeline,
    runtime_overrides: Dict[str, Any] | None = None,
    io_overrides: Dict[str, Any] | None = None,
):
    return pipeline.run(
        svc_type="sp-SVC",
        runtime_overrides=runtime_overrides or {},
        io_overrides=io_overrides or {},
    )


def sc_svc(
    *,
    pipeline: REVISEPipeline,
    runtime_overrides: Dict[str, Any] | None = None,
    io_overrides: Dict[str, Any] | None = None,
):
    return pipeline.run(
        svc_type="sc-SVC",
        runtime_overrides=runtime_overrides or {},
        io_overrides=io_overrides or {},
    )
