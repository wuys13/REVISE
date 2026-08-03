from __future__ import annotations

import tempfile
from pathlib import Path

from revise.backend.contracts import EvaluationPolicy
from revise.backend.contracts import InputValidationPolicy
from revise.backend.ops.tacco_runtime import require_tacco
from revise.config.runner_conf import REQUIRED_IO_BY_MODE_TASK
from revise.config.runner_conf import pm_on_cell_path_from_st_path
from revise.config.runner_conf import resolve_input_specs
from revise.config.runner_conf import resolved_input_path
from revise.io import REVISEInputService
from revise.utils import completed_artifact, input_identities, write_json


class ModeValidationPolicy(InputValidationPolicy):
    def validate(self, ctx) -> None:
        runtime = ctx.runtime
        mode = runtime.get("mode")
        if mode not in {"application", "benchmark"}:
            raise ValueError(f"Unsupported mode: {mode}")
        task = runtime.get("task")
        if task in (None, ""):
            raise ValueError("runtime.task is required")

        io_cfg = ctx.io
        data_root = io_cfg.get("data_root")
        if not data_root:
            raise ValueError("io.data_root is required")
        if not Path(data_root).exists():
            raise FileNotFoundError(f"io.data_root does not exist: {data_root}")

        sample_name = io_cfg.get("sample_name")
        if not sample_name:
            raise ValueError("io.sample_name is required")
        output_root = io_cfg.get("output_root")
        if not output_root:
            raise ValueError("io.output_root is required")

        key = (mode, task)
        if key not in REQUIRED_IO_BY_MODE_TASK:
            supported = sorted(f"{m}:{t}" for (m, t) in REQUIRED_IO_BY_MODE_TASK.keys())
            raise ValueError(
                f"Unsupported mode/task combination: {mode}:{task}. "
                f"Supported combinations: {supported}"
            )

        required = REQUIRED_IO_BY_MODE_TASK[key]
        missing = sorted(k for k in required if io_cfg.get(k) in (None, ""))
        if missing:
            raise ValueError(
                f"Missing required io keys for {mode}:{task}: {missing}"
            )

        if any(
            str(ctx.merged_config["ot"][phase]["solver"]) == "tacco"
            for phase in ("ga", "lr")
        ):
            require_tacco()

        specs = resolve_input_specs(runtime, io_cfg)
        ctx.input_specs = specs
        input_service = REVISEInputService.from_context(ctx)
        try:
            identities = input_identities(specs)
            ctx.pm_on_cell = None
            if str(task) == "sc_svc_sr":
                st_path = resolved_input_path(specs, "st", "")
                pm_path = pm_on_cell_path_from_st_path(st_path)
                ctx.pm_on_cell, pm_identity = input_service.snapshot_pm_on_cell(
                    pm_path
                )
                if pm_identity is not None:
                    identities.append(pm_identity)
            identities.sort(key=lambda identity: identity["role"])
            ctx.input_identities = identities
        except (OSError, RuntimeError, TypeError, ValueError) as exc:
            ctx.input_identities = []
            raise ValueError(
                "Invalid input identities: "
                f"{type(exc).__name__}: {exc}"
            ) from exc
        report = input_service.preflight(
            specs,
            runtime=runtime,
            columns=ctx.columns,
        )
        run_dir = Path(ctx.run_dir)
        run_dir.mkdir(parents=True, exist_ok=True)
        publication_dir = Path(output_root) / str(sample_name)
        publication_dir.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            dir=publication_dir,
            prefix=".preflight-write-",
        ):
            pass
        report["output"] = {
            "run_dir": str(run_dir),
            "publication_dir": str(publication_dir),
            "writable": True,
        }
        report_path = run_dir / "preflight.json"
        write_json(report_path, report)
        ctx.record_artifact(completed_artifact("preflight", report_path))

class ModeEvaluationPolicy(EvaluationPolicy):
    def should_evaluate(self, ctx) -> bool:
        runtime = ctx.runtime
        if runtime.get("mode") == "benchmark":
            return True
        return bool(ctx.merged_config.get("benchmark", {}).get("evaluate", False))
