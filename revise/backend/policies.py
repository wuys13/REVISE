from __future__ import annotations

import tempfile
from pathlib import Path

from revise.backend.contracts import EvaluationPolicy
from revise.backend.contracts import InputValidationPolicy
from revise.backend.ops.tacco_runtime import require_tacco
from revise.config.runner_conf import REQUIRED_IO_BY_MODE_TASK
from revise.config.runner_conf import resolve_input_specs
from revise.io import REVISEInputService
from revise.utils import completed_artifact, fingerprint_paths, write_json


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

        self._validate_solver_compatibility(ctx)
        if any(
            str(ctx.merged_config["ot"][phase]["solver"]) == "tacco"
            for phase in ("ga", "lr")
        ):
            require_tacco()

        specs = resolve_input_specs(runtime, io_cfg)
        ctx.input_specs = specs
        try:
            ctx.data_fingerprint = fingerprint_paths(specs)
        except (OSError, RuntimeError, TypeError, ValueError) as exc:
            ctx.data_fingerprint = None
            ctx.data_fingerprint_error = {
                "type": type(exc).__name__,
                "message": str(exc),
            }
            raise ValueError(
                "Invalid input content fingerprint: "
                f"{type(exc).__name__}: {exc}"
            ) from exc
        report = REVISEInputService.from_context(ctx).preflight(
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

    @staticmethod
    def _validate_solver_compatibility(ctx) -> None:
        conditioning = ctx.merged_config.get("posterior_conditioning", {}) or {}
        uses_reference = bool(conditioning.get("enabled", False)) and str(
            conditioning.get("mode", "off")
        ) == "reference"
        lr_solver = str(ctx.merged_config["ot"]["lr"]["solver"])
        if uses_reference and lr_solver == "tacco":
            raise ValueError(
                "TACCO local refinement is incompatible with reference "
                "posterior conditioning; use conditioning mode 'cost'/'off' "
                "or set ot.lr.solver=pot"
            )


class ModeEvaluationPolicy(EvaluationPolicy):
    def should_evaluate(self, ctx) -> bool:
        runtime = ctx.runtime
        if runtime.get("mode") == "benchmark":
            return True
        return bool(ctx.merged_config.get("benchmark", {}).get("evaluate", False))
