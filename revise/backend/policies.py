from __future__ import annotations

import tempfile
from pathlib import Path

from revise.backend.contracts import EvaluationPolicy
from revise.backend.contracts import InputValidationPolicy
from revise.backend.ops.tacco_runtime import require_tacco
from revise.config.runner_conf import REQUIRED_IO_BY_MODE_TASK
from revise.config.runner_conf import pm_on_cell_path_from_data_root
from revise.config.runner_conf import resolve_input_specs
from revise.io import REVISEInputService
from revise.io.input_service import PMOnCellSnapshotError
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
        direct_application = bool(
            mode == "application"
            and ctx.st_adata is not None
            and ctx.sc_ref_adata is not None
        )
        if not direct_application:
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

        required = (
            {"st_path", "sc_ref_path"}
            if direct_application
            else REQUIRED_IO_BY_MODE_TASK[key]
        )
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
        identities = []
        try:
            identities = input_identities(specs)
            ctx.pm_on_cell = None
            if (mode, task) == ("application", "sc_svc_super_resolution"):
                pm_path = io_cfg.get("pm_on_cell_path") or None
            elif (mode, task) == ("benchmark", "sc_svc_sr"):
                pm_path = (
                    pm_on_cell_path_from_data_root(str(data_root))
                    if data_root
                    else None
                )
            else:
                pm_path = None
            if pm_path is not None:
                if not Path(pm_path).exists():
                    raise FileNotFoundError(
                        f"io.pm_on_cell_path does not exist: {pm_path}"
                    )
                ctx.pm_on_cell, pm_identity = input_service.snapshot_pm_on_cell(
                    pm_path
                )
                if pm_identity is not None:
                    identities.append(pm_identity)
            identities.sort(key=lambda identity: identity["role"])
            ctx.input_identities = identities
        except (OSError, RuntimeError, TypeError, ValueError) as exc:
            if isinstance(exc, PMOnCellSnapshotError):
                identities.append(exc.identity)
            ctx.input_identities = sorted(
                identities,
                key=lambda identity: identity["role"],
            )
            raise ValueError(
                "Invalid input identities: "
                f"{type(exc).__name__}: {exc}"
            ) from exc
        direct_inputs = ctx.st_adata is not None and ctx.sc_ref_adata is not None
        if direct_inputs:
            if mode != "application":
                raise ValueError("preloaded AnnData is only supported for application runs")
            report = {"inputs": "preloaded"}
        else:
            validation_columns = dict(ctx.columns)
            if mode == "application" and str(task) == "sc_svc":
                validation_columns["select_cell_type"] = (
                    ctx.merged_config.get("sc", {}) or {}
                ).get("select_ct")
            reference_filter = ctx.application_config_metadata.get(
                "reference_filter", {}
            )
            report = input_service.preflight(
                specs,
                runtime=runtime,
                columns=validation_columns,
                reference_filter_column=reference_filter.get("column"),
                reference_filter_value=reference_filter.get("value"),
            )
        run_dir = Path(ctx.run_dir)
        run_dir.mkdir(parents=True, exist_ok=True)
        publication_dir = (
            Path(output_root)
            if mode == "application"
            else Path(output_root) / str(sample_name)
        )
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
        return bool(ctx.merged_config["benchmark"]["evaluate"])
