#!/usr/bin/env python3
"""The readable REVISE Application entrypoint.

The YAML is the baseline.  Explicit command-line/Python values are applied
after loading it, then the effective application configuration is mapped to
the shared engine.  Keeping that assembly here makes the complete user flow
visible in one file while Benchmark continues to use the same engine method.
"""

from __future__ import annotations

import argparse
from contextlib import redirect_stdout
import copy
from dataclasses import dataclass
import json
from pathlib import Path
import sys
from typing import Any, Mapping

from revise._version import __version__
from revise.application.config import (
    ApplicationConfig,
    ApplicationConfigError,
    compile_application_config,
    load_application_yaml,
)
from revise.framework import REVISEPipeline
from revise.utils import build_run_dir, completed_artifact


@dataclass(frozen=True)
class ApplicationExecution:
    status: str
    svc_type: str
    output_paths: dict[str, Path]
    pipeline: dict[str, Any]
    results: dict[str, Any] | None = None
    preflight: Path | None = None
    application_config: dict[str, Any] | None = None


_LEGACY_FLAGS = {
    "--sample-name": "output.name",
    "--st-file": "inputs.st.path",
    "--sc-ref-file": "inputs.reference.path",
    "--data-root": "inputs.*.path",
    "--output-root": "output.dir",
    "--patient-key": "prepare the reference before running",
}
_REMOVED_FLAGS = {
    "--set": "generic engine overrides are not supported",
    "--profile": "engine profiles are managed by the package",
    "--config-engine": "the engine config is managed by the package",
}


class _ApplicationArgumentParser(argparse.ArgumentParser):
    def parse_known_args(self, args=None, namespace=None):
        raw_args = sys.argv[1:] if args is None else list(args)
        for token in raw_args:
            flag = str(token).split("=", 1)[0]
            if flag in _LEGACY_FLAGS:
                self.error(f"{flag} was removed; use YAML field {_LEGACY_FLAGS[flag]}")
            if flag in _REMOVED_FLAGS:
                self.error(f"{flag} is not supported; {_REMOVED_FLAGS[flag]}")
        return super().parse_known_args(raw_args, namespace)


def build_parser() -> argparse.ArgumentParser:
    parser = _ApplicationArgumentParser(
        description="Reconstruct one SVC from an application YAML",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Choose a template, edit its inputs and output name, then run:\n"
            "  python reconstruct.py --config configs/application/Xenium_T.yaml --dry-run\n"
            "  python reconstruct.py --config configs/application/Xenium_T.yaml\n\n"
            "The installed command uses this same implementation:\n"
            "  revise-reconstruct --config configs/application/Xenium_T.yaml"
        ),
    )
    parser.add_argument("--version", action="version", version=f"revise-reconstruct {__version__}")

    application = parser.add_argument_group("Application")
    application.add_argument("--config", required=True, help="application YAML baseline")
    application.add_argument("--svc-type", choices=("sp-SVC", "sc-SVC", "sc-SVC-sr"), default=None,
                             help="override application.svc_type")

    inputs = parser.add_argument_group("Inputs")
    inputs.add_argument("--root-dir", default=None, help="override paths.root_dir")
    inputs.add_argument("--st-path", default=None, help="override inputs.st.path")
    inputs.add_argument("--st-format", choices=("h5ad", "spatialdata", "auto"), default=None)
    inputs.add_argument("--sc-ref-path", default=None, help="override inputs.reference.path")
    inputs.add_argument("--spatialdata-table", default=None)
    inputs.add_argument("--spatialdata-element", default=None)
    inputs.add_argument("--pm-on-cell-path", default=None, help="only for sc-SVC-sr")

    ot = parser.add_argument_group("Shared OT")
    ot.add_argument("--ot-method", choices=("pot", "tacco"), default=None,
                    help="set the GA and LR solver together")

    ga = parser.add_argument_group("Global Anchoring")
    ga.add_argument("--cell-type-col", default=None,
                    help="override global_anchoring.broad_column")

    lr = parser.add_argument_group("Local Refinement")
    lr.add_argument("--sub-cell-type-col", default=None, help="only for sc-SVC")
    lr.add_argument("--select-ct", default=None, help="only for sc-SVC")
    lr.add_argument("--local-refinement-strength", type=float, default=None,
                    help="only for sp-SVC and sc-SVC-sr")

    output = parser.add_argument_group("Output")
    output.add_argument("--output-dir", default=None)
    output.add_argument("--output-name", default=None)

    execution = parser.add_argument_group("Execution")
    execution.add_argument("--seed", type=int, default=None)
    execution.add_argument("--dry-run", action="store_true",
                           help="preflight inputs and route without writing H5AD")
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_parser().parse_args(argv)


def _set_value(document: dict[str, Any], section: str, key: str, value: Any) -> None:
    document.setdefault(section, {})[key] = value


def _set_nested(mapping: dict[str, Any], key: str, value: Any) -> None:
    mapping[key] = value


def _apply_overrides(document: dict[str, Any], values: Mapping[str, Any]) -> dict[str, Any]:
    """Apply only explicitly supplied high-level values to a YAML document."""
    overrides: dict[str, Any] = {}

    def set_if(name: str, section: str, key: str) -> None:
        value = values.get(name)
        if value is not None:
            _set_value(document, section, key, value)
            overrides[name] = value

    set_if("svc_type", "application", "svc_type")
    set_if("root_dir", "paths", "root_dir")
    set_if("ot_method", "algorithm", "ot_method")
    set_if("cell_type_col", "global_anchoring", "broad_column")
    set_if("sub_cell_type_col", "local_refinement", "subtype_column")
    set_if("select_cell_type", "local_refinement", "select_cell_type")
    set_if("local_refinement_strength", "local_refinement", "strength")
    set_if("output_dir", "output", "dir")
    set_if("output_name", "output", "name")
    set_if("seed", "execution", "seed")

    if values.get("st_path") is not None:
        _set_nested(document.setdefault("inputs", {}).setdefault("st", {}), "path", values["st_path"])
        overrides["st_path"] = values["st_path"]
    if values.get("st_format") is not None:
        _set_nested(document.setdefault("inputs", {}).setdefault("st", {}), "format", values["st_format"])
        overrides["st_format"] = values["st_format"]
    if values.get("sc_ref_path") is not None:
        _set_nested(document.setdefault("inputs", {}).setdefault("reference", {}), "path", values["sc_ref_path"])
        overrides["sc_ref_path"] = values["sc_ref_path"]
    if values.get("spatialdata_table") is not None:
        _set_nested(
            document.setdefault("inputs", {}).setdefault("st", {}).setdefault("spatialdata", {}),
            "table", values["spatialdata_table"],
        )
        overrides["spatialdata_table"] = values["spatialdata_table"]
    if values.get("spatialdata_element") is not None:
        _set_nested(
            document.setdefault("inputs", {}).setdefault("st", {}).setdefault("spatialdata", {}),
            "element", values["spatialdata_element"],
        )
        overrides["spatialdata_element"] = values["spatialdata_element"]
    if values.get("pm_on_cell_path") is not None:
        _set_nested(
            document.setdefault("inputs", {}),
            "pm_on_cell", {"path": values["pm_on_cell_path"]},
        )
        overrides["pm_on_cell_path"] = values["pm_on_cell_path"]
    return overrides


def _output_paths(
    config: ApplicationConfig,
    run_dir: str | Path,
) -> dict[str, Path]:
    output_dir = Path(run_dir)
    if config.svc_type == "sc-SVC":
        return {
            "spatial": output_dir / f"{config.output_name}_spatial.h5ad",
            "expression": output_dir / f"{config.output_name}_expr.h5ad",
        }
    return {"svc": output_dir / f"{config.output_name}.h5ad"}


def _engine_overrides(
    config: ApplicationConfig,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    runtime = {"seed": config.seed} if config.seed is not None else {}
    io = {
        "st_path": str(config.st_path),
        "sc_ref_path": str(config.reference_path),
        "pm_on_cell_path": str(config.pm_on_cell_path) if config.pm_on_cell_path else "",
        "output_root": str(config.output_dir),
        "sample_name": config.output_name,
        "patient_key": "",
        "save_outputs": False,
        "input_format": config.st_format,
        "data_root": "",
        "st_file": "",
        "sc_ref_file": "",
    }
    if config.st_format in {"spatialdata", "auto"}:
        io["spatialdata_path"] = str(config.st_path)
        if config.spatialdata_table is not None:
            io["spatialdata_table"] = config.spatialdata_table
        if config.spatialdata_element is not None:
            io["spatialdata_spatial_element"] = config.spatialdata_element

    algorithm: dict[str, Any] = {"columns": {"cell_type_col": config.broad_column}}
    if config.subtype_column is not None:
        algorithm["columns"]["sub_cell_type_col"] = config.subtype_column
    if config.select_cell_type is not None:
        algorithm["sc"] = {"select_ct": config.select_cell_type}
    if config.local_refinement_strength is not None:
        algorithm["local_refinement"] = {"strength": config.local_refinement_strength}
    if config.ot_method is not None:
        algorithm["ot"] = {
            "ga": {"solver": config.ot_method},
            "lr": {"solver": config.ot_method},
        }
    return runtime, io, algorithm


def _application_metadata(
    config: ApplicationConfig,
    *,
    cli_overrides: Mapping[str, Any],
    output_paths: Mapping[str, Path],
    dry_run: bool,
) -> dict[str, Any]:
    return {
        "source_path": config.source_path,
        "source_sha256": config.config_sha256,
        "cli_overrides": dict(cli_overrides),
        "declared_root": config.declared_root,
        "resolved_root": str(config.resolved_root),
        "cwd": str(config.cwd),
        "resolved_inputs": config.resolved_inputs,
        "output_name": config.output_name,
        "output_paths": {key: str(path) for key, path in output_paths.items()},
        "effective_action": "preflight" if dry_run else "run",
    }


def _pipeline_evidence(svc) -> dict[str, Any]:
    evidence = svc.summary()
    route = svc.provenance.get("route", {})
    evidence.update(
        profile=svc.provenance.get("profile"),
        task=route.get("task"),
        strategy=route.get("strategy"),
        route=route,
    )
    return evidence


def _cleanup_output_files(paths: list[Path], ctx) -> None:
    cleanup_errors = []
    for path in paths:
        try:
            path.unlink(missing_ok=True)
        except OSError as exc:
            cleanup_errors.append(
                {
                    "path": str(path),
                    "type": type(exc).__name__,
                    "message": str(exc),
                }
            )
    if cleanup_errors:
        ctx.provenance.setdefault("output_cleanup_errors", []).extend(
            cleanup_errors
        )
        ctx.logger.error(
            "[framework] failed to clean %d direct output(s)",
            len(cleanup_errors),
        )


def _write_outputs(
    config: ApplicationConfig,
    output_paths: Mapping[str, Path],
    ctx,
) -> dict[str, Any]:
    """Write final artifacts directly into this run's only output directory."""
    outputs = dict(ctx.svc.artifacts.get("outputs", {})) if ctx.svc else {}
    required = {
        "sp-SVC": {"svc": "sp_svc"},
        "sc-SVC": {"spatial": "sc_svc_spatial", "expression": "sc_svc_expr"},
        "sc-SVC-sr": {"svc": "sc_svc_dec"},
    }[config.svc_type]
    missing = [key for key in required.values() if key not in outputs]
    if missing:
        raise RuntimeError(
            f"{config.svc_type} did not produce required output(s): "
            f"{', '.join(missing)}"
        )

    run_dir = Path(ctx.run_dir).resolve()
    targets = {role: Path(path).resolve() for role, path in output_paths.items()}
    outside = [str(path) for path in targets.values() if path.parent != run_dir]
    if outside:
        raise RuntimeError(
            "Application outputs must be direct children of the run directory: "
            f"run_dir={run_dir}; invalid={outside}"
        )
    if len(set(targets.values())) != len(targets):
        raise RuntimeError("Application output roles must use distinct file paths")

    written: dict[str, Any] = {}
    created: list[Path] = []
    metadata = {
        "svc_type": config.svc_type,
        "output_name": config.output_name,
        "profile": ctx.profile,
        "run_manifest": str(Path(ctx.run_dir) / "provenance.json"),
        "selected_cell_type": config.select_cell_type,
        "ot": ctx.merged_config.get("ot"),
    }
    try:
        for role, artifact_key in required.items():
            target = targets[role]
            if target.exists():
                raise FileExistsError(
                    f"Refusing to overwrite an existing run output: {target}"
                )
            created.append(target)
            adata = outputs[artifact_key].copy()
            adata.uns["revise_reconstruction"] = dict(metadata, output_role=role)
            adata.write_h5ad(target)
            written[role] = adata
    except BaseException:
        _cleanup_output_files(created, ctx)
        raise
    return written


def _persist_run_outputs(
    config: ApplicationConfig,
    output_paths: Mapping[str, Path],
    ctx,
    returned_outputs: dict[str, Any],
) -> None:
    """Attach direct-write outputs to the run lifecycle for failure cleanup."""
    had_results = "results" in ctx.provenance
    previous_results = copy.deepcopy(ctx.provenance.get("results"))
    svc_had_results = bool(ctx.svc and "results" in ctx.svc.provenance)
    previous_svc_results = (
        copy.deepcopy(ctx.svc.provenance.get("results")) if ctx.svc else None
    )
    created: list[Path] = []
    output_artifacts = {
        (f"output:{role}", str(Path(path)))
        for role, path in output_paths.items()
    }

    def cleanup() -> None:
        _cleanup_output_files(created, ctx)
        returned_outputs.clear()
        ctx.artifact_records[:] = [
            record
            for record in ctx.artifact_records
            if (record.get("role"), record.get("path")) not in output_artifacts
        ]
        if had_results:
            ctx.provenance["results"] = previous_results
        else:
            ctx.provenance.pop("results", None)
        if ctx.svc is not None:
            if svc_had_results:
                ctx.svc.provenance["results"] = previous_svc_results
            else:
                ctx.svc.provenance.pop("results", None)

    ctx.register_output_failure_cleanup(cleanup)
    try:
        written = _write_outputs(config, output_paths, ctx)
        created.extend(Path(path) for path in output_paths.values())
        results = {
            role: {
                "path": str(path),
                "filename": Path(path).name,
                "svc_type": config.svc_type,
            }
            for role, path in output_paths.items()
        }
        ctx.provenance["results"] = results
        if ctx.svc is not None:
            ctx.svc.provenance["results"] = copy.deepcopy(results)
        for role, path in output_paths.items():
            ctx.record_artifact(completed_artifact(f"output:{role}", path))
        returned_outputs.update(written)
    except BaseException:
        ctx.cleanup_failed_outputs()
        raise


def run_application(
    config: str | Path,
    *,
    svc_type: str | None = None,
    root_dir: str | None = None,
    st_path: str | None = None,
    st_format: str | None = None,
    sc_ref_path: str | None = None,
    spatialdata_table: str | None = None,
    spatialdata_element: str | None = None,
    pm_on_cell_path: str | None = None,
    ot_method: str | None = None,
    cell_type_col: str | None = None,
    sub_cell_type_col: str | None = None,
    select_cell_type: str | None = None,
    local_refinement_strength: float | None = None,
    output_dir: str | None = None,
    output_name: str | None = None,
    seed: int | None = None,
    dry_run: bool = False,
) -> ApplicationExecution:
    source, document = load_application_yaml(config)
    cli_overrides = _apply_overrides(
        document,
        {
            "svc_type": svc_type,
            "root_dir": root_dir,
            "st_path": st_path,
            "st_format": st_format,
            "sc_ref_path": sc_ref_path,
            "spatialdata_table": spatialdata_table,
            "spatialdata_element": spatialdata_element,
            "pm_on_cell_path": pm_on_cell_path,
            "ot_method": ot_method,
            "cell_type_col": cell_type_col,
            "sub_cell_type_col": sub_cell_type_col,
            "select_cell_type": select_cell_type,
            "local_refinement_strength": local_refinement_strength,
            "output_dir": output_dir,
            "output_name": output_name,
            "seed": seed,
        },
    )
    effective = compile_application_config(document, source=source)
    runtime, io, algorithm = _engine_overrides(effective)
    run_dir = build_run_dir(
        output_root=str(effective.output_dir),
        sample_name=effective.output_name,
        route_key=f"application:{effective.svc_type}",
        io_cfg=io,
        mode="application",
    ).resolve()
    output_paths = _output_paths(effective, run_dir)
    metadata = _application_metadata(
        effective,
        cli_overrides=cli_overrides,
        output_paths=output_paths,
        dry_run=dry_run,
    )
    returned_outputs: dict[str, Any] = {}

    def finalize(ctx):
        if not dry_run:
            _persist_run_outputs(effective, output_paths, ctx, returned_outputs)

    svc = REVISEPipeline().run(
        svc_type=effective.svc_type,
        cf=None,
        runtime_overrides=runtime,
        io_overrides=io,
        algorithm_overrides=algorithm,
        dry_run=dry_run,
        finalize_callback=finalize,
        application_config_metadata=metadata,
        run_directory=run_dir,
    )
    pipeline = _pipeline_evidence(svc)
    if dry_run:
        return ApplicationExecution(
            status="preflight_passed",
            svc_type=effective.svc_type,
            output_paths=output_paths,
            pipeline=pipeline,
            preflight=Path(svc.provenance["run_dir"]) / "preflight.json",
            application_config=metadata,
        )
    return ApplicationExecution(
        status="succeeded",
        svc_type=effective.svc_type,
        output_paths=output_paths,
        pipeline=pipeline,
        results=returned_outputs,
        application_config=metadata,
    )


def _execution_payload(execution: ApplicationExecution) -> dict[str, Any]:
    payload = {
        "status": execution.status,
        "svc_type": execution.svc_type,
        "outputs": {key: str(path) for key, path in execution.output_paths.items()},
        "application_config": execution.application_config,
        "pipeline": execution.pipeline,
    }
    if execution.preflight is not None:
        payload["preflight"] = str(execution.preflight)
    if execution.results is not None:
        payload["shapes"] = {key: list(value.shape) for key, value in execution.results.items()}
    return payload


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        with redirect_stdout(sys.stderr):
            execution = run_application(
                args.config,
                svc_type=args.svc_type,
                root_dir=args.root_dir,
                st_path=args.st_path,
                st_format=args.st_format,
                sc_ref_path=args.sc_ref_path,
                spatialdata_table=args.spatialdata_table,
                spatialdata_element=args.spatialdata_element,
                pm_on_cell_path=args.pm_on_cell_path,
                ot_method=args.ot_method,
                cell_type_col=args.cell_type_col,
                sub_cell_type_col=args.sub_cell_type_col,
                select_cell_type=args.select_ct,
                local_refinement_strength=args.local_refinement_strength,
                output_dir=args.output_dir,
                output_name=args.output_name,
                seed=args.seed,
                dry_run=args.dry_run,
            )
    except ApplicationConfigError as exc:
        parser.error(str(exc))
    except KeyboardInterrupt:
        raise
    except Exception as exc:
        context = getattr(exc, "_revise_failure_context", None)
        detail = f"revise-reconstruct: error: {exc}"
        if context:
            detail += f"\nrun_dir: {context['run_dir']}\nmanifest: {context['manifest_path']}"
        parser.exit(1, detail + "\n")
    print(json.dumps(_execution_payload(execution), indent=2, ensure_ascii=False))


__all__ = ["ApplicationConfigError", "ApplicationExecution", "build_parser", "main", "parse_args", "run_application"]


if __name__ == "__main__":
    main()
