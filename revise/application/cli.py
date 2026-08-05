"""Canonical command-line interface for application reconstruction."""

from __future__ import annotations

import argparse
from contextlib import redirect_stdout
from dataclasses import dataclass
import json
from importlib import resources
from pathlib import Path
import sys

from revise._version import __version__

from .request import (
    ApplicationConfigError,
    load_application_request,
)


_OFFICIAL_TEMPLATES = frozenset(
    {
        "Xenium_T.yaml",
        "Xenium_Fib.yaml",
        "Xenium_Mono.yaml",
        "VisiumHD.yaml",
        "Visium.yaml",
    }
)
_MOVED_FLAGS = {
    "--svc-type": "application.svc_type",
    "--sample-name": "application.sample_name",
    "--st-file": "inputs.st.file",
    "--sc-ref-file": "inputs.reference.file",
    "--data-root": "inputs.data_root",
    "--output-root": "output.path",
    "--ot-method": "algorithm.ot_method",
    "--seed": "execution.seed",
    "--patient-key": "inputs.reference.patient_key",
    "--select-ct": "local_refinement.select_cell_type",
    "--cell-type-col": "global_anchoring.broad_column",
    "--sub-cell-type-col": "local_refinement.subtype_column",
    "--local-refinement-strength": "local_refinement.strength",
}
_LEGACY_LAYOUT_FLAGS = {"--st-file", "--sc-ref-file", "--data-root"}
_REMOVED_ASSIGNMENT_FLAGS = {
    "--local-refinement-guidance",
    "--local-refinement-compatibility-mode",
    "--posterior-mode",
    "--posterior-strict",
    "--posterior-key",
    "--posterior-beta",
    "--posterior-min-affinity",
    "--posterior-cost-strength",
}


@dataclass(frozen=True)
class _ApplicationConfigSource:
    label: str
    payload: bytes
    path: Path | None


def _local_source(path: Path) -> _ApplicationConfigSource:
    resolved = path.resolve()
    return _ApplicationConfigSource(
        label=str(resolved),
        payload=resolved.read_bytes(),
        path=resolved,
    )


def _canonical_template_name(value: str) -> str | None:
    path = Path(value)
    if path.is_absolute():
        return None
    parts = path.parts
    if len(parts) == 1 and parts[0] in _OFFICIAL_TEMPLATES:
        return parts[0]
    if len(parts) == 3 and parts[:2] == ("configs", "application"):
        if parts[2] in _OFFICIAL_TEMPLATES:
            return parts[2]
    return None


def _resolve_application_source(config_path: str | Path) -> _ApplicationConfigSource:
    """Resolve an explicit application file or an official packaged template."""
    value = str(config_path)
    path = Path(value).expanduser()
    if path.exists():
        if not path.is_file():
            raise ApplicationConfigError(
                f"Cannot read application config {path}: not a file"
            )
        try:
            return _local_source(path)
        except OSError as exc:
            raise ApplicationConfigError(
                f"Cannot read application config {path}: {exc}"
            ) from exc

    template_name = _canonical_template_name(value)
    if template_name is None:
        raise ApplicationConfigError(
            f"Cannot read application config {path}: file does not exist"
        )

    if Path(value).parts == (template_name,):
        mirror = Path.cwd() / "configs" / "application" / template_name
        if mirror.is_file():
            try:
                return _local_source(mirror)
            except OSError as exc:
                raise ApplicationConfigError(
                    f"Cannot read application config {mirror}: {exc}"
                ) from exc

    try:
        payload = (
            resources.files("revise.application")
            .joinpath("templates", template_name)
            .read_bytes()
        )
    except (FileNotFoundError, ModuleNotFoundError, OSError) as exc:
        raise ApplicationConfigError(
            f"Cannot read packaged application template {template_name}: {exc}"
        ) from exc
    return _ApplicationConfigSource(
        label=f"package:revise.application.templates/{template_name}",
        payload=payload,
        path=None,
    )


class _ApplicationArgumentParser(argparse.ArgumentParser):
    def parse_known_args(self, args=None, namespace=None):
        raw_args = sys.argv[1:] if args is None else list(args)
        for token in raw_args:
            flag = str(token).split("=", 1)[0]
            if flag in _MOVED_FLAGS:
                legacy_note = (
                    "; these fields require inputs.mode: legacy_layout"
                    if flag in _LEGACY_LAYOUT_FLAGS
                    else ""
                )
                self.error(
                    f"{flag} moved to application YAML field {_MOVED_FLAGS[flag]}"
                    f"{legacy_note}"
                )
            if flag == "--set":
                self.error(
                    "--set is not supported; the full engine config is internally "
                    "managed by the package"
                )
            if flag in _REMOVED_ASSIGNMENT_FLAGS:
                self.error(
                    f"{flag} is no longer supported; use local_refinement fields "
                    "in the application YAML"
                )
        return super().parse_known_args(raw_args, namespace)


def build_parser() -> argparse.ArgumentParser:
    parser = _ApplicationArgumentParser(
        description="Reconstruct one SVC from an application YAML",
        epilog=(
            "Choose the template that matches the ST data:\n"
            "  Xenium_T.yaml    segmented-cell T-cell data (sc-SVC)\n"
            "  Xenium_Fib.yaml  segmented-cell fibroblast data (sc-SVC)\n"
            "  Xenium_Mono.yaml segmented-cell mono/macro data (sc-SVC)\n"
            "  VisiumHD.yaml    high-resolution bins or pseudo-cells (sp-SVC)\n"
            "  Visium.yaml      multi-cell spots (sc-SVC-sr)\n\n"
            "Source checkout:\n"
            "  python reconstruct.py --config configs/application/VisiumHD.yaml --dry-run\n"
            "  python reconstruct.py --config configs/application/VisiumHD.yaml\n\n"
            "Installed command:\n"
            "  revise-reconstruct --config configs/application/VisiumHD.yaml --dry-run\n"
            "  revise-reconstruct --config configs/application/VisiumHD.yaml"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"revise-reconstruct {__version__}",
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to the selected application YAML",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Force preflight even when execution.action is run",
    )
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_parser().parse_args(argv)


def _execution_payload(execution) -> dict:
    """Convert the service result into the small JSON contract of the CLI."""
    if isinstance(execution, dict):
        return execution
    if execution.status == "ready":
        return {
            "status": execution.status,
            "svc_type": execution.svc_type,
            "preflight": str(execution.preflight),
            "pipeline": execution.pipeline,
        }
    if execution.svc_type == "sc-SVC":
        return {
            "svc_type": execution.svc_type,
            "outputs": {
                key: str(path) for key, path in execution.output_path.items()
            },
            "shapes": {
                key: list(adata.shape) for key, adata in execution.result.items()
            },
            "pipeline": execution.pipeline,
        }
    return {
        "svc_type": execution.svc_type,
        "output": str(execution.output_path),
        "shape": list(execution.result.shape),
        "pipeline": execution.pipeline,
    }


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        source = _resolve_application_source(args.config)
        request = load_application_request(
            source.path,
            dry_run=args.dry_run,
            _payload=source.payload,
            _source_path=source.label,
        )
    except ApplicationConfigError as exc:
        parser.error(str(exc))

    from .service import execute_application

    try:
        with redirect_stdout(sys.stderr):
            execution = execute_application(request)
    except KeyboardInterrupt:
        raise
    except Exception as exc:
        parser.exit(1, f"revise-reconstruct: error: {exc}\n")
    payload = _execution_payload(execution)
    print(json.dumps(payload, indent=2, ensure_ascii=False))


__all__ = ["ApplicationConfigError", "build_parser", "main", "parse_args"]


if __name__ == "__main__":
    main()
