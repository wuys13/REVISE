"""Command-line surface for the Application entrypoint."""

from __future__ import annotations

import argparse
from contextlib import redirect_stdout
import json
import sys

from revise._version import __version__

from .config import ApplicationConfigError


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
    parser.add_argument(
        "--version",
        action="version",
        version=f"revise-reconstruct {__version__}",
    )

    application = parser.add_argument_group("Application")
    application.add_argument("--config", required=True, help="application YAML baseline")
    application.add_argument(
        "--svc-type",
        choices=("sp-SVC", "sc-SVC", "sc-SVC-sr"),
        default=None,
        help="override application.svc_type",
    )

    inputs = parser.add_argument_group("Inputs")
    inputs.add_argument("--root-dir", default=None, help="override paths.root_dir")
    inputs.add_argument("--st-path", default=None, help="override inputs.st.path")
    inputs.add_argument(
        "--st-format", choices=("h5ad", "spatialdata", "auto"), default=None
    )
    inputs.add_argument("--sc-ref-path", default=None, help="override inputs.reference.path")
    inputs.add_argument("--spatialdata-table", default=None)
    inputs.add_argument("--spatialdata-element", default=None)
    inputs.add_argument("--pm-on-cell-path", default=None, help="only for sc-SVC-sr")

    ot = parser.add_argument_group("Shared OT")
    ot.add_argument(
        "--ot-method",
        choices=("pot", "tacco"),
        default=None,
        help="set the GA and LR solver together",
    )

    ga = parser.add_argument_group("Global Anchoring")
    ga.add_argument(
        "--cell-type-col",
        default=None,
        help="override global_anchoring.broad_column",
    )

    lr = parser.add_argument_group("Local Refinement")
    lr.add_argument("--sub-cell-type-col", default=None, help="only for sc-SVC")
    lr.add_argument("--select-ct", default=None, help="only for sc-SVC")
    lr.add_argument(
        "--local-refinement-strength",
        type=float,
        default=None,
        help="only for sp-SVC and sc-SVC-sr",
    )

    output = parser.add_argument_group("Output")
    output.add_argument("--output-dir", default=None)
    output.add_argument("--output-name", default=None)

    execution = parser.add_argument_group("Execution")
    execution.add_argument("--seed", type=int, default=None)
    execution.add_argument(
        "--dry-run",
        action="store_true",
        help="preflight inputs and route without writing H5AD",
    )
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_parser().parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    import reconstruct

    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        with redirect_stdout(sys.stderr):
            _, report = reconstruct._execute_application(
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
    print(json.dumps(report, indent=2, ensure_ascii=False))


__all__ = ["build_parser", "main", "parse_args"]
