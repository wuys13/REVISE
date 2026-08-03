"""Command-line interface for application reconstruction."""

from __future__ import annotations

import argparse
import json
import sys
from contextlib import redirect_stdout
from pathlib import Path

from revise._version import __version__

from .service import APPLICATION_ROUTES

_REMOVED_ASSIGNMENT_FLAGS = frozenset(
    {
        "--local-refinement-guidance",
        "--local-refinement-compatibility-mode",
        "--posterior-mode",
        "--posterior-strict",
        "--posterior-key",
        "--posterior-beta",
        "--posterior-min-affinity",
        "--posterior-cost-strength",
    }
)
_MIGRATION_ERROR = (
    "Assignment guidance options were removed; "
    "use --local-refinement-strength"
)


class _ApplicationArgumentParser(argparse.ArgumentParser):
    def parse_known_args(self, args=None, namespace=None):
        raw_args = sys.argv[1:] if args is None else list(args)
        if any(
            str(token).split("=", 1)[0] in _REMOVED_ASSIGNMENT_FLAGS
            for token in raw_args
        ):
            self.error(_MIGRATION_ERROR)
        return super().parse_known_args(raw_args, namespace)


def build_parser() -> argparse.ArgumentParser:
    parser = _ApplicationArgumentParser(
        description="Reconstruct one SVC through revise-svc",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"revise-reconstruct {__version__}",
    )
    parser.add_argument(
        "--svc-type",
        required=True,
        choices=tuple(APPLICATION_ROUTES),
    )
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--st-file", required=True)
    parser.add_argument("--sc-ref-file", required=True)
    parser.add_argument("--data-root", required=True)
    parser.add_argument("--output-root", default="output/reconstruct")
    parser.add_argument("--config", default="revise/revise.yaml")
    parser.add_argument(
        "--ot-method",
        choices=("pot", "tacco"),
        default=None,
        help=(
            "Set both Global Anchoring and Local Refinement OT solvers; "
            "standard sc-SVC defaults to TACCO, while 'pot' explicitly "
            "selects a different algorithm"
        ),
    )
    parser.add_argument(
        "--local-refinement-strength",
        type=float,
        default=None,
        help=(
            "Override posterior-conditioned local OT strength for sp-SVC "
            "or sc-SVC-sr"
        ),
    )
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--patient-key", default="Patient")
    parser.add_argument("--select-ct", default="all")
    parser.add_argument("--cell-type-col", default=None)
    parser.add_argument("--sub-cell-type-col", default=None)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate resolved inputs and dependencies without reconstruction",
    )
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_parser().parse_args(argv)


def main() -> None:
    args = parse_args()
    from .service import _run_pipeline, reconstruct

    if args.dry_run:
        with redirect_stdout(sys.stderr):
            profile, _, svc = _run_pipeline(args, dry_run=True)
        run_dir = Path(svc.provenance["run_dir"])
        summary = svc.summary()
        summary.update(profile=profile, route=svc.provenance.get("route_key"))
        payload = {
            "status": "ready",
            "svc_type": args.svc_type,
            "preflight": str(run_dir / "preflight.json"),
            "pipeline": summary,
        }
    else:
        with redirect_stdout(sys.stderr):
            result, output_path, pipeline_summary = reconstruct(args)
        if args.svc_type == "sc-SVC":
            payload = {
                "svc_type": args.svc_type,
                "outputs": {key: str(path) for key, path in output_path.items()},
                "shapes": {key: list(adata.shape) for key, adata in result.items()},
                "pipeline": pipeline_summary,
            }
        else:
            payload = {
                "svc_type": args.svc_type,
                "output": str(output_path),
                "shape": list(result.shape),
                "pipeline": pipeline_summary,
            }
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
