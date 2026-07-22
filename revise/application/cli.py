"""Command-line interface for application reconstruction."""

from __future__ import annotations

import argparse
import json
import sys
from contextlib import redirect_stdout
from pathlib import Path

from revise._version import __version__

from .service import APPLICATION_ROUTES


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
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
    parser.add_argument("--ot-method", choices=("pot", "tacco"), default=None)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--patient-key", default="Patient")
    parser.add_argument("--select-ct", default="all")
    parser.add_argument("--cell-type-col", default=None)
    parser.add_argument("--sub-cell-type-col", default=None)
    parser.add_argument("--sc-mapping", choices=("mean", "random"), default="mean")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate resolved inputs and dependencies without reconstruction",
    )
    parser.add_argument(
        "--set",
        dest="set_overrides",
        action="append",
        default=[],
        metavar="KEY=VALUE",
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
        payload = {
            "svc_type": args.svc_type,
            "output": str(output_path),
            "shape": list(result.shape),
            "pipeline": pipeline_summary,
        }
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
