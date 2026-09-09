"""Command-line adapter for Sim2Real-ST pseudo-spot preparation."""

from __future__ import annotations

import argparse
from typing import Sequence

from revise.preprocess.sim2real_pseudospot.workflow import (
    build_real_pseudospots,
    propose_regions,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Prepare real Xenium pseudo-spots for Sim2Real-ST."
    )
    stages = parser.add_subparsers(dest="stage", required=True)
    for stage in ("propose", "build"):
        command = stages.add_parser(stage)
        command.add_argument("--config", required=True)
        command.add_argument("--sample", required=True)
        if stage == "build":
            command.add_argument("--confirmation", required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.stage == "propose":
        result = propose_regions(args.config, args.sample)
        print(result.image_path)
        print(result.table_path)
        print(result.proposal_path)
        return 0
    result = build_real_pseudospots(args.config, args.sample, args.confirmation)
    print(result.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
