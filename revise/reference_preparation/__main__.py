"""Command-line entry for reference preparation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from .api import prepare_reference


def _host_runner(reference_path: Path, reconstruction_config: Path):
    from .host import run_global_anchoring

    return run_global_anchoring(reference_path, reconstruction_config)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Prepare a REVISE single-cell reference")
    parser.add_argument("--config", required=True, help="Preparation YAML path")
    args = parser.parse_args(argv)
    result = prepare_reference(args.config, ga_runner=_host_runner)
    print(json.dumps({
        "reference_path": str(result.reference_path),
        "reference_config_path": str(result.reference_config_path),
        "report_path": str(result.report_path),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
