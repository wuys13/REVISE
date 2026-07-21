#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path


GATE_IDS = (
    "unit",
    "integration",
    "installed_cli",
    "package",
    "docs",
    "repository",
    "tacco",
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--head-sha", required=True)
    parser.add_argument("--run-id", required=True, type=int)
    parser.add_argument("--run-attempt", required=True, type=int)
    parser.add_argument("--run-url", required=True)
    parser.add_argument("--gate", required=True, choices=GATE_IDS)
    parser.add_argument(
        "--tested-python",
        required=True,
        action="append",
        choices=("3.10", "3.11"),
    )
    parser.add_argument("--runner", default="ubuntu-latest")
    parser.add_argument("--os", default=os.environ.get("RUNNER_OS", "Linux"))
    parser.add_argument(
        "--architecture", default=os.environ.get("RUNNER_ARCH", "X64")
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": 1,
        "gate": args.gate,
        "conclusion": "success",
        "head_sha": args.head_sha,
        "run_id": args.run_id,
        "run_attempt": args.run_attempt,
        "run_url": args.run_url,
        "tested_python": args.tested_python,
        "runner": args.runner,
        "os": args.os,
        "architecture": args.architecture,
    }
    path = args.output_dir / f"gate-{args.gate}.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
