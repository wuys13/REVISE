#!/usr/bin/env python3
from __future__ import annotations

import argparse
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tools.release_manifest import sha256_file


def _artifact(path: Path, role: str) -> dict:
    return {
        "role": role,
        "filename": path.name,
        "size": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def _git(value: str) -> str:
    result = subprocess.run(
        ["git", "rev-parse", value],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode:
        raise RuntimeError(result.stderr.strip())
    return result.stdout.strip()


def _release_source() -> dict | None:
    expected_commit = os.environ.get("REVISE_SOURCE_COMMIT")
    if not expected_commit:
        return None
    actual_commit = _git("HEAD")
    if actual_commit != expected_commit:
        raise RuntimeError(
            f"source commit mismatch: expected {expected_commit}, got {actual_commit}"
        )
    status = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=True,
    ).stdout
    if status:
        raise RuntimeError("release source checkout must be clean")
    return {"commit": actual_commit, "tree": _git("HEAD^{tree}")}


def build(output_dir: Path, report_path: Path) -> dict:
    output_dir = output_dir.resolve()
    report_path = report_path.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    if any(output_dir.iterdir()):
        raise RuntimeError(f"output directory must be empty: {output_dir}")
    expected_tools = {"setuptools": "80.9.0", "wheel": "0.45.1"}
    actual_tools = {name: version(name) for name in expected_tools}
    if actual_tools != expected_tools:
        raise RuntimeError(
            f"build toolchain mismatch: expected {expected_tools}, got {actual_tools}"
        )
    source = _release_source()

    previous = Path.cwd()
    try:
        os.chdir(ROOT)
        from setuptools.build_meta import build_sdist, build_wheel

        sdist_name = build_sdist(str(output_dir))
        wheel_name = build_wheel(str(output_dir))
    finally:
        os.chdir(previous)

    sdist = output_dir / sdist_name
    wheel = output_dir / wheel_name
    if sorted(output_dir.iterdir()) != sorted([sdist, wheel]):
        raise RuntimeError("build must produce exactly one wheel and one sdist")
    payload = {
        "schema_version": 1,
        "build": {
            "python": {
                "implementation": platform.python_implementation(),
                "version": platform.python_version(),
            },
            "backend": "setuptools.build_meta",
            "tools": {
                "setuptools": actual_tools["setuptools"],
                "wheel": actual_tools["wheel"],
            },
            "isolation": False,
        },
        "artifacts": [_artifact(wheel, "wheel"), _artifact(sdist, "sdist")],
    }
    if source is not None:
        payload["source"] = source
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()
    build(args.output_dir, args.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
