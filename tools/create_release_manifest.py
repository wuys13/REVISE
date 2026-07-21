#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import subprocess
import sys

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - exercised by the Python 3.10 CI job
    import tomli as tomllib

import jsonschema


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tools.release_manifest import sha256_file, validate_declared_files


def _git(repo: Path, value: str) -> str:
    result = subprocess.run(
        ["git", "rev-parse", value],
        cwd=repo,
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode:
        raise RuntimeError(result.stderr.strip())
    return result.stdout.strip()


def _git_status(repo: Path) -> str:
    result = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=repo,
        text=True,
        capture_output=True,
        check=True,
    )
    return result.stdout


def _parse_gate(value: str) -> tuple[str, Path]:
    gate_id, separator, path = value.partition("=")
    if not separator or not gate_id or not path:
        raise argparse.ArgumentTypeError("gate must use ID=PATH")
    return gate_id, Path(path)


def create_manifest(args: argparse.Namespace) -> dict:
    source_root = args.source_root.resolve()
    assets_root = args.assets_dir.resolve()
    output = args.output.resolve()
    if output.is_relative_to(source_root):
        raise ValueError("detached release manifest must be written outside the source tree")
    if output.parent != assets_root or output.name != "release-manifest.json":
        raise ValueError(
            "detached release manifest must use release-manifest.json in the assets directory"
        )

    template = json.loads(args.template.read_text())
    schema = json.loads(args.schema.read_text())
    gate_schema = json.loads(args.gate_schema.read_text())
    build_report = json.loads(args.build_report.read_text())
    reserved_paths = {args.build_report.resolve()}
    reserved_paths.update(path.resolve() for _gate_id, path in args.gate)
    reserved_paths.update(
        assets_root / item["filename"] for item in build_report["artifacts"]
    )
    if output in reserved_paths:
        raise ValueError("release manifest output collides with declared evidence")
    commit = _git(source_root, f"{args.tag}^{{commit}}")
    tree = _git(source_root, f"{commit}^{{tree}}")
    if _git(source_root, "HEAD") != commit:
        raise ValueError("checked-out HEAD does not match the release tag commit")
    if _git_status(source_root):
        raise ValueError("release source checkout must be clean")
    if build_report.get("source") != {"commit": commit, "tree": tree}:
        raise ValueError("build report source does not match the release tag")
    if args.tag.removeprefix("v") != template["release"]["version"]:
        raise ValueError("release tag does not match template version")

    manifest = copy.deepcopy(template)
    pyproject = tomllib.loads((source_root / "pyproject.toml").read_text())
    if pyproject["build-system"]["build-backend"] != manifest["build"]["backend"]:
        raise ValueError("manifest build backend does not match pyproject")
    available_extras = sorted(pyproject["project"]["optional-dependencies"])
    if manifest["support"]["available_extras"] != available_extras:
        raise ValueError("manifest available extras do not match pyproject")
    if not set(manifest["support"]["tested_extras"]) <= set(available_extras):
        raise ValueError("tested extras must be available extras")
    if set(manifest["support"]["tested_extras"]) - {"tacco"}:
        raise ValueError("tested extras require a corresponding release gate")
    if manifest["docs"]["version"] != manifest["release"]["version"]:
        raise ValueError("docs version does not match release version")
    if manifest["docs"]["channel"] != manifest["release"]["channel"]:
        raise ValueError("docs channel does not match release channel")
    if len(manifest["support"]["tested_platforms"]) != 1:
        raise ValueError("tested platforms must be derived from the one CI matrix")
    manifest["release"]["tag"] = args.tag
    manifest["source"].update({"commit": commit, "tree": tree})
    manifest["artifacts"] = build_report["artifacts"]
    manifest["build"].update(build_report["build"])
    manifest["build"]["workflow_run"] = args.run_url
    manifest["build"]["constraints"] = [
        {
            "python": python,
            "path": f"constraints/python-{python}.txt",
            "sha256": sha256_file(
                source_root / "constraints" / f"python-{python}.txt"
            ),
        }
        for python in ("3.10", "3.11")
    ]
    gates = []
    gate_reports = []
    for gate_id, report_path in args.gate:
        report_path = report_path.resolve()
        try:
            relative_report = report_path.relative_to(assets_root).as_posix()
        except ValueError as exc:
            raise ValueError("gate reports must be inside the release assets directory") from exc
        report = json.loads(report_path.read_text())
        jsonschema.validate(
            report,
            gate_schema,
            format_checker=jsonschema.FormatChecker(),
        )
        if report["gate"] != gate_id:
            raise ValueError(f"gate report ID mismatch: {gate_id}")
        if report["head_sha"] != commit:
            raise ValueError(f"gate report head SHA mismatch: {gate_id}")
        if report["run_id"] != args.run_id:
            raise ValueError(f"gate report run ID mismatch: {gate_id}")
        if report["run_attempt"] != args.run_attempt:
            raise ValueError(f"gate report run attempt mismatch: {gate_id}")
        if report["run_url"] != args.run_url:
            raise ValueError(f"gate report run URL mismatch: {gate_id}")
        gate_reports.append(report)
        gates.append(
            {
                "id": report["gate"],
                "conclusion": report["conclusion"],
                "run_id": report["run_id"],
                "run_attempt": report["run_attempt"],
                "head_sha": report["head_sha"],
                "report": relative_report,
                "sha256": sha256_file(report_path),
            }
        )
    by_id = {report["gate"]: report for report in gate_reports}
    for gate_id in ("unit", "installed_cli"):
        if by_id.get(gate_id, {}).get("tested_python") != ["3.10", "3.11"]:
            raise ValueError(f"{gate_id} report must cover Python 3.10 and 3.11")
    platform_reports = [
        (report["runner"], report["os"], report["architecture"])
        for report in gate_reports
    ]
    expected_platform = manifest["support"]["tested_platforms"][0]
    expected_platform_tuple = (
        expected_platform["runner"],
        expected_platform["os"],
        expected_platform["architecture"],
    )
    if any(platform != expected_platform_tuple for platform in platform_reports):
        raise ValueError("gate report platform does not match tested support")
    manifest["ci"] = {"head_sha": commit, "gates": gates}
    manifest["docs"].update(
        {"source_commit": commit, "workflow_run": args.docs_run_url or args.run_url}
    )

    jsonschema.Draft202012Validator.check_schema(schema)
    jsonschema.validate(
        manifest,
        schema,
        format_checker=jsonschema.FormatChecker(),
    )
    validate_declared_files(
        manifest,
        source_root=source_root,
        assets_root=assets_root,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, default=ROOT)
    parser.add_argument(
        "--template",
        type=Path,
        default=ROOT / "release/0.1.0rc1/release-manifest.template.json",
    )
    parser.add_argument(
        "--schema",
        type=Path,
        default=ROOT / "release/0.1.0rc1/release-manifest.schema.json",
    )
    parser.add_argument(
        "--gate-schema",
        type=Path,
        default=ROOT / "release/0.1.0rc1/gate-report.schema.json",
    )
    parser.add_argument("--assets-dir", required=True, type=Path)
    parser.add_argument("--build-report", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--run-id", required=True, type=int)
    parser.add_argument("--run-attempt", required=True, type=int)
    parser.add_argument("--run-url", required=True)
    parser.add_argument("--docs-run-url")
    parser.add_argument("--gate", action="append", required=True, type=_parse_gate)
    args = parser.parse_args()
    create_manifest(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
