#!/usr/bin/env python3
"""Validate the committed boundary of the clean REVISE repository."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path, PurePosixPath


MAX_TRACKED_BYTES = 5 * 1024 * 1024
REQUIRED_PATHS = {
    ".github/workflows/ci.yml",
    ".gitignore",
    "LICENSE",
    "MANIFEST.in",
    "MIGRATION.md",
    "README.md",
    "legacy-assets.json",
    "pyproject.toml",
    "revise/__init__.py",
    "tools/check_repository.py",
}
ALLOWED_TOP_LEVEL = {
    ".github",
    ".gitignore",
    ".readthedocs.yaml",
    "LICENSE",
    "MANIFEST.in",
    "MIGRATION.md",
    "README.md",
    "application_sc_SVC_recon.py",
    "application_sc_SVC_recon.sh",
    "application_sp_SVC_recon.py",
    "application_sp_SVC_recon.sh",
    "benchmark_launcher.py",
    "benchmark_main.py",
    "benchmark_main.sh",
    "constraints",
    "docs",
    "legacy-assets.json",
    "png",
    "pyproject.toml",
    "reconstruct.py",
    "release",
    "revise",
    "scripts",
    "tests",
    "tools",
}
FORBIDDEN_PREFIXES = (
    "build/",
    "dist/",
    "docs/benchmark/",
    "docs/case/",
    "docs/plans/",
    "docs/superpowers/",
    "logo/",
    "release/legacy/",
    "reproduce/",
)
FORBIDDEN_EXACT = {
    "application_sc_SVC_analysis_case.ipynb",
    "benchmark_identity_free_ablation.py",
    "release/0.1.0rc1/gallery-manifest.source.json",
    "release/0.1.0rc1/status.json",
    "release/clean-export-allowlist.txt",
    "tests/test_gallery_contract.py",
    "tests/test_legacy_release_index.py",
    "tests/test_release_state.py",
    "tests/test_repository_export.py",
    "tools/check_release_state.py",
    "tools/export_clean_repository.py",
    "tools/gallery_bundle.py",
}


def _git(root: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["git", *args],
        cwd=root,
        capture_output=True,
        check=False,
    )


def _committed_tree(root: Path) -> tuple[list[tuple[str, str, str]], list[str]]:
    count = _git(root, "rev-list", "--count", "HEAD")
    if count.returncode:
        return [], ["repository must have one committed root"]
    errors = []
    shallow = _git(root, "rev-parse", "--is-shallow-repository")
    if shallow.returncode or shallow.stdout.strip() != b"false":
        errors.append("shallow repositories are forbidden")
    if count.stdout.strip() != b"1":
        errors.append("clean repository must have exactly one commit")
    commit = _git(root, "cat-file", "commit", "HEAD")
    headers = commit.stdout.partition(b"\n\n")[0].splitlines()
    if commit.returncode or any(line.startswith(b"parent ") for line in headers):
        errors.append("clean repository commit must have no parent")
    if _git(root, "diff", "--quiet", "HEAD", "--").returncode:
        errors.append("clean repository has uncommitted tracked changes")
    tree = _git(root, "ls-tree", "-r", "-z", "HEAD")
    if tree.returncode:
        return [], errors + ["repository root must be a Git working tree"]
    entries = []
    for record in tree.stdout.split(b"\0"):
        if not record:
            continue
        header, separator, raw_path = record.partition(b"\t")
        if not separator:
            return [], errors + ["invalid Git tree entry"]
        mode, object_type, object_id = header.decode("ascii").split()
        if object_type != "blob":
            return [], errors + ["clean repository may contain only blobs"]
        entries.append((raw_path.decode("utf-8"), mode, object_id))
    reachable = _git(root, "rev-list", "--objects", "HEAD")
    all_objects = _git(
        root,
        "cat-file",
        "--batch-all-objects",
        "--batch-check=%(objectname)",
    )
    reachable_ids = {line.split()[0] for line in reachable.stdout.splitlines()}
    all_ids = set(all_objects.stdout.splitlines())
    if (
        reachable.returncode
        or all_objects.returncode
        or reachable_ids != all_ids
    ):
        errors.append("clean repository contains objects outside the HEAD closure")
    return entries, errors


def _blob(root: Path, object_id: str) -> bytes:
    result = _git(root, "cat-file", "blob", object_id)
    if result.returncode:
        raise ValueError("tracked Git blob is unreadable")
    return result.stdout


def _validate_index(root: Path, content: bytes) -> list[str]:
    try:
        payload = json.loads(content.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError):
        return ["legacy asset index is invalid"]
    source = payload.get("source") if isinstance(payload, dict) else None
    if (
        set(payload) != {"schema_version", "source", "excluded_assets"}
        or payload["schema_version"] != 1
        or not isinstance(source, dict)
        or set(source) != {"commit", "repository", "repository_name"}
        or source["repository_name"] != "REVISE-legacy"
        or not re.fullmatch(r"[0-9a-f]{40}", source["commit"])
        or not re.fullmatch(
            r"https://github\.com/[^/]+/REVISE-legacy(?:\.git)?",
            source["repository"],
        )
        or not isinstance(payload["excluded_assets"], list)
    ):
        return ["legacy asset index is invalid"]
    seen = set()
    for item in payload["excluded_assets"]:
        path = item.get("path") if isinstance(item, dict) else None
        pure = PurePosixPath(path) if isinstance(path, str) else None
        expected_retrieval = (
            f"git -C REVISE-legacy show {source['commit']}:{path}"
            if pure is not None
            else None
        )
        if (
            not isinstance(item, dict)
            or set(item)
            != {"category", "kind", "path", "reason", "retrieval", "sha256", "size"}
            or not path
            or pure is None
            or pure.is_absolute()
            or pure.as_posix() != path
            or ".." in pure.parts
            or path in seen
            or item["kind"] not in {"file", "symlink"}
            or not isinstance(item["category"], str)
            or not item["category"]
            or not isinstance(item["reason"], str)
            or not item["reason"]
            or item["retrieval"] != expected_retrieval
            or not re.fullmatch(r"[0-9a-f]{64}", item["sha256"])
            or not isinstance(item["size"], int)
            or item["size"] < 0
        ):
            return ["legacy asset index is invalid"]
        seen.add(path)
    object_probe = _git(root, "cat-file", "-e", f"{source['commit']}^{{commit}}")
    if object_probe.returncode == 0:
        return ["clean repository unexpectedly contains the legacy source commit"]
    return []


def validate_repository(root: Path) -> list[str]:
    root = root.resolve()
    tracked, errors = _committed_tree(root)
    paths = {path for path, _mode, _object_id in tracked}
    for required in sorted(REQUIRED_PATHS - paths):
        errors.append(f"required tracked path is missing: {required}")
    contents = {}
    for path, mode, object_id in tracked:
        if mode not in {"100644", "100755"}:
            if mode == "120000":
                errors.append(f"tracked symlink is forbidden: {path}")
            else:
                errors.append(f"unsupported tracked mode {mode}: {path}")
            continue
        try:
            content = _blob(root, object_id)
        except ValueError:
            errors.append(f"tracked path is not readable: {path}")
            continue
        contents[path] = content
        if len(content) > MAX_TRACKED_BYTES:
            errors.append(f"tracked file exceeds 5 MiB: {path}")
        if (
            path.casefold().endswith(".ipynb")
            or path in FORBIDDEN_EXACT
            or path.startswith(FORBIDDEN_PREFIXES)
            or "/__pycache__/" in f"/{path}/"
            or path.endswith((".pyc", ".pyo"))
            or ".egg-info/" in path
        ):
            errors.append(f"forbidden tracked path: {path}")
        if path.split("/", 1)[0] not in ALLOWED_TOP_LEVEL:
            errors.append(f"unapproved top-level tracked path: {path}")
        if (
            path in {"MIGRATION.md", "README.md"}
            or path.startswith("docs/source/")
        ) and b"zenodo" in content.lower():
            errors.append(
                f"Zenodo migration is not part of the clean repository: {path}"
            )
    if "legacy-assets.json" in contents:
        errors.extend(_validate_index(root, contents["legacy-assets.json"]))
    else:
        errors.append("legacy asset index is invalid")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    args = parser.parse_args()
    errors = validate_repository(args.root)
    if errors:
        for error in errors:
            print(error, file=sys.stderr)
        return 1
    print(f"valid clean repository: {args.root.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
