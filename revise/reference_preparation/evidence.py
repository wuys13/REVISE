"""Resolve reference-only inputs with durable file identity evidence."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Any

from .config import ConfigError, load_reference_config


REPORT_KIND = "revise.reference_preparation.report"


def _state(path: Path) -> tuple[int, int, int, int]:
    stat = path.stat()
    return (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns)


def stable_file_identity(path: str | Path) -> dict[str, Any]:
    """Hash one regular file and reject changes observed while reading it."""
    source = Path(path).expanduser().resolve()
    digest = hashlib.sha256()
    with source.open("rb") as handle:
        before = os.fstat(handle.fileno())
        if not os.path.isfile(source):
            raise ValueError(f"Expected a regular file: {source}")
        while block := handle.read(1024 * 1024):
            digest.update(block)
        after = os.fstat(handle.fileno())
    observed = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
    completed = (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
    if observed != completed or _state(source) != completed:
        raise RuntimeError(f"File changed while hashing: {source}")
    return {
        "path": str(source),
        "sha256": digest.hexdigest(),
        "size_bytes": before.st_size,
    }


def _stable_json(path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    identity = stable_file_identity(path)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ConfigError(f"Invalid reference preparation report {path}: {error}") from error
    if stable_file_identity(path) != identity:
        raise RuntimeError(f"File changed while parsing: {path}")
    if not isinstance(payload, dict):
        raise ConfigError(f"Reference preparation report must be a mapping: {path}")
    return payload, identity


def _report_reference_path(value: object, report_path: Path) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise ConfigError("Preparation report selected_reference_path must be a non-empty string")
    path = Path(value).expanduser()
    return (path if path.is_absolute() else report_path.parent / path).resolve()


def resolve_reference_input(config_path: str | Path) -> tuple[Path, dict[str, Any]]:
    """Resolve a minimal reference YAML and verify its sibling report when present.

    The returned evidence is flat and deterministic so callers can embed it in
    run provenance without retaining parser or AnnData objects.
    """
    config = Path(config_path).expanduser().resolve()
    config_identity = stable_file_identity(config)
    resolved = load_reference_config(config).path.resolve()
    if stable_file_identity(config) != config_identity:
        raise RuntimeError(f"Reference configuration changed while resolving: {config}")
    reference_identity = stable_file_identity(resolved)
    report_path = config.parent / "report.json"

    evidence: dict[str, Any] = {
        "schema_version": 1,
        "origin": "reference_config_no_report",
        "verification": "unverified_no_report",
        "config_path": config_identity["path"],
        "config_sha256": config_identity["sha256"],
        "report_path": str(report_path.resolve()),
        "report_sha256": None,
        "reference_path": reference_identity["path"],
        "reference_sha256": reference_identity["sha256"],
        "reference_size_bytes": reference_identity["size_bytes"],
    }
    if not report_path.exists():
        return resolved, evidence

    report, report_identity = _stable_json(report_path)
    if report.get("report_kind") != REPORT_KIND or report.get("schema_version") != 1:
        raise ConfigError(f"Unrecognized reference preparation report: {report_path}")
    if report.get("status") != "selected":
        raise ConfigError(
            f"Reference preparation report is not successful: {report.get('status')!r}"
        )
    reported_path = _report_reference_path(report.get("selected_reference_path"), report_path)
    if reported_path != resolved:
        raise ConfigError(
            "Reference configuration path does not match report selected_reference_path"
        )
    reported_hash = report.get("selected_reference_sha256")
    if reported_hash != reference_identity["sha256"]:
        raise ConfigError("Reference file SHA-256 does not match preparation report")
    if report.get("mode") == "paired" and report.get("paired_output_sha256") != reported_hash:
        raise ConfigError("Paired output SHA-256 does not match selected reference SHA-256")

    evidence.update({
        "origin": REPORT_KIND,
        "verification": "verified_report",
        "report_sha256": report_identity["sha256"],
    })
    return resolved, evidence


def assert_reference_unchanged(evidence: dict[str, Any]) -> None:
    """Re-resolve an input and raise when any evidence field has changed."""
    required = {
        "config_path", "config_sha256", "report_path", "report_sha256",
        "reference_path", "reference_sha256", "reference_size_bytes",
    }
    missing = required - evidence.keys()
    if missing:
        raise ValueError(f"Reference evidence is missing fields: {sorted(missing)}")

    _, current = resolve_reference_input(evidence["config_path"])
    if current != evidence:
        changed = sorted(
            key for key in set(current) | set(evidence)
            if current.get(key) != evidence.get(key)
        )
        raise RuntimeError(
            f"Reference input changed after resolution; evidence fields: {changed}"
        )


__all__ = [
    "REPORT_KIND",
    "assert_reference_unchanged",
    "resolve_reference_input",
    "stable_file_identity",
]
