from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _safe_relative(root: Path, value: str) -> Path:
    relative = PurePosixPath(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"unsafe relative path: {value}")
    path = root.joinpath(*relative.parts)
    if not path.is_file():
        raise ValueError(f"declared file does not exist: {value}")
    return path


def _check_file(path: Path, item: dict, label: str) -> None:
    if path.stat().st_size != item["size"]:
        raise ValueError(f"{label} size mismatch: {path.name}")
    if sha256_file(path) != item["sha256"]:
        raise ValueError(f"{label} hash mismatch: {path.name}")


def validate_declared_files(
    manifest: dict,
    *,
    source_root: Path,
    assets_root: Path,
) -> None:
    artifact_roles = [item["role"] for item in manifest["artifacts"]]
    if sorted(artifact_roles) != ["sdist", "wheel"]:
        raise ValueError("manifest must declare exactly one wheel and one sdist")
    artifact_names = [item["filename"] for item in manifest["artifacts"]]
    if len(artifact_names) != len(set(artifact_names)):
        raise ValueError("artifact filenames must be unique")
    for item in manifest["artifacts"]:
        expected_suffix = ".whl" if item["role"] == "wheel" else ".tar.gz"
        if not item["filename"].endswith(expected_suffix):
            raise ValueError(f"artifact role/filename mismatch: {item['role']}")
        path = _safe_relative(assets_root, item["filename"])
        _check_file(path, item, "artifact")

    constraint_pythons = [item["python"] for item in manifest["build"]["constraints"]]
    if sorted(constraint_pythons) != ["3.10", "3.11"]:
        raise ValueError("constraints must cover Python 3.10 and 3.11 exactly once")
    for item in manifest["build"]["constraints"]:
        path = _safe_relative(source_root, item["path"])
        if sha256_file(path) != item["sha256"]:
            raise ValueError(f"constraint hash mismatch: {item['path']}")

    if manifest["source"]["commit"] != manifest["ci"]["head_sha"]:
        raise ValueError("CI head SHA does not match source commit")
    if manifest["source"]["commit"] != manifest["docs"]["source_commit"]:
        raise ValueError("docs source commit does not match source commit")
    gate_ids = []
    for gate in manifest["ci"]["gates"]:
        if gate["head_sha"] != manifest["source"]["commit"]:
            raise ValueError(f"gate head SHA mismatch: {gate['id']}")
        report = _safe_relative(assets_root, gate["report"])
        if sha256_file(report) != gate["sha256"]:
            raise ValueError(f"gate report hash mismatch: {gate['id']}")
        report_payload = json.loads(report.read_text())
        expected = {
            "gate": gate["id"],
            "conclusion": gate["conclusion"],
            "head_sha": gate["head_sha"],
            "run_id": gate["run_id"],
            "run_attempt": gate["run_attempt"],
        }
        actual = {key: report_payload.get(key) for key in expected}
        if actual != expected:
            raise ValueError(f"gate report content mismatch: {gate['id']}")
        if report_payload.get("run_url") != manifest["build"]["workflow_run"]:
            raise ValueError(f"gate report run URL mismatch: {gate['id']}")
        gate_ids.append(gate["id"])
    if len(gate_ids) != len(set(gate_ids)):
        raise ValueError("CI gate IDs must be unique")

    limitation_ids = [item["id"] for item in manifest["limitations"]]
    if len(limitation_ids) != len(set(limitation_ids)):
        raise ValueError("limitation IDs must be unique")
    for limitation in manifest["limitations"]:
        _safe_relative(source_root, limitation["source"])
