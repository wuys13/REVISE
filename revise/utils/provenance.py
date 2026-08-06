from __future__ import annotations

import hashlib
import json
import os
import platform
import stat
import tempfile
from contextlib import contextmanager
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as pkg_version
from pathlib import Path
from typing import Any, Dict, Iterable


def hash_jsonable(payload: Any) -> str:
    serialized = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()


def _with_relative_path(
    record: Dict[str, Any], relative_path: str | None
) -> Dict[str, Any]:
    if relative_path is not None:
        record["path"] = relative_path
    return record


def _file_record(path: Path, *, relative_path: str | None = None) -> Dict[str, Any]:
    record = {"kind": "file", "sha256": sha256_file(path)}
    return _with_relative_path(record, relative_path)


def _path_content_record(
    path: Path,
    *,
    relative_path: str | None = None,
    ancestor_directories: frozenset[tuple[int, int]] = frozenset(),
) -> Dict[str, Any]:
    try:
        path_stat = os.lstat(path)
    except FileNotFoundError:
        return _with_relative_path({"kind": "missing"}, relative_path)

    if stat.S_ISLNK(path_stat.st_mode):
        return _path_content_record(
            path.resolve(strict=True),
            relative_path=relative_path,
            ancestor_directories=ancestor_directories,
        )
    if stat.S_ISREG(path_stat.st_mode):
        return _file_record(path, relative_path=relative_path)
    if not stat.S_ISDIR(path_stat.st_mode):
        raise ValueError(
            "Input fingerprint does not support special files: "
            f"{relative_path or '<root>'}"
        )

    directory_identity = (int(path_stat.st_dev), int(path_stat.st_ino))
    if directory_identity in ancestor_directories:
        raise ValueError(
            "Input fingerprint encountered a symbolic-link directory cycle: "
            f"{relative_path or '<root>'}"
        )
    descendants = ancestor_directories | {directory_identity}
    with os.scandir(path) as entries:
        names = sorted(entry.name for entry in entries)
    members = []
    for name in names:
        logical_path = name if relative_path is None else f"{relative_path}/{name}"
        members.append(
            _path_content_record(
                path / name,
                relative_path=logical_path,
                ancestor_directories=descendants,
            )
        )
    return _with_relative_path(
        {"kind": "directory", "members": members}, relative_path
    )


def input_identities(specs: Iterable[Any]) -> list[Dict[str, str]]:
    """Record each resolved external input by role, path, and content SHA-256."""
    identities = []
    roles = set()
    for spec in specs:
        role = str(spec.role)
        if not role:
            raise ValueError("input role must be non-empty")
        if role in roles:
            raise ValueError(f"duplicate input role: {role}")
        roles.add(role)
        path = Path(spec.path)
        content = _path_content_record(path)
        content_sha = (
            str(content["sha256"])
            if content["kind"] == "file"
            else hash_jsonable(content)
        )
        identities.append(
            {"role": role, "path": str(path), "sha256": content_sha}
        )
    return identities


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    """Return a streamed SHA-256 identity for a file."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def completed_artifact(role: str, path: str | Path) -> Dict[str, Any]:
    """Describe a file only after its write has completed successfully."""
    artifact_path = Path(path)
    return {
        "role": str(role),
        "path": str(artifact_path),
        "status": "completed",
        "size": artifact_path.stat().st_size,
        "sha256": sha256_file(artifact_path),
    }


@contextmanager
def exclusive_run_directory(run_dir: str | Path):
    """Prevent concurrent writers and preserve an unfinished run envelope."""
    directory = Path(run_dir)
    directory.mkdir(parents=True, exist_ok=True)
    lock_path = directory / ".revise-run.lock"
    try:
        lock_path.mkdir()
    except FileExistsError as exc:
        raise RuntimeError(
            f"A REVISE run is already active in {directory}"
        ) from exc

    try:
        manifest_path = directory / "provenance.json"
        if manifest_path.exists():
            try:
                manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as exc:
                raise RuntimeError(
                    f"Cannot safely replace unreadable provenance {manifest_path}"
                ) from exc
            if manifest.get("run", {}).get("status") == "running":
                raise RuntimeError(
                    f"Existing running provenance must be inspected before "
                    f"reusing {directory}"
                )
        yield
    finally:
        lock_path.rmdir()


def collect_software_versions(merged_config: Dict[str, Any]) -> Dict[str, str]:
    packages = [
        ("REVISE", "revise-svc"),
        ("NumPy", "numpy"),
        ("SciPy", "scipy"),
        ("Pandas", "pandas"),
        ("AnnData", "anndata"),
        ("h5py", "h5py"),
    ]
    solver_packages = {"pot": ("POT", "POT"), "tacco": ("tacco", "tacco")}
    selected_solvers = dict.fromkeys(
        str(merged_config["ot"][phase]["solver"]).strip().lower()
        for phase in ("ga", "lr")
    )
    packages.extend(solver_packages[solver] for solver in selected_solvers)

    versions = {"Python": platform.python_version()}
    for label, distribution in packages:
        try:
            versions[label] = pkg_version(distribution)
        except PackageNotFoundError:
            versions[label] = "not-installed"
    return versions


def write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            json.dump(
                payload,
                handle,
                indent=2,
                ensure_ascii=False,
                default=str,
            )
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
    except BaseException:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
        raise
