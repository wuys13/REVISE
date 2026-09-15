"""Render reconstruction-impact reports from published aspect artifacts only.

The source notebook is a repository-owned template.  Rendering intentionally
fails when that template is absent instead of silently substituting a different
notebook, so an installed package must be used with the repository source tree
or an equivalent checked-out report template.
"""

from __future__ import annotations

from hashlib import sha256
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
from typing import Any, Iterable, Mapping


_SUCCESS_STATUSES = {"succeeded", "reused"}
_ASPECT_STATE = Path(".revise/analysis/reconstruction_impact.json")
_SOURCE_NOTEBOOK = (
    Path(__file__).resolve().parents[2]
    / "reproduce/case/reconstruction_impact/Reconstruction_Impact_Report.ipynb"
)


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot read JSON artifact {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"JSON artifact must be an object: {path}")
    return value


def _inside(path: Path, parent: Path) -> bool:
    try:
        path.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def _records(state: Mapping[str, Any]) -> list[dict[str, Any]]:
    raw = state.get("artifact_records", [])
    if isinstance(raw, Mapping):
        raw = list(raw.values())
    if not isinstance(raw, list):
        raise ValueError("reconstruction-impact artifact_records must be a list")
    result = []
    for record in raw:
        if not isinstance(record, Mapping):
            raise ValueError("reconstruction-impact artifact records must be mappings")
        result.append(dict(record))
    return result


def _task_identity(root: Path, state: Mapping[str, Any]) -> tuple[str, str, str | None]:
    metadata: dict[str, Any] = {}
    reconstruction = root / "reconstruction.json"
    if reconstruction.is_symlink():
        raise ValueError(f"Reconstruction metadata must not be a symlink: {reconstruction}")
    if reconstruction.is_file():
        metadata = _read_json(reconstruction)

    sample_value = state.get("sample_id") or metadata.get("sample_id")
    if not isinstance(sample_value, str) or not sample_value.strip():
        raise ValueError("Reconstruction-impact task identity requires sample_id metadata")
    sample_id = sample_value

    if "cell_type" in state:
        cell_value = state["cell_type"]
    elif "cell_type" in metadata:
        cell_value = metadata["cell_type"]
    else:
        raise ValueError("Reconstruction-impact task identity requires cell_type metadata")
    cell_type = None if cell_value in (None, "") else str(cell_value)
    task_id = f"{sample_id}/{cell_type}" if cell_type else sample_id
    return task_id, sample_id, cell_type


def _slug(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._")
    if not slug:
        raise ValueError(f"Cannot derive a safe report task name from {value!r}")
    return slug


def _safe_relative(value: Any, *, label: str) -> str:
    """Validate a manifest-owned path before joining it to a directory."""
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a nonempty relative path")
    path = Path(value)
    if (
        path.is_absolute()
        or any(part in {"", ".", ".."} for part in path.parts)
        or path.as_posix() != value
        or value == "manifest.json"
    ):
        raise ValueError(f"{label} must be a safe relative path: {value!r}")
    return value


def _validate_task(root_value: str | Path) -> dict[str, Any]:
    root_input = Path(root_value)
    if root_input.is_symlink():
        raise ValueError(f"Task root must not be a symlink: {root_input}")
    root = root_input.resolve()
    if not root.is_dir():
        raise ValueError(f"Task root is not a real directory: {root}")
    state_path = root / _ASPECT_STATE
    if not state_path.is_file() or state_path.is_symlink():
        raise ValueError(f"Missing reconstruction-impact state: {state_path}")
    state = _read_json(state_path)
    task_id, sample_id, cell_type = _task_identity(root, state)
    published_value = state.get("published_directory")
    published = None
    if published_value:
        published_input = Path(str(published_value))
        if published_input.is_symlink():
            raise ValueError(f"Published analysis directory must not be a symlink: {published_input}")
        published = published_input.resolve()
        if not _inside(published, root):
            raise ValueError(f"Published analysis directory escapes task root: {published}")

    records = _records(state)
    latest_status = str(state.get("status", "unknown"))
    if latest_status in _SUCCESS_STATUSES and not records:
        raise ValueError("Succeeded reconstruction-impact state has no artifact_records")
    artifacts: dict[str, dict[str, Any]] = {}
    expected: set[Path] = set()
    if records:
        if published is None or not published.is_dir():
            raise ValueError(f"Published analysis directory is missing: {published}")
        for record in records:
            raw_path = record.get("path")
            if not isinstance(raw_path, str) or not raw_path.strip():
                raise ValueError("Every reconstruction-impact artifact needs an absolute path")
            source_input = Path(raw_path)
            if source_input.is_symlink():
                raise ValueError(f"Published artifact must not be a symlink: {source_input}")
            source = source_input.resolve()
            if not source.is_absolute() or not _inside(source, published) or not source.is_file():
                raise ValueError(f"Artifact is outside or missing from published analysis: {source}")
            expected.add(source)
            expected_hash = record.get("sha256")
            if not isinstance(expected_hash, str) or _sha256(source) != expected_hash:
                raise ValueError(f"Published artifact hash mismatch: {source}")
            role = record.get("role")
            if not isinstance(role, str) or not role.strip():
                raise ValueError(f"Published artifact has no role: {source}")
            if role in artifacts:
                raise ValueError(f"Duplicate published artifact role: {role}")
            relative = source.relative_to(published).as_posix()
            artifacts[role] = {
                "role": role,
                "source_path": str(source),
                "relative_path": relative,
                "sha256": expected_hash,
                "description": record.get("description", role),
            }
        actual = {
            path.resolve()
            for path in published.rglob("*")
            if path.is_file() or path.is_symlink()
        }
        if any(path.is_symlink() for path in published.rglob("*")) or actual != expected:
            missing = sorted(str(path) for path in expected - actual)
            extra = sorted(str(path) for path in actual - expected)
            raise ValueError(
                f"Published analysis ownership mismatch: missing={missing}, extra={extra}"
            )
    elif published is not None and published.exists():
        actual = [path for path in published.rglob("*") if path.is_file() or path.is_symlink()]
        if actual:
            raise ValueError("Published analysis has files but no artifact_records")

    stale = bool(records) and latest_status not in _SUCCESS_STATUSES
    return {
        "task_id": task_id,
        "sample_id": sample_id,
        "cell_type": cell_type,
        "root": str(root),
        "latest_status": latest_status,
        "error": state.get("error") or state.get("reason"),
        "stale": stale,
        "published_directory": str(published) if published is not None else None,
        "artifacts": artifacts,
    }


def _existing_report(output_dir: Path) -> dict[str, Any] | None:
    if not output_dir.exists():
        return None
    if output_dir.is_symlink() or not output_dir.is_dir():
        raise ValueError(f"Report output is not a real directory: {output_dir}")
    manifest_path = output_dir / "manifest.json"
    if not manifest_path.is_file():
        if any(output_dir.iterdir()):
            raise ValueError(f"Refusing to write non-owned report directory: {output_dir}")
        return None
    previous = _read_json(manifest_path)
    raw_artifacts = previous.get("artifacts", [])
    if not isinstance(raw_artifacts, list):
        raise ValueError("Previous report manifest artifacts must be a list")
    owned: dict[str, str] = {}
    for item in raw_artifacts:
        if not isinstance(item, Mapping):
            raise ValueError("Previous report manifest artifact must be a mapping")
        relative = _safe_relative(item.get("path"), label="Previous report artifact path")
        digest = item.get("sha256")
        if not isinstance(digest, str) or len(digest) != 64:
            raise ValueError(f"Previous report artifact has an invalid hash: {relative}")
        try:
            int(digest, 16)
        except ValueError as exc:
            raise ValueError(f"Previous report artifact has an invalid hash: {relative}") from exc
        if relative in owned:
            raise ValueError(f"Previous report manifest repeats artifact path: {relative}")
        owned[relative] = digest
    actual: dict[str, Path] = {}
    owned_directories: set[str] = set()
    for relative in owned:
        parent = Path(relative).parent
        while parent != Path("."):
            owned_directories.add(parent.as_posix())
            parent = parent.parent
    for path in output_dir.rglob("*"):
        if path.is_symlink():
            raise ValueError(f"Refusing to write report directory containing symlink: {path}")
        relative = path.relative_to(output_dir).as_posix()
        if path.is_dir():
            if relative not in owned_directories:
                raise ValueError(f"Refusing to write report directory with unowned directory: {path}")
        elif path.is_file() and relative != "manifest.json":
            if relative not in owned:
                raise ValueError(f"Refusing to write report directory with unowned file: {path}")
            actual[relative] = path
    if set(actual) != set(owned):
        missing = sorted(set(owned) - set(actual))
        extra = sorted(set(actual) - set(owned))
        raise ValueError(f"Previous report ownership mismatch: missing={missing}, extra={extra}")
    for relative, path in actual.items():
        if _sha256(path) != owned[relative]:
            raise ValueError(f"Previous report artifact hash mismatch: {path}")
    return previous


def _copy_inputs(stage: Path, tasks: list[dict[str, Any]]) -> list[dict[str, Any]]:
    consumed: list[dict[str, Any]] = []
    used_slugs: dict[str, str] = {}
    for task in tasks:
        task_slug = _slug(task["task_id"])
        previous_task = used_slugs.get(task_slug)
        if previous_task is not None and previous_task != task["task_id"]:
            raise ValueError(
                f"Report task IDs produce a colliding input slug: {previous_task!r} and {task['task_id']!r}"
            )
        used_slugs[task_slug] = task["task_id"]
        local_artifacts: dict[str, dict[str, Any]] = {}
        for role, artifact in sorted(task["artifacts"].items()):
            source = Path(artifact["source_path"])
            relative = Path("inputs") / task_slug / artifact["relative_path"]
            target = stage / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            if _sha256(target) != artifact["sha256"]:
                raise ValueError(f"Copied report input hash mismatch: {target}")
            local = {
                "role": role,
                "path": relative.as_posix(),
                "sha256": artifact["sha256"],
                "source_path": artifact["source_path"],
            }
            local_artifacts[role] = local
            consumed.append(
                {
                    "task_id": task["task_id"],
                    "role": role,
                    "source_path": artifact["source_path"],
                    "path": relative.as_posix(),
                    "sha256": artifact["sha256"],
                }
            )
        task["artifacts"] = local_artifacts
    return consumed


def _write_inputs(stage: Path, tasks: list[dict[str, Any]], consumed: list[dict[str, Any]]) -> None:
    payload = {
        "schema_version": 1,
        "tasks": tasks,
        "consumed_artifacts": consumed,
    }
    (stage / "inputs.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _execute_notebook(path: Path, working_directory: Path) -> None:
    import nbformat
    from nbclient import NotebookClient

    notebook = nbformat.read(path, as_version=4)
    client = NotebookClient(
        notebook,
        timeout=600,
        kernel_name="python3",
        resources={"metadata": {"path": str(working_directory)}},
    )
    client.execute()
    nbformat.write(notebook, path)


def _code_cells(path: Path) -> list[str]:
    notebook = _read_json(path)
    sources = []
    for cell in notebook.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        source = cell.get("source", "")
        sources.append("".join(source) if isinstance(source, list) else source)
    return sources


def _artifact_list(stage: Path) -> list[dict[str, str]]:
    paths = sorted(path for path in stage.rglob("*") if path.is_file())
    expected_roots = {"source.ipynb", "report.ipynb", "inputs.json"}
    for path in paths:
        relative = path.relative_to(stage).as_posix()
        if (
            relative not in expected_roots
            and not relative.startswith("inputs/")
            and not relative.startswith("figures/")
        ):
            raise ValueError(f"Report notebook produced an unowned file: {relative}")
    if not (stage / "source.ipynb").is_file() or not (stage / "report.ipynb").is_file():
        raise ValueError("Report notebooks were not materialized")
    records = []
    for path in paths:
        relative = path.relative_to(stage).as_posix()
        if relative == "source.ipynb":
            role = "source_notebook"
        elif relative == "report.ipynb":
            role = "report_notebook"
        elif relative == "inputs.json":
            role = "inputs"
        elif relative.startswith("inputs/"):
            role = f"input:{relative[len('inputs/'):]}"
        else:
            role = f"figure:{Path(relative).stem}"
        records.append({"role": role, "path": relative, "sha256": _sha256(path)})
    return records


def _failure_manifest(
    output_dir: Path,
    tasks: list[dict[str, Any]],
    consumed: list[dict[str, Any]],
    error: str,
    previous: Mapping[str, Any] | None,
) -> dict[str, Any]:
    manifest = {
        "schema_version": 1,
        "status": "failed",
        "report_directory": str(output_dir),
        "tasks": tasks,
        "consumed_artifacts": consumed,
        "error": error,
        "execution": {"requested": True, "status": "failed"},
    }
    if previous:
        manifest["previous_manifest_status"] = previous.get("status")
        manifest["artifacts"] = previous.get("artifacts", [])
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return manifest


def _publish_stage(stage: Path, destination: Path) -> None:
    """Atomically replace a report directory while retaining rollback state."""
    backup: Path | None = None
    if destination.exists():
        backup = Path(tempfile.mkdtemp(prefix=".impact-report-old-", dir=destination.parent))
        backup.rmdir()
        os.replace(destination, backup)
    try:
        os.replace(stage, destination)
    except BaseException:
        if backup is not None and backup.exists():
            os.replace(backup, destination)
        raise
    if backup is not None and backup.exists():
        shutil.rmtree(backup, ignore_errors=True)


def render_impact_report(
    task_roots: Iterable[str | Path],
    output_dir: str | Path,
    *,
    execute: bool = True,
) -> dict[str, Any]:
    """Render a report from one or more published reconstruction-impact tasks."""
    roots = [Path(value) for value in task_roots]
    if not roots:
        raise ValueError("render_impact_report requires at least one task root")
    destination_input = Path(output_dir)
    if destination_input.is_symlink():
        raise ValueError(f"Report output is not a real directory: {destination_input}")
    destination = destination_input.resolve()
    previous = _existing_report(destination)
    tasks = [_validate_task(root) for root in roots]
    task_ids = [task["task_id"] for task in tasks]
    if len(set(task_ids)) != len(task_ids):
        raise ValueError("Report task IDs must be unique")
    if not _SOURCE_NOTEBOOK.is_file():
        raise ValueError(f"Missing report source notebook: {_SOURCE_NOTEBOOK}")
    consumed: list[dict[str, Any]] = []
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".impact-report-", dir=destination.parent))
    try:
        shutil.copy2(_SOURCE_NOTEBOOK, stage / "source.ipynb")
        shutil.copy2(_SOURCE_NOTEBOOK, stage / "report.ipynb")
        consumed = _copy_inputs(stage, tasks)
        _write_inputs(stage, tasks, consumed)
        if execute:
            _execute_notebook(stage / "report.ipynb", stage)
        if _code_cells(stage / "source.ipynb") != _code_cells(stage / "report.ipynb"):
            raise ValueError("Executed report notebook code differs from source notebook")
        records = _artifact_list(stage)
        status = "succeeded" if all(task["latest_status"] in _SUCCESS_STATUSES for task in tasks) else "failed"
        manifest = {
            "schema_version": 1,
            "status": status,
            "report_directory": str(destination),
            "source_notebook": {
                "path": "source.ipynb",
                "sha256": next(item["sha256"] for item in records if item["role"] == "source_notebook"),
                "repository_path": str(_SOURCE_NOTEBOOK),
            },
            "tasks": tasks,
            "consumed_artifacts": consumed,
            "artifacts": records,
            "execution": {"requested": bool(execute), "status": "succeeded" if execute else "skipped"},
        }
        if status == "failed":
            manifest["reason"] = "At least one latest task state is failed or unavailable; published artifacts are stale"
        (stage / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        _publish_stage(stage, destination)
        return manifest
    except Exception as exc:
        if execute:
            return _failure_manifest(
                destination,
                tasks,
                consumed,
                f"{type(exc).__name__}: {exc}",
                previous,
            )
        raise
    finally:
        shutil.rmtree(stage, ignore_errors=True)
