from __future__ import annotations

import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
from typing import Any


CONFOUNDINGS = (
    "segmentation",
    "bin2cell",
    "batch_effect",
    "spot_size",
    "gene_panel",
    "gene_dropout",
)
_RECEIVED_SIGNAL: int | None = None
TERMINATION_GRACE_SECONDS = 2.0
KILL_REAP_SECONDS = 2.0


def _capture_signal(signum, _frame) -> None:
    global _RECEIVED_SIGNAL
    if _RECEIVED_SIGNAL is None:
        _RECEIVED_SIGNAL = int(signum)


def _positive_int(raw: str, *, name: str) -> int:
    try:
        value = int(raw)
    except ValueError as exc:
        raise ValueError(f"{name} must be a positive integer, got {raw!r}") from exc
    if value <= 0:
        raise ValueError(f"{name} must be a positive integer, got {raw!r}")
    return value


def _write_json(path: Path, payload: dict[str, Any]) -> None:
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
            json.dump(payload, handle, indent=2, ensure_ascii=False)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
    except BaseException:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
        raise


def _public_task(task: dict[str, Any]) -> dict[str, Any]:
    return {
        key: task[key]
        for key in (
            "id",
            "sample_name",
            "confounding",
            "status",
            "pid",
            "exit_code",
            "log_path",
            "error",
        )
    }


def _status_payload(
    tasks: list[dict[str, Any]],
    *,
    status: str,
    max_jobs: int,
    launcher_error: dict[str, str] | None = None,
) -> dict[str, Any]:
    counts = {
        state: sum(task["status"] == state for task in tasks)
        for state in (
            "pending",
            "running",
            "succeeded",
            "failed",
            "interrupted",
            "skipped",
        )
    }
    return {
        "schema_version": 1,
        "status": status,
        "max_jobs": max_jobs,
        "summary": counts,
        "tasks": [_public_task(task) for task in tasks],
        "launcher_error": launcher_error,
    }


def _write_status(
    status_path: Path,
    tasks: list[dict[str, Any]],
    *,
    status: str,
    max_jobs: int,
    launcher_error: dict[str, str] | None = None,
) -> None:
    _write_json(
        status_path,
        _status_payload(
            tasks,
            status=status,
            max_jobs=max_jobs,
            launcher_error=launcher_error,
        ),
    )


def _build_tasks() -> tuple[list[dict[str, Any]], Path]:
    raw_data_path = os.environ.get("RAW_DATA_PATH") or "./raw_data/Sim2Real-ST"
    run_stamp = time.strftime("%Y%m%d_%H%M%S")
    save_path = Path(
        os.environ.get("SAVE_PATH")
        or f"results_unified/benchmark_runs/{run_stamp}"
    )
    sample_patient = os.environ.get("SAMPLE_PATIENT") or "P2CRC"
    sample_parts = (os.environ.get("SAMPLE_PARTS") or "part1").split()
    if not sample_parts:
        sample_parts = ["part1"]
    config_path = os.environ.get("CONFIG_PATH") or "revise/revise.yaml"
    child_python = os.environ.get("BENCHMARK_CHILD_PYTHON") or sys.executable
    benchmark_main = Path(__file__).resolve().parents[2] / "benchmark_main.py"

    tasks = []
    for sample_part in sample_parts:
        sample_name = f"{sample_patient}/cut_{sample_part}"
        record_dir = Path("0_records") / f"{sample_patient}_{sample_part}"
        for confounding in CONFOUNDINGS:
            log_path = record_dir / f"{confounding}.log"
            tasks.append(
                {
                    "id": f"{sample_name}|{confounding}",
                    "sample_name": sample_name,
                    "confounding": confounding,
                    "status": "pending",
                    "pid": None,
                    "exit_code": None,
                    "log_path": str(log_path),
                    "error": None,
                    "command": [
                        child_python,
                        "-u",
                        str(benchmark_main),
                        "--config",
                        config_path,
                        "--confounding",
                        confounding,
                        "--data-root",
                        raw_data_path,
                        "--dataset-task",
                        confounding,
                        "--sample-name",
                        sample_name,
                        "--output-root",
                        str(save_path),
                    ],
                }
            )
    return tasks, save_path / "launcher_status.json"


def _process_group_exists(pgid: int) -> bool:
    try:
        os.killpg(pgid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _signal_groups(children: list[dict[str, Any]], signum: int) -> None:
    for child in children:
        try:
            os.killpg(child["pgid"], signum)
        except ProcessLookupError:
            pass


def _groups_quiesced(children: list[dict[str, Any]]) -> bool:
    quiesced = True
    for child in children:
        child["process"].poll()
        if child["process"].returncode is None or _process_group_exists(
            child["pgid"]
        ):
            quiesced = False
    return quiesced


def _wait_for_groups(children: list[dict[str, Any]], timeout: float) -> bool:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if _groups_quiesced(children):
            return True
        time.sleep(0.02)
    return _groups_quiesced(children)


def _drain_groups(children: list[dict[str, Any]], signum: int) -> bool:
    if not children:
        return True
    _signal_groups(children, signum)
    if _wait_for_groups(children, TERMINATION_GRACE_SECONDS):
        return True
    _signal_groups(children, signal.SIGKILL)
    return _wait_for_groups(children, KILL_REAP_SECONDS)


def _finalize_interrupted(active: list[dict[str, Any]], signum: int) -> None:
    _drain_groups(active, signum)
    for child in active:
        process = child["process"]
        process.poll()
        task = child["task"]
        task["exit_code"] = process.returncode
        task["status"] = "interrupted"
        if process.returncode is None or _process_group_exists(child["pgid"]):
            task["error"] = "child process group did not quiesce after SIGKILL"
        child["log_handle"].close()
    active.clear()


def _reap_completed(
    active: list[dict[str, Any]],
) -> bool:
    changed = False
    for child in list(active):
        process = child["process"]
        returncode = process.poll()
        if returncode is None:
            continue
        task = child["task"]
        task["exit_code"] = int(returncode)
        if _process_group_exists(child["pgid"]):
            task["error"] = (
                "child process group remained active after its leader exited"
            )
            if not _drain_groups([child], signal.SIGTERM):
                task["error"] += " and did not quiesce after SIGKILL"
            task["status"] = "failed"
        elif returncode == 0:
            task["status"] = "succeeded"
        else:
            task["status"] = "failed"
        child["log_handle"].close()
        active.remove(child)
        changed = True
    return changed


def _cleanup_after_error(active: list[dict[str, Any]]) -> None:
    _finalize_interrupted(active, signal.SIGTERM)


def run_tasks(
    tasks: list[dict[str, Any]],
    *,
    max_jobs: int,
    status_path: Path,
) -> int:
    active: list[dict[str, Any]] = []
    next_task = 0
    forwarded_signal = None
    current_task = None
    _write_status(status_path, tasks, status="running", max_jobs=max_jobs)

    try:
        while next_task < len(tasks) or active:
            changed = _reap_completed(active)

            if _RECEIVED_SIGNAL is not None and forwarded_signal is None:
                forwarded_signal = _RECEIVED_SIGNAL
                for task in tasks[next_task:]:
                    task["status"] = "skipped"
                next_task = len(tasks)
                _finalize_interrupted(active, forwarded_signal)
                changed = True

            while (
                forwarded_signal is None
                and _RECEIVED_SIGNAL is None
                and next_task < len(tasks)
                and len(active) < max_jobs
            ):
                task = tasks[next_task]
                current_task = task
                next_task += 1
                log_path = Path(task["log_path"])
                log_path.parent.mkdir(parents=True, exist_ok=True)
                log_handle = log_path.open("w", encoding="utf-8")
                try:
                    process = subprocess.Popen(
                        task["command"],
                        stdout=log_handle,
                        stderr=subprocess.STDOUT,
                        start_new_session=True,
                    )
                except OSError as exc:
                    log_handle.close()
                    task["status"] = "failed"
                    task["exit_code"] = 127
                    task["error"] = f"{type(exc).__name__}: {exc}"
                else:
                    task["status"] = "running"
                    task["pid"] = int(process.pid)
                    active.append(
                        {
                            "process": process,
                            "pgid": int(process.pid),
                            "task": task,
                            "log_handle": log_handle,
                        }
                    )
                    print(
                        f"Start sample={task['sample_name']}; "
                        f"cf={task['confounding']}; pid={process.pid}"
                    )
                current_task = None
                changed = True

            if changed:
                _write_status(
                    status_path,
                    tasks,
                    status="running",
                    max_jobs=max_jobs,
                )
            if active and not changed:
                time.sleep(0.02)
    except Exception as exc:
        _cleanup_after_error(active)
        if current_task is not None and current_task["status"] == "pending":
            current_task["status"] = "failed"
            current_task["error"] = f"{type(exc).__name__}: {exc}"
        for task in tasks:
            if task["status"] == "pending":
                task["status"] = "skipped"
        launcher_error = {"type": type(exc).__name__, "message": str(exc)}
        try:
            _write_status(
                status_path,
                tasks,
                status="failed",
                max_jobs=max_jobs,
                launcher_error=launcher_error,
            )
        except Exception:
            pass
        print(
            f"Benchmark launcher failed: {type(exc).__name__}: {exc}",
            file=sys.stderr,
        )
        return 1

    if forwarded_signal is not None:
        final_status = "interrupted"
        exit_code = 128 + int(forwarded_signal)
    elif any(task["status"] == "failed" for task in tasks):
        final_status = "failed"
        exit_code = 1
    else:
        final_status = "succeeded"
        exit_code = 0
    _write_status(status_path, tasks, status=final_status, max_jobs=max_jobs)
    summary = _status_payload(
        tasks,
        status=final_status,
        max_jobs=max_jobs,
    )["summary"]
    print(
        "Benchmark launcher finished: "
        f"status={final_status}; succeeded={summary['succeeded']}; "
        f"failed={summary['failed']}; interrupted={summary['interrupted']}; "
        f"skipped={summary['skipped']}"
    )
    return exit_code


def main() -> int:
    try:
        max_jobs = _positive_int(
            os.environ.get("BENCHMARK_MAX_JOBS", "6"),
            name="BENCHMARK_MAX_JOBS",
        )
    except ValueError as exc:
        print(f"benchmark launcher: {exc}", file=sys.stderr)
        return 2
    tasks, status_path = _build_tasks()
    print(f"Benchmark output root: {status_path.parent}")
    signal.signal(signal.SIGINT, _capture_signal)
    signal.signal(signal.SIGTERM, _capture_signal)
    return run_tasks(tasks, max_jobs=max_jobs, status_path=status_path)


if __name__ == "__main__":
    raise SystemExit(main())
