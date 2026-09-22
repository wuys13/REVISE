#!/usr/bin/env python3
"""Run one full reconstruction with a persistent receipt and peak-resource record."""

from __future__ import annotations

import argparse
import copy
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import traceback

import yaml


def now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat()


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def hash_jsonable(value) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def scientific_config(document: dict) -> dict:
    normalized = copy.deepcopy(document)
    normalized["inputs"]["st"].pop("path")
    normalized["inputs"]["reference"].pop("path")
    normalized["output"].pop("dir")
    return normalized


def run_text(command: list[str], *, cwd: Path) -> str:
    return subprocess.run(
        command,
        cwd=cwd,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    ).stdout.strip()


def read_integer(path: Path) -> int | None:
    try:
        value = path.read_text().strip()
    except OSError:
        return None
    if value == "max":
        return None
    try:
        return int(value)
    except ValueError:
        return None


def read_key_values(path: Path) -> dict[str, int] | None:
    try:
        lines = path.read_text().splitlines()
        return {key: int(value) for key, value in (line.split() for line in lines)}
    except (OSError, ValueError):
        return None


def cgroup_snapshot() -> dict:
    memory_current = read_integer(Path("/sys/fs/cgroup/memory.current"))
    memory_stat = read_key_values(Path("/sys/fs/cgroup/memory.stat"))
    inactive_file = memory_stat.get("inactive_file") if memory_stat is not None else None
    working_set = (
        max(0, memory_current - inactive_file)
        if memory_current is not None and inactive_file is not None
        else None
    )
    return {
        "memory_limit_bytes": read_integer(Path("/sys/fs/cgroup/memory.max")),
        "memory_current_bytes": memory_current,
        "memory_working_set_bytes": working_set,
        "memory_peak_bytes": read_integer(Path("/sys/fs/cgroup/memory.peak")),
        "memory_stat": (
            {key: memory_stat.get(key) for key in ("anon", "file", "inactive_file")}
            if memory_stat is not None
            else None
        ),
        "memory_events": read_key_values(Path("/sys/fs/cgroup/memory.events")),
        "cpu_max": (
            Path("/sys/fs/cgroup/cpu.max").read_text().strip()
            if Path("/sys/fs/cgroup/cpu.max").is_file()
            else None
        ),
    }


def process_tree(root_pid: int) -> set[int]:
    pending = [root_pid]
    result: set[int] = set()
    while pending:
        pid = pending.pop()
        if pid in result:
            continue
        result.add(pid)
        children = Path(f"/proc/{pid}/task/{pid}/children")
        try:
            pending.extend(int(value) for value in children.read_text().split())
        except (OSError, ValueError):
            pass
    return result


def process_memory(root_pid: int) -> dict:
    rss_kib = swap_kib = 0
    pids = process_tree(root_pid)
    live = 0
    for pid in pids:
        try:
            lines = Path(f"/proc/{pid}/status").read_text().splitlines()
        except OSError:
            continue
        live += 1
        for line in lines:
            if line.startswith("VmRSS:"):
                rss_kib += int(line.split()[1])
            elif line.startswith("VmSwap:"):
                swap_kib += int(line.split()[1])
    if live == 0:
        return {"process_count": None, "rss_bytes": None, "swap_bytes": None}
    return {"process_count": live, "rss_bytes": rss_kib * 1024, "swap_bytes": swap_kib * 1024}


def update_peak(peak: dict, key: str, value: int | None) -> None:
    if value is not None:
        peak[key] = value if peak.get(key) is None else max(peak[key], value)


def atomic_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--expected-st-sha256", required=True)
    parser.add_argument("--expected-reference-sha256", required=True)
    parser.add_argument("--notebook-id", required=True)
    parser.add_argument("--resource-role", default="primary")
    parser.add_argument("--iteration", default="iter-001")
    parser.add_argument("--poll-seconds", type=float, default=2.0)
    parser.add_argument("--receipt-seconds", type=float, default=30.0)
    parser.add_argument("--memory-stop-fraction", type=float, default=0.90)
    parser.add_argument("--term-grace-seconds", type=float, default=30.0)
    args = parser.parse_args()
    if args.poll_seconds <= 0 or args.receipt_seconds <= 0 or args.term_grace_seconds <= 0:
        parser.error("poll, receipt, and TERM grace intervals must be positive")
    if not 0 < args.memory_stop_fraction < 1:
        parser.error("--memory-stop-fraction must be between 0 and 1")

    repo = Path(__file__).resolve().parents[1]
    config_path = (args.config if args.config.is_absolute() else repo / args.config).resolve()
    run_root = (args.run_root if args.run_root.is_absolute() else repo / args.run_root).resolve()
    receipt_path = run_root / "receipt.json"
    audit_path = run_root / "audit" / "delivery-audit.json"
    donor_sidecar = run_root / "audit" / "random-donor-assignments.tsv.gz"
    package_record = run_root / "environment" / "pip-freeze.txt"
    if receipt_path.exists():
        raise SystemExit(
            f"refusing to overwrite an existing run receipt; use a new iteration or an explicit recovery flow: {receipt_path}"
        )
    document = yaml.safe_load(config_path.read_text())
    base_config_path = repo / "configs/acceptance/full/P2CRC_Xenium-random.yaml"
    base_document = yaml.safe_load(base_config_path.read_text())
    normalized_science = scientific_config(document)
    if normalized_science != scientific_config(base_document):
        raise SystemExit("run-scoped config changes frozen scientific parameters")
    if document["paths"]["root_dir"] != ".":
        raise SystemExit("run launcher requires paths.root_dir: .")
    st_path = (repo / document["inputs"]["st"]["path"]).resolve()
    reference_path = (repo / document["inputs"]["reference"]["path"]).resolve()
    output_dir = (repo / document["output"]["dir"]).resolve()
    expected_output_dir = run_root / "P2CRC_Xenium" / "random"
    if run_root.name != args.iteration or output_dir != expected_output_dir:
        raise SystemExit("--iteration must match both the run-root and configured output path")
    if output_dir.exists():
        raise SystemExit(f"refusing to reuse existing run target: {output_dir}")
    git_dirty = run_text(["git", "status", "--short"], cwd=repo).splitlines()
    if git_dirty:
        raise SystemExit(f"refusing to launch from a dirty worktree: {git_dirty}")
    expected_python = (repo / ".venv/bin/python").resolve()
    if Path(sys.executable).resolve() != expected_python:
        raise SystemExit(f"launcher must use the isolated project interpreter: {expected_python}")
    if not package_record.is_file():
        raise SystemExit(f"environment package record is missing: {package_record}")
    run_root.mkdir(parents=True, exist_ok=True)

    input_identity = {
        "spatial": {"path": str(st_path), "size": st_path.stat().st_size, "sha256": digest(st_path)},
        "reference": {"path": str(reference_path), "size": reference_path.stat().st_size, "sha256": digest(reference_path)},
    }
    if input_identity["spatial"]["sha256"] != args.expected_st_sha256:
        raise SystemExit("spatial input hash does not match the launch gate")
    if input_identity["reference"]["sha256"] != args.expected_reference_sha256:
        raise SystemExit("reference input hash does not match the launch gate")

    command = [sys.executable, str(repo / "reconstruct.py"), "--config", str(config_path)]
    receipt = {
        "schema_version": 1,
        "task_id": "p2-xenium-impact-20260922",
        "iteration": args.iteration,
        "owner_project": "revise",
        "upstream": ["sc_annotation", "tumor"],
        "downstream": ["analysis"],
        "code_commit": run_text(["git", "rev-parse", "HEAD"], cwd=repo),
        "branch": run_text(["git", "branch", "--show-current"], cwd=repo),
        "git_dirty": git_dirty,
        "execution_location": "remote",
        "notebook": args.notebook_id,
        "resource_role": args.resource_role,
        "python_interpreter": sys.executable,
        "python_version": sys.version,
        "package_record": {
            "path": str(package_record),
            "sha256": digest(package_record) if package_record.is_file() else None,
        },
        "config": {
            "path": str(config_path),
            "sha256": digest(config_path),
            "base_path": str(base_config_path),
            "scientific_config_sha256": hash_jsonable(normalized_science),
        },
        "input_identity": input_identity,
        "command": command,
        "environment": {
            key: value
            for key, value in {
                "OMP_NUM_THREADS": "1",
                "OPENBLAS_NUM_THREADS": "1",
                "MKL_NUM_THREADS": "1",
                "NUMEXPR_NUM_THREADS": "1",
                "VECLIB_MAXIMUM_THREADS": "1",
                "NUMBA_NUM_THREADS": "1",
                "PYTHONHASHSEED": "42",
                "MPLBACKEND": "Agg",
            }.items()
        },
        "output_entry": str(output_dir),
        "evidence_entry": str(receipt_path),
        "audit_entry": str(audit_path),
        "started_at": now(),
        "observed_at": now(),
        "status": "running",
        "acceptance_level": "engineering",
        "scientific_acceptance": False,
        "resources": {"baseline": cgroup_snapshot()},
        "resource_guard": {
            "memory_stop_fraction": args.memory_stop_fraction,
            "term_grace_seconds": args.term_grace_seconds,
            "action": "SIGTERM only the launched reconstruction process group",
        },
        "resubmit_safety": "Inspect the recorded PID/PGID and output target before any retry.",
    }
    atomic_json(receipt_path, receipt)

    child = None
    interrupted = None

    def stop(signum, _frame):
        nonlocal interrupted
        interrupted = signum
        if child is not None and child.poll() is None:
            os.killpg(child.pid, signal.SIGTERM)

    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    env = os.environ.copy()
    env.update(receipt["environment"])
    peak = {
        "process_count": None,
        "rss_bytes": None,
        "swap_bytes": None,
        "cgroup_memory_current_bytes": None,
        "cgroup_memory_working_set_bytes": None,
        "cgroup_anon_bytes": None,
        "cgroup_file_bytes": None,
        "cgroup_inactive_file_bytes": None,
    }
    resource_guard_triggered = None
    resource_guard_triggered_monotonic = None
    resource_guard_kill_sent = False
    try:
        child = subprocess.Popen(command, cwd=repo, env=env, start_new_session=True)
        receipt["process"] = {"pid": child.pid, "pgid": os.getpgid(child.pid)}
        atomic_json(receipt_path, receipt)
        last_receipt = time.monotonic()
        while child.poll() is None:
            current = process_memory(child.pid)
            for key in ("process_count", "rss_bytes", "swap_bytes"):
                update_peak(peak, key, current[key])
            cgroup = cgroup_snapshot()
            cgroup_current = cgroup.get("memory_current_bytes")
            cgroup_working_set = cgroup.get("memory_working_set_bytes")
            memory_stat = cgroup.get("memory_stat") or {}
            update_peak(peak, "cgroup_memory_current_bytes", cgroup_current)
            update_peak(peak, "cgroup_memory_working_set_bytes", cgroup_working_set)
            update_peak(peak, "cgroup_anon_bytes", memory_stat.get("anon"))
            update_peak(peak, "cgroup_file_bytes", memory_stat.get("file"))
            update_peak(peak, "cgroup_inactive_file_bytes", memory_stat.get("inactive_file"))
            memory_limit = receipt["resources"]["baseline"].get("memory_limit_bytes")
            if (
                resource_guard_triggered is None
                and memory_limit is not None
                and cgroup_working_set is not None
                and cgroup_working_set >= memory_limit * args.memory_stop_fraction
            ):
                resource_guard_triggered = {
                    "observed_at": now(),
                    "memory_current_bytes": cgroup_current,
                    "memory_working_set_bytes": cgroup_working_set,
                    "memory_limit_bytes": memory_limit,
                    "working_set_fraction": cgroup_working_set / memory_limit,
                    "memory_stat": cgroup.get("memory_stat"),
                    "memory_events": cgroup.get("memory_events"),
                }
                resource_guard_triggered_monotonic = time.monotonic()
                receipt["resource_guard"]["triggered"] = resource_guard_triggered
                receipt["resources"]["observed_peak"] = dict(peak)
                receipt["observed_at"] = now()
                atomic_json(receipt_path, receipt)
                os.killpg(child.pid, signal.SIGTERM)
            if (
                resource_guard_triggered_monotonic is not None
                and not resource_guard_kill_sent
                and time.monotonic() - resource_guard_triggered_monotonic >= args.term_grace_seconds
                and child.poll() is None
            ):
                os.killpg(child.pid, signal.SIGKILL)
                resource_guard_kill_sent = True
                receipt["resource_guard"]["kill_escalation"] = {
                    "observed_at": now(),
                    "reason": "launched process group did not exit within TERM grace period",
                }
                atomic_json(receipt_path, receipt)
            if time.monotonic() - last_receipt >= args.receipt_seconds:
                receipt["resources"]["observed_peak"] = dict(peak)
                receipt["process"]["state"] = "alive"
                receipt["observed_at"] = now()
                atomic_json(receipt_path, receipt)
                last_receipt = time.monotonic()
            time.sleep(args.poll_seconds)
        returncode = child.returncode
        receipt["process"]["state"] = "exited"
        receipt["returncode"] = returncode
        receipt["resources"]["observed_peak"] = peak
        receipt["resources"]["final"] = cgroup_snapshot()
        if resource_guard_triggered is not None:
            receipt["status"] = "stopped_resource_guard"
        elif interrupted is not None:
            receipt["status"] = "stopped"
            receipt["signal"] = interrupted
        elif returncode != 0:
            receipt["status"] = "reconstruction_failed"
        else:
            receipt["status"] = "reconstruction_exited_zero_pending_artifact_validation"
        receipt["observed_at"] = now()
        atomic_json(receipt_path, receipt)
        if resource_guard_triggered is not None:
            raise SystemExit(returncode if returncode != 0 else 75)
        if interrupted is not None:
            raise SystemExit(returncode if returncode != 0 else 128 + int(interrupted))
        if returncode != 0:
            raise SystemExit(returncode)

        audit_command = [
            sys.executable,
            str(repo / "scripts" / "verify_full_ist_delivery.py"),
            "--repo", str(repo),
            "--config", str(config_path),
            "--expected-st-sha256", args.expected_st_sha256,
            "--expected-reference-sha256", args.expected_reference_sha256,
            "--audit", str(audit_path),
            "--donor-sidecar", str(donor_sidecar),
        ]
        audit_result = subprocess.run(audit_command, cwd=repo, env=env)
        receipt["audit_command"] = audit_command
        receipt["audit_returncode"] = audit_result.returncode
        receipt["status"] = "technical_verified" if audit_result.returncode == 0 else "artifact_validation_failed"
        receipt["observed_at"] = now()
        receipt["resources"]["final"] = cgroup_snapshot()
        atomic_json(receipt_path, receipt)
        if audit_result.returncode != 0:
            raise SystemExit(audit_result.returncode)
    except BaseException as exc:
        child_alive = child is not None and child.poll() is None
        receipt["child_state_at_launcher_exception"] = "alive" if child_alive else "exited_or_not_started"
        if child_alive:
            try:
                os.killpg(child.pid, signal.SIGTERM)
                try:
                    child.wait(timeout=10)
                    receipt["child_cleanup"] = "SIGTERM sent; launched process group exited"
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL)
                    child.wait(timeout=10)
                    receipt["child_cleanup"] = "SIGTERM timed out; SIGKILL sent to launched process group"
            except ProcessLookupError:
                receipt["child_cleanup"] = "process group already exited"
        if receipt.get("status") == "running":
            receipt["status"] = "launcher_failed"
        receipt["error"] = f"{type(exc).__name__}: {exc}"
        receipt["traceback"] = traceback.format_exc()
        receipt["observed_at"] = now()
        receipt.setdefault("resources", {})["final"] = cgroup_snapshot()
        atomic_json(receipt_path, receipt)
        raise


if __name__ == "__main__":
    main()
