from __future__ import annotations

import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import pytest


ROOT = Path(__file__).resolve().parents[2]
LAUNCHER = ROOT / "reproduce" / "benchmark_main.sh"
CONFOUNDINGS = {
    "segmentation",
    "bin2cell",
    "batch_effect",
    "spot_size",
    "gene_panel",
    "gene_dropout",
}

FAKE_PYTHON = r"""
import fcntl
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

args = sys.argv[1:]
cf = args[args.index("--confounding") + 1]
sample_name = args[args.index("--sample-name") + 1]
sample_part = sample_name.rsplit("/", 1)[-1].removeprefix("cut_")
marker = cf if sample_part == "part1" else f"{sample_part}__{cf}"
state_dir = Path(os.environ["FAKE_BENCHMARK_STATE"])
state_dir.mkdir(parents=True, exist_ok=True)
state_path = state_dir / "concurrency.json"
pid = os.getpid()
(state_dir / f"{marker}.pid").write_text(str(pid))
(state_dir / f"{marker}.args.json").write_text(json.dumps(args))


def update_current(delta):
    with state_path.open("a+") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        handle.seek(0)
        raw = handle.read()
        state = json.loads(raw) if raw else {"current": 0, "maximum": 0}
        state["current"] += delta
        state["maximum"] = max(state["maximum"], state["current"])
        handle.seek(0)
        handle.truncate()
        json.dump(state, handle)
        handle.flush()
        fcntl.flock(handle, fcntl.LOCK_UN)


def interrupted(signum, _frame):
    (state_dir / f"{marker}.signal").write_text(str(signum))
    if os.environ.get("FAKE_BENCHMARK_IGNORE_SIGNALS") == "1":
        return
    raise SystemExit(128 + signum)


signal.signal(signal.SIGTERM, interrupted)
signal.signal(signal.SIGINT, interrupted)
update_current(1)
(state_dir / f"{marker}.started").touch()
try:
    if cf == os.environ.get("FAKE_BENCHMARK_GRANDCHILD_CF"):
        grandchild = subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(30)"]
        )
        (state_dir / f"{marker}.grandchild.pid").write_text(str(grandchild.pid))
        raise SystemExit(0)
    if cf == os.environ.get("FAKE_BENCHMARK_ACTIVE_GRANDCHILD_CF"):
        grandchild = subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(30)"]
        )
        (state_dir / f"{marker}.grandchild.pid").write_text(str(grandchild.pid))
    if cf == os.environ.get("FAKE_BENCHMARK_FAIL"):
        raise SystemExit(int(os.environ.get("FAKE_BENCHMARK_FAIL_CODE", "7")))
    release_path = os.environ.get("FAKE_BENCHMARK_RELEASE")
    if release_path:
        while not Path(release_path).exists():
            time.sleep(0.01)
    else:
        time.sleep(float(os.environ.get("FAKE_BENCHMARK_SLEEP", "0.05")))
    (state_dir / f"{marker}.completed").touch()
finally:
    update_current(-1)
"""


def _launcher_env(
    tmp_path: Path,
    *,
    jobs: int,
    sleep: float = 0.05,
    barrier: bool = False,
) -> dict:
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake_python = fake_bin / "python"
    fake_python.write_text(f"#!{sys.executable}\n{FAKE_PYTHON}")
    fake_python.chmod(0o755)
    env = os.environ.copy()
    env.update(
        {
            "BENCHMARK_LAUNCHER_PYTHON": sys.executable,
            "BENCHMARK_CHILD_PYTHON": str(fake_python),
            "FAKE_BENCHMARK_STATE": str(tmp_path / "state"),
            "FAKE_BENCHMARK_SLEEP": str(sleep),
            "BENCHMARK_MAX_JOBS": str(jobs),
            "SAVE_PATH": str(tmp_path / "outputs"),
            "SAMPLE_PARTS": "part1",
        }
    )
    if barrier:
        env["FAKE_BENCHMARK_RELEASE"] = str(tmp_path / "release")
    return env


def _status_records(tmp_path: Path) -> dict[str, dict]:
    manifest = _launcher_status(tmp_path)
    return {
        record["confounding"]: record
        for record in manifest["tasks"]
    }


def _launcher_status(tmp_path: Path) -> dict:
    return json.loads(
        (tmp_path / "outputs" / "launcher_status.json").read_text()
    )


def _live_fake_pids(tmp_path: Path) -> list[int]:
    live = []
    for path in (tmp_path / "state").glob("*.pid"):
        pid = int(path.read_text())
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            continue
        live.append(pid)
    return live


def _pid_is_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    return True


def _cleanup_fake_pids(tmp_path: Path) -> None:
    live = _live_fake_pids(tmp_path)
    for pid in live:
        try:
            os.kill(pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    deadline = time.monotonic() + 0.2
    while _live_fake_pids(tmp_path) and time.monotonic() < deadline:
        time.sleep(0.01)
    for pid in _live_fake_pids(tmp_path):
        try:
            os.kill(pid, signal.SIGKILL)
        except ProcessLookupError:
            pass


def _wait_for_started(tmp_path: Path, count: int, timeout: float = 5.0) -> None:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if len(list((tmp_path / "state").glob("*.started"))) >= count:
            return
        time.sleep(0.02)
    raise AssertionError(f"fewer than {count} children started before timeout")


def _wait_for_path(path: Path, timeout: float = 5.0) -> None:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if path.exists():
            return
        time.sleep(0.02)
    raise AssertionError(f"path did not appear before timeout: {path}")


def test_launcher_waits_for_all_successful_children_and_records_status(tmp_path):
    env = _launcher_env(tmp_path, jobs=6, barrier=True)
    process = subprocess.Popen(
        ["bash", str(LAUNCHER)],
        cwd=tmp_path,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    try:
        _wait_for_started(tmp_path, 6)
        assert process.poll() is None
        (tmp_path / "release").touch()
        stdout, stderr = process.communicate(timeout=10)
        records = _status_records(tmp_path)
        assert process.returncode == 0, (stdout, stderr)
        assert set(records) == CONFOUNDINGS
        assert {record["status"] for record in records.values()} == {"succeeded"}
        assert {record["exit_code"] for record in records.values()} == {0}
        assert {
            cf: record["pid"] for cf, record in records.items()
        } == {
            cf: int((tmp_path / "state" / f"{cf}.pid").read_text())
            for cf in CONFOUNDINGS
        }
        assert _launcher_status(tmp_path)["status"] == "succeeded"
        assert set(path.stem for path in (tmp_path / "state").glob("*.completed")) == CONFOUNDINGS
        assert _live_fake_pids(tmp_path) == []
    finally:
        if process.poll() is None:
            process.kill()
            process.wait(timeout=5)
        _cleanup_fake_pids(tmp_path)


def test_launcher_returns_nonzero_and_preserves_failed_child_identity(tmp_path):
    env = _launcher_env(tmp_path, jobs=2)
    env["FAKE_BENCHMARK_FAIL"] = "bin2cell"
    env["FAKE_BENCHMARK_FAIL_CODE"] = "7"
    try:
        result = subprocess.run(
            ["bash", str(LAUNCHER)],
            cwd=tmp_path,
            env=env,
            text=True,
            capture_output=True,
            timeout=10,
            check=False,
        )
        records = _status_records(tmp_path)
        assert result.returncode != 0
        assert set(records) == CONFOUNDINGS
        assert records["bin2cell"]["status"] == "failed"
        assert records["bin2cell"]["exit_code"] == 7
        assert _launcher_status(tmp_path)["status"] == "failed"
        assert all(
            record["status"] == "succeeded"
            for cf, record in records.items()
            if cf != "bin2cell"
        )
        assert set(
            path.stem for path in (tmp_path / "state").glob("*.started")
        ) == CONFOUNDINGS
        assert set(
            path.stem for path in (tmp_path / "state").glob("*.completed")
        ) == CONFOUNDINGS - {"bin2cell"}
        assert _live_fake_pids(tmp_path) == []
    finally:
        _cleanup_fake_pids(tmp_path)


def test_launcher_never_exceeds_configured_concurrency(tmp_path):
    env = _launcher_env(tmp_path, jobs=2, barrier=True)
    process = subprocess.Popen(
        ["bash", str(LAUNCHER)],
        cwd=tmp_path,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    try:
        _wait_for_started(tmp_path, 2)
        assert process.poll() is None
        time.sleep(0.2)
        assert len(list((tmp_path / "state").glob("*.started"))) == 2
        (tmp_path / "release").touch()
        stdout, stderr = process.communicate(timeout=10)
        state = json.loads((tmp_path / "state" / "concurrency.json").read_text())
        assert process.returncode == 0, (stdout, stderr)
        assert state == {"current": 0, "maximum": 2}
        records = _status_records(tmp_path)
        assert set(records) == CONFOUNDINGS
        assert {record["status"] for record in records.values()} == {"succeeded"}
        assert set(
            path.stem for path in (tmp_path / "state").glob("*.started")
        ) == CONFOUNDINGS
        assert set(
            path.stem for path in (tmp_path / "state").glob("*.completed")
        ) == CONFOUNDINGS
        assert _live_fake_pids(tmp_path) == []
    finally:
        if process.poll() is None:
            process.kill()
            process.wait(timeout=5)
        _cleanup_fake_pids(tmp_path)


@pytest.mark.parametrize(
    ("handled_signal", "expected_exit"),
    [(signal.SIGTERM, 143), (signal.SIGINT, 130)],
)
def test_launcher_forwards_interruption_and_leaves_no_children(
    tmp_path,
    handled_signal,
    expected_exit,
):
    env = _launcher_env(tmp_path, jobs=2, sleep=30)
    process = subprocess.Popen(
        ["bash", str(LAUNCHER)],
        cwd=tmp_path,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    try:
        _wait_for_started(tmp_path, 2)
        assert process.poll() is None
        process.send_signal(handled_signal)
        stdout, stderr = process.communicate(timeout=5)
        assert process.returncode == expected_exit, (stdout, stderr)

        deadline = time.monotonic() + 3
        while _live_fake_pids(tmp_path) and time.monotonic() < deadline:
            time.sleep(0.02)
        assert _live_fake_pids(tmp_path) == []

        records = _status_records(tmp_path)
        interrupted = [
            record for record in records.values() if record["status"] == "interrupted"
        ]
        skipped = [record for record in records.values() if record["status"] == "skipped"]
        assert len(interrupted) == 2
        assert len(skipped) == 4
        assert {record["exit_code"] for record in interrupted} == {
            expected_exit
        }
        assert _launcher_status(tmp_path)["status"] == "interrupted"
        interrupted_cfs = {record["confounding"] for record in interrupted}
        started_cfs = {
            path.stem for path in (tmp_path / "state").glob("*.started")
        }
        signaled_cfs = {
            path.stem for path in (tmp_path / "state").glob("*.signal")
        }
        assert started_cfs == signaled_cfs == interrupted_cfs
        for record in interrupted:
            cf = record["confounding"]
            assert record["pid"] == int(
                (tmp_path / "state" / f"{cf}.pid").read_text()
            )
            assert int(
                (tmp_path / "state" / f"{cf}.signal").read_text()
            ) == handled_signal
        assert all(record["pid"] is None for record in skipped)
    finally:
        if process.poll() is None:
            process.kill()
            process.wait(timeout=5)
        _cleanup_fake_pids(tmp_path)


@pytest.mark.parametrize("invalid_jobs", ["0", "not-a-number"])
def test_launcher_rejects_invalid_concurrency_before_starting_children(
    tmp_path,
    invalid_jobs,
):
    env = _launcher_env(tmp_path, jobs=1)
    env["BENCHMARK_MAX_JOBS"] = invalid_jobs
    result = subprocess.run(
        ["bash", str(LAUNCHER)],
        cwd=tmp_path,
        env=env,
        text=True,
        capture_output=True,
        timeout=5,
        check=False,
    )

    assert result.returncode != 0
    assert "must be a positive integer" in result.stderr
    assert not (tmp_path / "state").exists()
    assert not (tmp_path / "outputs" / "launcher_status.json").exists()


def test_launcher_cleans_lingering_process_group_after_leader_exit(tmp_path):
    env = _launcher_env(tmp_path, jobs=6)
    env["FAKE_BENCHMARK_GRANDCHILD_CF"] = "gene_panel"
    try:
        result = subprocess.run(
            ["bash", str(LAUNCHER)],
            cwd=tmp_path,
            env=env,
            text=True,
            capture_output=True,
            timeout=10,
            check=False,
        )
        records = _status_records(tmp_path)
        assert result.returncode != 0
        assert records["gene_panel"]["status"] == "failed"
        assert "process group" in records["gene_panel"]["error"]
        assert _live_fake_pids(tmp_path) == []
    finally:
        _cleanup_fake_pids(tmp_path)


def test_launcher_escalates_when_children_ignore_forwarded_signal(tmp_path):
    env = _launcher_env(tmp_path, jobs=2, sleep=30)
    env["FAKE_BENCHMARK_IGNORE_SIGNALS"] = "1"
    process = subprocess.Popen(
        ["bash", str(LAUNCHER)],
        cwd=tmp_path,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    try:
        _wait_for_started(tmp_path, 2)
        process.send_signal(signal.SIGTERM)
        stdout, stderr = process.communicate(timeout=6)
        assert process.returncode == 143, (stdout, stderr)
        records = _status_records(tmp_path)
        interrupted = [
            record for record in records.values() if record["status"] == "interrupted"
        ]
        assert len(interrupted) == 2
        assert {record["exit_code"] for record in interrupted} == {-signal.SIGKILL}
        assert _live_fake_pids(tmp_path) == []
    finally:
        if process.poll() is None:
            process.kill()
            process.wait(timeout=5)
        _cleanup_fake_pids(tmp_path)


def test_launcher_error_persists_failed_terminal_status(tmp_path):
    env = _launcher_env(tmp_path, jobs=2)
    (tmp_path / "0_records").write_text("path collision")
    result = subprocess.run(
        ["bash", str(LAUNCHER)],
        cwd=tmp_path,
        env=env,
        text=True,
        capture_output=True,
        timeout=5,
        check=False,
    )

    assert result.returncode != 0
    manifest = _launcher_status(tmp_path)
    assert manifest["status"] == "failed"
    assert manifest["launcher_error"]["type"] == "NotADirectoryError"
    assert manifest["tasks"][0]["status"] == "failed"
    assert {task["status"] for task in manifest["tasks"][1:]} == {"skipped"}


def test_blank_sample_parts_preserves_default_part1_task_set(tmp_path):
    env = _launcher_env(tmp_path, jobs=6)
    env["SAMPLE_PARTS"] = "   "
    try:
        result = subprocess.run(
            ["bash", str(LAUNCHER)],
            cwd=tmp_path,
            env=env,
            text=True,
            capture_output=True,
            timeout=10,
            check=False,
        )
        manifest = _launcher_status(tmp_path)
        assert result.returncode == 0, result.stderr
        assert len(manifest["tasks"]) == 6
        assert {task["sample_name"] for task in manifest["tasks"]} == {
            "P2CRC/cut_part1"
        }
    finally:
        _cleanup_fake_pids(tmp_path)


def test_empty_legacy_environment_values_preserve_defaults(tmp_path):
    env = _launcher_env(tmp_path, jobs=6)
    env.update(
        {
            "RAW_DATA_PATH": "",
            "SAVE_PATH": "",
            "SAMPLE_PATIENT": "",
            "CONFIG_PATH": "",
        }
    )
    try:
        result = subprocess.run(
            ["bash", str(LAUNCHER)],
            cwd=tmp_path,
            env=env,
            text=True,
            capture_output=True,
            timeout=10,
            check=False,
        )
        manifests = list(
            (tmp_path / "results_unified" / "benchmark_runs").rglob(
                "launcher_status.json"
            )
        )
        assert result.returncode == 0, result.stderr
        assert len(manifests) == 1
        manifest = json.loads(manifests[0].read_text())
        assert {task["sample_name"] for task in manifest["tasks"]} == {
            "P2CRC/cut_part1"
        }
        args = json.loads(
            (tmp_path / "state" / "segmentation.args.json").read_text()
        )
        assert args[args.index("--data-root") + 1] == "./raw_data/Sim2Real-ST"
        assert args[args.index("--config") + 1] == "revise/revise.yaml"
        assert args[args.index("--output-root") + 1].startswith(
            "results_unified/benchmark_runs/"
        )
    finally:
        _cleanup_fake_pids(tmp_path)


def test_interruption_cleans_active_child_process_group(tmp_path):
    env = _launcher_env(tmp_path, jobs=2, sleep=30)
    env["FAKE_BENCHMARK_ACTIVE_GRANDCHILD_CF"] = "segmentation"
    process = subprocess.Popen(
        ["bash", str(LAUNCHER)],
        cwd=tmp_path,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    try:
        _wait_for_started(tmp_path, 2)
        _wait_for_path(tmp_path / "state" / "segmentation.grandchild.pid")
        process.send_signal(signal.SIGTERM)
        stdout, stderr = process.communicate(timeout=6)
        assert process.returncode == 143, (stdout, stderr)
        assert _live_fake_pids(tmp_path) == []
        records = _status_records(tmp_path)
        assert records["segmentation"]["status"] == "interrupted"
    finally:
        if process.poll() is None:
            process.kill()
            process.wait(timeout=5)
        _cleanup_fake_pids(tmp_path)


def test_launcher_error_cleans_already_active_children(tmp_path):
    env = _launcher_env(tmp_path, jobs=7, barrier=True)
    env["SAMPLE_PARTS"] = "part1 part2"
    (tmp_path / "0_records").mkdir()
    (tmp_path / "0_records" / "P2CRC_part2").write_text("path collision")
    try:
        result = subprocess.run(
            ["bash", str(LAUNCHER)],
            cwd=tmp_path,
            env=env,
            text=True,
            capture_output=True,
            timeout=8,
            check=False,
        )
        manifest = _launcher_status(tmp_path)
        part1 = [
            task
            for task in manifest["tasks"]
            if task["sample_name"] == "P2CRC/cut_part1"
        ]
        part2 = [
            task
            for task in manifest["tasks"]
            if task["sample_name"] == "P2CRC/cut_part2"
        ]
        assert result.returncode != 0
        assert manifest["status"] == "failed"
        assert {task["status"] for task in part1} == {"interrupted"}
        assert part2[0]["status"] == "failed"
        assert {task["status"] for task in part2[1:]} == {"skipped"}
        assert all(task["pid"] is not None for task in part1)
        assert all(not _pid_is_alive(task["pid"]) for task in part1)
        assert _live_fake_pids(tmp_path) == []
    finally:
        _cleanup_fake_pids(tmp_path)


def test_multiple_sample_parts_preserve_all_task_identities(tmp_path):
    env = _launcher_env(tmp_path, jobs=6)
    env["SAMPLE_PARTS"] = "part1 part2"
    try:
        result = subprocess.run(
            ["bash", str(LAUNCHER)],
            cwd=tmp_path,
            env=env,
            text=True,
            capture_output=True,
            timeout=10,
            check=False,
        )
        manifest = _launcher_status(tmp_path)
        assert result.returncode == 0, result.stderr
        assert len(manifest["tasks"]) == 12
        assert {task["id"] for task in manifest["tasks"]} == {
            f"P2CRC/cut_{part}|{cf}"
            for part in ("part1", "part2")
            for cf in CONFOUNDINGS
        }
        assert {task["status"] for task in manifest["tasks"]} == {"succeeded"}
        observed_identities = set()
        for path in (tmp_path / "state").glob("*.args.json"):
            args = json.loads(path.read_text())
            observed_identities.add(
                (
                    args[args.index("--sample-name") + 1],
                    args[args.index("--confounding") + 1],
                )
            )
        assert observed_identities == {
            (f"P2CRC/cut_{part}", cf)
            for part in ("part1", "part2")
            for cf in CONFOUNDINGS
        }
        assert len(list((tmp_path / "0_records" / "P2CRC_part1").glob("*.log"))) == 6
        assert len(list((tmp_path / "0_records" / "P2CRC_part2").glob("*.log"))) == 6
    finally:
        _cleanup_fake_pids(tmp_path)
