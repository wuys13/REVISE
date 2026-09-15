"""Contract tests for published reconstruction-impact reports."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from revise.analysis.impact_report import render_impact_report


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _published_task(tmp_path: Path, *, status: str = "succeeded", error: str | None = None) -> tuple[Path, Path]:
    task_root = tmp_path / "sample-1" / "T-cell"
    published = task_root / "analysis" / "reconstruction_impact"
    published.mkdir(parents=True)
    frames = {
        "comparisons": pd.DataFrame(
            [
                {
                    "comparison_id": "cmp",
                    "sample_id": "sample-1",
                    "task_cell_type": "T-cell",
                    "scope": "T-cell",
                    "comparison_edge": "spatial_context",
                }
            ]
        ),
        "partition_summary": pd.DataFrame(
            [{"comparison_id": "cmp::partition::edge", "comparison_edge": "edge", "matched_accuracy": 1.0}]
        ),
        "partition_contingency": pd.DataFrame(
            [{"comparison_id": "cmp::partition::edge", "comparison_edge": "edge", "raw_cluster": "a", "recon_cluster": "x", "normalization": "absolute", "value": 1}]
        ),
        "spatial_window_metrics": pd.DataFrame(
            [{"comparison_id": "cmp", "window_id": "w1", "window_x": 1.0, "window_y": 2.0, "metric": "neff", "baseline": "raw_leiden", "delta": 0.5, "status": "computed"}]
        ),
        "spatial_region_extent_by_anatomy": pd.DataFrame(
            [{"comparison_id": "cmp", "level1_region": "Tumor", "region_type": "state", "region_area_um2": 4.0}]
        ),
        "anatomy_summary": pd.DataFrame(
            [{"comparison_id": "cmp", "level1_region": "Tumor", "tissue_windows": 1}]
        ),
        "raw_level2_assignments": pd.DataFrame(
            [{"comparison_id": "cmp", "unit_id": "u1", "raw_level2": "L2a", "confidence": 0.9}]
        ),
        "raw_level2_posterior": pd.DataFrame(
            [{"comparison_id": "cmp", "unit_id": "u1", "raw_level2": "L2a", "posterior": 0.9}]
        ),
        "audit": pd.DataFrame(),
    }
    records = []
    roles = {}
    for role, frame in frames.items():
        relative = {
            "comparisons": "comparisons.csv",
            "partition_summary": "partition/summary.csv",
            "partition_contingency": "partition/contingency.csv",
            "spatial_window_metrics": "spatial/window_metrics.csv",
            "spatial_region_extent_by_anatomy": "spatial/region_extent_by_anatomy.csv",
            "anatomy_summary": "anatomy/summary.csv",
            "raw_level2_assignments": "raw_level2/assignments.csv",
            "raw_level2_posterior": "raw_level2/posterior.csv",
            "audit": "audit.csv",
        }[role]
        path = published / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(path, index=False)
        roles[role] = relative
        records.append({"path": str(path), "sha256": _sha256(path), "role": role, "description": role})
    state = {
        "status": status,
        "sample_id": "sample-1",
        "cell_type": "T-cell",
        "published_directory": str(published),
        "artifact_records": records,
        "roles": roles,
    }
    if error:
        state["error"] = error
    control = task_root / ".revise" / "analysis"
    control.mkdir(parents=True)
    (control / "reconstruction_impact.json").write_text(json.dumps(state, indent=2), encoding="utf-8")
    return task_root, published


def test_report_rejects_damaged_published_artifact(tmp_path: Path):
    task_root, published = _published_task(tmp_path)
    (published / "partition" / "summary.csv").write_text("damaged\n", encoding="utf-8")

    with pytest.raises(ValueError, match="hash|artifact|ownership"):
        render_impact_report([task_root], tmp_path / "report", execute=False)


def test_report_marks_latest_failure_and_old_result_stale(tmp_path: Path):
    task_root, _ = _published_task(
        tmp_path,
        status="failed",
        error="RuntimeError: dependency unavailable",
    )

    manifest = render_impact_report([task_root], tmp_path / "report", execute=False)

    assert manifest["status"] == "failed"
    assert manifest["tasks"][0]["latest_status"] == "failed"
    assert manifest["tasks"][0]["stale"] is True
    assert "dependency unavailable" in manifest["tasks"][0]["error"]
    assert json.loads((tmp_path / "report" / "manifest.json").read_text())["status"] == "failed"


def test_report_rejects_non_owned_report_directory(tmp_path: Path):
    task_root, _ = _published_task(tmp_path)
    report = tmp_path / "report"
    report.mkdir()
    (report / "user-file.txt").write_text("keep", encoding="utf-8")

    with pytest.raises(ValueError, match="owned|report directory"):
        render_impact_report([task_root], report, execute=False)


def test_report_rejects_tampered_previous_manifest_path(tmp_path: Path):
    task_root, _ = _published_task(tmp_path)
    report = tmp_path / "report"
    report.mkdir()
    outside = tmp_path / "outside.txt"
    outside.write_text("keep", encoding="utf-8")
    (report / "manifest.json").write_text(
        json.dumps({
            "status": "succeeded",
            "artifacts": [{"path": "../outside.txt", "sha256": _sha256(outside)}],
        }),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="path|manifest|owned"):
        render_impact_report([task_root], report, execute=False)
    assert outside.read_text(encoding="utf-8") == "keep"


@pytest.mark.parametrize("extra", ["nested/manifest.json", "empty-dir"])
def test_report_rejects_unowned_nested_files_and_directories(tmp_path: Path, extra: str):
    task_root, _ = _published_task(tmp_path)
    report = tmp_path / "report"
    report.mkdir()
    path = report / extra
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix:
        path.write_text("unowned", encoding="utf-8")
    else:
        path.mkdir()
    (report / "manifest.json").write_text(
        json.dumps({"status": "succeeded", "artifacts": []}),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="unowned"):
        render_impact_report([task_root], report, execute=False)


def test_report_rejects_destination_symlink(tmp_path: Path):
    task_root, _ = _published_task(tmp_path)
    real = tmp_path / "real-report"
    real.mkdir()
    report = tmp_path / "report"
    report.symlink_to(real, target_is_directory=True)

    with pytest.raises(ValueError, match="real directory|symlink"):
        render_impact_report([task_root], report, execute=False)


def test_report_rejects_copied_input_hash_mismatch(monkeypatch, tmp_path: Path):
    from revise.analysis import impact_report

    task_root, _ = _published_task(tmp_path)
    original_copy = impact_report.shutil.copy2

    def copy_bad(source, target):
        if "inputs" in Path(target).parts:
            Path(target).write_bytes(b"tampered")
        else:
            original_copy(source, target)

    monkeypatch.setattr(impact_report.shutil, "copy2", copy_bad)
    with pytest.raises(ValueError, match="copied|hash"):
        render_impact_report([task_root], tmp_path / "report", execute=False)


def test_report_rejects_colliding_task_slugs(tmp_path: Path):
    from revise.analysis import impact_report

    stage = tmp_path / "stage"
    stage.mkdir()
    tasks = [{"task_id": "a/b", "artifacts": {}}, {"task_id": "a_b", "artifacts": {}}]
    with pytest.raises(ValueError, match="slug"):
        impact_report._copy_inputs(stage, tasks)


def test_report_publish_rolls_back_previous_directory(monkeypatch, tmp_path: Path):
    from revise.analysis import impact_report

    destination = tmp_path / "report"
    destination.mkdir()
    (destination / "old.txt").write_text("old", encoding="utf-8")
    stage = tmp_path / "stage"
    stage.mkdir()
    (stage / "new.txt").write_text("new", encoding="utf-8")
    original_replace = impact_report.os.replace
    calls = 0

    def fail_during_publish(source, target):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise RuntimeError("publish interrupted")
        return original_replace(source, target)

    monkeypatch.setattr(impact_report.os, "replace", fail_during_publish)
    with pytest.raises(RuntimeError, match="publish interrupted"):
        impact_report._publish_stage(stage, destination)
    assert (destination / "old.txt").read_text(encoding="utf-8") == "old"
    assert not (destination / "new.txt").exists()


def test_report_identity_uses_reconstruction_metadata_for_global_task(tmp_path: Path):
    task_root, _ = _published_task(tmp_path)
    state_path = task_root / ".revise" / "analysis" / "reconstruction_impact.json"
    state = json.loads(state_path.read_text())
    state.pop("sample_id")
    state.pop("cell_type")
    state_path.write_text(json.dumps(state, indent=2), encoding="utf-8")
    (task_root / "reconstruction.json").write_text(
        json.dumps({"sample_id": "sample-1", "cell_type": None, "modality": "hST"}),
        encoding="utf-8",
    )

    manifest = render_impact_report([task_root], tmp_path / "report", execute=False)

    assert manifest["tasks"][0]["task_id"] == "sample-1"
    assert manifest["tasks"][0]["cell_type"] is None


def _two_scope_task(tmp_path: Path) -> tuple[Path, Path]:
    task_root, published = _published_task(tmp_path)
    tables = {
        "comparisons.csv": pd.DataFrame(
            [
                {"comparison_id": "cmp-a", "sample_id": "sample-1", "task_cell_type": "T-cell", "scope": "scope-a", "comparison_edge": "spatial_context"},
                {"comparison_id": "cmp-b", "sample_id": "sample-1", "task_cell_type": "T-cell", "scope": "scope-b", "comparison_edge": "spatial_context"},
            ]
        ),
        "partition/summary.csv": pd.DataFrame(
            [
                {"comparison_id": "cmp-a", "comparison_edge": "edge", "scope": "scope-a", "matched_accuracy": 0.1},
                {"comparison_id": "cmp-b", "comparison_edge": "edge", "scope": "scope-b", "matched_accuracy": 0.9},
            ]
        ),
        "partition/contingency.csv": pd.DataFrame(
            [
                {"comparison_id": "cmp-a", "comparison_edge": "edge", "scope": "scope-a", "raw_cluster": "0", "recon_cluster": "0", "normalization": "absolute", "value": 1},
                {"comparison_id": "cmp-a", "comparison_edge": "edge", "scope": "scope-a", "raw_cluster": "1", "recon_cluster": "1", "normalization": "absolute", "value": 2},
                {"comparison_id": "cmp-b", "comparison_edge": "edge", "scope": "scope-b", "raw_cluster": "0", "recon_cluster": "0", "normalization": "absolute", "value": 10},
                {"comparison_id": "cmp-b", "comparison_edge": "edge", "scope": "scope-b", "raw_cluster": "1", "recon_cluster": "1", "normalization": "absolute", "value": 20},
            ]
        ),
        "spatial/window_metrics.csv": pd.DataFrame(
            [
                {"comparison_id": "cmp-a", "metric": "neff", "baseline": "raw_leiden", "scale": 12.0, "status": "computed", "window_x": 0.0, "window_y": 0.0, "delta": 1.0},
                {"comparison_id": "cmp-b", "metric": "neff", "baseline": "raw_leiden", "scale": 12.0, "status": "computed", "window_x": 10.0, "window_y": 10.0, "delta": 9.0},
            ]
        ),
        "spatial/region_extent_by_anatomy.csv": pd.DataFrame(
            [
                {"comparison_id": "cmp-a", "level1_region": "Tumor", "region_type": "state", "region_area_um2": 1.0},
                {"comparison_id": "cmp-b", "level1_region": "Tumor", "region_type": "state", "region_area_um2": 9.0},
            ]
        ),
        "raw_level2/assignments.csv": pd.DataFrame(
            [
                {"comparison_id": "cmp-a", "unit_id": "u1", "raw_level2": "L2a"},
                {"comparison_id": "cmp-b", "unit_id": "u2", "raw_level2": "L2b"},
            ]
        ),
    }
    for relative, frame in tables.items():
        path = published / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(path, index=False)
    state_path = task_root / ".revise" / "analysis" / "reconstruction_impact.json"
    state = json.loads(state_path.read_text())
    for record in state["artifact_records"]:
        record["sha256"] = _sha256(Path(record["path"]))
    state_path.write_text(json.dumps(state, indent=2), encoding="utf-8")
    return task_root, published


def _figure_summary(report: Path) -> list[dict]:
    notebook = json.loads((report / "report.ipynb").read_text())
    for cell in notebook["cells"]:
        if cell.get("cell_type") != "code":
            continue
        for output in cell.get("outputs", []):
            if output.get("output_type") != "stream":
                continue
            text = output.get("text", "")
            text = "".join(text) if isinstance(text, list) else text
            for line in text.splitlines():
                if line.startswith("FIGURE_SUMMARY_JSON="):
                    return json.loads(line.split("=", 1)[1])
    raise AssertionError("Executed report did not expose figure summary")


def test_report_keeps_figures_separate_by_comparison_keys(tmp_path: Path):
    task_root, _ = _two_scope_task(tmp_path)

    render_impact_report([task_root], tmp_path / "report", execute=True)
    summary = _figure_summary(tmp_path / "report")

    partition = [item for item in summary if item["kind"] == "partition_contingency"]
    assert {item["comparison_id"]: item["values"] for item in partition} == {
        "cmp-a": ["3"],
        "cmp-b": ["30"],
    }
    windows = [item for item in summary if item["kind"] == "spatial_windows"]
    assert {item["comparison_id"]: item["values"] for item in windows} == {
        "cmp-a": ["1"],
        "cmp-b": ["9"],
    }
    extent = [item for item in summary if item["kind"] == "anatomy_extent"]
    assert {item["comparison_id"]: item["values"] for item in extent} == {
        "cmp-a": ["1"],
        "cmp-b": ["9"],
    }
    level2 = [item for item in summary if item["kind"] == "raw_level2"]
    assert {item["comparison_id"]: item["values"] for item in level2} == {
        "cmp-a": ["1"],
        "cmp-b": ["1"],
    }


def test_report_executes_from_copied_inputs_and_registers_hashes(tmp_path: Path):
    task_root, published = _published_task(tmp_path)
    science_before = {path: _sha256(path) for path in published.rglob("*") if path.is_file()}

    manifest = render_impact_report([task_root], tmp_path / "report", execute=True)
    report = tmp_path / "report"

    assert manifest["status"] == "succeeded"
    assert (report / "inputs.json").is_file()
    assert (report / "source.ipynb").is_file()
    assert (report / "report.ipynb").is_file()
    assert manifest["consumed_artifacts"]
    assert manifest["artifacts"]
    for item in manifest["artifacts"]:
        assert (report / item["path"]).is_file()
        assert _sha256(report / item["path"]) == item["sha256"]
    source = json.loads((report / "source.ipynb").read_text())
    executed = json.loads((report / "report.ipynb").read_text())
    assert [cell["source"] for cell in source["cells"] if cell["cell_type"] == "code"] == [
        cell["source"] for cell in executed["cells"] if cell["cell_type"] == "code"
    ]
    assert list((report / "figures").glob("*.png"))
    assert science_before == {path: _sha256(path) for path in published.rglob("*") if path.is_file()}


def test_report_execution_failure_is_persisted_as_failed(monkeypatch, tmp_path: Path):
    from revise.analysis import impact_report

    task_root, _ = _published_task(tmp_path)

    def fail(*_args, **_kwargs):
        raise RuntimeError("kernel unavailable")

    monkeypatch.setattr(impact_report, "_execute_notebook", fail)
    manifest = render_impact_report([task_root], tmp_path / "report", execute=True)

    assert manifest["status"] == "failed"
    assert "kernel unavailable" in manifest["error"]
    assert json.loads((tmp_path / "report" / "manifest.json").read_text())["status"] == "failed"


def test_report_execution_failure_keeps_previous_published_files(monkeypatch, tmp_path: Path):
    from revise.analysis import impact_report

    task_root, _ = _published_task(tmp_path)
    report = tmp_path / "report"
    first = render_impact_report([task_root], report, execute=False)
    previous_files = {
        item["path"]: _sha256(report / item["path"])
        for item in first["artifacts"]
    }

    def fail(*_args, **_kwargs):
        raise RuntimeError("kernel unavailable")

    monkeypatch.setattr(impact_report, "_execute_notebook", fail)
    manifest = render_impact_report([task_root], report, execute=True)

    assert manifest["status"] == "failed"
    assert all(_sha256(report / path) == digest for path, digest in previous_files.items())
