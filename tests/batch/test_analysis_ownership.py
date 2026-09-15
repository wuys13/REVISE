"""Ownership and publication boundaries for analysis artifacts."""
import json

import pytest

from test_runner import fake_solver as _fake_solver


fake_solver = _fake_solver


def test_staging_rejects_files_not_declared_by_the_adapter(tmp_path):
    from revise.batch import analysis

    stage = tmp_path / "stage"
    stage.mkdir()
    (stage / "declared.csv").write_text("value\n1\n")
    (stage / "extra.csv").write_text("value\n2\n")

    with pytest.raises(ValueError, match="undeclared|artifact"):
        analysis._stage_artifacts(
            stage,
            tmp_path / "published",
            {"declared": {"path": "declared.csv", "description": "declared"}},
        )


def test_reuse_rejects_an_extra_published_file(tmp_path):
    from revise.batch import analysis
    from revise.batch.sample import file_identity

    destination = tmp_path / "analysis" / "aspect"
    destination.mkdir(parents=True)
    artifact = destination / "result.csv"
    artifact.write_text("value\n1\n")
    previous = {
        "status": "succeeded",
        "fingerprint": "same",
        "published_directory": str(destination),
        "artifact_records": [file_identity(artifact)],
    }
    (destination / "user-file.txt").write_text("keep")

    assert not analysis._analysis_reusable(previous, "same", destination)


def test_reuse_rejects_an_extra_empty_published_directory(tmp_path):
    from revise.batch import analysis
    from revise.batch.sample import file_identity

    destination = tmp_path / "analysis" / "aspect"
    destination.mkdir(parents=True)
    artifact = destination / "result.csv"
    artifact.write_text("value\n1\n")
    previous = {
        "status": "succeeded",
        "fingerprint": "same",
        "published_directory": str(destination),
        "artifact_records": [file_identity(artifact)],
    }
    (destination / "user-empty-directory").mkdir()

    assert not analysis._analysis_reusable(previous, "same", destination)


def test_empty_analysis_does_not_create_a_task_placeholder(tmp_path, fake_solver):
    from revise.batch.analysis import run_analysis_batch
    from revise.batch.runner import run_batch
    from test_runner import batch_config, make_sample, result_root

    sample = tmp_path / "data" / "one"
    make_sample(sample, cell_types=["T"])
    config = batch_config(tmp_path)
    assert run_batch(config)["summary"]["succeeded"] == 1
    placeholder = result_root(sample) / "T" / "analysis" / "analysis.json"
    placeholder.unlink()

    result = run_analysis_batch(config)

    assert result["summary"] == {
        "succeeded": 0,
        "reused": 0,
        "failed": 0,
        "blocked": 0,
        "unavailable": 0,
        "not_implemented": 0,
        "inactive": 0,
    }
    assert not placeholder.exists()
    assert json.loads((tmp_path / "results" / "analysis_status.json").read_text())["tasks"] == []


def test_empty_analysis_does_not_verify_or_publish_a_broken_handoff(tmp_path, fake_solver):
    from revise.batch.analysis import run_analysis_batch
    from revise.batch.runner import run_batch
    from test_runner import batch_config, make_sample, result_root

    sample = tmp_path / "data" / "one"
    make_sample(sample, cell_types=["T"])
    config = batch_config(tmp_path)
    assert run_batch(config)["summary"]["succeeded"] == 1
    root = result_root(sample) / "T"
    (root / "spatial.h5ad").write_bytes(b"broken")
    placeholder = root / "analysis" / "analysis.json"
    placeholder.unlink()

    result = run_analysis_batch(config)

    assert result["tasks"] == []
    assert not placeholder.exists()
