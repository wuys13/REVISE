"""Focused tests for the single analysis task API."""
from contextlib import nullcontext
import json
from types import SimpleNamespace

import pytest


def structured_adapter(context):
    (context.output_dir / "result.csv").write_text("value\n1\n")
    return {
        "artifacts": {
            "result": {"path": "result.csv", "description": "test result"},
        },
        "calculation": {
            "input_view": "test",
            "parameters": {},
            "comparison_basis": "test comparison",
        },
    }


def resource_changes_during_run(context):
    context.resources["marker"].write_text("changed")
    return structured_adapter(context)


def adapter_for_code_change(context):
    return structured_adapter(context)


def _resolved(tmp_path, document):
    input_root = tmp_path / "data"
    sample_root = input_root / "one"
    sample_root.mkdir(parents=True)
    output_root = tmp_path / "results"
    output_root.mkdir()
    return SimpleNamespace(
        root=sample_root,
        sample_id="one",
        document=document,
        config_chain=[],
        project_config=tmp_path / "batch.yaml",
        input_root=input_root,
        output_root=output_root,
    ), input_root, output_root


def _handoff(output_root):
    return {
        "schema_version": 1,
        "status": "succeeded",
        "sample_id": "one",
        "cell_type": None,
        "directory": str(output_root / "one"),
        "fingerprint": "reconstruction",
        "inputs": {"sources": {}},
        "outputs": {},
        "pairing": {"status": "available"},
    }


def _patch_task(monkeypatch, tmp_path, document):
    from revise.batch import analysis

    resolved, input_root, output_root = _resolved(tmp_path, document)
    config = tmp_path / "batch.yaml"
    config.write_text("schema_version: 2\n")
    monkeypatch.setattr(analysis, "load_batch_config", lambda _: ({}, input_root, output_root))
    monkeypatch.setattr(analysis, "_selected_sample", lambda *_: resolved.root)
    monkeypatch.setattr(analysis, "resolve_sample", lambda *_: resolved)
    monkeypatch.setattr(analysis, "_batch_lock", lambda _: nullcontext())
    monkeypatch.setattr(analysis, "verify_reconstruction_task",
                        lambda current, cell_type: _handoff(output_root))
    return analysis, config, output_root


def test_single_analysis_does_not_require_batch_inventory(tmp_path, monkeypatch):
    analysis, config, output_root = _patch_task(monkeypatch, tmp_path, {
        "modality": "hST",
        "analysis": {
            "partition": {
                "entrypoint": "test_analysis_api:structured_adapter",
                "version": "1",
            },
        },
    })

    result = analysis.run_analysis_task(config, "one", "partition")

    assert result["status"] == "succeeded"
    assert result["calculation"] == {
        "input_view": "test",
        "parameters": {},
        "comparison_basis": "test comparison",
    }
    assert result["artifacts"]["result"]["path"] == "result.csv"
    assert result["artifact_records"][0]["path"].endswith("analysis/partition/result.csv")
    assert result["configuration_audit"]["chain"] == []
    assert result["log_path"].endswith(".revise/analysis/partition.log")
    assert not (output_root / "batch_status.json").exists()


def test_empty_analysis_does_not_create_placeholder_task(tmp_path, monkeypatch):
    analysis, config, _ = _patch_task(monkeypatch, tmp_path, {"modality": "hST", "analysis": {}})

    with pytest.raises(ValueError, match="configured"):
        analysis.run_analysis_task(config, "one", "partition")


def test_casefold_colliding_aspect_names_are_rejected():
    from revise.batch.analysis import _specifications

    with pytest.raises(ValueError, match="collide case-insensitively"):
        _specifications({"analysis": {"Partition": {}, "partition": {}}})


def test_artifact_description_is_required():
    from revise.batch.analysis import _normalise_adapter_result

    with pytest.raises(ValueError, match="description"):
        _normalise_adapter_result({
            "artifacts": {"result": {"path": "result.csv"}},
            "calculation": {
                "input_view": "test", "parameters": {}, "comparison_basis": "test",
            },
        })


def test_reuse_rejects_artifact_outside_current_aspect(tmp_path):
    from revise.batch.analysis import _analysis_reusable
    from revise.batch.sample import file_identity

    destination = tmp_path / "analysis" / "partition"
    destination.mkdir(parents=True)
    artifact = tmp_path / "outside.csv"
    artifact.write_text("outside")
    identity = file_identity(artifact)
    previous = {
        "status": "succeeded",
        "fingerprint": "fingerprint",
        "published_directory": str(destination),
        "artifact_records": [identity],
    }

    assert not _analysis_reusable(previous, "fingerprint", destination)


def test_disabled_analysis_is_inactive_without_running_adapter(tmp_path, monkeypatch):
    analysis, config, output_root = _patch_task(monkeypatch, tmp_path, {
        "modality": "hST",
        "analysis": {
            "partition": {"enabled": False},
        },
    })

    result = analysis.run_analysis_task(config, "one", "partition")

    assert result["status"] == "inactive"
    assert result["reason"] == "disabled"
    assert not (output_root / "one" / "analysis" / "partition").exists()


def test_single_summary_only_reports_requested_aspect(tmp_path, monkeypatch):
    document = {
        "modality": "hST",
        "analysis": {
            "first": {"entrypoint": "test_analysis_api:structured_adapter", "version": "1"},
            "second": {"entrypoint": "test_analysis_api:structured_adapter", "version": "1"},
        },
    }
    analysis, config, output_root = _patch_task(monkeypatch, tmp_path, document)

    assert analysis.run_analysis_task(config, "one", "first")["status"] == "succeeded"
    assert analysis.run_analysis_task(config, "one", "second")["status"] == "succeeded"

    summary = json.loads((output_root / "one" / "analysis" / "analysis.json").read_text())
    assert [item["aspect"] for item in summary["aspects"]] == ["second"]
    assert (output_root / "one" / "analysis" / "first" / "result.csv").is_file()


def test_resource_change_during_run_is_rejected_before_publish(tmp_path, monkeypatch):
    marker = tmp_path / "marker.txt"
    marker.write_text("before")
    analysis, config, output_root = _patch_task(monkeypatch, tmp_path, {
        "modality": "hST",
        "analysis": {
            "partition": {
                "entrypoint": "test_analysis_api:resource_changes_during_run",
                "version": "1",
                "resources": {"marker": str(marker)},
            },
        },
    })

    result = analysis.run_analysis_task(config, "one", "partition")

    assert result["status"] == "failed"
    assert "resources or code changed" in result["error"]
    assert not (output_root / "one" / "analysis" / "partition" / "result.csv").exists()
    assert marker.read_text() == "changed"


def test_code_change_during_run_is_rejected_before_publish(tmp_path, monkeypatch):
    analysis, config, output_root = _patch_task(monkeypatch, tmp_path, {
        "modality": "hST",
        "analysis": {
            "partition": {
                "entrypoint": "test_analysis_api:adapter_for_code_change",
                "version": "1",
            },
        },
    })
    identities = iter([{"sha256": "before"}, {"sha256": "after"}])
    monkeypatch.setattr(analysis, "_analysis_code_identity", lambda source: next(identities))

    result = analysis.run_analysis_task(config, "one", "partition")

    assert result["status"] == "failed"
    assert "resources or code changed" in result["error"]
    assert not (output_root / "one" / "analysis" / "partition" / "result.csv").exists()


def test_disabled_aspect_keeps_prior_files_and_invalidates_state(tmp_path, monkeypatch):
    document = {
        "modality": "hST",
        "analysis": {
            "partition": {
                "entrypoint": "test_analysis_api:structured_adapter", "version": "1",
            },
        },
    }
    analysis, config, output_root = _patch_task(monkeypatch, tmp_path, document)
    assert analysis.run_analysis_task(config, "one", "partition")["status"] == "succeeded"
    artifact = output_root / "one" / "analysis" / "partition" / "result.csv"
    document["analysis"]["partition"]["enabled"] = False

    result = analysis.run_analysis_task(config, "one", "partition")

    state = json.loads((output_root / "one" / ".revise" / "analysis" / "partition.json").read_text())
    assert result["status"] == "inactive"
    assert state["status"] == "inactive"
    assert artifact.read_text() == "value\n1\n"


def test_removed_aspect_keeps_files_but_marks_control_inactive(tmp_path):
    from revise.batch.analysis import _invalidate_removed_aspects

    root = tmp_path / "one"
    state = root / ".revise" / "analysis" / "partition.json"
    artifact = root / "analysis" / "partition" / "result.csv"
    state.parent.mkdir(parents=True)
    artifact.parent.mkdir(parents=True)
    state.write_text(json.dumps({"status": "succeeded", "published_directory": str(artifact.parent)}))
    artifact.write_text("keep")

    _invalidate_removed_aspects(root, {})

    assert json.loads(state.read_text())["status"] == "inactive"
    assert artifact.read_text() == "keep"
