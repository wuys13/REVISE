from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise.backend.ops.assignment import AssignmentState, one_hot_assignment
from revise.config.loader import ResolvedConfig
from revise.framework import REVISEPipeline
from revise.recon.context import PipelineContext
from revise.utils import provenance
from revise.utils.io import build_run_dir


def _small_adata(*, obs_names=("cell-1", "cell-2")) -> AnnData:
    return AnnData(
        X=np.ones((len(obs_names), 2)),
        obs=pd.DataFrame(index=pd.Index(obs_names)),
        var=pd.DataFrame(index=pd.Index(["g1", "g2"])),
    )


def _write_application_inputs(data_root: Path, sample_name: str) -> None:
    st = _small_adata(obs_names=("spot-1", "spot-2"))
    st.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 1.0]])
    sc_ref = _small_adata()
    sc_ref.obs["Level1"] = ["A", "B"]
    st.write_h5ad(data_root / f"{sample_name}_Xenium.h5ad")
    sc_ref.write_h5ad(data_root / "adata_sc_all_reanno.h5ad")


def _write_benchmark_seg_inputs(data_root: Path, sample_name: str) -> None:
    st = _small_adata()
    st.obs["seg_error"] = [0, 0]
    st.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 1.0]])
    sc_ref = _small_adata()
    sc_ref.obs["Level1"] = ["A", "B"]
    gt = _small_adata()
    st_path = data_root / sample_name / "seg_1" / "Xenium.h5ad"
    st_path.parent.mkdir(parents=True, exist_ok=True)
    st.write_h5ad(st_path)
    sc_ref.write_h5ad(data_root / sample_name / "adata_sc_all_reanno.h5ad")
    gt.write_h5ad(data_root / sample_name / "selected_xenium.h5ad")


def _guidance_context(tmp_path) -> PipelineContext:
    config = ResolvedConfig(
        {
            "io": {"save_outputs": False},
            "ot": {
                "ga": {"solver": "pot"},
                "lr": {"solver": "pot"},
            },
            "local_refinement": {
                "guidance": "prefer",
                "compatibility": {
                    "mode": "cost",
                    "beta": 1.0,
                    "min_affinity": 0.05,
                    "strength": 0.2,
                },
            },
        },
        request_evidence={
            "assignment_guidance": {
                "configured_guidance": None,
                "configured_compatibility_mode": None,
                "resolution_source": "route_default",
                "deprecations": [],
            }
        },
    )
    return PipelineContext(
        merged_config=config,
        raw_config={},
        config_path="revise/revise.yaml",
        profile="test",
        runtime={"mode": "application", "svc_kind": "sp"},
        route_key="sp_svc:bin2cell",
        run_dir=tmp_path,
        logger=logging.getLogger("test-assignment-guidance-manifest"),
    )


def test_assignment_guidance_callback_writes_canonical_request_and_events(tmp_path):
    ctx = _guidance_context(tmp_path)
    pipeline = REVISEPipeline()
    ctx.set_provenance_callback(pipeline._write_final_metadata)
    left_assignment = AssignmentState(
        values=np.array([[0.8, 0.2], [0.1, 0.9]]),
        observation_labels=("spot-1", "spot-2"),
        category_labels=("A", "B"),
        source="global_anchoring",
        level="Level1",
        value_semantics="soft",
        lineage=[{"operation": "aggregate"}],
    )
    right_assignment = one_hot_assignment(
        ("A", "B"),
        observation_labels=("ref-1", "ref-2"),
        category_labels=("A", "B"),
        source="reference_argmax",
        level="Level1",
        lineage=[{"operation": "project"}],
    )

    ctx.assignment_guidance_callback(
        "start",
        problem_key="spot-1",
        route=ctx.route_key,
        operator="local_ot",
        phase="local_refinement",
        mode="prefer",
        applicability="applicable",
        numerics={"strength": 0.2},
        solver="pot",
        left_assignment=left_assignment,
        right_assignment=None,
    )
    ctx.assignment_guidance_callback(
        "attempt",
        problem_key="spot-1",
        availability="available",
    )
    attempted = json.loads((tmp_path / "provenance.json").read_text())
    attempted_event = attempted["assignment_guidance"]["events"][0]
    assert attempted_event["attempted"] is True
    assert attempted_event["availability"] == "available"
    assert attempted_event["outcome"] == "not_started"
    ctx.assignment_guidance_callback(
        "terminal",
        problem_key="spot-1",
        outcome="applied",
        right_assignment=right_assignment,
    )

    manifest = json.loads((tmp_path / "provenance.json").read_text())
    assert manifest["run"]["status"] == "running"
    evidence = manifest["assignment_guidance"]
    assert evidence["schema_version"] == 1
    assert evidence["configured"] == {
        "guidance": None,
        "compatibility_mode": None,
        "source": "route_default",
        "deprecations": [],
    }
    assert evidence["resolved"] == {
        "guidance": "prefer",
        "compatibility_mode": "cost",
        "beta": 1.0,
        "min_affinity": 0.05,
        "operator_strength": 0.2,
    }
    assert evidence["events"][0]["outcome"] == "applied"
    assert evidence["events"][0]["left_assignment"]["value_semantics"] == "soft"
    assert evidence["events"][0]["left_assignment"]["lineage"] == [
        {"operation": "aggregate"}
    ]
    assert evidence["events"][0]["right_assignment"]["value_semantics"] == (
        "one_hot"
    )
    assert evidence["events"][0]["right_assignment"]["lineage"] == [
        {"operation": "project"}
    ]
    assert evidence["summary"] == "applied"
    assert manifest["sr_allocation"] == []

    ctx.skip_pending_stages("test_complete")
    ctx.mark_run_succeeded()
    completed = json.loads((tmp_path / "provenance.json").read_text())
    assert completed["run"]["status"] == "succeeded"
    assert completed["assignment_guidance"]["events"][0]["outcome"] == "applied"


@pytest.mark.parametrize(
    ("route", "operator"),
    [
        ("sp_svc:bin2cell", "neighbor_ot"),
        ("sim2real:segmentation", "replacement_ot"),
        ("sim2real:bin2cell", "replacement_ot"),
        ("sc_svc:segmentation", "graph_edge"),
        ("sc_svc_sr:spot_size", "virtual_cell_ot"),
        ("sim2real:batch_effect", "virtual_cell_ot"),
        ("sim2real:spot_size", "virtual_cell_ot"),
        ("sim2real:gene_panel", "imputation_ot"),
        ("sim2real:gene_dropout", "imputation_ot"),
    ],
)
def test_nine_public_routes_share_one_bilateral_event_contract(
    tmp_path,
    route,
    operator,
):
    ctx = _guidance_context(tmp_path)
    ctx.route_key = route
    state = AssignmentState(
        values=np.array([[0.8, 0.2], [0.2, 0.8]]),
        observation_labels=("left-1", "left-2"),
        category_labels=("A", "B"),
        source="route_assignment",
        level="Level2" if operator == "graph_edge" else "Level1",
        value_semantics="soft",
        lineage=[{"operation": "route_projection", "route": route}],
    )
    ctx.assignment_guidance_callback(
        "start",
        problem_key=f"{route}:candidate",
        route=route,
        operator=operator,
        phase="local_refinement",
        mode="prefer",
        applicability="applicable",
        numerics={"beta": 1.0, "operator_strength": 0.2},
        solver="pot",
    )
    ctx.assignment_guidance_callback(
        "attempt",
        problem_key=f"{route}:candidate",
        availability="available",
        left_assignment=state,
        right_assignment=state,
    )
    ctx.assignment_guidance_callback(
        "terminal",
        problem_key=f"{route}:candidate",
        outcome="applied",
    )

    [event] = ctx.assignment_guidance.manifest()["events"]
    assert event["route"] == route
    assert event["operator"] == operator
    assert event["mode"] == "prefer"
    assert event["solver"] == "pot"
    assert event["outcome"] == "applied"
    for side in ("left_assignment", "right_assignment"):
        assert event[side]["source"] == "route_assignment"
        assert event[side]["value_semantics"] == "soft"
        assert event[side]["lineage"] == [
            {"operation": "route_projection", "route": route}
        ]


def test_assignment_guidance_transition_rolls_back_when_manifest_write_fails(tmp_path):
    ctx = _guidance_context(tmp_path)
    ctx.assignment_guidance_callback(
        "start",
        problem_key="spot-1",
        route=ctx.route_key,
        operator="local_ot",
        phase="local_refinement",
        mode="prefer",
        applicability="applicable",
        numerics={},
        solver="pot",
        left_assignment=None,
        right_assignment=None,
    )

    def fail_write(_ctx):
        raise OSError("manifest write failed")

    ctx.set_provenance_callback(fail_write, notify=False)
    with pytest.raises(OSError, match="manifest write failed"):
        ctx.assignment_guidance_callback(
            "attempt",
            problem_key="spot-1",
            availability="available",
        )

    assert ctx.assignment_guidance.events[0]["attempted"] is False
    assert ctx.assignment_guidance.events[0]["availability"] == "not_checked"


def test_run_cannot_succeed_with_attempted_guidance_event_not_terminal(tmp_path):
    ctx = _guidance_context(tmp_path)
    ctx.assignment_guidance_callback(
        "start",
        problem_key="unfinished",
        route=ctx.route_key,
        operator="local_ot",
        phase="local_refinement",
        mode="prefer",
        applicability="applicable",
        numerics={},
        solver="pot",
        left_assignment=None,
        right_assignment=None,
    )
    ctx.assignment_guidance_callback(
        "attempt",
        problem_key="unfinished",
        availability="available",
    )
    ctx.skip_pending_stages("test")

    with pytest.raises(RuntimeError, match="unfinished assignment guidance"):
        ctx.mark_run_succeeded()

    assert ctx.run_status == "running"
    assert ctx.assignment_guidance.events[0]["outcome"] == "not_started"


def test_atomic_json_failure_preserves_previous_manifest(monkeypatch, tmp_path):
    path = tmp_path / "provenance.json"
    provenance.write_json(path, {"generation": 1})

    def broken_dump(payload, handle, *args, **kwargs):
        handle.write('{"generation":')
        handle.flush()
        raise OSError("simulated write failure")

    monkeypatch.setattr(provenance.json, "dump", broken_dump)
    with pytest.raises(OSError, match="simulated write failure"):
        provenance.write_json(path, {"generation": 2})

    assert json.loads(path.read_text(encoding="utf-8")) == {"generation": 1}
    assert list(tmp_path.iterdir()) == [path]


def test_completed_artifact_uses_streamed_content_identity(tmp_path):
    path = tmp_path / "result.bin"
    path.write_bytes(b"first-content")
    expected = hashlib.sha256(path.read_bytes()).hexdigest()

    record = provenance.completed_artifact("output:primary", path)

    assert record == {
        "role": "output:primary",
        "path": str(path),
        "status": "completed",
        "size": len(b"first-content"),
        "sha256": expected,
    }

    path.write_bytes(b"other-content")
    assert path.stat().st_size == len(b"first-content")
    assert provenance.sha256_file(path) != expected


@pytest.mark.parametrize(
    ("recorded", "process_alive", "effective"),
    [
        ("succeeded", False, "succeeded"),
        ("failed", False, "failed"),
        ("interrupted", False, "interrupted"),
        ("running", True, "running"),
        ("running", False, "incomplete"),
    ],
)
def test_effective_run_status_does_not_rewrite_recorded_history(
    recorded,
    process_alive,
    effective,
):
    manifest = {"run": {"status": recorded}}

    assert provenance.effective_run_status(
        manifest,
        process_alive=process_alive,
    ) == effective
    assert manifest == {"run": {"status": recorded}}


def test_application_run_directories_do_not_collide_when_created_together(tmp_path):
    first = build_run_dir(str(tmp_path), "sample", "sp_svc:bin2cell")
    second = build_run_dir(str(tmp_path), "sample", "sp_svc:bin2cell")

    assert first != second
    assert first.parent == second.parent


def test_fixed_benchmark_directory_refuses_running_envelope(tmp_path):
    output_root = tmp_path / "benchmark-output"
    run_dir = output_root / "sample" / "seg_1"
    manifest_path = run_dir / "provenance.json"
    provenance.write_json(manifest_path, {"run": {"status": "running"}})

    with pytest.raises(RuntimeError, match="running.*provenance"):
        REVISEPipeline().run(
            profile="benchmark_seg",
            io_overrides={
                "data_root": str(tmp_path),
                "output_root": str(output_root),
                "sample_name": "sample",
                "seg_method": "seg_1",
            },
            dry_run=True,
        )

    assert json.loads(manifest_path.read_text(encoding="utf-8")) == {
        "run": {"status": "running"}
    }


def test_run_directory_lock_rejects_concurrent_owner(tmp_path):
    run_dir = Path(tmp_path) / "shared-run"

    with provenance.exclusive_run_directory(run_dir):
        with pytest.raises(RuntimeError, match="already active"):
            with provenance.exclusive_run_directory(run_dir):
                pass


def test_run_directory_lock_releases_after_success_and_exception(tmp_path):
    run_dir = Path(tmp_path) / "reusable-run"
    lock_path = run_dir / ".revise-run.lock"

    with provenance.exclusive_run_directory(run_dir):
        assert lock_path.is_dir()
    assert not lock_path.exists()

    with pytest.raises(RuntimeError, match="inside run"):
        with provenance.exclusive_run_directory(run_dir):
            raise RuntimeError("inside run")
    assert not lock_path.exists()

    with provenance.exclusive_run_directory(run_dir):
        assert lock_path.is_dir()


def test_fixed_benchmark_directory_can_be_reused_after_terminal_run(tmp_path):
    output_root = tmp_path / "reusable-benchmark"
    _write_benchmark_seg_inputs(tmp_path, "sample")
    options = {
        "profile": "benchmark_seg",
        "io_overrides": {
            "data_root": str(tmp_path),
            "output_root": str(output_root),
            "sample_name": "sample",
            "seg_method": "seg_1",
        },
        "dry_run": True,
    }

    REVISEPipeline().run(**options)
    first = json.loads(
        (output_root / "sample" / "seg_1" / "provenance.json").read_text()
    )
    REVISEPipeline().run(**options)
    second = json.loads(
        (output_root / "sample" / "seg_1" / "provenance.json").read_text()
    )

    assert first["run"]["status"] == "succeeded"
    assert second["run"]["status"] == "succeeded"
    assert not (output_root / "sample" / "seg_1" / ".revise-run.lock").exists()


def test_application_manifest_and_log_share_each_unique_run_directory(tmp_path):
    output_root = tmp_path / "application-output"
    _write_application_inputs(tmp_path, "sample")

    for _ in range(2):
        REVISEPipeline().run(
            profile="application_sp",
            io_overrides={
                "data_root": str(tmp_path),
                "output_root": str(output_root),
                "sample_name": "sample",
            },
            dry_run=True,
        )

    route_dir = output_root / "sample" / "sp_svc__bin2cell"
    manifests = sorted(route_dir.glob("*/provenance.json"))
    assert len(manifests) == 2
    assert len({path.parent for path in manifests}) == 2
    assert all((path.parent / "run.log").is_file() for path in manifests)
    assert {path for path in route_dir.iterdir()} == {
        manifest.parent for manifest in manifests
    }
