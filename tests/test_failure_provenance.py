from __future__ import annotations

import copy
import json
import logging
import os
import signal
import subprocess
import sys
import textwrap
import types
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise.recon.context import PipelineContext
from revise.recon.pipeline import UnifiedReconstructionPipeline
from revise.svc import SVC
from revise.utils.provenance import completed_artifact, effective_run_status, write_json


STAGE_NAMES = [
    "validate_inputs",
    "global_anchoring",
    "local_refinement",
    "finalize",
    "evaluate",
]


class InjectedTermination(Exception):
    pass


def _raise(stage: str, selected: str, interrupt: bool) -> None:
    if stage != selected:
        return
    if interrupt:
        raise KeyboardInterrupt(f"interrupted at {stage}")
    raise InjectedTermination(f"failed at {stage}")


class FailingValidation:
    def __init__(self, selected: str, interrupt: bool) -> None:
        self.selected = selected
        self.interrupt = interrupt

    def validate(self, ctx) -> None:
        _raise("validate_inputs", self.selected, self.interrupt)


class FailingStrategy:
    strategy_id = "FailingStrategy"

    def __init__(self, selected: str, interrupt: bool, output: AnnData) -> None:
        self.selected = selected
        self.interrupt = interrupt
        self.output = output

    def prepare_context(self, ctx) -> None:
        _raise("validate_inputs", self.selected, self.interrupt)
        ctx.real_st_adata = self.output.copy()

    def global_anchoring(self, ctx) -> None:
        if self.selected == "global_anchoring":
            ctx.record_ot_event("ga", "pot", "attempted")
        _raise("global_anchoring", self.selected, self.interrupt)

    def prepare_local_units(self, ctx) -> None:
        return None

    def build_graph(self, ctx) -> None:
        return None

    def build_ot_problem(self, ctx) -> None:
        return None

    def solve_ot(self, ctx) -> None:
        if self.selected == "local_refinement":
            ctx.record_ot_event("lr", "pot", "attempted")
        _raise("local_refinement", self.selected, self.interrupt)

    def update_expression(self, ctx) -> None:
        return None

    def finalize_svc(self, ctx) -> SVC:
        _raise("finalize", self.selected, self.interrupt)
        return SVC(
            expr=None,
            spatial=None,
            svc_kind="sc",
            artifacts={"outputs": {"result": self.output.copy()}},
        )


class EvaluateAlways:
    def should_evaluate(self, ctx) -> bool:
        return True


def _context(tmp_path) -> PipelineContext:
    return PipelineContext(
        merged_config={
            "io": {"save_outputs": False},
            "ot": {
                "ga": {"solver": "pot"},
                "lr": {"solver": "pot"},
            },
        },
        raw_config={},
        config_path="revise/revise.yaml",
        profile="test",
        runtime={"mode": "benchmark", "svc_kind": "sc"},
        route_key="test:test",
        run_dir=tmp_path,
        logger=logging.getLogger("test-failure-provenance"),
    )


def _install_metrics(monkeypatch, selected: str, interrupt: bool) -> None:
    metrics = types.ModuleType("revise.analysis.metrics")

    def compute_metric(*args, **kwargs):
        _raise("evaluate", selected, interrupt)
        return pd.DataFrame({"MSE": [0.0]}, index=["g1"])

    metrics.compute_metric = compute_metric
    monkeypatch.setitem(sys.modules, "revise.analysis.metrics", metrics)


def _pipeline(selected: str, interrupt: bool, output: AnnData):
    return UnifiedReconstructionPipeline(
        strategy=FailingStrategy(selected, interrupt, output),
        validation_policy=FailingValidation(selected, interrupt),
        evaluation_policy=EvaluateAlways(),
    )


def _attach_manifest_writer(ctx: PipelineContext) -> None:
    def persist(current: PipelineContext) -> None:
        write_json(
            current.run_dir / "provenance.json",
            {
                "run": {
                    "status": current.run_status,
                    "error": copy.deepcopy(current.run_error),
                },
                "stages": copy.deepcopy(current.stage_records),
                "ot_events": copy.deepcopy(current.ot_events),
                "artifacts": copy.deepcopy(current.artifact_records),
            },
        )

    ctx.set_provenance_callback(persist)


@pytest.mark.parametrize("selected", STAGE_NAMES)
def test_each_stage_failure_persists_terminal_truth(monkeypatch, tmp_path, selected):
    output = AnnData(
        X=np.array([[1.0]]),
        obs=pd.DataFrame(index=["cell-1"]),
        var=pd.DataFrame(index=["g1"]),
    )
    _install_metrics(monkeypatch, selected, interrupt=False)
    ctx = _context(tmp_path)
    _attach_manifest_writer(ctx)

    with pytest.raises(InjectedTermination, match=f"failed at {selected}"):
        _pipeline(selected, False, output).run(ctx)

    manifest = json.loads((tmp_path / "provenance.json").read_text())
    selected_index = STAGE_NAMES.index(selected)
    statuses = [stage["status"] for stage in manifest["stages"]]
    assert statuses[:selected_index] == ["succeeded"] * selected_index
    assert statuses[selected_index] == "failed"
    assert statuses[selected_index + 1 :] == ["skipped"] * (
        len(STAGE_NAMES) - selected_index - 1
    )
    assert all(
        stage["reason"] == "upstream_failure"
        for stage in manifest["stages"][selected_index + 1 :]
    )
    assert manifest["run"]["status"] == "failed"
    assert manifest["stages"][selected_index]["error"]["type"] == (
        "InjectedTermination"
    )
    assert not any(event["status"] == "completed" for event in manifest["ot_events"])


class _Registry:
    def __init__(self, strategy) -> None:
        self.strategy = strategy

    def get(self, strategy_id):
        return self.strategy


class _FrameworkStrategy:
    strategy_id = "FrameworkStrategy"

    def __init__(self, action: str) -> None:
        self.action = action

    def prepare_context(self, ctx) -> None:
        if self.action == "fail":
            raise RuntimeError("framework validation failed")

    def global_anchoring(self, ctx) -> None:
        if self.action == "sigterm":
            os.kill(os.getpid(), signal.SIGTERM)

    def prepare_local_units(self, ctx) -> None:
        return None

    def build_graph(self, ctx) -> None:
        return None

    def build_ot_problem(self, ctx) -> None:
        return None

    def solve_ot(self, ctx) -> None:
        return None

    def update_expression(self, ctx) -> None:
        return None

    def finalize_svc(self, ctx) -> SVC:
        return SVC(expr=None, spatial=None, svc_kind="sp", artifacts={"outputs": {}})


def _framework_manifest(output_root: Path, sample_name: str) -> dict:
    paths = list((output_root / sample_name).rglob("provenance.json"))
    assert len(paths) == 1, paths
    return json.loads(paths[0].read_text(encoding="utf-8"))


def _write_framework_inputs(data_root: Path, sample_name: str) -> None:
    st = AnnData(
        X=np.ones((2, 1)),
        obs=pd.DataFrame(index=["spot-1", "spot-2"]),
        var=pd.DataFrame(index=["g1"]),
    )
    st.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 1.0]])
    sc_ref = AnnData(
        X=np.ones((2, 1)),
        obs=pd.DataFrame({"Level1": ["A", "B"]}, index=["cell-1", "cell-2"]),
        var=pd.DataFrame(index=["g1"]),
    )
    st.write_h5ad(data_root / f"{sample_name}_Xenium.h5ad")
    sc_ref.write_h5ad(data_root / "adata_sc_all_reanno.h5ad")


def test_finalize_callback_failure_marks_finalize_and_run_failed(tmp_path):
    from revise.framework import REVISEPipeline

    output_root = tmp_path / "callback-failure"
    _write_framework_inputs(tmp_path, "callback-case")
    pipeline = REVISEPipeline()
    pipeline.registry = _Registry(_FrameworkStrategy("ok"))

    def fail_publication(ctx):
        raise OSError("public result write failed")

    with pytest.raises(OSError, match="public result write failed"):
        pipeline.run(
            profile="application_sp",
            io_overrides={
                "data_root": str(tmp_path),
                "output_root": str(output_root),
                "sample_name": "callback-case",
            },
            finalize_callback=fail_publication,
        )

    manifest = _framework_manifest(output_root, "callback-case")
    assert manifest["run"]["status"] == "failed"
    assert manifest["stages"][3]["status"] == "failed"
    assert not any(
        artifact["role"] == "public_result"
        and artifact["status"] == "completed"
        for artifact in manifest["artifacts"]
    )


def test_finalize_callback_can_register_public_result_before_success(tmp_path):
    from revise.framework import REVISEPipeline

    output_root = tmp_path / "callback-success"
    _write_framework_inputs(tmp_path, "callback-case")
    pipeline = REVISEPipeline()
    pipeline.registry = _Registry(_FrameworkStrategy("ok"))

    def publish(ctx):
        path = output_root / "callback-case" / "public.bin"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"published")
        ctx.record_artifact(completed_artifact("public_result", path))

    pipeline.run(
        profile="application_sp",
        io_overrides={
            "data_root": str(tmp_path),
            "output_root": str(output_root),
            "sample_name": "callback-case",
        },
        finalize_callback=publish,
    )

    manifest = _framework_manifest(output_root, "callback-case")
    assert manifest["run"]["status"] == "succeeded"
    assert manifest["stages"][3]["status"] == "succeeded"
    assert any(
        artifact["role"] == "public_result"
        and artifact["status"] == "completed"
        for artifact in manifest["artifacts"]
    )


@pytest.mark.parametrize("selected", STAGE_NAMES)
def test_framework_each_stage_failure_keeps_terminal_manifest(
    monkeypatch,
    tmp_path,
    selected,
):
    import revise.framework as framework

    output = AnnData(
        X=np.array([[1.0]]),
        obs=pd.DataFrame(index=["cell-1"]),
        var=pd.DataFrame(index=["g1"]),
    )
    _install_metrics(monkeypatch, selected, interrupt=False)
    monkeypatch.setattr(framework, "ModeEvaluationPolicy", EvaluateAlways)
    output_root = tmp_path / f"failure-{selected}"
    _write_framework_inputs(tmp_path, "failure-case")
    pipeline = framework.REVISEPipeline()
    pipeline.registry = _Registry(FailingStrategy(selected, False, output))

    with pytest.raises(InjectedTermination, match=f"failed at {selected}"):
        pipeline.run(
            profile="application_sp",
            io_overrides={
                "data_root": str(tmp_path),
                "output_root": str(output_root),
                "sample_name": "failure-case",
            },
        )

    manifest = _framework_manifest(output_root, "failure-case")
    selected_index = STAGE_NAMES.index(selected)
    statuses = [stage["status"] for stage in manifest["stages"]]
    assert statuses[:selected_index] == ["succeeded"] * selected_index
    assert statuses[selected_index] == "failed"
    assert statuses[selected_index + 1 :] == ["skipped"] * (
        len(STAGE_NAMES) - selected_index - 1
    )
    assert manifest["run"]["status"] == "failed"
    assert manifest["stages"][selected_index]["error"]["type"] == (
        "InjectedTermination"
    )


def test_framework_retries_terminal_manifest_without_masking_science_error(
    tmp_path,
):
    from revise.framework import REVISEPipeline

    class FailFirstTerminalWrite(REVISEPipeline):
        def __init__(self):
            super().__init__()
            self.failed_terminal_write = False

        def _write_final_metadata(self, ctx):
            if ctx.run_status == "failed" and not self.failed_terminal_write:
                self.failed_terminal_write = True
                raise OSError("transient provenance failure")
            return super()._write_final_metadata(ctx)

    output_root = tmp_path / "retry-output"
    _write_framework_inputs(tmp_path, "retry-case")
    pipeline = FailFirstTerminalWrite()
    pipeline.registry = _Registry(_FrameworkStrategy("fail"))

    with pytest.raises(RuntimeError, match="framework validation failed") as exc_info:
        pipeline.run(
            profile="application_sp",
            io_overrides={
                "data_root": str(tmp_path),
                "output_root": str(output_root),
                "sample_name": "retry-case",
            },
        )

    assert isinstance(exc_info.value.__cause__, OSError)
    manifest = _framework_manifest(output_root, "retry-case")
    assert manifest["run"]["status"] == "failed"
    assert manifest["run"]["error"]["message"] == "framework validation failed"


def test_persistent_terminal_write_failure_preserves_running_disk_truth(
    tmp_path,
):
    from revise.framework import REVISEPipeline

    class RefuseTerminalWrites(REVISEPipeline):
        def _write_final_metadata(self, ctx):
            if ctx.run_status == "failed":
                raise OSError("persistent provenance failure")
            return super()._write_final_metadata(ctx)

    output_root = tmp_path / "persistent-output"
    _write_framework_inputs(tmp_path, "persistent-case")
    pipeline = RefuseTerminalWrites()
    pipeline.registry = _Registry(_FrameworkStrategy("fail"))

    with pytest.raises(RuntimeError, match="framework validation failed") as exc_info:
        pipeline.run(
            profile="application_sp",
            io_overrides={
                "data_root": str(tmp_path),
                "output_root": str(output_root),
                "sample_name": "persistent-case",
            },
        )

    assert isinstance(exc_info.value.__cause__, OSError)
    manifest = _framework_manifest(output_root, "persistent-case")
    assert manifest["run"]["status"] == "running"
    assert effective_run_status(manifest, process_alive=False) == "incomplete"


def test_terminal_persistence_failure_rolls_back_in_memory_before_retry(tmp_path):
    ctx = _context(tmp_path)
    manifest_path = tmp_path / "provenance.json"
    failed_once = False

    def persist(current):
        nonlocal failed_once
        if current.run_status == "failed" and not failed_once:
            failed_once = True
            raise OSError("transient provenance failure")
        write_json(
            manifest_path,
            {"run": current.run_record, "stages": current.stage_records},
        )

    ctx.set_provenance_callback(persist)
    ctx.start_stage("validate_inputs")

    with pytest.raises(OSError, match="transient provenance failure"):
        ctx.terminate_stage("validate_inputs", RuntimeError("science failed"))

    assert ctx.run_status == "running"
    assert ctx.run_ended_at is None
    assert ctx.run_error is None
    assert ctx.stage_records[0]["status"] == "running"
    assert ctx.stage_records[0]["error"] is None
    assert all(stage["status"] == "pending" for stage in ctx.stage_records[1:])
    assert json.loads(manifest_path.read_text())["run"]["status"] == "running"

    ctx.terminate_stage("validate_inputs", RuntimeError("science failed"))
    assert ctx.run_status == "failed"
    assert json.loads(manifest_path.read_text())["run"]["status"] == "failed"


def test_framework_sigterm_is_persisted_and_host_handler_is_restored(tmp_path):
    from revise.framework import REVISEPipeline

    output_root = tmp_path / "sigterm-output"
    _write_framework_inputs(tmp_path, "sigterm-case")
    previous_handler = signal.getsignal(signal.SIGTERM)
    pipeline = REVISEPipeline()
    pipeline.registry = _Registry(_FrameworkStrategy("sigterm"))

    with pytest.raises(KeyboardInterrupt, match="received SIGTERM"):
        pipeline.run(
            profile="application_sp",
            io_overrides={
                "data_root": str(tmp_path),
                "output_root": str(output_root),
                "sample_name": "sigterm-case",
            },
        )

    assert signal.getsignal(signal.SIGTERM) == previous_handler
    manifest = _framework_manifest(output_root, "sigterm-case")
    assert manifest["run"]["status"] == "interrupted"
    assert manifest["stages"][1]["status"] == "interrupted"
    assert all(
        stage["status"] == "skipped" and stage["reason"] == "run_interrupted"
        for stage in manifest["stages"][2:]
    )


def test_framework_sigterm_chains_custom_handler_after_terminal_persistence(
    tmp_path,
):
    from revise.framework import REVISEPipeline

    output_root = tmp_path / "custom-sigterm-output"
    _write_framework_inputs(tmp_path, "custom-sigterm-case")
    original_handler = signal.getsignal(signal.SIGTERM)
    observed_statuses = []

    def custom_handler(_signum, _frame):
        assert signal.getsignal(signal.SIGTERM) is custom_handler
        manifest = _framework_manifest(output_root, "custom-sigterm-case")
        observed_statuses.append(manifest["run"]["status"])

    signal.signal(signal.SIGTERM, custom_handler)
    try:
        pipeline = REVISEPipeline()
        pipeline.registry = _Registry(_FrameworkStrategy("sigterm"))

        with pytest.raises(KeyboardInterrupt, match="received SIGTERM"):
            pipeline.run(
                profile="application_sp",
                io_overrides={
                    "data_root": str(tmp_path),
                    "output_root": str(output_root),
                    "sample_name": "custom-sigterm-case",
                },
            )

        assert observed_statuses == ["interrupted"]
        assert signal.getsignal(signal.SIGTERM) is custom_handler
    finally:
        signal.signal(signal.SIGTERM, original_handler)


def test_temporary_sigterm_handler_respects_ignored_host_signal():
    from revise.framework import _temporary_sigterm_handler

    original_handler = signal.getsignal(signal.SIGTERM)
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    try:
        with _temporary_sigterm_handler():
            assert signal.getsignal(signal.SIGTERM) is signal.SIG_IGN
        assert signal.getsignal(signal.SIGTERM) is signal.SIG_IGN
    finally:
        signal.signal(signal.SIGTERM, original_handler)


def test_abrupt_death_leaves_running_manifest_for_incomplete_inspection(tmp_path):
    output_root = tmp_path / "abrupt-output"
    _write_framework_inputs(tmp_path, "abrupt-case")
    code = f"""
import os

from revise.framework import REVISEPipeline

class AbruptStrategy:
    strategy_id = "AbruptStrategy"

    def prepare_context(self, ctx):
        os._exit(17)

class Registry:
    def get(self, strategy_id):
        return AbruptStrategy()

pipeline = REVISEPipeline()
pipeline.registry = Registry()
pipeline.run(
    profile="application_sp",
    io_overrides={{
        "data_root": {str(tmp_path)!r},
        "output_root": {str(output_root)!r},
        "sample_name": "abrupt-case",
    }},
)
"""
    result = subprocess.run(
        [sys.executable, "-c", textwrap.dedent(code)],
        cwd=Path(__file__).resolve().parents[1],
        check=False,
    )

    assert result.returncode == 17
    manifests = list(output_root.rglob("provenance.json"))
    assert len(manifests) == 1
    manifest_path = manifests[0]
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["run"]["status"] == "running"
    assert manifest["stages"][0]["status"] == "running"
    assert manifest["stages"][1]["status"] == "pending"
    assert effective_run_status(manifest, process_alive=False) == "incomplete"
    assert (manifest_path.parent / ".revise-run.lock").is_dir()


class _WrittenOutput:
    def __init__(self, payload: bytes, *, fail: bool = False) -> None:
        self.payload = payload
        self.fail = fail

    def write_h5ad(self, path) -> None:
        Path(path).write_bytes(self.payload)
        if self.fail:
            raise OSError("second artifact failed")


class _ArtifactStrategy(_FrameworkStrategy):
    def __init__(self) -> None:
        super().__init__("success")

    def finalize_svc(self, ctx) -> SVC:
        return SVC(
            expr=None,
            spatial=None,
            svc_kind="sp",
            artifacts={
                "outputs": {
                    "first": _WrittenOutput(b"complete"),
                    "second": _WrittenOutput(b"partial", fail=True),
                }
            },
        )


def test_only_completed_artifacts_are_registered_before_finalize_failure(tmp_path):
    from revise.framework import REVISEPipeline

    output_root = tmp_path / "artifact-output"
    _write_framework_inputs(tmp_path, "artifact-case")
    pipeline = REVISEPipeline()
    pipeline.registry = _Registry(_ArtifactStrategy())

    with pytest.raises(OSError, match="second artifact failed"):
        pipeline.run(
            profile="application_sp",
            io_overrides={
                "data_root": str(tmp_path),
                "output_root": str(output_root),
                "sample_name": "artifact-case",
            },
        )

    manifest = _framework_manifest(output_root, "artifact-case")
    assert manifest["run"]["status"] == "failed"
    assert [artifact["role"] for artifact in manifest["artifacts"]] == [
        "preflight",
        "output:first",
    ]
    assert all(
        artifact["status"] == "completed" for artifact in manifest["artifacts"]
    )
    assert all(len(artifact["sha256"]) == 64 for artifact in manifest["artifacts"])


@pytest.mark.parametrize("selected", STAGE_NAMES)
def test_each_stage_interrupt_persists_interrupted_truth(
    monkeypatch,
    tmp_path,
    selected,
):
    output = AnnData(
        X=np.array([[1.0]]),
        obs=pd.DataFrame(index=["cell-1"]),
        var=pd.DataFrame(index=["g1"]),
    )
    _install_metrics(monkeypatch, selected, interrupt=True)
    ctx = _context(tmp_path)
    _attach_manifest_writer(ctx)

    with pytest.raises(KeyboardInterrupt, match=f"interrupted at {selected}"):
        _pipeline(selected, True, output).run(ctx)

    manifest = json.loads((tmp_path / "provenance.json").read_text())
    selected_index = STAGE_NAMES.index(selected)
    assert manifest["run"]["status"] == "interrupted"
    assert manifest["stages"][selected_index]["status"] == "interrupted"
    assert all(
        stage["status"] == "skipped" and stage["reason"] == "run_interrupted"
        for stage in manifest["stages"][selected_index + 1 :]
    )
    assert not any(event["status"] == "completed" for event in manifest["ot_events"])
