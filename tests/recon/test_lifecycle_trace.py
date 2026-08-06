from __future__ import annotations

import logging
import json
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise.recon.context import PipelineContext
from revise.recon.pipeline import UnifiedReconstructionPipeline
from revise.svc import SVC


STAGES = [
    ("validate_inputs", "pipeline.validate_inputs"),
    ("global_anchoring", "strategy.global_anchoring"),
    ("local_refinement", "strategy.local_refinement"),
    ("finalize", "pipeline.finalize"),
    ("evaluate", "pipeline.evaluate"),
]


class AcceptInputs:
    def validate(self, ctx) -> None:
        return None


class EvaluationPolicy:
    def __init__(self, enabled: bool) -> None:
        self.enabled = enabled

    def should_evaluate(self, ctx) -> bool:
        return self.enabled


class RecordingStrategy:
    strategy_id = "RecordingStrategy"

    def __init__(self, *, outputs=None, ground_truth=None) -> None:
        self.calls: list[str] = []
        self.outputs = outputs
        self.ground_truth = ground_truth

    def prepare_context(self, ctx) -> None:
        self.calls.append("prepare_context")
        ctx.real_st_adata = self.ground_truth

    def global_anchoring(self, ctx) -> None:
        self.calls.append("global_anchoring")

    def prepare_local_units(self, ctx) -> None:
        self.calls.append("prepare_local_units")

    def build_graph(self, ctx) -> None:
        self.calls.append("build_graph")

    def build_ot_problem(self, ctx) -> None:
        self.calls.append("build_ot_problem")

    def solve_ot(self, ctx) -> None:
        self.calls.append("solve_ot")

    def update_expression(self, ctx) -> None:
        self.calls.append("update_expression")

    def finalize_svc(self, ctx) -> SVC:
        self.calls.append("finalize_svc")
        outputs = {} if self.outputs is None else self.outputs
        return SVC(
            expr=None,
            spatial=None,
            svc_kind="sc",
            artifacts={"outputs": outputs},
        )


def _adata(obs_names: list[str]) -> AnnData:
    return AnnData(
        X=np.arange(1, len(obs_names) + 1, dtype=float).reshape(-1, 1),
        obs=pd.DataFrame(index=pd.Index(obs_names)),
        var=pd.DataFrame(index=pd.Index(["g1"])),
    )


def _write_application_inputs(data_root: Path, sample_name: str) -> None:
    st = _adata(["spot-1", "spot-2"])
    st.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 1.0]])
    sc_ref = _adata(["cell-1", "cell-2"])
    sc_ref.obs["Level1"] = ["A", "B"]
    st.write_h5ad(data_root / f"{sample_name}_Xenium.h5ad")
    sc_ref.write_h5ad(data_root / "adata_sc_all_reanno.h5ad")


def _context(tmp_path: Path, *, dry_run: bool = False) -> PipelineContext:
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
        runtime={
            "mode": "application",
            "svc_kind": "sc",
            "strategy": "RecordingStrategy",
        },
        route_key="test:test",
        run_dir=tmp_path,
        logger=logging.getLogger("test-lifecycle-trace"),
        dry_run=dry_run,
    )


def _pipeline(strategy: RecordingStrategy, *, evaluate: bool):
    return UnifiedReconstructionPipeline(
        strategy=strategy,
        validation_policy=AcceptInputs(),
        evaluation_policy=EvaluationPolicy(evaluate),
    )


def _stage(ctx: PipelineContext, name: str) -> dict:
    return next(stage for stage in ctx.stage_records if stage["name"] == name)


def _install_metrics(monkeypatch) -> None:
    metrics = types.ModuleType("revise.analysis.metrics")
    metrics.compute_metric = lambda *args, **kwargs: pd.DataFrame(
        {"g1": [0.0]}, index=["MSE"]
    )
    monkeypatch.setitem(sys.modules, "revise.analysis.metrics", metrics)


def test_success_records_five_owned_terminal_stages(monkeypatch, tmp_path):
    _install_metrics(monkeypatch)
    data = _adata(["c1", "c2"])
    strategy = RecordingStrategy(outputs={"result": data}, ground_truth=data.copy())
    ctx = _context(tmp_path)

    _pipeline(strategy, evaluate=True).run(ctx)

    assert ctx.run_status == "succeeded"
    assert [(stage["name"], stage["owner"]) for stage in ctx.stage_records] == STAGES
    assert [stage["status"] for stage in ctx.stage_records] == ["succeeded"] * 5
    for stage in ctx.stage_records:
        assert set(stage) == {
            "name",
            "owner",
            "status",
            "started_at",
            "duration_seconds",
            "reason",
            "error",
        }
        assert stage["started_at"] is not None
        assert stage["duration_seconds"] is not None
        assert stage["duration_seconds"] >= 0
        assert stage["reason"] is None
        assert stage["error"] is None


def test_local_refinement_umbrella_runs_all_five_internal_operations_once(tmp_path):
    strategy = RecordingStrategy()
    ctx = _context(tmp_path)

    _pipeline(strategy, evaluate=False).run(ctx)

    assert strategy.calls == [
        "prepare_context",
        "global_anchoring",
        "prepare_local_units",
        "build_graph",
        "build_ot_problem",
        "solve_ot",
        "update_expression",
        "finalize_svc",
    ]
    assert [
        stage["name"]
        for stage in ctx.stage_records
        if stage["name"] == "local_refinement"
    ] == ["local_refinement"]
    assert _stage(ctx, "local_refinement")["status"] == "succeeded"


def test_dry_run_succeeds_validation_and_skips_remaining_stages(tmp_path):
    strategy = RecordingStrategy()
    ctx = _context(tmp_path, dry_run=True)

    _pipeline(strategy, evaluate=True).run(ctx)

    assert ctx.run_status == "succeeded"
    assert _stage(ctx, "validate_inputs")["status"] == "succeeded"
    assert strategy.calls == ["prepare_context"]
    for name, _ in STAGES[1:]:
        stage = _stage(ctx, name)
        assert stage["status"] == "skipped"
        assert stage["reason"] == "dry_run"
        assert stage["started_at"] is None
        assert stage["duration_seconds"] is None


def test_framework_dry_run_persists_the_same_five_stage_truth(
    monkeypatch,
    tmp_path,
):
    import revise.framework as framework

    monkeypatch.setattr(
        framework,
        "build_default_registry",
        lambda: (_ for _ in ()).throw(
            AssertionError("dry-run imported the scientific strategy registry")
        ),
    )
    output_root = tmp_path / "dry-run-output"
    _write_application_inputs(tmp_path, "dry-run-case")

    svc = framework.REVISEPipeline().run(
        svc_type="sp-SVC",
        io_overrides={
            "data_root": str(tmp_path),
            "output_root": str(output_root),
            "sample_name": "dry-run-case",
        },
        dry_run=True,
    )

    paths = list((output_root / "dry-run-case").rglob("provenance.json"))
    assert len(paths) == 1, paths
    manifest = json.loads(paths[0].read_text(encoding="utf-8"))
    assert manifest["run"]["status"] == "succeeded"
    assert manifest["run"]["dry_run"] is True
    assert manifest["stages"][0]["status"] == "succeeded"
    assert all(
        stage["status"] == "skipped" and stage["reason"] == "dry_run"
        for stage in manifest["stages"][1:]
    )
    assert "stage_trace" not in manifest
    assert svc.provenance["stages"] == manifest["stages"]


def test_disabled_evaluation_is_an_explicit_policy_skip(tmp_path):
    strategy = RecordingStrategy()
    ctx = _context(tmp_path)

    _pipeline(strategy, evaluate=False).run(ctx)

    stage = _stage(ctx, "evaluate")
    assert ctx.run_status == "succeeded"
    assert stage["status"] == "skipped"
    assert stage["reason"] == "policy_disabled"
    assert stage["started_at"] is None
    assert stage["duration_seconds"] is None


@pytest.mark.parametrize(
    ("outputs", "ground_truth", "reason"),
    [
        ({}, _adata(["c1"]), "no_outputs"),
        ({"result": _adata(["c1"])}, None, "ground_truth_unavailable"),
        ({"result": _adata(["pred"])}, _adata(["truth"]), "no_aligned_outputs"),
    ],
)
def test_unavailable_evaluation_inputs_are_explicit_skips(
    monkeypatch, tmp_path, outputs, ground_truth, reason
):
    _install_metrics(monkeypatch)
    strategy = RecordingStrategy(outputs=outputs, ground_truth=ground_truth)
    ctx = _context(tmp_path)

    _pipeline(strategy, evaluate=True).run(ctx)

    stage = _stage(ctx, "evaluate")
    assert ctx.run_status == "succeeded"
    assert stage["status"] == "skipped"
    assert stage["reason"] == reason
    assert stage["started_at"] is None
    assert stage["duration_seconds"] is None


def test_evaluation_decision_failure_is_owned_by_evaluate_stage(tmp_path):
    class BrokenEvaluationPolicy:
        def should_evaluate(self, ctx) -> bool:
            raise ValueError("evaluation policy failed")

    strategy = RecordingStrategy()
    ctx = _context(tmp_path)
    pipeline = UnifiedReconstructionPipeline(
        strategy=strategy,
        validation_policy=AcceptInputs(),
        evaluation_policy=BrokenEvaluationPolicy(),
    )

    with pytest.raises(ValueError, match="evaluation policy failed"):
        pipeline.run(ctx)

    stage = _stage(ctx, "evaluate")
    assert ctx.run_status == "failed"
    assert stage["status"] == "failed"
    assert stage["error"]["type"] == "ValueError"
