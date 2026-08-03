from pathlib import Path
from types import SimpleNamespace

from revise.backend.contracts import EvaluationPolicy, InputValidationPolicy, LocalRefinementStrategy
from revise.recon.context import PipelineContext
from revise.recon.pipeline import UnifiedReconstructionPipeline


class _Validation(InputValidationPolicy):
    def validate(self, ctx):
        return None


class _Evaluation(EvaluationPolicy):
    def should_evaluate(self, ctx):
        return False


class _SelectionGateStrategy(LocalRefinementStrategy):
    def prepare_context(self, ctx):
        return None

    def global_anchoring(self, ctx):
        ctx.artifacts["selection_review_required"] = True
        ctx.artifacts["selection_assessment_path"] = str(ctx.run_dir / "selection_assessment.json")

    def local_refinement(self, ctx):
        raise AssertionError("local refinement must not run before selection review")

    def solve_ot(self, ctx):
        raise AssertionError("local refinement must not run before selection review")

    def finalize_svc(self, ctx):
        raise AssertionError("finalization must not run before selection review")


def test_selection_review_gate_succeeds_without_publishing_or_refining(tmp_path):
    ctx = PipelineContext(
        merged_config={},
        raw_config={},
        config_path="revise/revise.yaml",
        profile="application_sc",
        runtime={"svc_kind": "sc", "platform": "sc_svc", "confounding": "segmentation"},
        route_key="sc_svc:segmentation",
        run_dir=Path(tmp_path),
        logger=SimpleNamespace(info=lambda *args, **kwargs: None),
    )

    svc = UnifiedReconstructionPipeline(
        strategy=_SelectionGateStrategy(),
        validation_policy=_Validation(),
        evaluation_policy=_Evaluation(),
    ).run(ctx)

    assert ctx.run_status == "succeeded"
    assert ctx.local_refinement_record["applied"] is False
    assert [record["status"] for record in ctx.stage_records] == [
        "succeeded",
        "succeeded",
        "skipped",
        "skipped",
        "skipped",
    ]
    assert svc.provenance["selection_review_required"] is True
    assert svc.provenance["selection_assessment"] == str(
        tmp_path / "selection_assessment.json"
    )
