from __future__ import annotations

import logging

from revise.recon.context import PipelineContext


def _context(tmp_path, *, task: str, local_refinement=None):
    config = {
        "ot": {
            "ga": {"solver": "pot"},
            "lr": {"solver": "pot"},
        }
    }
    if local_refinement is not None:
        config["local_refinement"] = local_refinement
    return PipelineContext(
        merged_config=config,
        raw_config={},
        config_path="revise/revise.yaml",
        profile="test",
        runtime={"task": task},
        route_key=f"test:{task}",
        run_dir=tmp_path,
        logger=logging.getLogger("test-local-refinement-record"),
    )


def test_conditioned_route_starts_with_one_minimal_control_record(tmp_path):
    ctx = _context(
        tmp_path,
        task="sp_svc",
        local_refinement={"strength": 0.2},
    )

    assert ctx.local_refinement_record == {
        "route": "test:sp_svc",
        "applied": False,
        "strength": 0.2,
    }

    ctx.record_local_refinement(True)
    assert ctx.local_refinement_record["applied"] is True

    ctx.record_local_refinement(False)
    assert ctx.local_refinement_record["applied"] is True


def test_sc_route_has_no_inactive_strength_identity_value(tmp_path):
    ctx = _context(tmp_path, task="sc_svc")

    assert ctx.local_refinement_record == {
        "route": "test:sc_svc",
        "applied": False,
        "strength": None,
    }
