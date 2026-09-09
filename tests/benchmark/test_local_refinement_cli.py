from __future__ import annotations

import argparse

import pytest

from revise.benchmark import cli


def test_cli_exposes_only_local_refinement_strength():
    help_text = cli.get_parser().format_help()

    assert "--local-refinement-strength" in help_text
    for removed in cli.REMOVED_ASSIGNMENT_GUIDANCE_FLAGS:
        assert removed not in help_text


def test_strength_builds_the_only_local_refinement_override():
    args = argparse.Namespace(
        local_refinement_strength=0.35,
        sr_refinement_preset=None,
        route="segmentation",
        evaluate=True,
    )

    assert cli._build_algorithm_overrides(args) == {
        "local_refinement": {"strength": 0.35},
        "benchmark": {"evaluate": True},
    }


@pytest.mark.parametrize("removed", sorted(cli.REMOVED_ASSIGNMENT_GUIDANCE_FLAGS))
def test_removed_cli_flags_fail_with_one_migration_message(removed, monkeypatch):
    parser = cli.get_parser()
    monkeypatch.setattr(
        parser,
        "error",
        lambda message: (_ for _ in ()).throw(ValueError(message)),
    )
    argv = [
        "--config",
        "segmentation.yaml",
        "--data-root",
        "data",
        "--sample-name",
        "sample",
        removed,
    ]
    if removed != "--posterior-strict":
        argv.append("value")

    with pytest.raises(
        ValueError,
        match=(
            "Assignment guidance options were removed; "
            "use --local-refinement-strength"
        ),
    ):
        parser.parse_args(argv)


def test_leaf_aggregation_does_not_rebuild_event_summaries():
    results = [
        {"local_refinement": {"route": "sim2real:segmentation", "applied": True, "strength": 0.2}},
        {"local_refinement": {"route": "sim2real:segmentation", "applied": False, "strength": 0.2}},
    ]

    assert cli._aggregate_local_refinement(results) == [
        {"route": "sim2real:segmentation", "applied": True, "strength": 0.2},
        {"route": "sim2real:segmentation", "applied": False, "strength": 0.2},
    ]
