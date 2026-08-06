from __future__ import annotations

import argparse
import sys

import pytest


def _args(**overrides):
    values = {
        "confounding": "segmentation",
        "local_refinement_strength": None,
        "sr_refinement_preset": None,
    }
    values.update(overrides)
    return argparse.Namespace(**values)


def test_benchmark_cli_does_not_expose_generic_set_overrides(monkeypatch):
    from revise.benchmark import cli

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_main.py",
            "--confounding",
            "segmentation",
            "--data-root",
            "data",
            "--sample-name",
            "sample",
        ],
    )

    assert not hasattr(cli.get_args(), "set_overrides")


def test_benchmark_cli_rejects_removed_set_option(monkeypatch):
    from revise.benchmark import cli

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_main.py",
            "--confounding",
            "segmentation",
            "--data-root",
            "data",
            "--sample-name",
            "sample",
            "--set",
            "graph.method=pca",
        ],
    )

    with pytest.raises(SystemExit) as exc_info:
        cli.get_args()

    assert exc_info.value.code == 2


def test_sr_refinement_preset_builds_nested_sc_configuration():
    from revise.benchmark.cli import _build_algorithm_overrides

    overrides = _build_algorithm_overrides(
        _args(
            confounding="spot_size",
            sr_refinement_preset="confidence_anchor",
        )
    )

    assert overrides["sc"]["sr_graph_agg_enabled"] is True
    assert overrides["sc"]["sr_graph_agg_low_conf_only"] is True
    assert overrides["sc"]["sr_graph_agg_anchor_only"] is True


def test_run_case_uses_only_the_benchmark_selector():
    from revise.benchmark.cli import _run_case

    captured = {}

    class Pipeline:
        def run(self, **kwargs):
            captured.update(kwargs)
            return type(
                "SVC",
                (),
                {
                    "provenance": {
                        "run_dir": "run",
                        "profile": "benchmark_seg",
                        "local_refinement": {"applied": True},
                    },
                    "summary": staticmethod(lambda: {"ok": True}),
                },
            )()

    result = _run_case(
        Pipeline(),
        confounding="segmentation",
        io_overrides={"sample_name": "sample"},
        runtime_seed=None,
        algorithm_overrides={"local_refinement": {"strength": 0.0}},
    )

    assert result["ok"] is True
    assert result["profile"] == "benchmark_seg"
    assert captured["svc_type"] is None
    assert captured["cf"] == "segmentation"
    assert captured["runtime_overrides"] == {"seed": None}


def test_seed_scope_preserves_fixed_run_seed_and_distinct_process_seeds():
    from revise.benchmark.cli import _runtime_seed_supplier

    args = argparse.Namespace(seed=42)
    fixed = _runtime_seed_supplier(args, "run")
    streamed = _runtime_seed_supplier(args, "process")

    assert fixed() == fixed() == 42
    assert streamed() != streamed()
