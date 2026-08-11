from __future__ import annotations

import argparse
import sys

import pytest


def _args(**overrides):
    values = {
        "route": "segmentation",
        "local_refinement_strength": None,
        "sr_refinement_preset": None,
    }
    values.update(overrides)
    return argparse.Namespace(**values)


def test_benchmark_yaml_is_the_only_route_selector():
    from revise.benchmark.cli import get_parser

    help_text = get_parser().format_help()
    assert "--config" in help_text
    assert "--confounding" not in help_text


def test_benchmark_cli_does_not_expose_generic_set_overrides(monkeypatch):
    from revise.benchmark import cli

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_main.py",
            "--config",
            "segmentation.yaml",
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
            "--config",
            "segmentation.yaml",
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
            route="spot_size",
            sr_refinement_preset="confidence_anchor",
        )
    )

    assert overrides["sc"]["sr_graph_agg_enabled"] is True
    assert overrides["sc"]["sr_graph_agg_low_conf_only"] is True
    assert overrides["sc"]["sr_graph_agg_anchor_only"] is True


def test_run_case_uses_only_the_benchmark_selector():
    from revise.benchmark.cli import _run_case
    from revise.utils.provenance import hash_jsonable

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
        route="segmentation",
        io_overrides={"sample_name": "sample"},
        runtime_seed=None,
        algorithm_overrides={"local_refinement": {"strength": 0.0}},
        benchmark_config_metadata={
            "source_sha256": "abc",
            "effective_request": {"route": "segmentation"},
            "cli_overrides": {"st_file": "override.h5ad"},
        },
    )

    assert result["ok"] is True
    assert result["profile"] == "benchmark_seg"
    assert captured["svc_type"] is None
    assert captured["cf"] == "segmentation"
    assert captured["runtime_overrides"] == {"seed": 42}
    expected_metadata = {
        "source_sha256": "abc",
        "effective_request": {
            "route": "segmentation",
            "io": {"sample_name": "sample"},
            "algorithm": {"local_refinement": {"strength": 0.0}},
            "runtime": {"seed": 42},
        },
        "cli_overrides": {"st_file": "override.h5ad"},
    }
    expected_metadata["effective_request_hash"] = hash_jsonable(
        expected_metadata["effective_request"]
    )
    assert captured["benchmark_config_metadata"] == expected_metadata


def test_seed_scope_preserves_fixed_run_seed_and_distinct_process_seeds():
    from revise.benchmark.cli import _runtime_seed_supplier

    args = argparse.Namespace(seed=42)
    fixed = _runtime_seed_supplier(args, "run")
    streamed = _runtime_seed_supplier(args, "process")

    assert fixed() == fixed() == 42
    assert streamed() != streamed()
