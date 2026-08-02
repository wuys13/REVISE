from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from revise.config import load_raw_config


CONFIG_PATH = Path(__file__).parents[2] / "revise" / "revise.yaml"


def _args(**overrides):
    values = {
        "confounding": "segmentation",
        "local_refinement_guidance": None,
        "local_refinement_compatibility_mode": None,
        "posterior_mode": None,
        "posterior_key": None,
        "posterior_beta": None,
        "posterior_min_affinity": None,
        "posterior_cost_strength": None,
        "posterior_strict": False,
        "sr_refinement_preset": None,
    }
    values.update(overrides)
    return argparse.Namespace(**values)


def _guidance_manifest(*outcomes):
    return {
        "schema_version": 2,
        "configured": {
            "guidance": None,
            "compatibility_mode": None,
            "source": "route_default",
            "deprecations": [],
        },
        "resolved": {
            "guidance": "prefer",
            "compatibility_mode": "cost",
            "beta": 1.0,
            "min_affinity": 0.05,
            "operator_strength": 1.0,
        },
        "events": [
            {
                "ordinal": index,
                "problem_key": "shared-problem",
                "outcome": outcome,
            }
            for index, outcome in enumerate(outcomes, start=1)
        ],
        "summary": outcomes[-1] if outcomes else "not_started",
    }


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

    args = cli.get_args()

    assert not hasattr(args, "set_overrides")


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


def test_legacy_posterior_key_help_states_rejection_and_replacement():
    from revise.benchmark import cli

    help_text = cli.get_parser().format_help()
    assert "--posterior-key" in help_text
    assert "Deprecated" in help_text
    assert "non-empty values are rejected" in help_text
    assert "route-provided Assignment State" in help_text


@pytest.mark.parametrize(
    ("profile", "confounding", "overrides", "expected"),
    [
        (
            "benchmark_seg",
            "segmentation",
            {},
            {"posterior_mode": "cost", "posterior_strict": False},
        ),
        (
            "benchmark_sr_spot_size",
            "spot_size",
            {},
            {"posterior_mode": "off", "posterior_strict": False},
        ),
        (
            "benchmark_seg",
            "segmentation",
            {
                "local_refinement": {
                    "guidance": "require",
                    "compatibility": {"mode": "reference"},
                }
            },
            {"posterior_mode": "reference", "posterior_strict": True},
        ),
    ],
)
def test_report_aliases_derive_from_actual_resolved_route_config(
    profile,
    confounding,
    overrides,
    expected,
):
    from revise.benchmark.cli import _resolved_report_aliases

    aliases = _resolved_report_aliases(
        raw_config=load_raw_config(CONFIG_PATH),
        profile=profile,
        platform="sim2real",
        confounding=confounding,
        algorithm_overrides=overrides,
    )

    assert {
        "posterior_mode": aliases["posterior_mode"],
        "posterior_strict": aliases["posterior_strict"],
    } == expected
    guidance = aliases["assignment_guidance"]
    assert guidance["schema_version"] == 2
    assert guidance["events"] == []
    assert guidance["summary"] == "not_started"


@pytest.mark.parametrize(
    ("leaf_outcomes", "expected_summary"),
    [
        ([("applied",), ("applied",)], "applied"),
        ([("applied",), ("fallback",)], "mixed"),
        ([("applied",), ("failed",)], "failed"),
    ],
)
def test_assignment_guidance_report_aggregates_leaf_events_with_case_identity(
    leaf_outcomes,
    expected_summary,
):
    from revise.benchmark.cli import _aggregate_assignment_guidance

    request = _guidance_manifest()
    results = [
        {
            "tag": f"case-{index}",
            "assignment_guidance": _guidance_manifest(*outcomes),
        }
        for index, outcomes in enumerate(leaf_outcomes, start=1)
    ]

    aggregate = _aggregate_assignment_guidance(
        results,
        request_manifest=request,
    )

    assert aggregate["schema_version"] == 2
    assert aggregate["configured"] == request["configured"]
    assert aggregate["resolved"] == request["resolved"]
    assert aggregate["summary"] == expected_summary
    assert [event["ordinal"] for event in aggregate["events"]] == [1, 2]
    assert [event["case_ordinal"] for event in aggregate["events"]] == [1, 2]
    assert [event["case_tag"] for event in aggregate["events"]] == [
        "case-1",
        "case-2",
    ]
    assert [event["leaf_ordinal"] for event in aggregate["events"]] == [1, 1]
    assert [event["problem_key"] for event in aggregate["events"]] == [
        "case-1::shared-problem",
        "case-2::shared-problem",
    ]


def test_assignment_guidance_report_zero_results_keeps_request_not_started():
    from revise.benchmark.cli import _aggregate_assignment_guidance

    request = _guidance_manifest()

    assert _aggregate_assignment_guidance(
        [],
        request_manifest=request,
    ) == request


def test_assignment_guidance_report_rejects_inconsistent_leaf_resolution():
    from revise.benchmark.cli import _aggregate_assignment_guidance

    request = _guidance_manifest()
    inconsistent = _guidance_manifest("applied")
    inconsistent["resolved"]["guidance"] = "require"

    with pytest.raises(
        ValueError,
        match="inconsistent assignment-guidance resolved",
    ):
        _aggregate_assignment_guidance(
            [
                {
                    "tag": "case-1",
                    "assignment_guidance": _guidance_manifest("applied"),
                },
                {
                    "tag": "case-2",
                    "assignment_guidance": inconsistent,
                },
            ],
            request_manifest=request,
        )


@pytest.mark.parametrize(
    ("confounding", "guidance", "mode", "expected_mode", "expected_strict"),
    [
        ("segmentation", None, None, "cost", False),
        ("spot_size", None, None, "off", False),
        ("segmentation", "require", "reference", "reference", True),
    ],
)
def test_main_reports_resolved_config_aliases(
    monkeypatch,
    confounding,
    guidance,
    mode,
    expected_mode,
    expected_strict,
):
    from revise.benchmark import cli

    printed = []

    class Pipeline:
        def __init__(self, config_path):
            self.config_path = config_path
            self.raw_config = load_raw_config(CONFIG_PATH)

    def fake_run_case(_pipeline, **kwargs):
        guidance_evidence = cli._resolved_report_aliases(
            raw_config=_pipeline.raw_config,
            profile=kwargs["profile"],
            platform=kwargs["platform"],
            confounding=kwargs["confounding"],
            algorithm_overrides=kwargs["algorithm_overrides"],
        )["assignment_guidance"]
        outcome = (
            "off"
            if guidance_evidence["resolved"]["guidance"] == "off"
            else "applied"
        )
        guidance_evidence["events"] = [
            {
                "ordinal": 1,
                "problem_key": "shared-problem",
                "outcome": outcome,
            }
        ]
        guidance_evidence["summary"] = outcome
        return {
            "ok": True,
            "profile": kwargs["profile"],
            "seed": kwargs["runtime_seed"],
            "run_dir": "run",
            "assignment_guidance": guidance_evidence,
            "summary": {},
            "error": None,
        }

    monkeypatch.setattr(cli, "REVISEPipeline", Pipeline)
    monkeypatch.setattr(cli, "_run_case", fake_run_case)
    monkeypatch.setattr(
        "builtins.print",
        lambda *args, **_kwargs: printed.append(" ".join(map(str, args))),
    )
    args = argparse.Namespace(
        platform="sim2real",
        confounding=confounding,
        data_root="data",
        dataset_task=None,
        sample_name="sample",
        st_file=None,
        gt_svc_file=None,
        sc_ref_file=None,
        output_root="output",
        sample_size=None,
        config=str(CONFIG_PATH),
        seed=42,
        seed_scope="run",
        local_refinement_guidance=guidance,
        local_refinement_compatibility_mode=mode,
        posterior_mode=None,
        posterior_key=None,
        posterior_beta=None,
        posterior_min_affinity=None,
        posterior_cost_strength=None,
        posterior_strict=False,
        sr_refinement_preset=None,
    )

    cli.main(args)

    report = json.loads(printed[-1])
    assert report["posterior_mode"] == expected_mode
    assert report["posterior_strict"] is expected_strict
    guidance = report["assignment_guidance"]
    assert guidance["schema_version"] == 2
    assert guidance["resolved"]["guidance"] == (
        "require" if expected_strict else (
            "off" if expected_mode == "off" else "prefer"
        )
    )
    assert guidance["resolved"]["compatibility_mode"] == (
        None if expected_mode == "off" else (mode or expected_mode)
    )
    expected_summary = "off" if expected_mode == "off" else "applied"
    assert guidance["summary"] == expected_summary
    assert guidance["events"]
    assert all("case_tag" in event for event in guidance["events"])
    assert len({event["problem_key"] for event in guidance["events"]}) == len(
        guidance["events"]
    )


def test_absent_guidance_flags_build_no_algorithm_override(monkeypatch):
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

    overrides = cli._build_algorithm_overrides(cli.get_args())

    # This consciously replaces U1's expected-to-change characterization:
    # omission now preserves the route-specific resolver default.
    assert overrides == {}


@pytest.mark.parametrize(
    ("guidance", "mode"),
    [
        ("off", None),
        ("prefer", "cost"),
        ("require", "cost"),
        ("prefer", "reference"),
        ("require", "reference"),
    ],
)
def test_new_guidance_options_build_independent_typed_overrides(guidance, mode):
    from revise.benchmark.cli import _build_algorithm_overrides

    overrides = _build_algorithm_overrides(
        _args(
            local_refinement_guidance=guidance,
            local_refinement_compatibility_mode=mode,
        )
    )

    assert overrides == {
        "local_refinement": {
            "guidance": guidance,
            **({"compatibility": {"mode": mode}} if mode is not None else {}),
        }
    }


def test_legacy_benchmark_options_translate_through_deprecated_raw_section():
    from revise.benchmark.cli import _build_algorithm_overrides

    overrides = _build_algorithm_overrides(
        _args(
            posterior_mode="cost",
            posterior_beta=2.5,
            posterior_min_affinity=0.1,
            posterior_cost_strength=3.0,
            posterior_strict=True,
        )
    )

    assert overrides == {
        "posterior_conditioning": {
            "mode": "cost",
            "beta": 2.5,
            "min_affinity": 0.1,
            "cost_strength": 3.0,
            "strict": True,
        }
    }


def test_legacy_posterior_key_is_rejected_by_cli_builder():
    from revise.benchmark.cli import _build_algorithm_overrides

    with pytest.raises(ValueError, match=r"posterior-key.*Assignment State"):
        _build_algorithm_overrides(_args(posterior_key="Level1"))


@pytest.mark.parametrize(
    "args",
    [
        _args(posterior_mode="off", posterior_strict=True),
        _args(
            local_refinement_guidance="off",
            posterior_mode="cost",
        ),
        _args(
            local_refinement_compatibility_mode="cost",
            posterior_mode="reference",
        ),
    ],
)
def test_conflicting_new_and_legacy_cli_flags_are_rejected(args):
    from revise.benchmark.cli import _build_algorithm_overrides

    with pytest.raises(ValueError, match="conflict"):
        _build_algorithm_overrides(args)


@pytest.mark.parametrize("posterior_mode", ["cost", "reference"])
def test_legacy_sr_prefer_mode_with_graph_preset_none_reaches_preflight(
    posterior_mode,
):
    from revise.benchmark.cli import _build_algorithm_overrides

    overrides = _build_algorithm_overrides(
        _args(
            confounding="spot_size",
            posterior_mode=posterior_mode,
            sr_refinement_preset="none",
        )
    )

    assert overrides == {
        "posterior_conditioning": {"mode": posterior_mode},
        "sc": {
            "sr_graph_agg_enabled": False,
            "sr_graph_agg_low_conf_only": False,
            "sr_graph_agg_anchor_only": False,
            "sr_graph_agg_conf_weighted_alpha": False,
        },
    }


def test_sr_refinement_preset_builds_nested_sc_configuration():
    from revise.benchmark.cli import _build_algorithm_overrides

    overrides = _build_algorithm_overrides(
        _args(
            confounding="spot_size",
            sr_refinement_preset="confidence_anchor",
        )
    )

    assert overrides["sc"] == {
        "sr_graph_agg_enabled": True,
        "sr_graph_agg_low_conf_only": True,
        "sr_graph_agg_low_conf_quantile": 0.1,
        "sr_graph_agg_anchor_only": True,
        "sr_graph_agg_anchor_high_conf_quantile": 0.9,
        "sr_graph_agg_confidence_mode": "auto",
        "sr_graph_agg_conf_weighted_alpha": True,
        "sr_graph_agg_conf_alpha_min": 0.0,
        "sr_graph_agg_conf_alpha_max": 0.25,
        "sr_graph_agg_conf_alpha_power": 1.0,
    }


def test_run_case_passes_none_seed_and_structured_algorithm_overrides():
    from revise.benchmark.cli import _run_case

    captured = {}

    class Pipeline:
        def _run_with_algorithm_overrides(self, **kwargs):
            captured.update(kwargs)

            class SVC:
                provenance = {
                    "run_dir": "run",
                    "assignment_guidance": {
                        "schema_version": 2,
                        "summary": "applied",
                    },
                }

                @staticmethod
                def summary():
                    return {"ok": True}

            return SVC()

    overrides = {"posterior_conditioning": {"enabled": False}}
    result = _run_case(
        Pipeline(),
        platform="sim2real",
        profile="benchmark_seg",
        confounding="segmentation",
        io_overrides={"sample_name": "sample"},
        runtime_seed=None,
        algorithm_overrides=overrides,
    )

    assert result["ok"] is True
    assert captured["runtime_overrides"]["seed"] is None
    assert captured["algorithm_overrides"] == overrides
    assert "set_overrides" not in captured
    assert result["manifest_path"] == "run/provenance.json"
    assert result["assignment_guidance"]["summary"] == "applied"


def test_pre_context_failure_does_not_reuse_stale_succeeded_manifest(tmp_path):
    from revise.benchmark.cli import _run_case

    output_root = tmp_path / "output"
    manifest = output_root / "sample" / "seg_1" / "provenance.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text(
        json.dumps(
            {
                "run": {"status": "succeeded"},
                "assignment_guidance": {
                    "schema_version": 2,
                    "summary": "applied",
                },
            }
        ),
        encoding="utf-8",
    )

    class Pipeline:
        raw_config = load_raw_config(CONFIG_PATH)

        def _run_with_algorithm_overrides(self, **kwargs):
            raise RuntimeError("pre-context failure")

    result = _run_case(
        Pipeline(),
        platform="sim2real",
        profile="benchmark_seg",
        confounding="segmentation",
        io_overrides={
            "output_root": str(output_root),
            "sample_name": "sample",
            "seg_method": "seg_1",
        },
        runtime_seed=17,
        algorithm_overrides={},
    )

    assert result["ok"] is False
    assert result["run_dir"] is None
    assert result["manifest_path"] is None
    assert result["assignment_guidance"] is None


def test_concurrent_lock_failure_does_not_reuse_stale_manifest(
    tmp_path,
    monkeypatch,
):
    from contextlib import contextmanager

    from revise import framework
    from revise.benchmark.cli import _run_case
    from revise.framework import REVISEPipeline

    output_root = tmp_path / "output"
    manifest = output_root / "sample" / "seg_1" / "provenance.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text(
        json.dumps(
            {
                "run": {"status": "succeeded"},
                "assignment_guidance": {
                    "schema_version": 2,
                    "summary": "applied",
                },
            }
        ),
        encoding="utf-8",
    )

    @contextmanager
    def locked(_run_dir):
        raise RuntimeError("run directory is locked")
        yield

    monkeypatch.setattr(framework, "exclusive_run_directory", locked)

    result = _run_case(
        REVISEPipeline(config_path=CONFIG_PATH),
        platform="sim2real",
        profile="benchmark_seg",
        confounding="segmentation",
        io_overrides={
            "output_root": str(output_root),
            "sample_name": "sample",
            "seg_method": "seg_1",
        },
        runtime_seed=17,
        algorithm_overrides={},
    )

    assert result["ok"] is False
    assert result["run_dir"] is None
    assert result["manifest_path"] is None
    assert result["assignment_guidance"] is None


def test_post_context_failure_keeps_current_invocation_manifest(
    tmp_path,
    monkeypatch,
):
    from revise.benchmark.cli import _run_case
    from revise.framework import REVISEPipeline

    output_root = tmp_path / "output"
    pipeline = REVISEPipeline(config_path=CONFIG_PATH)

    def fail_in_directory(**kwargs):
        manifest = kwargs["run_dir"] / "provenance.json"
        manifest.write_text(
            json.dumps(
                {
                    "run": {"status": "failed"},
                    "assignment_guidance": {
                        "schema_version": 2,
                        "summary": "failed",
                    },
                }
            ),
            encoding="utf-8",
        )
        raise RuntimeError("post-context failure")

    monkeypatch.setattr(pipeline, "_run_in_directory", fail_in_directory)

    result = _run_case(
        pipeline,
        platform="sim2real",
        profile="benchmark_seg",
        confounding="segmentation",
        io_overrides={
            "output_root": str(output_root),
            "sample_name": "sample",
            "seg_method": "seg_1",
        },
        runtime_seed=17,
        algorithm_overrides={},
    )

    assert result["ok"] is False
    assert result["run_dir"] == str(output_root / "sample" / "seg_1")
    assert result["manifest_path"] == str(
        output_root / "sample" / "seg_1" / "provenance.json"
    )
    assert result["assignment_guidance"]["summary"] == "failed"


@pytest.mark.parametrize(
    ("requested_scope", "reported_scope"),
    [("process", "process"), (None, "process"), ("run", "run")],
)
def test_seed_scope_is_explicit_per_case_and_report_hides_internal_overrides(
    monkeypatch,
    requested_scope,
    reported_scope,
):
    from revise.benchmark import cli

    observed_seeds = []
    printed = []

    class Pipeline:
        def __init__(self, config_path):
            self.config_path = config_path
            self.raw_config = load_raw_config(CONFIG_PATH)

    def fake_run_case(_pipeline, **kwargs):
        observed_seeds.append(kwargs["runtime_seed"])
        return {
            "ok": True,
            "profile": kwargs["profile"],
            "seed": kwargs["runtime_seed"],
            "run_dir": f"run-{len(observed_seeds)}",
            "summary": {},
            "error": None,
        }

    monkeypatch.setattr(cli, "REVISEPipeline", Pipeline)
    monkeypatch.setattr(cli, "_run_case", fake_run_case)
    monkeypatch.setattr(
        "builtins.print",
        lambda *args, **_kwargs: printed.append(" ".join(map(str, args))),
    )

    args = argparse.Namespace(
        platform="sim2real",
        confounding="spot_size",
        data_root="data",
        dataset_task=None,
        sample_name="sample",
        st_file=None,
        gt_svc_file=None,
        sc_ref_file=None,
        output_root="output",
        sample_size=None,
        config="revise/revise.yaml",
        seed=731,
        seed_scope=requested_scope,
        posterior_mode="off",
        posterior_key=None,
        posterior_beta=None,
        posterior_min_affinity=None,
        posterior_cost_strength=None,
        posterior_strict=False,
        sr_refinement_preset=None,
    )

    cli.main(args)

    if reported_scope == "process":
        expected_rng = np.random.RandomState(731)
        expected_seeds = [
            int(expected_rng.randint(0, np.iinfo(np.int32).max))
            for _ in cli.SPOT_SIZES
        ]
    else:
        expected_seeds = [731] * len(cli.SPOT_SIZES)
    assert observed_seeds == expected_seeds
    assert all(seed is not None for seed in observed_seeds)

    report = json.loads(printed[-1])
    assert report["seed"] == 731
    assert report["seed_scope"] == reported_scope
    assert [item["seed"] for item in report["results"]] == expected_seeds
    assert "algorithm_overrides" not in report
