from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import pytest


def _args(**overrides):
    values = {
        "confounding": "segmentation",
        "posterior_mode": "cost",
        "posterior_key": None,
        "posterior_beta": None,
        "posterior_min_affinity": None,
        "posterior_cost_strength": None,
        "posterior_strict": False,
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


def test_benchmark_options_build_typed_algorithm_overrides():
    from revise.benchmark.cli import _build_algorithm_overrides

    overrides = _build_algorithm_overrides(
        _args(
            posterior_key="major_type",
            posterior_beta=2.5,
            posterior_min_affinity=0.1,
            posterior_cost_strength=3.0,
            posterior_strict=True,
        )
    )

    assert overrides == {
        "posterior_conditioning": {
            "enabled": True,
            "mode": "cost",
            "posterior_key": "major_type",
            "beta": 2.5,
            "min_affinity": 0.1,
            "cost_strength": 3.0,
            "strict": True,
        }
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
                provenance = {"run_dir": "run"}

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
