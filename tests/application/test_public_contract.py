"""Application public vocabulary and entrypoint contract.

Covers: the 2.0 platform selector, route mapping, and shared root/package CLI.
Proof limit: does not execute scientific reconstruction or validate real datasets.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest


ROOT = Path(__file__).resolve().parents[2]


def _required_args(platform: str) -> list[str]:
    return [
        "--platform",
        platform,
        "--sample-name",
        "sample",
        "--st-file",
        "st.h5ad",
        "--sc-ref-file",
        "sc.h5ad",
        "--data-root",
        "data",
    ]


def test_parser_exposes_only_the_three_2x_platforms():
    from revise.application.cli import build_parser

    parser = build_parser()
    actions = {option: action for action in parser._actions for option in action.option_strings}

    assert "--svc-type" not in actions
    assert actions["--platform"].choices == ("hST", "iST", "sST")


@pytest.mark.parametrize("platform", ("hST", "iST", "sST"))
def test_parser_accepts_each_2x_platform(platform):
    from revise.application.cli import parse_args

    args = parse_args(_required_args(platform))

    assert args.platform == platform
    assert not hasattr(args, "svc_type")


def test_application_platforms_map_to_stable_internal_ids():
    from revise.application.service import APPLICATION_ROUTES

    assert set(APPLICATION_ROUTES) == {"hST", "iST", "sST"}
    assert {
        platform: route.route_id for platform, route in APPLICATION_ROUTES.items()
    } == {
        "hST": "sp_svc",
        "iST": "sc_svc",
        "sST": "sc_svc_sr",
    }
@pytest.mark.parametrize(
    ("platform", "profile", "route_id", "confounding", "output_key"),
    (
        ("hST", "application_sp", "sp_svc", "bin2cell", "sp_svc"),
        ("iST", "application_sc", "sc_svc", "segmentation", None),
        (
            "sST",
            "application_sc_sr",
            "sc_svc_sr",
            "spot_size",
            "sc_svc_dec",
        ),
    ),
)
def test_pipeline_receives_the_internal_route_for_each_public_type(
    monkeypatch,
    tmp_path,
    platform,
    profile,
    route_id,
    confounding,
    output_key,
):
    from revise.application import service

    captured = {}

    class Pipeline:
        def __init__(self, config_path):
            captured["config_path"] = config_path

        def run(self, **kwargs):
            captured.update(kwargs)
            return SimpleNamespace()

    monkeypatch.setattr(service, "REVISEPipeline", Pipeline)
    args = argparse.Namespace(
        platform=platform,
        config="revise/revise.yaml",
        seed=17,
        data_root=str(tmp_path),
        output_root=str(tmp_path / "out"),
        sample_name="sample",
        st_file="st.h5ad",
        sc_ref_file="sc.h5ad",
        patient_key="Patient",
        ot_method="pot",
        select_ct="all",
        cell_type_col="Level1",
        sub_cell_type_col="Level2",
        set_overrides=[],
    )

    actual_profile, actual_output_key, _ = service._run_pipeline(args)

    assert actual_profile == profile
    assert actual_output_key == output_key
    assert captured["runtime_overrides"] == {
        "platform": route_id,
        "confounding": confounding,
        "seed": 17,
    }


def test_root_application_wrapper_delegates_to_package_cli():
    import application_reconstruct
    from revise.application import cli

    assert application_reconstruct.main is cli.main
    assert not hasattr(application_reconstruct, "APPLICATION_ROUTES")


def test_root_application_help_uses_2x_vocabulary():
    result = subprocess.run(
        [sys.executable, "application_reconstruct.py", "--help"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "--platform {hST,iST,sST}" in result.stdout
    assert "--svc-type" not in result.stdout
