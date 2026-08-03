"""Application public vocabulary and entrypoint contract.

Covers: the 2.0 reconstruction selector, route mapping, and shared root/package CLI.
Proof limit: does not execute scientific reconstruction or validate real datasets.
"""

from __future__ import annotations

import argparse
from contextlib import redirect_stdout
from io import StringIO
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest


ROOT = Path(__file__).resolve().parents[2]


def _required_args(svc_type: str) -> list[str]:
    args = [
        "--svc-type",
        svc_type,
        "--sample-name",
        "sample",
        "--st-file",
        "st.h5ad",
        "--sc-ref-file",
        "sc.h5ad",
        "--data-root",
        "data",
    ]
    if svc_type == "sc-SVC":
        args.extend(["--select-ct", "T"])
    return args


def test_parser_exposes_only_the_three_2x_svc_types():
    from revise.application.cli import build_parser

    parser = build_parser()
    actions = {option: action for action in parser._actions for option in action.option_strings}

    assert "--platform" not in actions
    assert "--set" not in actions
    assert actions["--svc-type"].choices == ("hST-SVC", "iST-SVC", "sST-SVC")
    assert actions["--ist-mapping"].choices == ("mean", "random")
    assert isinstance(actions["--select-ct"], argparse._AppendAction)


@pytest.mark.parametrize("svc_type", ("hST-SVC", "iST-SVC", "sST-SVC"))
def test_parser_accepts_each_2x_svc_type(svc_type):
    from revise.application.cli import parse_args

    args = parse_args(_required_args(svc_type))

    assert args.svc_type == svc_type
    assert not hasattr(args, "platform")
    assert args.ist_mapping == ("mean" if svc_type == "iST-SVC" else None)


def test_parser_rejects_removed_sc_mapping_option():
    from revise.application.cli import parse_args

    with pytest.raises(SystemExit) as exc_info:
        parse_args(_required_args("iST-SVC") + ["--sc-mapping", "mean"])

    assert exc_info.value.code == 2


@pytest.mark.parametrize("svc_type", ("sp-SVC", "sc-SVC", "sc-SVC-sr"))
def test_parser_rejects_each_legacy_svc_type(svc_type):
    from revise.application.cli import parse_args

    with pytest.raises(SystemExit) as exc_info:
        parse_args(_required_args(svc_type))

    assert exc_info.value.code == 2


def test_parser_collects_repeated_iST_cell_types():
    from revise.application.cli import parse_args

    args = parse_args(
        _required_args("iST-SVC")
        + ["--select-ct", "T", "--select-ct", "B", "--select-ct", "T"]
    )

    assert args.select_ct == ["T", "B", "T"]


@pytest.mark.parametrize("svc_type", ("hST-SVC", "sST-SVC"))
def test_parser_rejects_iST_cell_type_selection_outside_iST(svc_type):
    from revise.application.cli import parse_args

    with pytest.raises(SystemExit) as exc_info:
        parse_args(_required_args(svc_type) + ["--select-ct", "T"])

    assert exc_info.value.code == 2


@pytest.mark.parametrize("svc_type", ("hST-SVC", "sST-SVC"))
def test_ist_mapping_is_rejected_outside_ist(svc_type):
    from revise.application.cli import parse_args

    with pytest.raises(SystemExit) as exc_info:
        parse_args(_required_args(svc_type) + ["--ist-mapping", "random"])

    assert exc_info.value.code == 2


@pytest.mark.parametrize("mapping", ("mean", "random"))
def test_ist_mapping_accepts_both_public_modes(mapping):
    from revise.application.cli import parse_args

    args = parse_args(_required_args("iST-SVC") + ["--ist-mapping", mapping])

    assert args.ist_mapping == mapping


def test_ot_method_help_explains_tacco_default_and_explicit_pot():
    from revise.application.cli import build_parser

    help_text = " ".join(build_parser().format_help().split())

    assert "iST-SVC defaults to TACCO" in help_text
    assert "'pot' explicitly selects a different algorithm" in help_text


@pytest.mark.parametrize(
    ("svc_type", "route"),
    (
        ("hST-SVC", "sp_svc:bin2cell"),
        ("iST-SVC", "sc_svc:segmentation"),
        ("sST-SVC", "sc_svc_sr:spot_size"),
    ),
)
def test_cli_reports_the_unified_single_output(
    monkeypatch,
    tmp_path,
    svc_type,
    route,
):
    from revise.application import cli, service

    monkeypatch.setattr(
        sys,
        "argv",
        ["revise-reconstruct", *_required_args(svc_type)],
    )
    monkeypatch.setattr(
        service,
        "reconstruct",
        lambda args: (
            SimpleNamespace(shape=(3, 5)),
            tmp_path / "SVC.h5ad",
            {"route": route},
        ),
    )

    stdout = StringIO()
    with redirect_stdout(stdout):
        cli.main()

    payload = json.loads(stdout.getvalue())
    assert payload == {
        "status": "succeeded",
        "svc_type": svc_type,
        "output": str(tmp_path / "SVC.h5ad"),
        "shape": [3, 5],
        "pipeline": {"route": route},
    }


def test_parser_rejects_removed_set_option():
    from revise.application.cli import parse_args

    with pytest.raises(SystemExit) as exc_info:
        parse_args(_required_args("hST-SVC") + ["--set", "graph.method=pca"])

    assert exc_info.value.code == 2


def test_omitted_cli_config_values_leave_the_yaml_in_control():
    from revise.application.cli import parse_args
    from revise.application.service import _build_algorithm_overrides

    args = parse_args(_required_args("hST-SVC"))

    assert args.seed is None
    assert args.cell_type_col is None
    assert args.sub_cell_type_col is None
    assert _build_algorithm_overrides(args) == {}


def test_application_routes_map_2x_vocabulary_to_unchanged_internal_ids():
    from revise.application.service import APPLICATION_ROUTES

    assert set(APPLICATION_ROUTES) == {"hST-SVC", "iST-SVC", "sST-SVC"}
    assert {
        svc_type: route.route_id for svc_type, route in APPLICATION_ROUTES.items()
    } == {
        "hST-SVC": "sp_svc",
        "iST-SVC": "sc_svc",
        "sST-SVC": "sc_svc_sr",
    }
    serialized = repr(APPLICATION_ROUTES)
    assert "application_sp" in serialized
    assert "application_sc" in serialized
    assert "application_sc_sr" in serialized


@pytest.mark.parametrize(
    ("svc_type", "profile", "route_id", "confounding", "output_key"),
    (
        ("hST-SVC", "application_sp", "sp_svc", "bin2cell", "sp_svc"),
        ("iST-SVC", "application_sc", "sc_svc", "segmentation", None),
        (
            "sST-SVC",
            "application_sc_sr",
            "sc_svc_sr",
            "spot_size",
            "sc_svc_dec",
        ),
    ),
)
@pytest.mark.parametrize("seed", (17, None))
def test_pipeline_receives_the_internal_route_for_each_public_type(
    monkeypatch,
    tmp_path,
    svc_type,
    profile,
    route_id,
    confounding,
    output_key,
    seed,
):
    from revise.application import service

    captured = {}

    class Pipeline:
        def __init__(self, config_path):
            captured["config_path"] = config_path

        def _run_with_algorithm_overrides(self, **kwargs):
            captured.update(kwargs)
            return SimpleNamespace()

    monkeypatch.setattr(service, "REVISEPipeline", Pipeline)
    args = argparse.Namespace(
        svc_type=svc_type,
        config="revise/revise.yaml",
        seed=seed,
        data_root=str(tmp_path),
        output_root=str(tmp_path / "out"),
        sample_name="sample",
        st_file="st.h5ad",
        sc_ref_file="sc.h5ad",
        patient_key="Patient",
        ot_method="pot",
        ist_mapping="mean" if svc_type == "iST-SVC" else None,
        select_ct=None,
        cell_type_col="Level1",
        sub_cell_type_col="Level2",
    )

    actual_profile, actual_output_key, _ = service._run_pipeline(args)

    assert actual_profile == profile
    assert actual_output_key == output_key
    expected_runtime = {
        "platform": route_id,
        "confounding": confounding,
    }
    if seed is not None:
        expected_runtime["seed"] = seed
    assert captured["runtime_overrides"] == expected_runtime
    assert captured["io_overrides"]["save_outputs"] is False
    expected_overrides = {
        "ot": {"ga": {"solver": "pot"}, "lr": {"solver": "pot"}},
        "columns": {
            "cell_type_col": "Level1",
            "sub_cell_type_col": "Level2",
        },
    }
    if svc_type == "iST-SVC":
        expected_overrides["sc"] = {"selection_review_gate": True}
    assert captured["algorithm_overrides"] == expected_overrides


def test_sst_forwards_configured_cell_type_columns():
    from revise.application.service import _build_algorithm_overrides

    args = SimpleNamespace(
        svc_type="sST-SVC",
        ot_method=None,
        cell_type_col="major_type",
        sub_cell_type_col="minor_type",
    )

    overrides = _build_algorithm_overrides(args)

    assert overrides == {
        "columns": {
            "cell_type_col": "major_type",
            "sub_cell_type_col": "minor_type",
        }
    }


def test_ist_cell_type_override_deduplicates_in_first_seen_order():
    from revise.application.service import _build_algorithm_overrides

    args = SimpleNamespace(
        svc_type="iST-SVC",
        ot_method=None,
        local_refinement_strength=None,
        cell_type_col=None,
        sub_cell_type_col=None,
        select_ct=["T", "B", "T"],
    )

    assert _build_algorithm_overrides(args) == {
        "sc": {
            "selection_review_gate": True,
            "select_ct": ["T", "B"],
        }
    }


def test_hst_forwards_configured_annotation_columns():
    from revise.application.service import _build_algorithm_overrides

    args = SimpleNamespace(
        svc_type="hST-SVC",
        ot_method=None,
        cell_type_col="major_type",
        sub_cell_type_col="minor_type",
    )

    overrides = _build_algorithm_overrides(args)

    assert overrides == {
        "columns": {
            "cell_type_col": "major_type",
            "sub_cell_type_col": "minor_type",
        }
    }


def test_root_application_wrapper_delegates_to_package_cli():
    import reconstruct
    from revise.application import cli

    assert reconstruct.main is cli.main
    assert not hasattr(reconstruct, "APPLICATION_ROUTES")


def test_root_application_help_uses_2x_vocabulary():
    result = subprocess.run(
        [sys.executable, "reconstruct.py", "--help"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "--svc-type {hST-SVC,iST-SVC,sST-SVC}" in result.stdout
    assert "--ist-mapping {mean,random}" in result.stdout
    assert "--select-ct" in result.stdout
    assert "--platform" not in result.stdout
