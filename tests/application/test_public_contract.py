"""Application public vocabulary and entrypoint contract.

Covers: the 1.x reconstruction selector, route mapping, and shared root/package CLI.
Proof limit: does not execute scientific reconstruction or validate real datasets.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest


ROOT = Path(__file__).resolve().parents[2]


def _required_args(svc_type: str) -> list[str]:
    return [
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


def test_parser_exposes_only_the_three_1x_svc_types():
    from revise.application.cli import build_parser

    parser = build_parser()
    actions = {option: action for action in parser._actions for option in action.option_strings}

    assert "--platform" not in actions
    assert "--set" not in actions
    assert actions["--svc-type"].choices == ("sp-SVC", "sc-SVC", "sc-SVC-sr")


@pytest.mark.parametrize("svc_type", ("sp-SVC", "sc-SVC", "sc-SVC-sr"))
def test_parser_accepts_each_1x_svc_type(svc_type):
    from revise.application.cli import parse_args

    args = parse_args(_required_args(svc_type))

    assert args.svc_type == svc_type
    assert not hasattr(args, "platform")


def test_parser_rejects_removed_sc_mapping_option():
    from revise.application.cli import parse_args

    with pytest.raises(SystemExit) as exc_info:
        parse_args(_required_args("sc-SVC") + ["--sc-mapping", "mean"])

    assert exc_info.value.code == 2


def test_ot_method_help_explains_tacco_default_and_explicit_pot():
    from revise.application.cli import build_parser

    help_text = build_parser().format_help()

    assert "standard sc-SVC defaults to TACCO" in help_text
    assert "'pot' explicitly selects a different algorithm" in help_text


def test_sc_svc_cli_reports_plural_outputs(monkeypatch, capsys, tmp_path):
    from revise.application import cli, service

    monkeypatch.setattr(
        sys,
        "argv",
        ["revise-reconstruct", *_required_args("sc-SVC")],
    )
    monkeypatch.setattr(
        service,
        "reconstruct",
        lambda args: (
            {
                "spatial": SimpleNamespace(shape=(3, 2)),
                "expression": SimpleNamespace(shape=(4, 5)),
            },
            {
                "spatial": tmp_path / "sc_SVC_spatial.h5ad",
                "expression": tmp_path / "sc_SVC_expr.h5ad",
            },
            {"route": "sc_svc:segmentation"},
        ),
    )

    cli.main()

    payload = json.loads(capsys.readouterr().out)
    assert payload["outputs"] == {
        "spatial": str(tmp_path / "sc_SVC_spatial.h5ad"),
        "expression": str(tmp_path / "sc_SVC_expr.h5ad"),
    }
    assert payload["shapes"] == {
        "spatial": [3, 2],
        "expression": [4, 5],
    }
    assert "output" not in payload
    assert "shape" not in payload


def test_parser_rejects_removed_set_option():
    from revise.application.cli import parse_args

    with pytest.raises(SystemExit) as exc_info:
        parse_args(_required_args("sp-SVC") + ["--set", "graph.method=pca"])

    assert exc_info.value.code == 2


def test_omitted_cli_config_values_leave_the_yaml_in_control():
    from revise.application.cli import parse_args
    from revise.application.service import _build_algorithm_overrides

    args = parse_args(_required_args("sp-SVC"))

    assert args.seed is None
    assert args.cell_type_col is None
    assert args.sub_cell_type_col is None
    assert _build_algorithm_overrides(args) == {}


def test_application_routes_use_internal_ids_without_2x_vocabulary():
    from revise.application.service import APPLICATION_ROUTES

    assert set(APPLICATION_ROUTES) == {"sp-SVC", "sc-SVC", "sc-SVC-sr"}
    assert {
        svc_type: route.route_id for svc_type, route in APPLICATION_ROUTES.items()
    } == {
        "sp-SVC": "sp_svc",
        "sc-SVC": "sc_svc",
        "sc-SVC-sr": "sc_svc_sr",
    }
    serialized = repr(APPLICATION_ROUTES)
    assert "hST" not in serialized
    assert "iST" not in serialized
    assert "sST" not in serialized


@pytest.mark.parametrize(
    ("svc_type", "profile", "route_id", "confounding", "output_key"),
    (
        ("sp-SVC", "application_sp", "sp_svc", "bin2cell", "sp_svc"),
        ("sc-SVC", "application_sc", "sc_svc", "segmentation", None),
        (
            "sc-SVC-sr",
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
        select_ct="all",
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
    assert captured["algorithm_overrides"] == {
        "ot": {"ga": {"solver": "pot"}, "lr": {"solver": "pot"}},
        "columns": {
            "cell_type_col": "Level1",
            "sub_cell_type_col": "Level2",
        },
        **({"sc": {"select_ct": "all"}} if svc_type == "sc-SVC" else {}),
    }


def test_sc_svc_sr_forwards_configured_cell_type_columns_without_sc_only_selection():
    from revise.application.service import _build_algorithm_overrides

    args = SimpleNamespace(
        svc_type="sc-SVC-sr",
        ot_method=None,
        select_ct="T",
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


def test_sp_svc_forwards_configured_annotation_columns():
    from revise.application.service import _build_algorithm_overrides

    args = SimpleNamespace(
        svc_type="sp-SVC",
        ot_method=None,
        select_ct="all",
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


def test_root_application_help_uses_1x_vocabulary():
    result = subprocess.run(
        [sys.executable, "reconstruct.py", "--help"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "--svc-type {sp-SVC,sc-SVC,sc-SVC-sr}" in result.stdout
    assert "--platform" not in result.stdout
