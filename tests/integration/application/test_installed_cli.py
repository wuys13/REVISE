from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

import pytest
from packaging.requirements import Requirement

from revise import __version__


ROOT = Path(__file__).resolve().parents[3]


def _run(command, *, cwd, env):
    return subprocess.run(
        [str(item) for item in command],
        cwd=cwd,
        env=env,
        text=True,
        capture_output=True,
        check=False,
        timeout=180,
    )


def _copy_installed_template(installed_cli, filename: str, destination: Path) -> None:
    script = (
        "from importlib.resources import files; from pathlib import Path; "
        f"Path({str(destination)!r}).write_bytes("
        f"files('revise.application').joinpath('templates', {filename!r}).read_bytes())"
    )
    copy = _run(
        [installed_cli["python"], "-c", script],
        cwd=destination.parent,
        env=installed_cli["env"],
    )
    assert copy.returncode == 0, copy.stderr


@pytest.fixture(scope="module")
def installed_cli(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("installed-cli")
    integration_python = Path(
        os.environ.get("REVISE_INTEGRATION_PYTHON", sys.executable)
    )
    wheel_dir = tmp_path / "wheel"
    wheel_dir.mkdir()
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env["NUMBA_CACHE_DIR"] = str(tmp_path / "numba-cache")
    env["MPLCONFIGDIR"] = str(tmp_path / "mpl-cache")
    supplied_wheel = os.environ.get("REVISE_WHEEL")
    if supplied_wheel:
        wheel = Path(supplied_wheel).resolve()
        assert wheel.is_file()
    else:
        source = tmp_path / "source"
        shutil.copytree(
            ROOT,
            source,
            ignore=shutil.ignore_patterns(
                ".git",
                ".agents",
                ".codegraph",
                ".codex",
                "build",
                "dist",
                "output",
                "raw_data",
                "results",
                "*.egg-info",
                "__pycache__",
                ".pytest_cache",
            ),
        )
        build = _run(
            [
                integration_python,
                "-m",
                "pip",
                "wheel",
                ".",
                "--no-deps",
                "--no-build-isolation",
                "--wheel-dir",
                str(wheel_dir),
            ],
            cwd=source,
            env=env,
        )
        assert build.returncode == 0, build.stderr
        wheel = next(wheel_dir.glob("revise_svc-*.whl"))

    venv = tmp_path / "venv"
    create = _run(
        [integration_python, "-m", "venv", "--system-site-packages", venv],
        cwd=tmp_path,
        env=env,
    )
    assert create.returncode == 0, create.stderr
    python = venv / "bin" / "python"
    install = _run(
        [
            str(python),
            "-m",
            "pip",
            "install",
            "--force-reinstall",
            "--no-deps",
            str(wheel),
        ],
        cwd=tmp_path,
        env=env,
    )
    assert install.returncode == 0, install.stderr

    probe = _run(
        [python, "-c", "import revise; print(revise.__file__)"],
        cwd=tmp_path,
        env=env,
    )
    assert probe.returncode == 0, probe.stderr
    installed_package = Path(probe.stdout.strip())
    assert str(venv) in str(installed_package)
    return {
        "root": tmp_path,
        "python": python,
        "source_python": integration_python,
        "command": venv / "bin" / "revise-reconstruct",
        "env": env,
        "wheel": wheel,
    }


def test_built_wheel_has_canonical_metadata_and_contents(installed_cli):
    with ZipFile(installed_cli["wheel"]) as archive:
        names = archive.namelist()
        metadata_name = next(name for name in names if name.endswith(".dist-info/METADATA"))
        entry_points_name = next(
            name for name in names if name.endswith(".dist-info/entry_points.txt")
        )
        metadata = archive.read(metadata_name).decode("utf-8")
        entry_points = archive.read(entry_points_name).decode("utf-8")

    assert f"Version: {__version__}" in metadata
    assert "Name: revise-svc" in metadata
    assert "Requires-Python: <3.12,>=3.10" in metadata
    assert "reconstruct.py" in names
    assert "revise/application/config.py" in names
    assert "revise/benchmark/cli.py" in names
    assert "revise/benchmark/launcher.py" in names
    assert "revise/revise.yaml" not in names
    application_templates = {
        "revise/application/templates/VisiumHD.yaml",
        "revise/application/templates/Xenium.yaml",
        "revise/application/templates/Visium.yaml",
    }
    benchmark_templates = {
        f"revise/benchmark/templates/{route}.yaml"
        for route in (
            "segmentation",
            "bin2cell",
            "batch_effect",
            "spot_size",
            "gene_panel",
            "gene_dropout",
        )
    }
    assert {
        name for name in names if name.startswith("revise/application/templates/")
    } == application_templates
    assert {
        name for name in names if name.startswith("revise/benchmark/templates/")
    } == benchmark_templates
    assert any(".dist-info/" in name and name.endswith("/LICENSE") for name in names)
    assert not any(name.startswith("tests/") for name in names)
    assert "revise-reconstruct = reconstruct:main" in entry_points

    requirements = [
        Requirement(line.removeprefix("Requires-Dist: "))
        for line in metadata.splitlines()
        if line.startswith("Requires-Dist: ")
    ]
    base = {requirement.name.lower() for requirement in requirements if requirement.marker is None}
    assert {"pot", "leidenalg"} <= base
    assert not base & {
        "tacco",
        "gseapy",
        "networkx",
        "omicverse",
        "cellphonedb",
        "spatialdata",
    }
    omicverse_stack = {
        "omicverse",
        "torch-geometric",
        "torch",
        "setuptools",
        "transformers",
    }
    expected_extras = {
        "tacco": {"tacco"},
        "pathway": {"gseapy", "networkx", *omicverse_stack},
        "cci": {"cellphonedb", *omicverse_stack},
        "trajectory": omicverse_stack,
        "spatialdata": {"spatialdata"},
    }
    for extra, expected in expected_extras.items():
        selected = {
            requirement.name.lower()
            for requirement in requirements
            if requirement.marker is not None
            and requirement.marker.evaluate({"extra": extra})
        }
        assert selected == expected


def test_built_wheel_installs_console_help(installed_cli):
    command = installed_cli["command"]
    env = installed_cli["env"]
    cwd = installed_cli["root"]

    help_result = _run([command, "--help"], cwd=cwd, env=env)
    assert help_result.returncode == 0, help_result.stderr
    assert "--config CONFIG" in help_result.stdout
    assert "--select-ct SELECT_CT" in help_result.stdout
    assert "--svc-type" not in help_result.stdout
    assert "--set" not in help_result.stdout


def test_installed_cli_normalizes_cluster_override_before_reconstruction(installed_cli):
    root = installed_cli["root"] / "xenium-override"
    root.mkdir()
    config_path = root / "Xenium.yaml"
    _copy_installed_template(installed_cli, "Xenium.yaml", config_path)
    probe = (
        "import json; "
        "from revise.application.config import (compile_application_config, "
        "load_application_yaml, override_select_cell_type); "
        "from revise.application.publication import output_paths; "
        f"source, document = load_application_yaml({str(config_path)!r}); "
        "config = compile_application_config(document, source=source); "
        "config = override_select_cell_type(config, 'Mono/Macro'); "
        "print(json.dumps({'select_cell_type': config.select_cell_type, "
        "'output_dir': str(config.output_dir), "
        "'outputs': {key: str(path) for key, path in output_paths(config).items()}}))"
    )
    result = _run(
        [installed_cli["python"], "-c", probe],
        cwd=root,
        env=installed_cli["env"],
    )

    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    expected_dir = root / "results" / "sc_SVC_case" / "P2CRC_Xenium" / "Mono_Macro"
    assert payload == {
        "select_cell_type": "Mono_Macro",
        "output_dir": str(expected_dir),
        "outputs": {
            "spatial": str(expected_dir / "spatial.h5ad"),
            "expression": str(expected_dir / "expr.h5ad"),
        },
    }



def test_installed_wheel_benchmark_module_and_refinement_option(installed_cli):
    python = installed_cli["python"]
    root = installed_cli["root"]
    env = installed_cli["env"]

    help_result = _run(
        [python, "-m", "revise.benchmark.cli", "--help"],
        cwd=root,
        env=env,
    )
    assert help_result.returncode == 0, help_result.stderr
    assert "--local-refinement-strength" in help_result.stdout

    probe = _run(
        [
            python,
            "-c",
            (
                "import json, revise.benchmark.cli as cli;"
                "print(json.dumps({'module': cli.__file__}))"
            ),
        ],
        cwd=root,
        env=env,
    )
    assert probe.returncode == 0, probe.stderr
    payload = json.loads(probe.stdout)
    assert str(installed_cli["root"] / "venv") in payload["module"]
