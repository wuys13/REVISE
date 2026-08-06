from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
from zipfile import ZipFile

import numpy as np
import pandas as pd
import pytest
import yaml
from anndata import AnnData
from packaging.requirements import Requirement
from scipy import sparse

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


def _write_inputs(st_path: Path, reference_path: Path) -> None:
    st_path.parent.mkdir(parents=True, exist_ok=True)
    reference_path.parent.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(17)
    st = AnnData(
        X=sparse.csr_matrix(rng.poisson(3, size=(52, 52)) + 1),
        obs=pd.DataFrame(index=[f"spot-{index}" for index in range(52)]),
        var=pd.DataFrame(index=[f"g{index}" for index in range(52)]),
    )
    st.obsm["spatial"] = np.column_stack(
        [np.arange(52, dtype=float), np.arange(52, dtype=float) % 7]
    )
    sc_ref = AnnData(
        X=sparse.csr_matrix(rng.poisson(3, size=(52, 52)) + 1),
        obs=pd.DataFrame(
            {"Level1": ["A"] * 52, "Level2": ["A1"] * 52},
            index=[f"cell-{index}" for index in range(52)],
        ),
        var=pd.DataFrame(index=st.var_names.copy()),
    )
    st.write_h5ad(st_path)
    sc_ref.write_h5ad(reference_path)


def _write_application_config(
    path: Path,
    *,
    st_path: str,
    reference_path: str,
    output_path: str,
) -> None:
    path.write_text(
        yaml.safe_dump(
            {
                "schema_version": 1,
                "application": {"svc_type": "sp-SVC"},
                "paths": {"root_dir": "."},
                "algorithm": {"ot_method": "pot"},
                "inputs": {
                    "st": {"path": st_path, "format": "h5ad"},
                    "reference": {
                        "path": reference_path,
                        "format": "h5ad",
                    },
                },
                "global_anchoring": {"broad_column": "Level1"},
                "local_refinement": {"strength": 0.2},
                "output": {"dir": output_path, "name": "sample_sp-SVC"},
                "execution": {"seed": 42},
            },
            sort_keys=False,
        ),
        encoding="utf-8",
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
            cwd=ROOT,
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
    assert "revise/revise.yaml" in names
    assert {
        "revise/application/templates/Xenium_T.yaml",
        "revise/application/templates/Xenium_Fib.yaml",
        "revise/application/templates/Xenium_Mono.yaml",
        "revise/application/templates/VisiumHD.yaml",
        "revise/application/templates/Visium.yaml",
    } <= set(names)
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
    expected_extras = {
        "tacco": {"tacco"},
        "pathway": {"gseapy", "networkx", "omicverse"},
        "cci": {"cellphonedb", "omicverse"},
        "trajectory": {"omicverse"},
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


def test_built_wheel_installs_console_help_and_version(installed_cli):
    command = installed_cli["command"]
    env = installed_cli["env"]
    cwd = installed_cli["root"]

    help_result = _run([command, "--help"], cwd=cwd, env=env)
    assert help_result.returncode == 0, help_result.stderr
    assert "--config CONFIG" in help_result.stdout
    assert "--dry-run" in help_result.stdout
    assert "Xenium_T.yaml" in help_result.stdout
    assert "Xenium_Fib.yaml" in help_result.stdout
    assert "Xenium_Mono.yaml" in help_result.stdout
    assert "VisiumHD.yaml" in help_result.stdout
    assert "Visium.yaml" in help_result.stdout
    assert "--svc-type" not in help_result.stdout
    assert "--set" not in help_result.stdout

    removed_set = _run(
        [
            str(command),
            "--config",
            "run.yaml",
            "--set",
            "graph.method=pca",
        ],
        cwd=installed_cli["root"],
        env=env,
    )
    assert removed_set.returncode == 2
    assert "internally managed" in removed_set.stderr

    version = _run(
        [str(command), "--version"],
        cwd=cwd,
        env=env,
    )

    assert version.returncode == 0, version.stderr
    assert version.stdout.strip() == f"revise-reconstruct {__version__}"


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


def test_installed_cli_preflight_runs_outside_checkout(installed_cli):
    root = installed_cli["root"] / "visium-hd-preflight"
    root.mkdir()
    _write_inputs(
        root / "raw_data/Real_application/P1CRC_HD.h5ad",
        root / "raw_data/Real_application/adata_sc_all_reanno.h5ad",
    )
    result = _run(
        [
            installed_cli["command"],
            "--config",
            "configs/application/VisiumHD.yaml",
            "--dry-run",
        ],
        cwd=root,
        env=installed_cli["env"],
    )

    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["status"] == "preflight_passed"
    assert payload["pipeline"]["profile"] == "application_sp"
    assert payload["pipeline"]["route"]["application_route"] == "sp-SVC"
    assert "confounding" not in payload["pipeline"]["route"]
    provenance = json.loads(
        (Path(payload["preflight"]).parent / "provenance.json").read_text()
    )
    assert provenance["application_config"]["source_path"] == (
        "package:revise.application.templates/VisiumHD.yaml"
    )
    assert not list((root / "output").rglob("*.h5ad"))


def test_source_and_installed_preflights_match(installed_cli):
    source_python = installed_cli["source_python"]
    runs = []
    for name, prefix in (
        (
            "source",
            [source_python, ROOT / "reconstruct.py"],
        ),
        ("installed", [installed_cli["command"]]),
    ):
        root = installed_cli["root"] / f"{name}-preflight"
        root.mkdir()
        _write_inputs(root / "data/sample_st.h5ad", root / "data/sc.h5ad")
        config_path = root / "application.yaml"
        _write_application_config(
            config_path,
            st_path="data/sample_st.h5ad",
            reference_path="data/sc.h5ad",
            output_path="output",
        )
        result = _run(
            [*prefix, "--config", config_path, "--dry-run"],
            cwd=root,
            env=installed_cli["env"],
        )
        assert result.returncode == 0, result.stderr
        payload = json.loads(result.stdout)
        preflight_path = Path(payload["preflight"])
        assert preflight_path.is_file()
        manifest_path = preflight_path.with_name("provenance.json")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        assert manifest["run"]["status"] == "succeeded"
        assert manifest["application_config"]["source_path"] == str(
            config_path.resolve()
        )
        assert manifest["application_config"]["resolved_root"] == str(root.resolve())
        assert manifest["application_config"]["resolved_inputs"]["output_dir"] == str(
            (root / "output").resolve()
        )
        run_dir = preflight_path.parent
        assert manifest["application_config"]["output_paths"]["svc"] == str(
            (run_dir / "sample_sp-SVC.h5ad").resolve()
        )
        assert payload["outputs"]["svc"] == str(
            (run_dir / "sample_sp-SVC.h5ad").resolve()
        )
        assert not list((root / "output").rglob("*.h5ad"))
        runs.append((payload, manifest))

    assert runs[0][0]["status"] == runs[1][0]["status"] == "preflight_passed"
    assert runs[0][0]["pipeline"]["profile"] == runs[1][0]["pipeline"]["profile"]
    assert runs[0][0]["pipeline"]["route"] == runs[1][0]["pipeline"]["route"]


def test_installed_xenium_template_requires_tacco_extra_without_fallback(installed_cli):
    root = installed_cli["root"] / "xenium-tacco-gate"
    root.mkdir()
    config_path = root / "Xenium_T.yaml"
    _copy_installed_template(installed_cli, "Xenium_T.yaml", config_path)
    document = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    document["local_refinement"]["select_cell_type"] = "A"
    config_path.write_text(yaml.safe_dump(document, sort_keys=False), encoding="utf-8")
    _write_inputs(
        root / "raw_data/Real_application/P1CRC_Xenium.h5ad",
        root / "raw_data/Real_application/adata_sc_all_reanno.h5ad",
    )

    blocker = root / "dependency-blocker"
    blocker.mkdir()
    (blocker / "sitecustomize.py").write_text(
        "import importlib.abc, sys\n"
        "class BlockTacco(importlib.abc.MetaPathFinder):\n"
        "    def find_spec(self, fullname, path=None, target=None):\n"
        "        if fullname == 'tacco' or fullname.startswith('tacco.'):\n"
        "            raise ModuleNotFoundError('blocked optional dependency', name='tacco')\n"
        "sys.meta_path.insert(0, BlockTacco())\n",
        encoding="utf-8",
    )
    env = installed_cli["env"].copy()
    env["PYTHONPATH"] = str(blocker)
    result = _run(
        [installed_cli["command"], "--config", config_path, "--dry-run"],
        cwd=root,
        env=env,
    )

    assert result.returncode != 0
    assert 'pip install "revise-svc[tacco]"' in result.stderr
    assert "REVISE does not fall back automatically" in result.stderr
    assert document["algorithm"] == {"ot_method": "tacco"}
