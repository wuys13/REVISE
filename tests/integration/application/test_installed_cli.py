from __future__ import annotations

import hashlib
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
from anndata import read_h5ad
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


def _write_inputs(data_root: Path) -> None:
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
            {
                "Level1": ["A", "B"] * 26,
                "Level2": ["A1", "B1"] * 26,
            },
            index=[f"cell-{index}" for index in range(52)],
        ),
        var=pd.DataFrame(index=st.var_names.copy()),
    )
    st.write_h5ad(data_root / "sample_st.h5ad")
    sc_ref.write_h5ad(data_root / "sc.h5ad")


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
    assert str(venv) in probe.stdout.strip()
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

    assert __version__ == "2.0.0rc1"
    assert "Version: 2.0.0rc1" in metadata
    assert "Name: revise-svc" in metadata
    assert "Requires-Python: <3.12,>=3.10" in metadata
    assert "revise/application/cli.py" in names
    assert "revise/benchmark/cli.py" in names
    assert "revise/benchmark/launcher.py" in names
    assert "revise/revise.yaml" in names
    assert any(".dist-info/" in name and name.endswith("/LICENSE") for name in names)
    assert not any(name.startswith("tests/") for name in names)
    assert "revise-reconstruct = revise.application.cli:main" in entry_points

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
    assert "--ot-method" in help_result.stdout
    assert "--dry-run" in help_result.stdout
    assert "--svc-type {hST-SVC,iST-SVC,sST-SVC}" in help_result.stdout
    assert "--ist-mapping {mean,random}" in help_result.stdout
    assert "--set" not in help_result.stdout
    assert "--select-ct" not in help_result.stdout
    assert "--sc-mapping" not in help_result.stdout

    for removed_option in (
        ["--set", "graph.method=pca"],
        ["--select-ct", "T"],
        ["--sc-mapping", "mean"],
    ):
        removed = _run(
            [
                str(command),
                "--svc-type",
                "iST-SVC",
                "--sample-name",
                "sample",
                "--st-file",
                "st.h5ad",
                "--sc-ref-file",
                "sc.h5ad",
                "--data-root",
                "data",
                *removed_option,
            ],
            cwd=installed_cli["root"],
            env=env,
        )
        assert removed.returncode == 2
        assert f"unrecognized arguments: {removed_option[0]}" in removed.stderr

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


@pytest.mark.parametrize(
    ("svc_type", "profile", "route", "mapping_args"),
    (
        ("hST-SVC", "application_sp", "sp_svc:bin2cell", []),
        (
            "iST-SVC",
            "application_sc",
            "sc_svc:segmentation",
            ["--ist-mapping", "mean"],
        ),
        ("sST-SVC", "application_sc_sr", "sc_svc_sr:spot_size", []),
    ),
)
def test_installed_cli_preflight_runs_outside_checkout(
    installed_cli,
    svc_type,
    profile,
    route,
    mapping_args,
):
    root = installed_cli["root"]
    data_root = root / f"dry-data-{svc_type}"
    data_root.mkdir()
    _write_inputs(data_root)
    output_root = root / f"dry-output-{svc_type}"
    result = _run(
        [
            installed_cli["command"],
            "--svc-type",
            svc_type,
            "--sample-name",
            "sample",
            "--st-file",
            "st.h5ad",
            "--sc-ref-file",
            "sc.h5ad",
            "--data-root",
            data_root,
            "--output-root",
            output_root,
            "--ot-method",
            "pot",
            *mapping_args,
            "--dry-run",
        ],
        cwd=root,
        env=installed_cli["env"],
    )

    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["status"] == "ready"
    assert payload["svc_type"] == svc_type
    assert payload["pipeline"]["profile"] == profile
    assert payload["pipeline"]["route"] == route
    assert not list(output_root.rglob("*.h5ad"))


def test_source_and_installed_minimal_pot_runs_match(installed_cli):
    source_python = installed_cli["source_python"]
    availability = _run(
        [source_python, "-c", "import scanpy, ot, squidpy"],
        cwd=installed_cli["root"],
        env=installed_cli["env"],
    )
    assert availability.returncode == 0, (
        "release integration interpreter lacks mandatory base dependencies: "
        f"{availability.stderr}"
    )

    data_root = installed_cli["root"] / "pot-data"
    data_root.mkdir()
    _write_inputs(data_root)
    config = yaml.safe_load((ROOT / "revise/revise.yaml").read_text(encoding="utf-8"))
    config["defaults"]["preprocess"].update(
        st_min_counts=1,
        st_min_cells=1,
        sc_min_counts=1,
        sc_min_cells=1,
    )
    config["defaults"]["graph"].update(
        method="pca",
        n_neighbors=5,
        exp_neighbors=5,
    )
    config["defaults"]["plot"]["enabled"] = False
    config_path = installed_cli["root"] / "minimal-pot.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    common = [
        "--svc-type",
        "hST-SVC",
        "--sample-name",
        "sample",
        "--st-file",
        "st.h5ad",
        "--sc-ref-file",
        "sc.h5ad",
        "--data-root",
        str(data_root),
        "--ot-method",
        "pot",
        "--config",
        str(config_path),
    ]
    runs = []
    for name, prefix in (
        (
            "source",
            [source_python, ROOT / "reconstruct.py"],
        ),
        ("installed", [installed_cli["command"]]),
    ):
        output_root = installed_cli["root"] / f"{name}-output"
        result = _run(
            [*prefix, *common, "--output-root", output_root],
            cwd=installed_cli["root"],
            env=installed_cli["env"],
        )
        assert result.returncode == 0, result.stderr
        payload = json.loads(result.stdout)
        public = output_root / "sample" / "hST-SVC" / "SVC.h5ad"
        assert public.is_file()
        assert Path(payload["output"]) == public
        manifest_path = next(output_root.rglob("provenance.json"))
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        assert manifest["run"]["status"] == "succeeded"
        assert manifest["result"] == {
            "filename": "SVC.h5ad",
            "type": "hST-SVC",
        }
        assert "result_hash" not in manifest
        assert manifest["stages"][3]["status"] == "succeeded"
        # Each A/B fixture group has 26 observations, so the runner's existing
        # <=50 rule skips local OT while still exercising full publication.
        assert manifest["local_refinement"] == {
            "route": "sp_svc:bin2cell",
            "applied": False,
            "strength": 0.2,
        }
        public_artifacts = [
            artifact
            for artifact in manifest["artifacts"]
            if artifact["role"] == "public_result"
            and artifact["status"] == "completed"
        ]
        assert len(public_artifacts) == 1
        assert Path(public_artifacts[0]["path"]) == public
        assert public_artifacts[0]["sha256"] == hashlib.sha256(
            public.read_bytes()
        ).hexdigest()
        published = read_h5ad(public)
        assert published.uns["revise_reconstruction"] == {
            "schema_version": 2,
            "svc_type": "hST-SVC",
        }
        assert "ot_events" not in manifest
        runs.append((payload, manifest))

    assert runs[0][0]["shape"] == runs[1][0]["shape"] == [52, 52]
    assert runs[0][0]["pipeline"]["profile"] == runs[1][0]["pipeline"]["profile"]
    assert runs[0][0]["pipeline"]["route"] == runs[1][0]["pipeline"]["route"]
