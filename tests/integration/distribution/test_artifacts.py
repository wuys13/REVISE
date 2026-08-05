from __future__ import annotations

import hashlib
import base64
import csv
import io
import os
from pathlib import Path
import subprocess
import sys
import tarfile
from zipfile import ZipFile

import pytest


ROOT = Path(__file__).resolve().parents[3]


def _run(command, *, cwd, env):
    return subprocess.run(
        [str(item) for item in command],
        cwd=cwd,
        env=env,
        text=True,
        capture_output=True,
        check=False,
        timeout=240,
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


@pytest.fixture(scope="module")
def built_distributions(tmp_path_factory):
    root = tmp_path_factory.mktemp("distribution-artifacts")
    supplied_dist = os.environ.get("REVISE_DIST_DIR")
    dist = Path(supplied_dist).resolve() if supplied_dist else root / "dist"
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env["PYTHONNOUSERSITE"] = "1"
    if not supplied_dist:
        command = [
            os.environ.get("REVISE_INTEGRATION_PYTHON", sys.executable),
            "-m",
            "build",
            "--outdir",
            dist,
            ROOT,
        ]
        result = _run(command, cwd=root, env=env)
        assert result.returncode == 0, result.stderr
    wheels = list(dist.glob("*.whl"))
    sdists = list(dist.glob("*.tar.gz"))
    assert len(wheels) == 1
    assert len(sdists) == 1
    artifacts = {"wheel": wheels[0], "sdist": sdists[0]}
    return root, artifacts, env


def test_distribution_contents_match_runtime_and_source_contract(built_distributions):
    _, artifacts, _ = built_distributions
    with ZipFile(artifacts["wheel"]) as archive:
        wheel_names = archive.namelist()
        entry_points = archive.read(
            next(name for name in wheel_names if name.endswith("entry_points.txt"))
        ).decode()
        metadata = archive.read(
            next(name for name in wheel_names if name.endswith("/METADATA"))
        ).decode()
        record_name = next(name for name in wheel_names if name.endswith("/RECORD"))
        records = {
            row[0]: row[1:]
            for row in csv.reader(io.StringIO(archive.read(record_name).decode()))
        }
        assert set(records) == set(wheel_names)
        for name in wheel_names:
            if name == record_name:
                assert records[name] == ["", ""]
                continue
            encoded_hash, size = records[name]
            digest = base64.urlsafe_b64encode(
                hashlib.sha256(archive.read(name)).digest()
            ).rstrip(b"=").decode()
            assert encoded_hash == f"sha256={digest}"
            assert int(size) == len(archive.read(name))
    with tarfile.open(artifacts["sdist"], "r:gz") as archive:
        members = {
            member.name.split("/", 1)[-1]: member
            for member in archive.getmembers()
            if member.isfile()
        }
        sdist_names = list(members)

    assert "revise/revise.yaml" in wheel_names
    packaged_templates = {
        f"revise/application/templates/{filename}"
        for filename in (
            "Xenium_T.yaml",
            "Xenium_Fib.yaml",
            "Xenium_Mono.yaml",
            "VisiumHD.yaml",
            "Visium.yaml",
        )
    }
    assert packaged_templates <= set(wheel_names)
    assert "reconstruct.py" not in wheel_names
    assert "revise-reconstruct = revise.application.cli:main" in entry_points
    assert (
        "revise-build-histology-priors = revise.preprocess.cli:main"
        in entry_points
    )
    assert (
        "revise-compute-biological-metrics = revise.analysis.cli:main"
        in entry_points
    )
    assert "Version: 0.1.0rc1" in metadata
    assert "Requires-Python: <3.12,>=3.10" in metadata
    assert any(name.endswith("/LICENSE") for name in wheel_names)
    assert not any(name.startswith("tests/") for name in wheel_names)
    assert not any(name.startswith("constraints/") for name in wheel_names)
    assert not any(name.startswith("release/") for name in wheel_names)

    expected_sdist = {
        "LICENSE",
        "README.md",
        "pyproject.toml",
        "MANIFEST.in",
        "reconstruct.py",
        "reproduce/benchmark_main.py",
        "reproduce/benchmark_main.sh",
    }
    for pattern in [
        "revise/**/*.py",
        "revise/**/*.yaml",
        "configs/application/*.yaml",
        "constraints/*.txt",
    ]:
        expected_sdist.update(
            path.relative_to(ROOT).as_posix() for path in ROOT.glob(pattern)
        )
    assert expected_sdist <= set(sdist_names)
    generated_sdist = {
        "PKG-INFO",
        "setup.cfg",
    }
    unexpected_sdist = set(sdist_names) - expected_sdist
    assert all(
        name in generated_sdist or name.startswith("revise_svc.egg-info/")
        for name in unexpected_sdist
    )
    assert not any(name.endswith("release-manifest.json") for name in sdist_names)
    assert not any(name.startswith("tests/") for name in sdist_names)
    with tarfile.open(artifacts["sdist"], "r:gz") as archive:
        for source_path in expected_sdist:
            archived = archive.extractfile(members[source_path])
            assert archived is not None
            assert archived.read() == (ROOT / source_path).read_bytes()

    expected_wheel_sources = {
        path.relative_to(ROOT).as_posix()
        for pattern in ["revise/**/*.py", "revise/**/*.yaml"]
        for path in ROOT.glob(pattern)
    }
    with ZipFile(artifacts["wheel"]) as archive:
        actual_wheel_sources = {
            name for name in wheel_names if name.startswith("revise/")
        }
        assert actual_wheel_sources == expected_wheel_sources
        for source_path in expected_wheel_sources:
            assert archive.read(source_path) == (ROOT / source_path).read_bytes()
        for package_path in packaged_templates:
            filename = Path(package_path).name
            assert archive.read(package_path) == (
                ROOT / "configs" / "application" / filename
            ).read_bytes()
        license_name = next(
            name for name in wheel_names if name.endswith("/licenses/LICENSE")
        )
        assert archive.read(license_name) == (ROOT / "LICENSE").read_bytes()

    with tarfile.open(artifacts["sdist"], "r:gz") as archive:
        for package_path in packaged_templates:
            filename = Path(package_path).name
            packaged = archive.extractfile(members[package_path])
            mirrored = archive.extractfile(members[f"configs/application/{filename}"])
            assert packaged is not None
            assert mirrored is not None
            assert packaged.read() == mirrored.read()


@pytest.mark.parametrize("role", ["wheel", "sdist"])
def test_each_distribution_installs_outside_checkout(built_distributions, role):
    root, artifacts, env = built_distributions
    artifact_hashes = {
        artifact_role: _sha256(path) for artifact_role, path in artifacts.items()
    }
    venv = root / f"{role}-venv"
    source_python = os.environ.get("REVISE_INTEGRATION_PYTHON", sys.executable)
    clean_install = os.environ.get("REVISE_CLEAN_INSTALL") == "1"
    venv_args = [source_python, "-m", "venv"]
    if not clean_install:
        venv_args.extend(["--without-pip", "--system-site-packages"])
    venv_args.append(venv)
    create = _run(
        venv_args,
        cwd=root,
        env=env,
    )
    assert create.returncode == 0, create.stderr
    python = venv / "bin" / "python"
    if clean_install:
        version_probe = _run(
            [python, "-c", "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')"],
            cwd=root,
            env=env,
        )
        assert version_probe.returncode == 0, version_probe.stderr
        python_minor = version_probe.stdout.strip()
        constraint = ROOT / "constraints" / f"python-{python_minor}.txt"
        toolchain = _run(
            [
                python,
                "-m",
                "pip",
                "install",
                "-c",
                constraint,
                "setuptools==80.9.0",
                "wheel==0.45.1",
            ],
            cwd=root,
            env=env,
        )
        assert toolchain.returncode == 0, toolchain.stderr
    install = _run(
        [
            python,
            "-m",
            "pip",
            "install",
            *([] if clean_install else ["--no-deps"]),
            *(
                [
                    "-c",
                    ROOT / "constraints" / f"python-{python_minor}.txt",
                ]
                if clean_install
                else []
            ),
            "--no-build-isolation",
            artifacts[role],
        ],
        cwd=root,
        env=env,
    )
    assert install.returncode == 0, install.stderr
    probe = _run(
        [
            python,
            "-c",
            "import importlib.metadata as m, revise; "
            "print(m.version('revise-svc')); print(revise.__file__)",
        ],
        cwd=root,
        env=env,
    )
    assert probe.returncode == 0, probe.stderr
    assert "0.1.0rc1" in probe.stdout
    assert str(venv) in probe.stdout
    assert str(ROOT) not in probe.stdout
    version = _run(
        [venv / "bin" / "revise-reconstruct", "--version"],
        cwd=root,
        env=env,
    )
    assert version.returncode == 0, version.stderr
    assert "0.1.0rc1" in version.stdout
    histology_help = _run(
        [venv / "bin" / "revise-build-histology-priors", "--help"],
        cwd=root,
        env=env,
    )
    assert histology_help.returncode == 0, histology_help.stderr
    metrics_help = _run(
        [venv / "bin" / "revise-compute-biological-metrics", "--help"],
        cwd=root,
        env=env,
    )
    assert metrics_help.returncode == 0, metrics_help.stderr
    assert "--st-h5ad" in histology_help.stdout
    for artifact_role, artifact in artifacts.items():
        assert _sha256(artifact) == artifact_hashes[artifact_role]
