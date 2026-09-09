from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10
    import tomli as tomllib


ROOT = Path(__file__).resolve().parents[2]


def _project():
    with (ROOT / "pyproject.toml").open("rb") as handle:
        return tomllib.load(handle)["project"]


def test_base_metadata_keeps_reconstruction_essentials_and_excludes_optional_domains():
    project = _project()
    dependencies = project["dependencies"]
    normalized = [item.lower() for item in dependencies]

    assert project["requires-python"] == ">=3.10,<3.12"
    assert any(item.startswith("pot") for item in normalized)
    assert any(item.startswith("leidenalg") for item in normalized)
    assert not any(item.startswith("tacco") for item in normalized)
    assert not any(item.startswith("gseapy") for item in normalized)
    assert not any(item.startswith("networkx") for item in normalized)
    assert not any(item.startswith("omicverse") for item in normalized)
    assert not any(item.startswith("scvelo") for item in normalized)


def test_optional_extras_match_supported_domain_boundaries():
    extras = _project()["optional-dependencies"]

    assert set(extras) >= {
        "tacco",
        "pathway",
        "cci",
        "trajectory",
        "spatialdata",
        "dev",
    }
    assert extras["tacco"] == ["tacco==0.5.0"]
    assert any(item.lower().startswith("gseapy") for item in extras["pathway"])
    assert any(item.lower().startswith("omicverse") for item in extras["pathway"])
    for extra in ("pathway", "cci", "trajectory"):
        assert "omicverse==1.7.5" in extras[extra]
        assert "torch-geometric" in extras[extra]
        assert "torch" in extras[extra]
        assert "setuptools<81" in extras[extra]
        assert "transformers<5" in extras[extra]
    assert "cellphonedb>=5,<6" in extras["cci"]
    assert extras["spatialdata"] == ["spatialdata>=0.2"]


def test_base_import_does_not_import_optional_domain_packages():
    code = """
import sys
import revise

blocked = (
    'tacco', 'gseapy', 'networkx', 'omicverse', 'cellphonedb',
    'torch_geometric', 'scvelo', 'spatialdata'
)
loaded = [name for name in blocked if name in sys.modules]
assert loaded == [], loaded
"""
    env = os.environ.copy()
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr


def test_pathway_missing_dependency_names_the_extra(monkeypatch):
    import builtins
    import importlib
    import types

    scanpy = types.ModuleType("scanpy")
    monkeypatch.setitem(sys.modules, "scanpy", scanpy)
    analysis = types.ModuleType("revise.analysis")
    analysis.__path__ = [str(ROOT / "revise" / "analysis")]
    monkeypatch.setitem(sys.modules, "revise.analysis", analysis)
    monkeypatch.delitem(sys.modules, "revise.analysis.bio", raising=False)
    bio = importlib.import_module("revise.analysis.bio")

    real_import = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name == "gseapy":
            raise ModuleNotFoundError("missing", name="gseapy")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    try:
        bio._require_gseapy()
    except ImportError as exc:
        assert "revise-svc[pathway]" in str(exc)
    else:
        raise AssertionError("missing pathway dependency did not fail")
