from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

import pytest


ROOT = Path(__file__).resolve().parents[2]
HOST_ENV_KEYS = ("NUMBA_DISABLE_JIT", "NUMBA_CACHE_DIR", "MPLCONFIGDIR")
HOST_ENV_VALUES = {
    "NUMBA_DISABLE_JIT": "0",
    "NUMBA_CACHE_DIR": str(Path(tempfile.gettempdir()) / "revise-u6-numba"),
    "MPLCONFIGDIR": str(Path(tempfile.gettempdir()) / "revise-u6-mpl"),
}


def _run_python(code: str, *, env: dict[str, str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-c", code],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )


@pytest.mark.parametrize("preexisting", [False, True])
def test_import_revise_does_not_mutate_host_environment(preexisting):
    env = os.environ.copy()
    expected = {}
    for key in HOST_ENV_KEYS:
        if preexisting:
            expected[key] = env[key] = HOST_ENV_VALUES[key]
        else:
            env.pop(key, None)

    code = f"""
import os
import revise

keys = {HOST_ENV_KEYS!r}
expected = {expected!r}
actual = {{key: os.environ[key] for key in keys if key in os.environ}}
assert actual == expected, (expected, actual)
"""
    result = _run_python(code, env=env)

    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("preexisting", [False, True])
def test_set_global_seed_preserves_pythonhashseed_and_seeds_rngs(preexisting):
    env = os.environ.copy()
    if preexisting:
        env["PYTHONHASHSEED"] = "98765"
    else:
        env.pop("PYTHONHASHSEED", None)
    expected_hash_seed = env.get("PYTHONHASHSEED")

    code = f"""
import os
import random
import numpy as np
from revise.utils import set_global_seed

set_global_seed(123)
observed = (random.random(), np.random.random())
random.seed(123)
np.random.seed(123)
expected = (random.random(), np.random.random())
assert observed == expected, (expected, observed)
assert os.environ.get("PYTHONHASHSEED") == {expected_hash_seed!r}
"""
    result = _run_python(code, env=env)

    assert result.returncode == 0, result.stderr


def test_import_revise_utils_exposes_no_gseapy_stub_or_fake_module():
    code = """
import sys
import revise.utils

assert not hasattr(revise.utils, "ensure_gseapy_stub")
assert "gseapy" not in sys.modules, sys.modules["gseapy"]
"""
    env = os.environ.copy()
    env["MPLCONFIGDIR"] = str(Path(tempfile.gettempdir()) / "revise-u6-mpl")
    result = _run_python(code, env=env)

    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("operation", ["enrichment", "network"])
@pytest.mark.parametrize("failure", ["missing", "broken"])
def test_nonempty_gseapy_features_fail_closed_with_install_guidance(
    operation, failure
):
    code = f"""
import builtins
from pathlib import Path
import sys
import types

import pandas as pd
import revise

analysis_package = types.ModuleType("revise.analysis")
analysis_package.__path__ = [str(Path.cwd() / "revise" / "analysis")]
sys.modules["revise.analysis"] = analysis_package
sys.modules["scanpy"] = types.ModuleType("scanpy")

from revise.analysis.bio import get_enrichment, pathway_network_plot

real_import = builtins.__import__
failure = {failure!r}

def blocked_import(name, *args, **kwargs):
    if name == "gseapy" or name.startswith("gseapy."):
        if failure == "missing":
            raise ModuleNotFoundError("No module named 'gseapy'", name="gseapy")
        raise ImportError("broken gseapy installation")
    return real_import(name, *args, **kwargs)

builtins.__import__ = blocked_import
try:
    if {operation!r} == "enrichment":
        get_enrichment(["G1"], "gene_sets.gmt")
    else:
        pathway_network_plot(pd.DataFrame({{"Term": ["T1"]}}))
except Exception as exc:
    assert "revise-svc[pathway]" in str(exc).lower(), str(exc)
else:
    raise AssertionError("nonempty gseapy feature did not fail closed")
"""
    env = os.environ.copy()
    env["MPLCONFIGDIR"] = str(Path(tempfile.gettempdir()) / "revise-u6-mpl")
    result = _run_python(code, env=env)

    assert result.returncode == 0, result.stderr


def test_lazy_gseapy_features_use_an_available_module():
    code = """
from pathlib import Path
from types import ModuleType, SimpleNamespace
import sys

import pandas as pd
import revise

analysis_package = ModuleType("revise.analysis")
analysis_package.__path__ = [str(Path.cwd() / "revise" / "analysis")]
sys.modules["revise.analysis"] = analysis_package
sys.modules["scanpy"] = ModuleType("scanpy")

from revise.analysis import bio

calls = []
gseapy = ModuleType("gseapy")

def enrichr(**kwargs):
    calls.append(("enrichr", kwargs))
    return SimpleNamespace(results=pd.DataFrame({"Term": ["T1"]}))

def enrichment_map(pathway, top_term):
    calls.append(("enrichment_map", top_term))
    nodes = pd.DataFrame(
        {"Genes": ["G1,G2"], "p_inv": [1.0], "Term": ["T1"]}
    )
    edges = pd.DataFrame(
        {
            "src_idx": [0],
            "targ_idx": [0],
            "jaccard_coef": [0.5],
            "overlap_coef": [0.5],
            "overlap_genes": ["G1"],
        }
    )
    return nodes, edges

gseapy.enrichr = enrichr
gseapy.enrichment_map = enrichment_map
sys.modules["gseapy"] = gseapy
bio.plt.show = lambda: None

result = bio.get_enrichment(["G1"], "gene_sets.gmt")
bio.pathway_network_plot(pd.DataFrame({"Term": ["T1"]}), top_term=3)

assert result["Term"].tolist() == ["T1"]
assert calls[0][0] == "enrichr"
assert calls[0][1]["gene_list"] == ["G1"]
assert calls[1] == ("enrichment_map", 3)
"""
    env = os.environ.copy()
    env["MPLCONFIGDIR"] = str(Path(tempfile.gettempdir()) / "revise-u6-mpl")
    result = _run_python(code, env=env)

    assert result.returncode == 0, result.stderr


def test_empty_gseapy_features_do_not_import_gseapy():
    code = """
import builtins
from pathlib import Path
import sys
import types

import pandas as pd
import revise

analysis_package = types.ModuleType("revise.analysis")
analysis_package.__path__ = [str(Path.cwd() / "revise" / "analysis")]
sys.modules["revise.analysis"] = analysis_package
sys.modules["scanpy"] = types.ModuleType("scanpy")

from revise.analysis.bio import get_enrichment, pathway_network_plot

real_import = builtins.__import__

def blocked_import(name, *args, **kwargs):
    if name == "gseapy" or name.startswith("gseapy."):
        raise AssertionError("empty input attempted to import gseapy")
    return real_import(name, *args, **kwargs)

builtins.__import__ = blocked_import
assert get_enrichment([], "gene_sets.gmt").empty
assert pathway_network_plot(pd.DataFrame()) is None
assert "gseapy" not in sys.modules
"""
    env = os.environ.copy()
    env["MPLCONFIGDIR"] = str(Path(tempfile.gettempdir()) / "revise-u6-mpl")
    result = _run_python(code, env=env)

    assert result.returncode == 0, result.stderr
