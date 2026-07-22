from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

import pytest


ROOT = Path(__file__).resolve().parents[3]


@pytest.mark.parametrize("import_order", ["revise_then_tacco", "tacco_then_revise"])
def test_real_tacco_050_completes_global_and_local_smoke(import_order):
    code = f"""
import importlib.metadata
import logging
import os
from types import SimpleNamespace

import_order = {import_order!r}
host_keys = ("NUMBA_DISABLE_JIT", "NUMBA_CACHE_DIR", "MPLCONFIGDIR")

def host_snapshot():
    return {{key: os.environ.get(key) for key in host_keys}}

if import_order == "revise_then_tacco":
    before_revise = host_snapshot()
    import revise
    assert host_snapshot() == before_revise
    import tacco
else:
    import tacco
    before_revise = host_snapshot()
    import revise
    assert host_snapshot() == before_revise

import numba
import numpy as np
import pandas as pd
from anndata import AnnData

from revise.backend.kernels.global_anchoring import GlobalAnchoringKernel
from revise.backend.ops.local_ot import solve_local_ot
from revise.backend.ops.tacco_runtime import require_tacco

require_tacco()
assert importlib.metadata.version("tacco") == "0.5.0"

target = AnnData(
    X=np.array([[3.0, 1.0], [1.0, 3.0]]),
    obs=pd.DataFrame(index=["spot1", "spot2"]),
    var=pd.DataFrame(index=["g1", "g2"]),
)
reference = AnnData(
    X=np.array([[4.0, 0.0], [0.0, 4.0], [3.0, 1.0], [1.0, 3.0]]),
    obs=pd.DataFrame(
        {{"Level1": ["A", "B", "A", "B"]}},
        index=["cell1", "cell2", "cell3", "cell4"],
    ),
    var=pd.DataFrame(index=["g1", "g2"]),
)

global_events = []
config = SimpleNamespace(
    annotate_mode="tacco",
    cell_type_col="Level1",
    confidence_col="Confidence",
    unknown_key="Unknown",
    ot_event_callback=lambda *event: global_events.append(event),
)
annotated = GlobalAnchoringKernel(config, logging.getLogger("tacco-smoke")).run(
    target, reference
)

local_events = []
coupling = solve_local_ot(
    [0.5, 0.5],
    [0.5, 0.5],
    [[0.0, 1.0], [1.0, 0.0]],
    method="tacco",
    event_callback=lambda *event: local_events.append(event),
)

@numba.njit
def add_one(value):
    return value + 1

assert add_one(1) == 2
assert add_one.signatures
assert annotated.obsm["Level1"].shape == (2, 2)
assert coupling.shape == (2, 2)
assert global_events == [
    ("ga", "tacco", "attempted"),
    ("ga", "tacco", "completed"),
]
assert local_events == [
    ("lr", "tacco", "attempted"),
    ("lr", "tacco", "completed"),
]
"""
    env = os.environ.copy()
    env.pop("NUMBA_DISABLE_JIT", None)
    env.pop("NUMBA_CACHE_DIR", None)
    env["MPLCONFIGDIR"] = str(Path(tempfile.gettempdir()) / "revise-u6-tacco-mpl")
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
