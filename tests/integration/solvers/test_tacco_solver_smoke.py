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
from revise.backend.kernels.local_anchoring import LocalAnchoringKernel
from revise.backend.ops.local_ot import solve_local_ot
from revise.backend.ops.tacco_runtime import require_tacco
from revise.config.runner_conf import ApplicationScConf

require_tacco()
assert importlib.metadata.version("tacco") == "0.5.0"

target = AnnData(
    X=np.array(
        [
            [8.0, 1.0, 1.0],
            [7.0, 2.0, 1.0],
            [1.0, 8.0, 1.0],
            [2.0, 7.0, 1.0],
        ]
    ),
    obs=pd.DataFrame(index=["spot1", "spot2", "spot3", "spot4"]),
    var=pd.DataFrame(index=["g1", "g2", "g3"]),
)
reference = AnnData(
    X=np.array(
        [
            [9.0, 1.0, 1.0],
            [8.0, 2.0, 1.0],
            [1.0, 9.0, 1.0],
            [2.0, 8.0, 1.0],
            [7.0, 3.0, 1.0],
            [3.0, 7.0, 1.0],
        ]
    ),
    obs=pd.DataFrame(
        {{
            "Level1": ["A", "A", "B", "B", "A", "B"],
            "Level2": ["A1", "A2", "B1", "B2", "A1", "B1"],
        }},
        index=["cell1", "cell2", "cell3", "cell4", "cell5", "cell6"],
    ),
    var=pd.DataFrame(index=["g1", "g2", "g3"]),
)

config = ApplicationScConf(
    sample_name="sample",
    raw_data_path="data",
    result_root_path="output",
    st_file="sp.h5ad",
    sc_ref_file="sc.h5ad",
    annotate_mode="tacco",
    rec_ot_method="tacco",
    cell_type_col="Level1",
    confidence_col="Confidence",
    unknown_key="Unknown",
    tacco_annotate_multi_center=1,
    tacco_annotate_lamb=0.001,
)
annotated = GlobalAnchoringKernel(config, logging.getLogger("tacco-smoke")).run(
    target, reference
)
local = LocalAnchoringKernel(config, logging.getLogger("tacco-smoke"))
annotated_level2 = local.run(annotated, reference, cell_type_col="Level2")
annotated_level2.obs["SVC_cluster"] = pd.Categorical(["0", "0", "1", "1"])
annotated_reference = local.run(
    reference,
    annotated_level2,
    cell_type_col="SVC_cluster",
)

coupling = solve_local_ot(
    [0.5, 0.5],
    [0.5, 0.5],
    [[0.0, 1.0], [1.0, 0.0]],
    method="tacco",
)

@numba.njit
def add_one(value):
    return value + 1

assert add_one(1) == 2
assert add_one.signatures
assert annotated.obsm["Level1"].shape == (4, 2)
assert annotated.obsm["Level1"].index.equals(annotated.obs_names)
assert annotated.obsm["Level1"].columns.tolist() == ["A", "B"]
np.testing.assert_allclose(
    annotated.obsm["Level1"].sum(axis=1).to_numpy(),
    np.ones(annotated.n_obs),
    rtol=0,
    atol=1e-6,
)
assert annotated.obs["Level1"].tolist() == annotated.obsm["Level1"].idxmax(axis=1).tolist()
assert annotated_level2.obsm["Level2"].shape == (4, 4)
assert annotated_reference.obsm["SVC_cluster"].shape == (6, 2)
assert coupling.shape == (2, 2)
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
