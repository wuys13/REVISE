from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
RTOL = 1e-10
ATOL = 1e-12

PROBE = r"""
import json
import logging
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd

from revise.backend.kernels.spot_sr import SpotSrKernel
from revise.backend.kernels.ot import OTKernel
from revise.utils.deterministic import canonical_config_projection
from revise.utils.provenance import hash_jsonable

seed = int(sys.argv[1])
input_root = Path(sys.argv[2])
output_root = sys.argv[3]
input_root.mkdir(parents=True)
st_path = input_root / "st.bin"
sc_path = input_root / "sc.bin"
st_path.write_bytes(b"st-content")
sc_path.write_bytes(b"sc-content")

config = {
    "runtime": {
        "seed": seed,
        "deterministic": True,
        "mode": "application",
        "application_route": "sc-SVC",
        "application_mode": "sr",
        "task": "sc_svc_super_resolution",
        "svc_kind": "sc",
        "strategy": "ScSvcSuperResolutionApplicationStrategy",
    },
    "io": {
        "data_root": str(input_root),
        "output_root": output_root,
        "sample_name": "sample",
        "st_file": st_path.name,
        "sc_ref_file": sc_path.name,
        "gt_svc_file": "unused.h5ad",
        "spatialdata_path": None,
        "spot_size": 50,
        "save_outputs": True,
    },
    "graph": {"alpha": 0.2},
    "ot": {
        "ga": {"solver": "pot", "pot": {"reg": 0.1}},
        "lr": {"solver": "pot", "pot": {"reg": 0.1}},
    },
}
svc_obs = pd.DataFrame(
    [
        {"spot_name": spot, "cell_id": f"{spot}-{index}"}
        for spot in ("spot-a", "spot-b")
        for index in range(6)
    ]
)
quota = pd.DataFrame(
    [[1, 2, 3], [2, 2, 2]],
    index=["spot-a", "spot-b"],
    columns=["A", "B", "C"],
)
kernel = SpotSrKernel(
    SimpleNamespace(
        pm_on_cell_file=str(input_root / "missing-pm.csv"),
        svc_completeness=True,
        sr_assignment_seed=seed,
    ),
    logging.getLogger("cross-process-probe"),
)
assigned = kernel.assign_cell_types_random(svc_obs, quota)
ordered = assigned.sort_values(["spot_name", "cell_id"], kind="stable")

coupling = OTKernel.couple(
    [0.4, 0.6],
    [0.5, 0.5],
    [[0.0, 1.0], [1.0, 0.0]],
    method="pot",
    pot_reg=0.1,
    pot_reg_m=0.0,
)
observed_quota = {
    spot: {
        cell_type: int(count)
        for cell_type, count in group["cell_type"].value_counts().sort_index().items()
    }
    for spot, group in assigned.groupby("spot_name", sort=True)
}
print(json.dumps({
    "input_paths": [str(st_path), str(sc_path)],
    "ordered_cell_ids": ordered["cell_id"].tolist(),
    "ordered_labels": ordered["cell_type"].tolist(),
    "quota_counts": observed_quota,
    "coupling": np.asarray(coupling).tolist(),
    "config_hash": hash_jsonable(canonical_config_projection(config)),
}))
"""


def _probe(tmp_path: Path, *, seed: int, hash_seed: int) -> dict:
    name = f"seed-{seed}-hash-{hash_seed}"
    env = os.environ.copy()
    env["PYTHONHASHSEED"] = str(hash_seed)
    result = subprocess.run(
        [
            os.environ.get("REVISE_INTEGRATION_PYTHON", sys.executable),
            "-c",
            PROBE,
            str(seed),
            str(tmp_path / name / "inputs"),
            str(tmp_path / name / "outputs"),
        ],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def test_fresh_processes_match_normalized_semantic_output(tmp_path):
    first = _probe(tmp_path, seed=17, hash_seed=1)
    second = _probe(tmp_path, seed=17, hash_seed=987654)

    assert first["input_paths"] != second["input_paths"]
    assert first["ordered_cell_ids"] == second["ordered_cell_ids"]
    assert first["ordered_labels"] == second["ordered_labels"]
    assert first["quota_counts"] == second["quota_counts"]
    first_coupling = np.asarray(first["coupling"], dtype=float)
    second_coupling = np.asarray(second["coupling"], dtype=float)
    assert first_coupling.shape == second_coupling.shape == (2, 2)
    assert np.isfinite(first_coupling).all()
    np.testing.assert_allclose(
        first_coupling,
        second_coupling,
        rtol=RTOL,
        atol=ATOL,
    )
    assert first["config_hash"] == second["config_hash"]


def test_fresh_process_seed_change_changes_positions_and_preserves_quota(tmp_path):
    first = _probe(tmp_path, seed=17, hash_seed=1)
    second = _probe(tmp_path, seed=18, hash_seed=1)

    assert first["ordered_labels"] != second["ordered_labels"]
    assert first["quota_counts"] == second["quota_counts"] == {
        "spot-a": {"A": 1, "B": 2, "C": 3},
        "spot-b": {"A": 2, "B": 2, "C": 2},
    }
    assert first["config_hash"] != second["config_hash"]
