from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData


ROOT = Path(__file__).resolve().parents[1]


def _write_inputs(data_root: Path) -> None:
    st = AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame(index=["spot-1", "spot-2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    st.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 1.0]])
    sc = AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame(
            {"Level1": ["A", "B"], "Level2": ["A1", "B1"]},
            index=["cell-1", "cell-2"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    st.write_h5ad(data_root / "sample_st.h5ad")
    sc.write_h5ad(data_root / "sc.h5ad")


def test_application_cli_dry_run_performs_preflight_without_reconstruction(tmp_path):
    _write_inputs(tmp_path)
    output_root = tmp_path / "output"

    result = subprocess.run(
        [
            sys.executable,
            "application_reconstruct.py",
            "--svc-type",
            "sp-SVC",
            "--sample-name",
            "sample",
            "--st-file",
            "st.h5ad",
            "--sc-ref-file",
            "sc.h5ad",
            "--data-root",
            str(tmp_path),
            "--output-root",
            str(output_root),
            "--dry-run",
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["status"] == "ready"
    assert payload["svc_type"] == "sp-SVC"
    assert payload["pipeline"]["route"] == "sp_svc:bin2cell"
    assert "platform" not in payload
    assert Path(payload["preflight"]).is_file()
    assert not list(output_root.rglob("*.h5ad"))
