from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData
import yaml


ROOT = Path(__file__).resolve().parents[2]


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
    config_path = tmp_path / "VisiumHD.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "schema_version": 1,
                "application": {"svc_type": "sp-SVC", "sample_name": "sample"},
                "paths": {"root_dir": str(tmp_path)},
                "algorithm": {},
                "inputs": {
                    "mode": "direct",
                    "st": {"path": "sample_st.h5ad", "format": "h5ad"},
                    "reference": {
                        "path": "sc.h5ad",
                        "format": "h5ad",
                        "patient_key": "Patient",
                    },
                },
                "global_anchoring": {"broad_column": "Level1"},
                "local_refinement": {"strength": 0.2},
                "output": {"path": "output"},
                "execution": {"action": "run", "seed": 42},
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    result = subprocess.run(
        [
            sys.executable,
            "reconstruct.py",
            "--config",
            str(config_path),
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
    assert payload["pipeline"]["route"]["application_route"] == "sp-SVC"
    assert "confounding" not in payload["pipeline"]["route"]
    assert "platform" not in payload
    assert Path(payload["preflight"]).is_file()
    assert not list(output_root.rglob("*.h5ad"))
    provenance = json.loads(
        (Path(payload["preflight"]).parent / "provenance.json").read_text()
    )
    assert provenance["application_config"] == {
        "source_path": str(config_path.resolve()),
        "source_sha256": __import__("hashlib").sha256(config_path.read_bytes()).hexdigest(),
        "declared_root": str(tmp_path),
        "resolved_root": str(tmp_path.resolve()),
        "cwd": str(ROOT.resolve()),
        "resolved_paths": {
            "st": str((tmp_path / "sample_st.h5ad").resolve()),
            "reference": str((tmp_path / "sc.h5ad").resolve()),
            "output": str(output_root.resolve()),
        },
        "declared_action": "run",
        "effective_action": "preflight",
        "dry_run_override": True,
    }
