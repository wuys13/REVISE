from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise.config.runner_conf import InputSpec
from revise.io.input_service import REVISEInputService


def _write(path: Path, labels, *, patient=None):
    obs = pd.DataFrame({"Level1": labels, "Level2": [f"{x}_sub" for x in labels]})
    if patient is not None:
        obs["Patient"] = patient
    adata = AnnData(
        X=np.ones((len(obs), 2), dtype=float),
        obs=obs,
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    adata.obsm["spatial"] = np.zeros((len(obs), 2), dtype=float)
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path)


def test_preflight_rejects_selected_cell_type_missing_after_patient_filter(tmp_path):
    st_path = tmp_path / "st.h5ad"
    reference_path = tmp_path / "reference.h5ad"
    _write(st_path, ["T", "T"], patient=["P1", "P1"])
    _write(reference_path, ["T", "Fibroblast"], patient=["P1", "P2"])

    service = REVISEInputService(
        {
            "sample_name": "P1",
            "patient_key": "Patient",
        }
    )
    with pytest.raises(ValueError, match="select_cell_type.*Fibroblast.*P1"):
        service.preflight(
            (InputSpec("st", st_path), InputSpec("sc_ref", reference_path)),
            runtime={"mode": "application", "task": "sc_svc"},
            columns={
                "cell_type_col": "Level1",
                "sub_cell_type_col": "Level2",
                "select_cell_type": "Fibroblast",
            },
        )


def test_preflight_accepts_slash_normalized_selected_cell_type(tmp_path):
    st_path = tmp_path / "st.h5ad"
    reference_path = tmp_path / "reference.h5ad"
    _write(st_path, ["Mono/Macro", "Mono/Macro"])
    _write(reference_path, ["Mono/Macro", "T"])

    report = REVISEInputService().preflight(
        (InputSpec("st", st_path), InputSpec("sc_ref", reference_path)),
        runtime={"mode": "application", "task": "sc_svc"},
        columns={
            "cell_type_col": "Level1",
            "sub_cell_type_col": "Level2",
            "select_cell_type": "Mono_Macro",
        },
    )

    assert report["status"] == "ready"
