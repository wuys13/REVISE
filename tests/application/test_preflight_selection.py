from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise.config.runner_conf import InputSpec
from revise.io.input_service import REVISEInputService


def test_cluster_label_normalization_is_available_without_application_config():
    from revise.utils.labels import normalize_cell_type_label

    assert normalize_cell_type_label(" Mono/Macro ") == "Mono_Macro"


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


def test_preflight_rejects_selected_cell_type_missing_from_full_reference(tmp_path):
    st_path = tmp_path / "st.h5ad"
    reference_path = tmp_path / "reference.h5ad"
    _write(st_path, ["T", "T"], patient=["P1", "P1"])
    _write(reference_path, ["T", "Fibroblast"], patient=["P1", "P2"])

    service = REVISEInputService({"sample_name": "output"})
    with pytest.raises(ValueError, match="select_cell_type.*Mono_Macro"):
        service.preflight(
            (InputSpec("st", st_path), InputSpec("sc_ref", reference_path)),
            runtime={"mode": "application", "task": "sc_svc"},
            columns={
                "cell_type_col": "Level1",
                "sub_cell_type_col": "Level2",
                "select_cell_type": "Mono_Macro",
            },
        )


@pytest.mark.parametrize("selected", ["Mono/Macro", "Mono_Macro"])
def test_preflight_accepts_config_normalized_selected_cell_type(tmp_path, selected):
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
            "select_cell_type": selected.replace("/", "_"),
        },
    )

    assert report["status"] == "ready"


def test_preflight_accepts_reference_labels_normalized_like_selected_cell_type(tmp_path):
    st_path = tmp_path / "st.h5ad"
    reference_path = tmp_path / "reference.h5ad"
    _write(st_path, ["Mono/Macro", "Mono/Macro"])
    _write(reference_path, [" Mono/Macro ", "T"])

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


@pytest.mark.parametrize(
    ("filter_column", "filter_value", "message"),
    [
        ("Donor", "P2", "filter_column.*Donor"),
        ("Patient", "P3", "filter_value.*P3.*matched no rows"),
    ],
)
def test_preflight_rejects_unusable_reference_filter(
    tmp_path,
    filter_column,
    filter_value,
    message,
):
    st_path = tmp_path / "st.h5ad"
    reference_path = tmp_path / "reference.h5ad"
    _write(st_path, ["T", "T"])
    _write(reference_path, ["T", "Fibroblast"], patient=["P1", "P2"])

    with pytest.raises(ValueError, match=message):
        REVISEInputService().preflight(
            (InputSpec("st", st_path), InputSpec("sc_ref", reference_path)),
            runtime={"mode": "application", "task": "sc_svc"},
            columns={
                "cell_type_col": "Level1",
                "sub_cell_type_col": "Level2",
                "select_cell_type": "T",
            },
            reference_filter_column=filter_column,
            reference_filter_value=filter_value,
        )
