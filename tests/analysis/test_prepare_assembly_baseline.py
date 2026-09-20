from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, read_h5ad

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts/prepare_assembly_baseline.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("prepare_assembly_baseline_test", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_carriers(tmp_path: Path, *, duplicate: bool = False):
    paths = {}
    for offset, cell_type in enumerate(("T", "Mono_Macro", "Fibroblast")):
        ids = [f"{cell_type}-{i}" for i in range(5)]
        if duplicate and cell_type == "Fibroblast":
            ids[0] = "T-0"
        if cell_type == "T":
            labels = ["T", "T", pd.NA, "T", "T"]
        elif cell_type == "Mono_Macro":
            labels = ["Mono/Macro"] * 5
        else:
            labels = ["Fibroblast"] * 5
        obs = pd.DataFrame(
            {
                "Level1": pd.Series(labels, index=ids, dtype="object"),
                "SVC_cluster": pd.Series(
                    [str(i % 2) for i in range(5)], index=ids, dtype="object"
                ),
            },
            index=ids,
        )
        adata = AnnData(
            np.arange(25, dtype=np.float32).reshape(5, 5),
            obs=obs,
            obsm={"spatial": np.column_stack([np.arange(5) + offset, np.arange(5)]).astype(float)},
        )
        path = tmp_path / f"{cell_type}.h5ad"
        adata.write_h5ad(path)
        paths[cell_type] = path
    return paths


def test_assemble_baseline_reads_obs_spatial_and_writes_zero_gene_axis(tmp_path):
    module = _load_module()
    carriers = _write_carriers(tmp_path)
    output = tmp_path / "baseline.h5ad"

    summary = module.assemble_baseline(carriers, output)

    assert summary["n_obs"] == 15
    assert summary["n_vars"] == 0
    baseline = read_h5ad(output)
    assert baseline.shape == (15, 0)
    assert baseline.obs_names.tolist() == [
        f"{cell_type}-{i}"
        for cell_type in ("T", "Mono_Macro", "Fibroblast")
        for i in range(5)
    ]
    assert baseline.obs["Level1"].tolist()[5:10] == ["Mono_Macro"] * 5
    assert pd.isna(baseline.obs["Level1"].iloc[2])
    np.testing.assert_array_equal(
        baseline.obsm["spatial"],
        np.concatenate(
            [
                np.column_stack([np.arange(5) + offset, np.arange(5)]).astype(float)
                for offset in range(3)
            ]
        ),
    )
    assert baseline.obs["SVC_cluster"].tolist() == [
        str(i % 2) for _ in range(3) for i in range(5)
    ]
    provenance = baseline.uns["assembly_baseline"]
    assert bool(provenance["zero_gene_axis"]) is True
    assert bool(provenance["expression_loaded"]) is False
    assert provenance["coordinate_provenance"]["unit"] == "pixel"
    assert provenance["coordinate_provenance"]["microns_per_coordinate"] == 0.2125
    assert set(provenance["source_sha256"]) == {"T", "Mono_Macro", "Fibroblast"}
    for cell_type, path in carriers.items():
        assert provenance["source_sha256"][cell_type] == hashlib.sha256(
            path.read_bytes()
        ).hexdigest()


def test_assemble_baseline_rejects_cross_carrier_duplicate_ids(tmp_path):
    module = _load_module()
    carriers = _write_carriers(tmp_path, duplicate=True)
    with pytest.raises(ValueError, match="cross-carrier duplicate"):
        module.assemble_baseline(carriers, tmp_path / "baseline.h5ad")


def test_baseline_script_uses_backed_reader_and_never_references_source_x():
    source = SCRIPT.read_text(encoding="utf-8")
    assert 'read_h5ad(source_path, backed="r")' in source
    assert "native.X" not in source
