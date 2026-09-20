from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml


def _module():
    path = Path(__file__).parents[2] / "scripts" / "verify_mini_delivery.py"
    spec = spec_from_file_location("verify_mini_delivery", path)
    module = module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _write_delivery(tmp_path, *, zero_support=False):
    module = _module()
    module.ROOT = tmp_path
    run = tmp_path / "run"
    delivery = run / "delivery" / "P2CRC_Visium_mini" / "default"
    delivery.mkdir(parents=True)
    original = ad.AnnData(
        np.array([[1.0, 3.0], [2.0, 2.0]]),
        obs=pd.DataFrame(index=["spot-1", "spot-2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    original.obsm["spatial"] = np.array([[1.0, 2.0], [3.0, 4.0]])
    input_path = tmp_path / "inputs" / "P2CRC_Visium_mini.h5ad"
    input_path.parent.mkdir()
    original.write_h5ad(input_path)
    original.write_h5ad(delivery / "raw.h5ad")
    first_gene = [0.0, 0.0, 5000.0] if zero_support else [1000.0, 1500.0, 5000.0]
    svc = ad.AnnData(
        np.column_stack([first_gene, [3000.0, 4500.0, 5000.0]]),
        obs=pd.DataFrame(
            {"spot_name": ["spot-1", "spot-1", "spot-2"]},
            index=["cell-1", "cell-2", "cell-3"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    svc.obsm["spatial"] = np.array([[1.0, 2.0], [1.0, 2.0], [3.0, 4.0]])
    svc.write_h5ad(delivery / "SVC.h5ad")
    (delivery / "sample.yaml").write_text(yaml.safe_dump({
        "sample_id": "P2CRC_Visium_mini",
        "files": {"raw": "raw.h5ad", "svc": "SVC.h5ad"},
        "expression": {"raw": {}, "svc": {}},
        "spatial": {"key": "spatial", "unit": "pixel"},
    }))
    (run / "input_manifest.json").write_text(json.dumps({"samples": {
        "P2CRC_Visium": {
            "mini_input": "inputs/P2CRC_Visium_mini.h5ad",
            "mini_sha256": module.digest(input_path),
        }
    }}))
    return module, run


def test_sst_delivery_reports_exact_parent_gene_conservation(tmp_path):
    module, run = _write_delivery(tmp_path)

    record = module.verify(run, "P2CRC_Visium", "default")

    conservation = record["parent_gene_conservation"]
    assert record["status"] == conservation["status"] == "passed"
    assert conservation["zero_support_entries"] == 0
    assert conservation["zero_support_target_mass"] == 0.0
    assert conservation["zero_support_target_mass_fraction"] == 0.0
    assert conservation["zero_support_max_target"] == 0.0
    assert conservation["supported_nonconserving_entries"] == 0
    assert conservation["cell_total"] == {"min": 4000.0, "median": 6000.0, "max": 10000.0}


def test_sst_delivery_records_zero_support_as_partial_not_conserved(tmp_path):
    module, run = _write_delivery(tmp_path, zero_support=True)

    record = module.verify(run, "P2CRC_Visium", "default")

    conservation = record["parent_gene_conservation"]
    assert record["status"] == conservation["status"] == "partial"
    assert conservation["zero_support_entries"] == 1
    assert conservation["zero_support_target_mass"] == 2500.0
    assert conservation["zero_support_target_mass_fraction"] == 0.125
    assert conservation["zero_support_max_target"] == 2500.0
    assert conservation["supported_nonconserving_entries"] == 0
    assert conservation["zero_support_examples"] == [{
        "parent_id": "spot-1",
        "gene": "g1",
        "normalized_parent_target": 2500.0,
        "aggregated_svc": 0.0,
    }]


def test_sst_delivery_does_not_hide_supported_numeric_mismatch(tmp_path):
    module, run = _write_delivery(tmp_path)
    path = run / "delivery/P2CRC_Visium_mini/default/SVC.h5ad"
    svc = ad.read_h5ad(path)
    svc.X[0, 0] *= 0.5
    svc.write_h5ad(path)

    record = module.verify(run, "P2CRC_Visium", "default")

    assert record["status"] == "failed"
    assert record["parent_gene_conservation"]["supported_nonconserving_entries"] == 1
