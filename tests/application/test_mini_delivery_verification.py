from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path
import sys

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml


def _module():
    path = Path(__file__).parents[2] / "scripts" / "verify_mini_delivery.py"
    spec = spec_from_file_location("verify_mini_delivery", path)
    module = module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _write_delivery(tmp_path, *, zero_support=False, diagnostic=True):
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
    if diagnostic:
        true_zero_pairs = np.array([[0, 0]], dtype=np.int32) if zero_support else np.empty((0, 2), dtype=np.int32)
        svc.uns["sst_parent_gene_correction"] = {
            "schema_version": 1,
            "operator": "share_then_target_no_epsilon",
            "target_definition": "internally normalized parent expression on the SVC gene axis",
            "current_sum_dtype": "float64",
            "positive_target_entries": 4,
            "positive_support_entries": 3 if zero_support else 4,
            "true_zero_support_entries": int(zero_support),
            "true_zero_support_target_mass": 2500.0 if zero_support else 0.0,
            "true_zero_support_target_max": 2500.0 if zero_support else 0.0,
            "min_positive_current_sum": 1.0,
            "true_zero_support_parent_ids": np.array(["spot-1"]) if zero_support else np.array([], dtype=str),
            "true_zero_support_pairs": true_zero_pairs,
        }
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
    assert conservation["pre_correction_support_status"] == "observed"
    assert conservation["true_zero_support_entries"] == 0
    assert conservation["true_zero_support_target_mass"] == 0.0
    assert conservation["true_zero_support_target_mass_fraction"] == 0.0
    assert conservation["true_zero_support_target_max"] == 0.0
    assert conservation["positive_support_nonconserving_entries"] == 0
    assert conservation["cell_total"] == {"min": 4000.0, "median": 6000.0, "max": 10000.0}


def test_sst_delivery_records_zero_support_as_partial_not_conserved(tmp_path):
    module, run = _write_delivery(tmp_path, zero_support=True)

    record = module.verify(run, "P2CRC_Visium", "default")

    conservation = record["parent_gene_conservation"]
    assert record["status"] == conservation["status"] == "partial"
    assert conservation["true_zero_support_entries"] == 1
    assert conservation["true_zero_support_target_mass"] == 2500.0
    assert conservation["true_zero_support_target_mass_fraction"] == 0.125
    assert conservation["true_zero_support_target_max"] == 2500.0
    assert conservation["parent_count_with_true_zero_support"] == 1
    assert conservation["parent_max_missing_mass_fraction"] == 0.25
    assert conservation["positive_support_nonconserving_entries"] == 0
    assert conservation["true_zero_support_examples"] == [{
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
    assert record["parent_gene_conservation"]["positive_support_nonconserving_entries"] == 1


def test_sst_delivery_rejects_tiny_positive_output_for_true_zero_support(tmp_path):
    module, run = _write_delivery(tmp_path, zero_support=True)
    path = run / "delivery/P2CRC_Visium_mini/default/SVC.h5ad"
    svc = ad.read_h5ad(path)
    svc.X[0, 0] = np.nextafter(0.0, 1.0)
    svc.write_h5ad(path)

    record = module.verify(run, "P2CRC_Visium", "default")

    assert record["status"] == "failed"
    assert record["parent_gene_conservation"]["true_zero_nonzero_output_entries"] == 1


def test_sst_delivery_without_pre_correction_diagnostic_is_unobserved(tmp_path):
    module, run = _write_delivery(tmp_path, zero_support=True, diagnostic=False)

    record = module.verify(run, "P2CRC_Visium", "default")

    conservation = record["parent_gene_conservation"]
    assert record["status"] == conservation["status"] == "unobserved"
    assert conservation["pre_correction_support_status"] == "unobserved"
    assert conservation["effective_zero_entries"] == 1
    assert "true_zero_support_entries" not in conservation


def test_sst_delivery_cli_exits_two_when_pre_correction_support_is_unobserved(
    tmp_path,
    monkeypatch,
):
    module, run = _write_delivery(tmp_path, zero_support=True, diagnostic=False)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "verify_mini_delivery.py",
            "--run",
            str(run),
            "--sample",
            "P2CRC_Visium",
        ],
    )

    with pytest.raises(SystemExit, match="2"):
        module.main()


def test_sst_correction_diagnostic_round_trips_and_its_sparse_mask_is_verified(tmp_path):
    module, run = _write_delivery(tmp_path, zero_support=True)
    path = run / "delivery/P2CRC_Visium_mini/default/SVC.h5ad"
    stored = ad.read_h5ad(path)
    diagnostic = stored.uns["sst_parent_gene_correction"]
    np.testing.assert_array_equal(diagnostic["true_zero_support_parent_ids"], ["spot-1"])
    np.testing.assert_array_equal(diagnostic["true_zero_support_pairs"], [[0, 0]])

    record = module.verify(run, "P2CRC_Visium", "default")

    assert record["parent_gene_conservation"]["true_zero_support_entries"] == 1
