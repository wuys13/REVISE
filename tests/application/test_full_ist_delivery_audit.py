from copy import deepcopy
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
import yaml

from revise.utils.provenance import hash_jsonable
from scripts.verify_full_ist_delivery import digest, verify


ROOT = Path(__file__).resolve().parents[2]


def _case(tmp_path, *, second_patient="P2CRC"):
    base_document = yaml.safe_load(
        (ROOT / "configs/acceptance/full/P2CRC_Xenium-random.yaml").read_text()
    )
    base_path = tmp_path / "configs/acceptance/full/P2CRC_Xenium-random.yaml"
    base_path.parent.mkdir(parents=True)
    base_path.write_text(yaml.safe_dump(base_document, sort_keys=False))

    run_document = deepcopy(base_document)
    run_document["inputs"]["st"]["path"] = "inputs/st.h5ad"
    run_document["inputs"]["reference"]["path"] = "inputs/ref.h5ad"
    run_document["output"]["dir"] = "results/run/random"
    config_path = tmp_path / "configs/acceptance/runs/iter-001/P2CRC_Xenium-random.yaml"
    config_path.parent.mkdir(parents=True)
    config_path.write_text(yaml.safe_dump(run_document, sort_keys=False))

    inputs = tmp_path / "inputs"
    inputs.mkdir()
    raw_obs = pd.DataFrame(
        {"Level1": ["T", "B", "B"], "Level2": ["T1", "B1", "B2"]},
        index=["s1", "s2", "s3"],
    )
    original = ad.AnnData(
        sparse.csr_matrix([[1.0], [2.0], [3.0]]),
        obs=raw_obs,
        var=pd.DataFrame(index=["g1"]),
    )
    original.obsm["spatial"] = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
    st_path = inputs / "st.h5ad"
    original.write_h5ad(st_path)

    reference = ad.AnnData(
        sparse.csr_matrix([[11.0, 12.0], [21.0, 22.0], [31.0, 32.0]]),
        obs=pd.DataFrame(
            {
                "Patient": ["P2CRC", second_patient, "P2CRC"],
                "Level1": ["T", "B", "B"],
                "Level2": ["T1", "B1", "B2"],
            },
            index=["d1", "d2", "d3"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    reference_path = inputs / "ref.h5ad"
    reference.write_h5ad(reference_path)

    output = tmp_path / "results/run/random"
    output.mkdir(parents=True)
    raw = original.copy()
    raw.obs["revise_Level1"] = pd.Categorical(["T", "B", "B"])
    raw.obs["revise_Level2"] = pd.Categorical(["T1", "B1", "B2"])
    raw.write_h5ad(output / "raw.h5ad")

    donor_ids = ["d1", "d2"]
    donor_clusters = ['["T","0"]', '["B","0"]']
    svc = ad.AnnData(
        reference[donor_ids].X.copy(),
        obs=pd.DataFrame(index=["s1", "s2"]),
        var=reference.var.copy(),
    )
    svc.obsm["spatial"] = original.obsm["spatial"][:2].copy()
    svc.obs["SVC_cluster"] = pd.Categorical(donor_clusters)
    svc.obs["revise_ist_donor_id"] = donor_ids
    svc.obs["revise_ist_donor_cluster"] = donor_clusters
    svc.obs["revise_Level1"] = pd.Categorical(["T", "B"])
    svc.obs["revise_Level2"] = pd.Categorical(["T1", "B1"])
    svc.uns["revise_reconstruction"] = {
        "ist_mapping": "random",
        "effective_seed": 42,
        "reference_filter": {"column": "Patient", "value": "P2CRC"},
        "processed_cell_types": ["B", "T"],
        "skipped_cell_types": {},
        "donor_column": "revise_ist_donor_id",
        "donor_cluster_column": "revise_ist_donor_cluster",
        "donor_sha256": hash_jsonable(donor_ids),
        "donor_cluster_sha256": hash_jsonable(donor_clusters),
    }
    svc.write_h5ad(output / "SVC.h5ad")
    sample = {
        "schema_version": 1,
        "sample_id": "P2CRC_Xenium",
        "files": {"raw": "raw.h5ad", "svc": "SVC.h5ad"},
        "columns": {
            "broad": "revise_Level1",
            "subtype": "revise_Level2",
            "reconstruction": "SVC_cluster",
        },
        "expression": {
            "raw": {"identity": "raw_counts"},
            "svc": {"identity": "reference_expression_random"},
        },
        "spatial": {"key": "spatial", "unit": "pixel", "microns_per_coordinate": 0.2125},
    }
    (output / "sample.yaml").write_text(yaml.safe_dump(sample, sort_keys=False))
    return {
        "repo": tmp_path,
        "config_path": config_path,
        "expected_st_sha256": digest(st_path),
        "expected_reference_sha256": digest(reference_path),
        "audit_path": tmp_path / "audit.json",
        "donor_sidecar": tmp_path / "donors.tsv.gz",
    }


def test_full_random_audit_allows_svc_reference_genes_absent_from_raw(tmp_path):
    record = verify(**_case(tmp_path))

    assert record["status"] == "technical_verified"
    assert record["gene_intersections"] == {
        "raw_svc": 1,
        "reference_svc": 2,
        "raw_reference": 1,
    }
    assert record["random_donor_audit"]["source_reference_rows_verified"] == 2


def test_full_random_audit_rejects_donor_outside_p2_reference_scope(tmp_path):
    with pytest.raises(AssertionError, match="outside Patient=P2CRC"):
        verify(**_case(tmp_path, second_patient="P1CRC"))


def test_full_random_audit_checks_every_reused_donor_assignment(tmp_path):
    arguments = _case(tmp_path)
    svc_path = tmp_path / "results/run/random/SVC.h5ad"
    svc = ad.read_h5ad(svc_path)
    donor_ids = ["d1", "d1"]
    donor_clusters = ['["T","0"]', '["T","0"]']
    svc.X = sparse.csr_matrix([[11.0, 12.0], [999.0, 999.0]])
    svc.obs["SVC_cluster"] = pd.Categorical(donor_clusters)
    svc.obs["revise_ist_donor_id"] = donor_ids
    svc.obs["revise_ist_donor_cluster"] = donor_clusters
    svc.obs["revise_Level1"] = pd.Categorical(["T", "T"])
    svc.uns["revise_reconstruction"]["donor_sha256"] = hash_jsonable(donor_ids)
    svc.uns["revise_reconstruction"]["donor_cluster_sha256"] = hash_jsonable(donor_clusters)
    svc.write_h5ad(svc_path)

    with pytest.raises(AssertionError, match="donor expression differs"):
        verify(**arguments)
