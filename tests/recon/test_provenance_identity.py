from __future__ import annotations

import copy
from hashlib import sha256
import json
import logging
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise.config.runner_conf import resolve_input_specs
from revise.config import (
    load_raw_config,
    merge_unified_config,
)
from revise.framework import REVISEPipeline
from revise.recon.context import PipelineContext
from revise.utils.deterministic import canonical_config_projection
from revise.utils.provenance import hash_jsonable, input_identities


CONFIG_PATH = Path(__file__).parents[2] / "revise" / "revise.yaml"


def _resolved_config(profile, algorithm_overrides=None):
    return merge_unified_config(
        raw_config=load_raw_config(CONFIG_PATH),
        profile=profile,
        runtime_overrides={},
        io_overrides={},
        algorithm_overrides=algorithm_overrides or {},
    )


def _config(tmp_path: Path) -> dict:
    return {
        "runtime": {
            "seed": 17,
            "deterministic": True,
            "platform": "sp_svc",
            "confounding": "bin2cell",
            "mode": "application",
            "task": "sp_svc",
            "svc_kind": "sp",
            "strategy": "SpSvcApplicationStrategy",
        },
        "io": {
            "data_root": str(tmp_path / "data"),
            "output_root": str(tmp_path / "output"),
            "sample_name": "sample",
            "st_file": "st.h5ad",
            "sc_ref_file": "sc.h5ad",
            "gt_svc_file": "gt.h5ad",
            "spatialdata_path": str(tmp_path / "input.zarr"),
            "seg_method": "seg_1",
            "spot_size": 50,
            "patient_key": "Patient",
            "sample_size": None,
            "save_outputs": True,
            "input_format": "h5ad",
            "spatialdata_reader": "zarr",
            "spatialdata_table": None,
            "spatialdata_spatial_element": None,
            "spatialdata_coordinate_system": "global",
        },
        "graph": {"method": "joint", "alpha": 0.5, "n_neighbors": 10},
        "ot": {
            "ga": {"solver": "pot", "pot": {"reg": 0.1}},
            "lr": {"solver": "pot", "pot": {"reg": 0.05}},
        },
        "plot": {"enabled": False},
    }


def _write_inputs(root: Path) -> None:
    root.mkdir(parents=True)
    st = AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame(index=["spot-1", "spot-2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    st.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 1.0]])
    sc_ref = AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame(
            {"Level1": ["A", "B"]},
            index=["cell-1", "cell-2"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    st.write_h5ad(root / "sample_Xenium.h5ad")
    sc_ref.write_h5ad(root / "adata_sc_all_reanno.h5ad")


@pytest.mark.parametrize("payload", [{"value": float("nan")}, {"value": object()}])
def test_hash_jsonable_rejects_non_json_semantic_values(payload):
    with pytest.raises((TypeError, ValueError)):
        hash_jsonable(payload)


def test_config_hash_excludes_only_input_and_output_locators(tmp_path):
    baseline = _config(tmp_path)
    projected = canonical_config_projection(baseline)
    assert set(projected["io"]).isdisjoint(
        {
            "data_root",
            "output_root",
            "st_file",
            "sc_ref_file",
            "gt_svc_file",
            "spatialdata_path",
        }
    )
    baseline_hash = hash_jsonable(projected)

    for key in (
        "data_root",
        "output_root",
        "st_file",
        "sc_ref_file",
        "gt_svc_file",
        "spatialdata_path",
    ):
        changed = copy.deepcopy(baseline)
        changed["io"][key] = f"different-{key}"
        assert hash_jsonable(canonical_config_projection(changed)) == baseline_hash


@pytest.mark.parametrize(
    ("section", "key", "value"),
    [
        ("runtime", "seed", 18),
        ("runtime", "confounding", "segmentation"),
        ("io", "sample_name", "other-sample"),
        ("io", "save_outputs", False),
        ("io", "spatialdata_table", "other-table"),
        ("graph", "alpha", 0.25),
    ],
)
def test_semantic_config_changes_change_hash(tmp_path, section, key, value):
    baseline = _config(tmp_path)
    changed = copy.deepcopy(baseline)
    changed[section][key] = value

    assert hash_jsonable(canonical_config_projection(changed)) != hash_jsonable(
        canonical_config_projection(baseline)
    )


@pytest.mark.parametrize("phase", ["ga", "lr"])
def test_solver_changes_change_config_hash(tmp_path, phase):
    baseline = _config(tmp_path)
    changed = copy.deepcopy(baseline)
    changed["ot"][phase]["solver"] = "tacco"

    assert hash_jsonable(canonical_config_projection(changed)) != hash_jsonable(
        canonical_config_projection(baseline)
    )


def test_effective_local_refinement_strength_changes_config_hash():
    baseline = _resolved_config("application_sp")
    changed = _resolved_config(
        "application_sp",
        {"local_refinement": {"strength": 0.5}},
    )

    assert hash_jsonable(canonical_config_projection(baseline)) != hash_jsonable(
        canonical_config_projection(changed)
    )


def test_route_default_and_explicit_equal_strength_share_hash():
    default = _resolved_config("application_sp")
    explicit = _resolved_config(
        "application_sp",
        {"local_refinement": {"strength": 0.2}},
    )

    assert hash_jsonable(canonical_config_projection(default)) == hash_jsonable(
        canonical_config_projection(explicit)
    )


def test_sc_route_omits_inactive_local_refinement_identity():
    resolved = _resolved_config("application_sc")

    assert "local_refinement" not in resolved
    assert "local_refinement" not in canonical_config_projection(resolved)


def test_pipeline_manifest_records_one_identity_per_input_role(tmp_path):
    left_data = tmp_path / "left-data"
    right_data = tmp_path / "right-data"
    _write_inputs(left_data)
    shutil.copytree(left_data, right_data, copy_function=shutil.copy2)

    manifests = []
    configs = []
    for name, data_root in (("left", left_data), ("right", right_data)):
        svc = REVISEPipeline().run(
            profile="application_sp",
            io_overrides={
                "data_root": str(data_root),
                "output_root": str(tmp_path / f"{name}-output"),
                "sample_name": "sample",
            },
            dry_run=True,
        )
        run_dir = Path(svc.provenance["run_dir"])
        manifests.append(json.loads((run_dir / "provenance.json").read_text()))
        configs.append(json.loads((run_dir / "merged_config.json").read_text()))

    assert manifests[0]["schema_version"] == manifests[1]["schema_version"] == 2
    assert manifests[0]["run_dir"] != manifests[1]["run_dir"]
    assert manifests[0]["run"]["started_at"] != manifests[1]["run"]["started_at"]
    assert configs[0]["io"]["data_root"] != configs[1]["io"]["data_root"]
    assert configs[0]["io"]["output_root"] != configs[1]["io"]["output_root"]
    assert manifests[0]["config_hash"] == manifests[1]["config_hash"]
    assert "data_fingerprint" not in manifests[0]
    assert "data_fingerprint" not in manifests[1]
    assert [record["role"] for record in manifests[0]["input_identities"]] == [
        "sc_ref",
        "st",
    ]
    assert [record["sha256"] for record in manifests[0]["input_identities"]] == [
        record["sha256"] for record in manifests[1]["input_identities"]
    ]
    assert [record["path"] for record in manifests[0]["input_identities"]] != [
        record["path"] for record in manifests[1]["input_identities"]
    ]
    for manifest, config in zip(manifests, configs):
        expected_specs = resolve_input_specs(config["runtime"], config["io"])
        assert manifest["runtime_seed"] == config["runtime"]["seed"]
        assert manifest["config_hash"] == hash_jsonable(
            canonical_config_projection(config)
        )
        assert manifest["input_identities"] == sorted(
            input_identities(expected_specs),
            key=lambda identity: identity["role"],
        )


def test_pipeline_computes_input_identities_once(
    monkeypatch,
    tmp_path,
):
    import revise.backend.policies as policies

    data_root = tmp_path / "data"
    _write_inputs(data_root)
    real_identities = policies.input_identities
    calls = []

    def counted(specs):
        calls.append(tuple(specs))
        return real_identities(specs)

    monkeypatch.setattr(policies, "input_identities", counted)
    REVISEPipeline().run(
        profile="application_sp",
        io_overrides={
            "data_root": str(data_root),
            "output_root": str(tmp_path / "output"),
            "sample_name": "sample",
        },
        dry_run=True,
    )

    assert len(calls) == 1


def test_sc_sr_manifest_adds_optional_pm_identity_and_isolates_pm_changes(
    tmp_path,
):
    data_root = tmp_path / "data"
    _write_inputs(data_root)
    pm_path = data_root / "sample_Xenium_PM_on_cell.csv"
    pm_path.write_text(",A,B\nc1,1,0\nc2,0,1\n", encoding="utf-8")

    manifests = []
    for name, replacement in (
        ("before", None),
        ("after", ",A,B\nc1,0,1\nc2,1,0\n"),
    ):
        if replacement is not None:
            pm_path.write_text(replacement, encoding="utf-8")
        svc = REVISEPipeline().run(
            profile="application_sc_sr",
            io_overrides={
                "data_root": str(data_root),
                "output_root": str(tmp_path / name),
                "sample_name": "sample",
            },
            dry_run=True,
        )
        manifests.append(
            json.loads(
                (Path(svc.provenance["run_dir"]) / "provenance.json").read_text()
            )
        )

    before = {item["role"]: item for item in manifests[0]["input_identities"]}
    after = {item["role"]: item for item in manifests[1]["input_identities"]}
    assert [item["role"] for item in manifests[0]["input_identities"]] == [
        "pm_on_cell",
        "sc_ref",
        "st",
    ]
    assert set(before) == set(after) == {"st", "sc_ref", "pm_on_cell"}
    assert before["st"]["sha256"] == after["st"]["sha256"]
    assert before["sc_ref"]["sha256"] == after["sc_ref"]["sha256"]
    assert before["pm_on_cell"]["sha256"] != after["pm_on_cell"]["sha256"]
    assert all("data_fingerprint" not in manifest for manifest in manifests)


def test_invalid_pm_preserves_all_read_input_identities(tmp_path):
    data_root = tmp_path / "data"
    output_root = tmp_path / "output"
    _write_inputs(data_root)
    pm_path = data_root / "sample_Xenium_PM_on_cell.csv"
    payload = b",A,B\nc1,invalid,0\nc2,0,1\n"
    pm_path.write_bytes(payload)

    with pytest.raises(ValueError, match="pm_on_cell"):
        REVISEPipeline().run(
            profile="application_sc_sr",
            io_overrides={
                "data_root": str(data_root),
                "output_root": str(output_root),
                "sample_name": "sample",
            },
            dry_run=True,
        )

    manifest_path = next(output_root.rglob("provenance.json"))
    manifest = json.loads(manifest_path.read_text())
    identities = {item["role"]: item for item in manifest["input_identities"]}
    assert set(identities) == {"pm_on_cell", "sc_ref", "st"}
    assert identities["pm_on_cell"]["sha256"] == sha256(payload).hexdigest()


def test_input_identity_failure_persists_terminal_manifest(monkeypatch, tmp_path):
    import revise.backend.policies as policies

    data_root = tmp_path / "data"
    output_root = tmp_path / "output"
    _write_inputs(data_root)

    def fail_identities(_specs):
        raise PermissionError("synthetic identity denial")

    monkeypatch.setattr(policies, "input_identities", fail_identities)
    with pytest.raises(ValueError, match="input identities"):
        REVISEPipeline().run(
            profile="application_sp",
            io_overrides={
                "data_root": str(data_root),
                "output_root": str(output_root),
                "sample_name": "sample",
            },
            dry_run=True,
        )

    manifests = list(output_root.rglob("provenance.json"))
    assert len(manifests) == 1
    manifest = json.loads(manifests[0].read_text())
    assert manifest["run"]["status"] == "failed"
    assert manifest["input_identities"] == []
    assert "input_identities_error" not in manifest
    assert manifest["run"]["error"]["type"] == "ValueError"


def test_invalid_semantic_config_fails_before_run_envelope(tmp_path):
    output_root = tmp_path / "output"
    with pytest.raises(ValueError, match="Out of range float values"):
        REVISEPipeline()._run_with_algorithm_overrides(
            profile="application_sp",
            io_overrides={
                "data_root": str(tmp_path / "data"),
                "output_root": str(output_root),
                "sample_name": "sample",
            },
            algorithm_overrides={"graph": {"alpha": float("nan")}},
            dry_run=True,
        )

    assert not output_root.exists()


def test_manifest_marks_unresolved_inputs_with_null_fingerprint(tmp_path):
    config = _config(tmp_path)
    ctx = PipelineContext(
        merged_config=config,
        raw_config={},
        config_path="revise/revise.yaml",
        profile="application_sp",
        runtime=config["runtime"],
        route_key="sp_svc:bin2cell",
        run_dir=tmp_path / "run",
        logger=logging.getLogger("test-unresolved-provenance"),
    )
    ctx.provenance["result"] = {
        "filename": "SVC.h5ad",
        "type": "sp-SVC",
    }
    ctx.provenance["results"] = {
        "spatial": {
            "filename": "sc_SVC_spatial.h5ad",
            "type": "sc-SVC",
        },
        "expression": {
            "filename": "sc_SVC_expr.h5ad",
            "type": "sc-SVC",
        },
    }
    REVISEPipeline.__new__(REVISEPipeline)._write_final_metadata(ctx)

    manifest = json.loads((ctx.run_dir / "provenance.json").read_text())
    assert manifest["schema_version"] == 2
    assert manifest["result"] == {
        "filename": "SVC.h5ad",
        "type": "sp-SVC",
    }
    assert manifest["results"] == {
        "spatial": {
            "filename": "sc_SVC_spatial.h5ad",
            "type": "sc-SVC",
        },
        "expression": {
            "filename": "sc_SVC_expr.h5ad",
            "type": "sc-SVC",
        },
    }
    assert ctx.provenance["results"] == manifest["results"]
    assert manifest["input_identities"] == []
    assert "data_fingerprint" not in manifest
