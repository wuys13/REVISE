from __future__ import annotations

import copy
import json
import logging
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise.config.runner_conf import InputSpec
from revise.config.runner_conf import resolve_input_specs
from revise.framework import REVISEPipeline
from revise.recon.context import PipelineContext
from revise.utils.deterministic import canonical_config_projection
from revise.utils.provenance import fingerprint_paths, hash_jsonable


def _config(tmp_path: Path) -> dict:
    return {
        "runtime": {
            "seed": 17,
            "deterministic": True,
            "platform": "hST",
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


def test_file_fingerprint_binds_roles_and_bytes_not_location_or_metadata(tmp_path):
    left = tmp_path / "left"
    right = tmp_path / "right"
    left.mkdir()
    right.mkdir()
    left_st = left / "st.bin"
    left_sc = left / "sc.bin"
    left_st.write_bytes(b"st-content")
    left_sc.write_bytes(b"sc-content")
    right_st = right / "renamed-st.bin"
    right_sc = right / "renamed-sc.bin"
    shutil.copy2(left_st, right_st)
    shutil.copy2(left_sc, right_sc)
    os.utime(right_st, ns=(1_000_000_000, 1_000_000_000))

    baseline = fingerprint_paths(
        (InputSpec("st", str(left_st)), InputSpec("sc_ref", str(left_sc)))
    )
    moved = fingerprint_paths(
        (InputSpec("sc_ref", str(right_sc)), InputSpec("st", str(right_st)))
    )
    swapped = fingerprint_paths(
        (InputSpec("st", str(right_sc)), InputSpec("sc_ref", str(right_st)))
    )

    assert baseline == moved
    assert baseline != swapped


def test_equal_size_equal_mtime_byte_change_changes_fingerprint(tmp_path):
    path = tmp_path / "input.bin"
    path.write_bytes(b"first-value")
    before_stat = path.stat()
    before = fingerprint_paths((InputSpec("st", str(path)),))

    path.write_bytes(b"other-value")
    os.utime(path, ns=(before_stat.st_atime_ns, before_stat.st_mtime_ns))
    after_stat = path.stat()
    after = fingerprint_paths((InputSpec("st", str(path)),))

    assert after_stat.st_size == before_stat.st_size
    assert after_stat.st_mtime_ns == before_stat.st_mtime_ns
    assert after != before


def test_missing_empty_file_and_empty_directory_have_distinct_identities(tmp_path):
    empty_file = tmp_path / "empty.bin"
    empty_file.touch()
    empty_directory = tmp_path / "empty.zarr"
    empty_directory.mkdir()

    missing_a = fingerprint_paths((InputSpec("st", str(tmp_path / "missing-a")),))
    missing_b = fingerprint_paths((InputSpec("st", str(tmp_path / "missing-b")),))
    missing_other_role = fingerprint_paths(
        (InputSpec("sc_ref", str(tmp_path / "missing-a")),)
    )
    empty_file_id = fingerprint_paths((InputSpec("st", str(empty_file)),))
    empty_directory_id = fingerprint_paths((InputSpec("st", str(empty_directory)),))

    assert missing_a == missing_b
    assert len({missing_a, missing_other_role, empty_file_id, empty_directory_id}) == 4


def test_directory_fingerprint_tracks_relative_structure_and_nested_bytes(tmp_path):
    left = tmp_path / "left.zarr"
    right = tmp_path / "right.zarr"
    (left / "tables").mkdir(parents=True)
    (left / "empty-group").mkdir()
    (left / ".zgroup").write_bytes(b"group-v1")
    nested = left / "tables" / "chunk-0"
    nested.write_bytes(b"chunk-v1")
    shutil.copytree(left, right)

    baseline = fingerprint_paths((InputSpec("st", str(left)),))
    assert baseline == fingerprint_paths((InputSpec("st", str(right)),))

    before_stat = (right / "tables" / "chunk-0").stat()
    (right / "tables" / "chunk-0").write_bytes(b"chunk-v2")
    os.utime(
        right / "tables" / "chunk-0",
        ns=(before_stat.st_atime_ns, before_stat.st_mtime_ns),
    )
    assert baseline != fingerprint_paths((InputSpec("st", str(right)),))

    shutil.rmtree(right)
    shutil.copytree(left, right)
    (right / "tables" / "chunk-0").rename(right / "tables" / "chunk-renamed")
    assert baseline != fingerprint_paths((InputSpec("st", str(right)),))


def test_fingerprint_rejects_duplicate_roles_and_resolves_root_symbolic_links(
    tmp_path,
):
    path = tmp_path / "input.bin"
    path.write_bytes(b"content")
    duplicate = (InputSpec("st", str(path)), InputSpec("st", str(path)))
    with pytest.raises(ValueError, match="duplicate input role"):
        fingerprint_paths(duplicate)

    link = tmp_path / "link.bin"
    link.symlink_to(path)
    assert fingerprint_paths((InputSpec("st", str(link)),)) == fingerprint_paths(
        (InputSpec("st", str(path)),)
    )


def test_fingerprint_does_not_treat_unenumerable_directory_as_empty(
    monkeypatch,
    tmp_path,
):
    import revise.utils.provenance as provenance

    empty = tmp_path / "empty.zarr"
    blocked = tmp_path / "blocked.zarr"
    empty.mkdir()
    blocked.mkdir()
    real_scandir = provenance.os.scandir

    def guarded_scandir(path):
        if Path(path) == blocked:
            raise PermissionError("blocked fixture")
        return real_scandir(path)

    monkeypatch.setattr(provenance.os, "scandir", guarded_scandir)
    empty_id = fingerprint_paths((InputSpec("st", str(empty)),))
    with pytest.raises(PermissionError, match="blocked fixture"):
        fingerprint_paths((InputSpec("st", str(blocked)),))
    assert empty_id


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


def test_pipeline_manifest_uses_location_independent_identities(tmp_path):
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
    assert manifests[0]["data_fingerprint"] == manifests[1]["data_fingerprint"]
    assert manifests[0]["data_fingerprint"] is not None
    for manifest, config in zip(manifests, configs):
        expected_specs = resolve_input_specs(config["runtime"], config["io"])
        assert manifest["config_hash"] == hash_jsonable(
            canonical_config_projection(config)
        )
        assert manifest["data_fingerprint"] == fingerprint_paths(expected_specs)


def test_pipeline_computes_content_fingerprint_once(
    monkeypatch,
    tmp_path,
):
    import revise.backend.policies as policies

    data_root = tmp_path / "data"
    _write_inputs(data_root)
    real_fingerprint = policies.fingerprint_paths
    calls = []

    def counted(specs):
        calls.append(tuple(specs))
        return real_fingerprint(specs)

    monkeypatch.setattr(policies, "fingerprint_paths", counted)
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


def test_fingerprint_failure_persists_terminal_manifest(monkeypatch, tmp_path):
    import revise.backend.policies as policies

    data_root = tmp_path / "data"
    output_root = tmp_path / "output"
    _write_inputs(data_root)

    def fail_fingerprint(_specs):
        raise PermissionError("synthetic fingerprint denial")

    monkeypatch.setattr(policies, "fingerprint_paths", fail_fingerprint)
    with pytest.raises(ValueError, match="content fingerprint"):
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
    assert manifest["data_fingerprint"] is None
    assert manifest["data_fingerprint_error"]["type"] == "PermissionError"


def test_invalid_semantic_config_fails_before_run_envelope(tmp_path):
    output_root = tmp_path / "output"
    with pytest.raises(ValueError, match="Out of range float values"):
        REVISEPipeline().run(
            profile="application_sp",
            io_overrides={
                "data_root": str(tmp_path / "data"),
                "output_root": str(output_root),
                "sample_name": "sample",
            },
            set_overrides=["graph.alpha=.nan"],
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
        route_key="hST:bin2cell",
        run_dir=tmp_path / "run",
        logger=logging.getLogger("test-unresolved-provenance"),
    )
    ctx.provenance["result"] = {
        "filename": "SVC.h5ad",
        "type": "hST",
    }
    REVISEPipeline.__new__(REVISEPipeline)._write_final_metadata(ctx)

    manifest = json.loads((ctx.run_dir / "provenance.json").read_text())
    assert manifest["schema_version"] == 2
    assert manifest["result"] == {
        "filename": "SVC.h5ad",
        "type": "hST",
    }
    assert manifest["data_fingerprint"] is None
