from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise.framework import REVISEPipeline
from revise.utils import provenance
from revise.utils.io import build_run_dir


def _small_adata(*, obs_names=("cell-1", "cell-2")) -> AnnData:
    return AnnData(
        X=np.ones((len(obs_names), 2)),
        obs=pd.DataFrame(index=pd.Index(obs_names)),
        var=pd.DataFrame(index=pd.Index(["g1", "g2"])),
    )


def _write_application_inputs(data_root: Path, sample_name: str) -> None:
    st = _small_adata(obs_names=("spot-1", "spot-2"))
    st.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 1.0]])
    sc_ref = _small_adata()
    sc_ref.obs["Level1"] = ["A", "B"]
    st.write_h5ad(data_root / f"{sample_name}_Xenium.h5ad")
    sc_ref.write_h5ad(data_root / "adata_sc_all_reanno.h5ad")


def _write_benchmark_seg_inputs(data_root: Path, sample_name: str) -> None:
    st = _small_adata()
    st.obs["seg_error"] = [0, 0]
    st.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 1.0]])
    sc_ref = _small_adata()
    sc_ref.obs["Level1"] = ["A", "B"]
    gt = _small_adata()
    st_path = data_root / sample_name / "seg_1" / "Xenium.h5ad"
    st_path.parent.mkdir(parents=True, exist_ok=True)
    st.write_h5ad(st_path)
    sc_ref.write_h5ad(data_root / sample_name / "adata_sc_all_reanno.h5ad")
    gt.write_h5ad(data_root / sample_name / "selected_xenium.h5ad")


def test_atomic_json_failure_preserves_previous_manifest(monkeypatch, tmp_path):
    path = tmp_path / "provenance.json"
    provenance.write_json(path, {"generation": 1})

    def broken_dump(payload, handle, *args, **kwargs):
        handle.write('{"generation":')
        handle.flush()
        raise OSError("simulated write failure")

    monkeypatch.setattr(provenance.json, "dump", broken_dump)
    with pytest.raises(OSError, match="simulated write failure"):
        provenance.write_json(path, {"generation": 2})

    assert json.loads(path.read_text(encoding="utf-8")) == {"generation": 1}
    assert list(tmp_path.iterdir()) == [path]


def test_completed_artifact_uses_streamed_content_identity(tmp_path):
    path = tmp_path / "result.bin"
    path.write_bytes(b"first-content")
    expected = hashlib.sha256(path.read_bytes()).hexdigest()

    record = provenance.completed_artifact("output:primary", path)

    assert record == {
        "role": "output:primary",
        "path": str(path),
        "status": "completed",
        "size": len(b"first-content"),
        "sha256": expected,
    }

    path.write_bytes(b"other-content")
    assert path.stat().st_size == len(b"first-content")
    assert provenance.sha256_file(path) != expected


@pytest.mark.parametrize(
    ("recorded", "process_alive", "effective"),
    [
        ("succeeded", False, "succeeded"),
        ("failed", False, "failed"),
        ("interrupted", False, "interrupted"),
        ("running", True, "running"),
        ("running", False, "incomplete"),
    ],
)
def test_effective_run_status_does_not_rewrite_recorded_history(
    recorded,
    process_alive,
    effective,
):
    manifest = {"run": {"status": recorded}}

    assert provenance.effective_run_status(
        manifest,
        process_alive=process_alive,
    ) == effective
    assert manifest == {"run": {"status": recorded}}


def test_application_run_directories_do_not_collide_when_created_together(tmp_path):
    first = build_run_dir(str(tmp_path), "sample", "sp_svc:bin2cell")
    second = build_run_dir(str(tmp_path), "sample", "sp_svc:bin2cell")

    assert first != second
    assert first.parent == second.parent


def test_fixed_benchmark_directory_refuses_running_envelope(tmp_path):
    output_root = tmp_path / "benchmark-output"
    run_dir = output_root / "sample" / "seg_1"
    manifest_path = run_dir / "provenance.json"
    provenance.write_json(manifest_path, {"run": {"status": "running"}})

    with pytest.raises(RuntimeError, match="running.*provenance"):
        REVISEPipeline().run(
            profile="benchmark_seg",
            io_overrides={
                "data_root": str(tmp_path),
                "output_root": str(output_root),
                "sample_name": "sample",
                "seg_method": "seg_1",
            },
            dry_run=True,
        )

    assert json.loads(manifest_path.read_text(encoding="utf-8")) == {
        "run": {"status": "running"}
    }


def test_run_directory_lock_rejects_concurrent_owner(tmp_path):
    run_dir = Path(tmp_path) / "shared-run"

    with provenance.exclusive_run_directory(run_dir):
        with pytest.raises(RuntimeError, match="already active"):
            with provenance.exclusive_run_directory(run_dir):
                pass


def test_run_directory_lock_releases_after_success_and_exception(tmp_path):
    run_dir = Path(tmp_path) / "reusable-run"
    lock_path = run_dir / ".revise-run.lock"

    with provenance.exclusive_run_directory(run_dir):
        assert lock_path.is_dir()
    assert not lock_path.exists()

    with pytest.raises(RuntimeError, match="inside run"):
        with provenance.exclusive_run_directory(run_dir):
            raise RuntimeError("inside run")
    assert not lock_path.exists()

    with provenance.exclusive_run_directory(run_dir):
        assert lock_path.is_dir()


def test_fixed_benchmark_directory_can_be_reused_after_terminal_run(tmp_path):
    output_root = tmp_path / "reusable-benchmark"
    _write_benchmark_seg_inputs(tmp_path, "sample")
    options = {
        "profile": "benchmark_seg",
        "io_overrides": {
            "data_root": str(tmp_path),
            "output_root": str(output_root),
            "sample_name": "sample",
            "seg_method": "seg_1",
        },
        "dry_run": True,
    }

    REVISEPipeline().run(**options)
    first = json.loads(
        (output_root / "sample" / "seg_1" / "provenance.json").read_text()
    )
    REVISEPipeline().run(**options)
    second = json.loads(
        (output_root / "sample" / "seg_1" / "provenance.json").read_text()
    )

    assert first["run"]["status"] == "succeeded"
    assert second["run"]["status"] == "succeeded"
    assert not (output_root / "sample" / "seg_1" / ".revise-run.lock").exists()


def test_application_manifest_and_log_share_each_unique_run_directory(tmp_path):
    output_root = tmp_path / "application-output"
    _write_application_inputs(tmp_path, "sample")

    for _ in range(2):
        REVISEPipeline().run(
            profile="application_sp",
            io_overrides={
                "data_root": str(tmp_path),
                "output_root": str(output_root),
                "sample_name": "sample",
            },
            dry_run=True,
        )

    route_dir = output_root / "sample" / "sp_svc__bin2cell"
    manifests = sorted(route_dir.glob("*/provenance.json"))
    assert len(manifests) == 2
    assert len({path.parent for path in manifests}) == 2
    assert all((path.parent / "run.log").is_file() for path in manifests)
    assert {path for path in route_dir.iterdir()} == {
        manifest.parent for manifest in manifests
    }
