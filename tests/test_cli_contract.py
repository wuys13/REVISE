from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, read_h5ad
from revise import __version__
from revise.svc import SVC


ROOT = Path(__file__).resolve().parents[1]


def test_package_cli_exposes_version_without_required_run_arguments():
    result = subprocess.run(
        [sys.executable, "-m", "revise.cli", "--version"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == f"revise-reconstruct {__version__}"


def test_canonical_application_cli_rejects_benchmark_only_spot_size(monkeypatch):
    from revise import cli

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "revise-reconstruct",
            "--platform",
            "sST",
            "--sample-name",
            "sample",
            "--st-file",
            "st.h5ad",
            "--sc-ref-file",
            "sc.h5ad",
            "--data-root",
            "data",
            "--spot-size",
            "100",
        ],
    )

    with pytest.raises(SystemExit, match="2"):
        cli.get_args()


def test_root_wrapper_delegates_only_to_the_canonical_main():
    import reconstruct
    from revise import cli

    assert reconstruct.main is cli.main
    assert not hasattr(reconstruct, "ROUTES")


def test_canonical_cli_passes_publication_into_pipeline_finalize(monkeypatch, tmp_path):
    from revise import cli

    captured = {}

    class Pipeline:
        def __init__(self, config_path):
            captured["config_path"] = config_path

        def run(self, **kwargs):
            captured.update(kwargs)
            return SimpleNamespace()

    monkeypatch.setattr(cli, "REVISEPipeline", Pipeline)
    callback = lambda ctx: None
    args = SimpleNamespace(
        platform="hST",
        config="revise/revise.yaml",
        seed=17,
        data_root=str(tmp_path),
        output_root=str(tmp_path / "out"),
        sample_name="sample",
        st_file="st.h5ad",
        sc_ref_file="sc.h5ad",
        patient_key="Patient",
        spot_size=50,
        ot_method="pot",
        select_ct="all",
        cell_type_col="Level1",
        sub_cell_type_col="Level2",
        set_overrides=[],
    )

    cli._run_pipeline(args, finalize_callback=callback)

    assert captured["finalize_callback"] is callback
    assert captured["profile"] == "application_sp"
    assert captured["runtime_overrides"]["confounding"] == "bin2cell"
    assert "spot_size" not in captured["io_overrides"]


def test_public_result_links_to_manifest_and_registers_artifact(tmp_path):
    from revise import cli

    output = AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame(index=["cell-1", "cell-2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    records = []
    run_dir = tmp_path / "run"
    ctx = SimpleNamespace(
        svc=SVC(
            expr=output,
            spatial=output,
            svc_kind="sp",
            artifacts={"outputs": {"sp_svc": output}},
        ),
        run_dir=run_dir,
        merged_config={"ot": {"ga": {"solver": "pot"}, "lr": {"solver": "pot"}}},
        ot_events=[
            {"phase": "ga", "solver": "pot", "status": "completed", "call": 1}
        ],
        record_artifact=records.append,
    )
    args = SimpleNamespace(
        platform="hST",
        output_root=str(tmp_path / "out"),
        sample_name="sample",
        seed=17,
        ot_method="pot",
    )

    _, path = cli._build_public_result(
        args,
        "application_sp",
        "sp_svc",
        ctx,
    )

    published = read_h5ad(path)
    provenance = published.uns["revise_reconstruction"]
    expected_run_dir = Path(os.path.relpath(run_dir, start=path.parent)).as_posix()
    assert provenance["run_dir"] == expected_run_dir
    assert provenance["run_manifest"] == f"{expected_run_dir}/provenance.json"
    assert json.loads(provenance["ot_events"])[0]["status"] == "completed"
    assert records[0]["role"] == "public_result"
    assert records[0]["status"] == "completed"


def test_public_result_write_failure_preserves_previous_result(
    monkeypatch,
    tmp_path,
):
    from revise import cli

    output = AnnData(
        X=np.ones((1, 1)),
        obs=pd.DataFrame(index=["cell-1"]),
        var=pd.DataFrame(index=["g1"]),
    )
    output_path = tmp_path / "out" / "sample" / "hST-SVC.h5ad"
    output_path.parent.mkdir(parents=True)
    output_path.write_bytes(b"previous-valid-result")
    ctx = SimpleNamespace(
        svc=SVC(
            expr=output,
            spatial=output,
            svc_kind="sp",
            artifacts={"outputs": {"sp_svc": output}},
        ),
        run_dir=tmp_path / "run",
        merged_config={"ot": {"ga": {"solver": "pot"}, "lr": {"solver": "pot"}}},
        ot_events=[],
        record_artifact=lambda artifact: None,
    )
    args = SimpleNamespace(
        platform="hST",
        output_root=str(tmp_path / "out"),
        sample_name="sample",
        seed=17,
        ot_method="pot",
    )

    def partial_write(self, path, *args, **kwargs):
        Path(path).write_bytes(b"partial")
        raise OSError("simulated H5AD failure")

    monkeypatch.setattr(AnnData, "write_h5ad", partial_write)

    with pytest.raises(OSError, match="simulated H5AD failure"):
        cli._build_public_result(args, "application_sp", "sp_svc", ctx)

    assert output_path.read_bytes() == b"previous-valid-result"
    assert list(output_path.parent.iterdir()) == [output_path]
