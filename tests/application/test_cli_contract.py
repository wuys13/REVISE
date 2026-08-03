from __future__ import annotations

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


class PublicationContext(SimpleNamespace):
    def set_pending_publication(self, *, commit, rollback):
        self._publication_commit = commit
        self._publication_rollback = rollback

    def rollback_pending_publication(self):
        rollback = getattr(self, "_publication_rollback", None)
        self._publication_commit = None
        self._publication_rollback = None
        if rollback is not None:
            rollback()


ROOT = Path(__file__).resolve().parents[2]


def test_package_cli_exposes_version_without_required_run_arguments():
    result = subprocess.run(
        [sys.executable, "-m", "revise.application.cli", "--version"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == f"revise-reconstruct {__version__}"


def test_canonical_application_cli_rejects_benchmark_only_spot_size(monkeypatch):
    from revise.application import cli

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "revise-reconstruct",
            "--svc-type",
            "sST-SVC",
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
        cli.parse_args()


def test_canonical_application_cli_exposes_only_local_refinement_strength():
    from revise.application import cli

    args = cli.parse_args(
        [
            "--svc-type",
            "hST-SVC",
            "--sample-name",
            "sample",
            "--st-file",
            "st.h5ad",
            "--sc-ref-file",
            "sc.h5ad",
            "--data-root",
            "data",
            "--local-refinement-strength",
            "0.4",
        ]
    )

    assert args.local_refinement_strength == pytest.approx(0.4)
    help_text = cli.build_parser().format_help()
    assert "--local-refinement-strength" in help_text
    assert "--local-refinement-guidance" not in help_text
    assert "--posterior-conditioning-mode" not in help_text


@pytest.mark.parametrize(
    "select_ct",
    [None, "", " ", "all", "*", "__all__", "all_cell_types"],
)
def test_sc_svc_cli_requires_one_concrete_cell_type(select_ct, capsys):
    from revise.application import cli

    argv = [
        "--svc-type",
        "sc-SVC",
        "--sample-name",
        "sample",
        "--st-file",
        "st.h5ad",
        "--sc-ref-file",
        "sc.h5ad",
        "--data-root",
        "data",
    ]
    if select_ct is not None:
        argv.extend(["--select-ct", select_ct])

    with pytest.raises(SystemExit, match="2"):
        cli.parse_args(argv)

    assert "--select-ct must name one concrete cell type" in capsys.readouterr().err


def test_canonical_application_cli_removed_guidance_flag_has_migration_error(monkeypatch):
    from revise.application import cli

    parser = cli.build_parser()

    def fail(message):
        raise ValueError(message)

    monkeypatch.setattr(parser, "error", fail)
    with pytest.raises(
        ValueError,
        match=(
            "Assignment guidance options were removed; "
            "use --local-refinement-strength"
        ),
    ):
        parser.parse_args(["--posterior-mode", "cost"])


def test_root_wrapper_delegates_only_to_the_canonical_main():
    import reconstruct
    from revise.application import cli

    assert reconstruct.main is cli.main
    assert not hasattr(reconstruct, "APPLICATION_ROUTES")


def test_canonical_cli_passes_publication_into_pipeline_finalize(monkeypatch, tmp_path):
    from revise.application import service

    captured = {}

    class Pipeline:
        def __init__(self, config_path):
            captured["config_path"] = config_path

        def _run_with_algorithm_overrides(self, **kwargs):
            captured.update(kwargs)
            return SimpleNamespace()

    monkeypatch.setattr(service, "REVISEPipeline", Pipeline)
    def callback(ctx):
        return None

    args = SimpleNamespace(
        svc_type="hST-SVC",
        ist_mapping=None,
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
        local_refinement_strength=0.4,
        cell_type_col="Level1",
        sub_cell_type_col="Level2",
    )

    service._run_pipeline(args, finalize_callback=callback)

    assert captured["finalize_callback"] is callback
    assert captured["profile"] == "application_sp"
    assert captured["runtime_overrides"]["confounding"] == "bin2cell"
    assert "spot_size" not in captured["io_overrides"]
    assert captured["algorithm_overrides"]["local_refinement"] == {
        "strength": 0.4,
    }


@pytest.mark.parametrize(
    ("svc_type", "profile", "svc_kind", "output_key", "expected_type"),
    [
        ("hST-SVC", "application_sp", "sp", "sp_svc", "hST-SVC"),
        (
            "sST-SVC",
            "application_sc_sr",
            "sc",
            "sc_svc_dec",
            "sST-SVC",
        ),
    ],
)
def test_public_result_links_to_manifest_and_registers_artifact(
    tmp_path,
    svc_type,
    profile,
    svc_kind,
    output_key,
    expected_type,
):
    from revise.application import service

    output = AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame(index=["cell-1", "cell-2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    records = []
    ctx = PublicationContext(
        svc=SVC(
            expr=output,
            spatial=output,
            svc_kind=svc_kind,
            artifacts={"outputs": {output_key: output}},
        ),
        run_dir=tmp_path / "run",
        merged_config={"ot": {"ga": {"solver": "pot"}, "lr": {"solver": "pot"}}},
        provenance={},
        record_artifact=records.append,
    )
    args = SimpleNamespace(
        svc_type=svc_type,
        output_root=str(tmp_path / "out"),
        sample_name="sample",
        seed=17,
        ot_method="pot",
    )

    result, path = service._build_public_result(
        args,
        profile,
        output_key,
        ctx,
    )

    published = read_h5ad(path)
    assert result.uns["revise_reconstruction"] == {
        "schema_version": 2,
        "svc_type": expected_type,
        "ist_mapping": None,
        "effective_seed": None,
        "expression_source": None,
        "donor_column": None,
        "donor_sha256": None,
    }
    assert published.uns["revise_reconstruction"] == {
        "schema_version": 2,
        "svc_type": expected_type,
    }
    assert path == tmp_path / "out" / "sample" / expected_type / "SVC.h5ad"
    assert ctx.provenance["result"] == {
        "filename": "SVC.h5ad",
        "type": expected_type,
    }
    assert records[0]["role"] == "public_result"
    assert records[0]["status"] == "completed"


def test_public_result_write_failure_preserves_previous_result(
    monkeypatch,
    tmp_path,
):
    from revise.application import service

    output = AnnData(
        X=np.ones((1, 1)),
        obs=pd.DataFrame(index=["cell-1"]),
        var=pd.DataFrame(index=["g1"]),
    )
    output_path = tmp_path / "out" / "sample" / "hST-SVC" / "SVC.h5ad"
    output_path.parent.mkdir(parents=True)
    output_path.write_bytes(b"previous-valid-result")
    ctx = PublicationContext(
        svc=SVC(
            expr=output,
            spatial=output,
            svc_kind="sp",
            artifacts={"outputs": {"sp_svc": output}},
        ),
        run_dir=tmp_path / "run",
        merged_config={"ot": {"ga": {"solver": "pot"}, "lr": {"solver": "pot"}}},
        provenance={},
        record_artifact=lambda artifact: None,
    )
    args = SimpleNamespace(
        svc_type="hST-SVC",
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
        service._build_public_result(args, "application_sp", "sp_svc", ctx)

    assert output_path.read_bytes() == b"previous-valid-result"
    assert list(output_path.parent.iterdir()) == [output_path]


def test_public_result_manifest_failure_restores_previous_result(tmp_path):
    from revise.application import service

    output = AnnData(
        X=np.ones((1, 1)),
        obs=pd.DataFrame(index=["cell-1"]),
        var=pd.DataFrame(index=["g1"]),
    )
    output_path = tmp_path / "out" / "sample" / "hST-SVC" / "SVC.h5ad"
    output_path.parent.mkdir(parents=True)
    output_path.write_bytes(b"previous-valid-result")
    ctx = PublicationContext(
        svc=SVC(
            expr=output,
            spatial=output,
            svc_kind="sp",
            artifacts={"outputs": {"sp_svc": output}},
        ),
        run_dir=tmp_path / "run",
        merged_config={"ot": {"ga": {"solver": "pot"}, "lr": {"solver": "pot"}}},
        provenance={},
        artifact_records=[],
    )

    def fail_record_artifact(artifact):
        ctx.artifact_records.append(artifact)
        raise OSError("simulated manifest failure")

    ctx.record_artifact = fail_record_artifact
    args = SimpleNamespace(
        svc_type="hST-SVC",
        output_root=str(tmp_path / "out"),
        sample_name="sample",
        seed=17,
        ot_method="pot",
    )

    with pytest.raises(OSError, match="simulated manifest failure"):
        service._build_public_result(args, "application_sp", "sp_svc", ctx)

    assert output_path.read_bytes() == b"previous-valid-result"
    assert "result" not in ctx.provenance
    assert ctx.artifact_records == []
    assert list(output_path.parent.iterdir()) == [output_path]


def test_public_result_rejects_route_type_mismatch_before_publishing(tmp_path):
    from revise.application import service

    output = AnnData(
        X=np.ones((1, 1)),
        obs=pd.DataFrame(index=["cell-1"]),
        var=pd.DataFrame(index=["g1"]),
    )
    output_root = tmp_path / "out"
    ctx = PublicationContext(
        svc=SVC(
            expr=output,
            spatial=output,
            svc_kind="sc",
            artifacts={"outputs": {"sp_svc": output}},
        ),
        run_dir=tmp_path / "run",
        merged_config={"ot": {"ga": {"solver": "pot"}, "lr": {"solver": "pot"}}},
        provenance={},
        record_artifact=lambda artifact: None,
    )
    args = SimpleNamespace(
        svc_type="hST-SVC",
        output_root=str(output_root),
        sample_name="sample",
        seed=17,
        ot_method="pot",
    )

    with pytest.raises(
        ValueError,
        match="SVC type 'hST-SVC' requires internal kind 'sp'",
    ):
        service._build_public_result(args, "application_sp", "sp_svc", ctx)

    assert not output_root.exists()
