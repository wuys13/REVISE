from __future__ import annotations

import argparse
import builtins
import json
from types import SimpleNamespace


def test_legacy_sc_wrapper_emits_json_summary(monkeypatch, tmp_path):
    import application_sc_SVC_recon as legacy

    class Pipeline:
        def __init__(self, config_path):
            pass

        def run(self, **kwargs):
            return SimpleNamespace(
                summary=lambda: {"status": "ok"},
                provenance={"run_dir": str(tmp_path / "run")},
            )

    monkeypatch.setattr(legacy, "REVISEPipeline", Pipeline)
    monkeypatch.setattr(
        legacy,
        "_publish_notebook_outputs",
        lambda args, run_dir: {"sc_SVC_expr": "expr.h5ad"},
    )
    printed = []
    monkeypatch.setattr(builtins, "print", lambda value: printed.append(value))
    args = argparse.Namespace(
        config="revise/revise.yaml",
        platform="iST",
        confounding=None,
        profile=None,
        sample_name="sample",
        st_file="st.h5ad",
        data_root=str(tmp_path),
        output_root=str(tmp_path / "out"),
        sc_ref_file="sc.h5ad",
        patient_key="Patient",
        select_ct="all",
        cell_type_col="Level1",
        sub_cell_type_col="Level2",
        compatibility_mode=False,
    )

    legacy.main(args)

    summary = json.loads(printed[0])
    assert summary["status"] == "ok"
    assert summary["notebook_outputs"]["sc_SVC_expr"] == "expr.h5ad"
