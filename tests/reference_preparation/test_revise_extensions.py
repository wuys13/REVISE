from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml

from revise.reference_preparation import (
    assert_reference_unchanged,
    prepare_reference,
    resolve_reference_input,
)
from revise.reference_preparation.config import Candidate, ConfigError
from revise.reference_preparation.evidence import REPORT_KIND
from revise.reference_preparation.ga_contract import (
    GAResponse,
    GlobalAnchoringResult,
    validate_global_anchoring_result,
)
from revise.reference_preparation.screening import screen_candidates
from revise.reference_preparation.__main__ import main as cli_main


def _write_yaml(path: Path, payload) -> Path:
    path.write_text(yaml.safe_dump(payload), encoding="utf-8")
    return path


def _adata(path: Path) -> Path:
    ad.AnnData(
        X=np.array([[1.0, 2.0], [3.0, 4.0]]),
        obs=pd.DataFrame({"pair": ["A", "B"]}, index=["c1", "c2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    ).write_h5ad(path)
    return path


def _response(
    probability: float = 0.8,
    *,
    ids=("u1", "u2"),
    labels=("A", "B"),
    metadata=None,
) -> GAResponse:
    matrix = np.array([
        [probability, 1 - probability],
        [probability, 1 - probability],
    ])
    return GAResponse(
        GlobalAnchoringResult(matrix, list(ids), list(labels), metadata or {}),
        list(ids),
    )


def _candidate(tmp_path: Path, name: str) -> Candidate:
    path = tmp_path / f"{name}.h5ad"
    path.write_bytes(name.encode("utf-8"))
    return Candidate(name, path)


def test_callback_metadata_is_sanitized_to_controlled_json_fields():
    runtime = object()
    original = {
        "backend": "test",
        "effective_seed": 7,
        "cell_type_labels": ["A", "B"],
        "private_runtime": runtime,
    }
    validated = validate_global_anchoring_result(
        _response(metadata=original).result, ["u1", "u2"]
    )
    assert validated.metadata == {
        "backend": "test",
        "effective_seed": 7,
        "cell_type_labels": ["A", "B"],
        "normalization": "adapter_declared_probability",
    }
    assert original["private_runtime"] is runtime
    with pytest.raises(ValueError, match="effective_seed"):
        validate_global_anchoring_result(
            _response(metadata={"effective_seed": {"nested": True}}).result,
            ["u1", "u2"],
        )


def test_screen_requires_one_st_axis_and_one_execution_context(tmp_path):
    candidates = [_candidate(tmp_path, name) for name in ("baseline", "axis", "input")]

    def runner(path, _):
        metadata = {
            "reconstruction_config_sha256": "config",
            "st_input_sha256": "st" if path.stem != "input" else "other-st",
            "effective_parameters_sha256": "parameters",
            "effective_seed": 3,
            "effective_solver": "solver",
        }
        ids = ("u1", "u2") if path.stem != "axis" else ("u2", "u1")
        return _response(ids=ids, metadata=metadata)

    report = screen_candidates(candidates, tmp_path / "base.yaml", runner)
    assert report["candidates"][0]["status"] == "success"
    assert report["candidates"][1]["status"] == "failed"
    assert "first valid candidate" in report["candidates"][1]["failure_reason"]
    assert report["candidates"][2]["status"] == "failed"
    assert "st_input_sha256" in report["candidates"][2]["failure_reason"]


def test_screen_records_allowed_label_and_gene_differences(tmp_path):
    candidates = [_candidate(tmp_path, name) for name in ("a", "b")]

    def runner(path, _):
        if path.stem == "a":
            return _response(labels=("A", "B"), metadata={
                "scoring_genes_sha256": "genes-a", "scoring_n_vars": 10,
            })
        return _response(0.9, labels=("A", "C"), metadata={
            "scoring_genes_sha256": "genes-b", "scoring_n_vars": 11,
        })

    report = screen_candidates(candidates, tmp_path / "base.yaml", runner)
    second = report["candidates"][1]
    assert second["status"] == "success"
    assert second["comparability"] == {
        "st_unit_ids": "equal",
        "cell_type_labels": "different",
        "scoring_genes": "different",
        "execution_context": "equal",
    }


def test_selected_reference_hash_must_match_callback_metadata(tmp_path):
    candidate = _candidate(tmp_path, "chosen")
    with pytest.raises(RuntimeError, match="does not match"):
        screen_candidates(
            [candidate], tmp_path / "base.yaml",
            lambda *_: _response(metadata={"reference_sha256": "wrong"}),
        )


def test_paired_report_has_source_output_and_selected_digests(tmp_path):
    source = _adata(tmp_path / "pool.h5ad")
    config = _write_yaml(tmp_path / "prepare.yaml", {
        "schema_version": 1,
        "mode": "paired",
        "source": "pool.h5ad",
        "pair_column": "pair",
        "pair_key": "A",
        "output_dir": "run",
    })
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    result = prepare_reference(config)
    report = json.loads(result.report_path.read_text(encoding="utf-8"))
    output_hash = hashlib.sha256(result.reference_path.read_bytes()).hexdigest()
    assert report["report_kind"] == REPORT_KIND
    assert report["source_sha256"] == source_hash
    assert report["paired_output_sha256"] == output_hash
    assert report["selected_reference_sha256"] == output_hash


def test_resolve_verified_prepared_reference_and_detect_later_change(tmp_path):
    _adata(tmp_path / "pool.h5ad")
    config = _write_yaml(tmp_path / "prepare.yaml", {
        "schema_version": 1, "mode": "paired", "source": "pool.h5ad",
        "pair_column": "pair", "pair_key": "A", "output_dir": "run",
    })
    result = prepare_reference(config)
    resolved, evidence = resolve_reference_input(result.reference_config_path)
    assert resolved == result.reference_path
    assert evidence == {
        "schema_version": 1,
        "origin": REPORT_KIND,
        "verification": "verified_report",
        "config_path": str(result.reference_config_path),
        "config_sha256": hashlib.sha256(result.reference_config_path.read_bytes()).hexdigest(),
        "report_path": str(result.report_path),
        "report_sha256": hashlib.sha256(result.report_path.read_bytes()).hexdigest(),
        "reference_path": str(result.reference_path),
        "reference_sha256": hashlib.sha256(result.reference_path.read_bytes()).hexdigest(),
        "reference_size_bytes": result.reference_path.stat().st_size,
    }
    assert_reference_unchanged(evidence)
    result.reference_path.write_bytes(result.reference_path.read_bytes() + b"changed")
    with pytest.raises(ConfigError, match="SHA-256"):
        assert_reference_unchanged(evidence)


def test_manual_reference_without_report_is_explicitly_unverified(tmp_path):
    reference = tmp_path / "manual.h5ad"
    reference.write_bytes(b"manual")
    config = _write_yaml(tmp_path / "reference.yaml", {
        "schema_version": 1,
        "reference": {"path": "manual.h5ad", "format": "h5ad"},
    })
    resolved, evidence = resolve_reference_input(config)
    assert resolved == reference
    assert evidence["origin"] == "reference_config_no_report"
    assert evidence["verification"] == "unverified_no_report"
    assert evidence["report_sha256"] is None
    assert_reference_unchanged(evidence)


def test_present_malformed_or_failed_report_is_not_treated_as_manual(tmp_path):
    reference = tmp_path / "manual.h5ad"
    reference.write_bytes(b"manual")
    config = _write_yaml(tmp_path / "reference.yaml", {
        "schema_version": 1,
        "reference": {"path": "manual.h5ad", "format": "h5ad"},
    })
    (tmp_path / "report.json").write_text("{}", encoding="utf-8")
    with pytest.raises(ConfigError, match="Unrecognized"):
        resolve_reference_input(config)
    (tmp_path / "report.json").write_text(json.dumps({
        "schema_version": 1,
        "report_kind": REPORT_KIND,
        "status": "failed",
    }), encoding="utf-8")
    with pytest.raises(ConfigError, match="not successful"):
        resolve_reference_input(config)


def test_cli_paired_mode_does_not_import_host_backend(tmp_path, monkeypatch, capsys):
    _adata(tmp_path / "pool.h5ad")
    config = _write_yaml(tmp_path / "prepare.yaml", {
        "schema_version": 1, "mode": "paired", "source": "pool.h5ad",
        "pair_column": "pair", "pair_key": "A", "output_dir": "run",
    })
    monkeypatch.setitem(sys.modules, "revise.reference_preparation.host", None)
    assert cli_main(["--config", str(config)]) == 0
    result = json.loads(capsys.readouterr().out)
    assert Path(result["reference_path"]).is_file()
    assert Path(result["reference_config_path"]).is_file()
    assert Path(result["report_path"]).is_file()
