"""Fixture contracts for the batch reconstruction-impact adapter."""

from dataclasses import dataclass
import json
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from revise.batch.analysis import AnalysisContext


def test_hst_parent_cannot_overwrite_global_scope():
    from revise.analysis.reconstruction_impact_adapter import _parent_specs

    with pytest.raises(ValueError, match="reserved"):
        _parent_specs({"modality": "hST"}, {"parent_labels": {"global": "T"}})


def _write_carrier(path: Path, ids, labels, *, clusters=None, offset=0):
    obs = pd.DataFrame({"Level1": labels}, index=ids)
    if clusters is not None:
        obs["SVC_cluster"] = clusters
    data = ad.AnnData(
        X=np.asarray([[1 + offset, 2 + offset, 3 + offset]] * len(ids), dtype=float),
        obs=obs,
        var=pd.DataFrame(index=["g0", "g1", "g2"]),
    )
    data.obsm["spatial"] = np.asarray([[float(i), float(i)] for i in range(len(ids))])
    data.write_h5ad(path)
    return path


def _reference(path: Path):
    reference = ad.AnnData(
        X=np.ones((4, 3), dtype=float),
        obs=pd.DataFrame(
            {
                "Level1": ["Fibroblast", "Fibroblast", "Mono/Macro", "Mono/Macro"],
                "Level2": ["F1", "F2", "M1", "M2"],
            },
            index=["r0", "r1", "r2", "r3"],
        ),
        var=pd.DataFrame(index=["g0", "g1", "g2"]),
    )
    reference.write_h5ad(path)
    return path


def _context(tmp_path: Path, *, modality: str, mapping: str | None, raw, spatial,
             parameters: dict, resources: dict[str, Path]):
    outputs = {"svc": {"path": str(spatial)}} if modality == "hST" else {
        "spatial": {"path": str(spatial)},
    }
    reconstruction = {
        "status": "succeeded",
        "sample_id": "fixture",
        "cell_type": "Mono_Macro" if modality == "iST" else None,
        "modality": modality,
        "ist_mapping": mapping,
        "inputs": {"sources": {"spatial": {"path": str(raw)}}},
        "outputs": outputs,
        "pairing": {"status": "available"},
        "coordinates": {"unit": "pixel", "microns_per_coordinate": 1.0},
        "fingerprint": "fixture-reconstruction",
    }
    return AnalysisContext(
        reconstruction,
        tmp_path / "stage",
        parameters,
        resources,
    )


def _parameters(*, route_kind, parent_labels, sample_n_units=None, parent_sample_n_units=None,
                mode="fixed_within_level1"):
    return {
        "route_kind": route_kind,
        "level1_column": "Level1",
        "parent_labels": parent_labels,
        "partition_change": {
            "mode": mode,
            "level1_resolution_candidates": [0.3],
            "within_level1_resolution": 0.5,
            "random_state": 13,
            "n_top_genes": 2,
            "raw_qc_min_genes": 1,
            "raw_qc_min_cells": 1,
            "sample_n_units": sample_n_units,
            "parent_sample_n_units": parent_sample_n_units,
        },
        "spatial_region": {
            "microns_per_coordinate": 1.0,
            "candidate_window_sides_um": [10.0],
            "min_parent_units": 1,
            "rarefaction_draws": 1,
            "threshold_bootstraps": 1,
            "cell_equivalent_um": 1.0,
        },
        "raw_level2_mapping": {
            "method": "pot",
            "level1_column": "Level1",
            "level2_column": "Level2",
            "pot": {"reg": 0.1, "reg_m": 0.0, "reg_type": "entropy"},
        },
    }


@dataclass
class _FakeMapping:
    labels: pd.Series
    posterior: pd.DataFrame
    assignments: pd.DataFrame
    audit: dict


def _fake_partition(calls):
    def run(raw, reconstructed, **kwargs):
        ids = raw.obs_names
        edge = "raw_to_final_svc" if kwargs["route_kind"] == "sc_svc" else "raw_to_recon_expression"
        assignments = pd.DataFrame(
            {
                "unit_id": ids.astype(str),
                "raw_cluster": ["r0"] * len(ids),
                "matched_raw_cluster": ["r0"] * len(ids),
                "recon_cluster": ["c0"] * len(ids),
                "unit_changed": [False] * len(ids),
                "comparison_edge": edge,
            },
            index=ids,
        )
        comparison = SimpleNamespace(assignments=assignments)
        calls.append({"raw": raw.copy(), "reconstructed": reconstructed.copy(), "kwargs": kwargs})
        return SimpleNamespace(
            comparisons={edge: comparison},
            feature_names=["g0", "g1"],
            audit={"n_units": int(len(ids)), "excluded_raw_qc_units": 0},
            matched_cluster_status="ok",
        )

    return run


def _fake_mapping(calls):
    def map_labels(raw_parent, reference, **kwargs):
        calls.append({"raw": raw_parent.copy(), "reference": reference.copy(), "kwargs": kwargs})
        labels = pd.Series(["L2"] * raw_parent.n_obs, index=raw_parent.obs_names, name="raw_level2")
        posterior = pd.DataFrame({"L2": [1.0] * raw_parent.n_obs}, index=raw_parent.obs_names)
        assignments = pd.DataFrame({"unit_id": labels.index, "raw_level2": labels}, index=labels.index)
        return _FakeMapping(labels, posterior, assignments, {"method": kwargs["method"]})

    return map_labels


def test_hst_adapter_samples_global_and_parents_from_covered_ids(monkeypatch, tmp_path):
    from revise.analysis import reconstruction_impact_adapter as adapter

    ids = [f"u{i}" for i in range(6)]
    labels = ["Fibroblast"] * 3 + ["T"] * 3
    raw = _write_carrier(tmp_path / "raw.h5ad", ids, labels)
    spatial = _write_carrier(tmp_path / "spatial.h5ad", ids, ["stale"] * 6)
    reference = _reference(tmp_path / "reference.h5ad")
    partition_calls = []
    mapping_calls = []
    writer_calls = []
    monkeypatch.setattr(adapter, "run_partition_analysis", _fake_partition(partition_calls))
    monkeypatch.setattr(adapter, "map_raw_level2_labels", _fake_mapping(mapping_calls))
    monkeypatch.setattr(
        adapter,
        "_write_outputs",
        lambda output_dir, bundle: writer_calls.append(bundle) or {
            "audit": {"path": "audit.json", "description": "fixture"}
        },
    )

    context = _context(
        tmp_path,
        modality="hST",
        mapping=None,
        raw=raw,
        spatial=spatial,
        parameters=_parameters(
            route_kind="sp_svc",
            parent_labels={"Fibroblast": "Fibroblast", "T": "T"},
            sample_n_units=4,
            parent_sample_n_units=2,
            mode="level1_ari",
        ),
        resources={"raw_level2_reference": reference},
    )

    result = adapter.run(context)

    assert result["artifacts"]["audit"]["path"] == "audit.json"
    assert list(writer_calls[0]["scopes"]) == ["global", "Fibroblast", "T"]
    assert len(partition_calls) == 3
    assert len(mapping_calls) == 2
    assert partition_calls[0]["raw"].n_obs == 4
    assert partition_calls[1]["raw"].n_obs == partition_calls[2]["raw"].n_obs == 2
    assert partition_calls[1]["raw"].obs_names.tolist() != partition_calls[2]["raw"].obs_names.tolist()
    for call in partition_calls:
        assert call["reconstructed"].obs["Level1"].tolist() == call["raw"].obs["Level1"].tolist()
        assert call["kwargs"]["random_state"] == 13
    assert partition_calls[0]["kwargs"]["resolution_mode"] == "level1_ari"
    assert partition_calls[1]["kwargs"]["resolution_mode"] == "fixed_within_level1"
    assert writer_calls[0]["anatomy"].scale_audit["origin_x_um"] == 0.0
    assert writer_calls[0]["scopes"]["global"]["audit"]["sampling_seed"] == 13
    assert writer_calls[0]["scopes"]["global"]["spatial"] is None
    assert writer_calls[0]["scopes"]["global"]["raw_level2"] is None
    assert writer_calls[0]["scopes"]["Fibroblast"]["audit"]["sampling_seed"] != writer_calls[0]["scopes"]["T"]["audit"]["sampling_seed"]
    assert writer_calls[0]["resources"]["raw_level2_reference_identity"]["sha256"]


def test_hst_global_keeps_non_parent_covered_labels_out_of_parent_level2_mapping(monkeypatch, tmp_path):
    from revise.analysis import reconstruction_impact_adapter as adapter

    ids = [f"u{i}" for i in range(5)]
    raw = _write_carrier(
        tmp_path / "raw.h5ad",
        ids,
        ["Fibroblast", "Fibroblast", "Fibroblast", "Tumor", "Intestinal Epithelial"],
    )
    spatial = _write_carrier(tmp_path / "spatial.h5ad", ids, ["stale"] * len(ids))
    reference = _reference(tmp_path / "reference.h5ad")
    partition_calls = []
    mapping_calls = []
    writer_calls = []
    monkeypatch.setattr(adapter, "run_partition_analysis", _fake_partition(partition_calls))
    monkeypatch.setattr(adapter, "map_raw_level2_labels", _fake_mapping(mapping_calls))
    monkeypatch.setattr(adapter, "compute_anatomy_regions", lambda **_kwargs: SimpleNamespace(scale_audit={}))
    monkeypatch.setattr(adapter, "compute_spatial_impact", lambda **_kwargs: SimpleNamespace())
    monkeypatch.setattr(adapter, "_write_outputs", lambda _output_dir, bundle: writer_calls.append(bundle) or {"audit": {"path": "audit.json"}})

    adapter.run(
        _context(
            tmp_path,
            modality="hST",
            mapping=None,
            raw=raw,
            spatial=spatial,
            parameters=_parameters(
                route_kind="sp_svc",
                parent_labels={"Fibroblast": "Fibroblast"},
            ),
            resources={"raw_level2_reference": reference},
        )
    )

    assert [call["raw"].obs_names.tolist() for call in partition_calls] == [ids, ids[:3]]
    assert [call["raw"].obs_names.tolist() for call in mapping_calls] == [ids[:3]]
    global_scope = writer_calls[0]["scopes"]["global"]
    assert global_scope["spatial"] is None
    assert global_scope["raw_level2"] is None
    assert global_scope["audit"]["spatial_status"] == "not_applicable_parent_baseline"
    assert global_scope["audit"]["raw_level2_status"] == "not_applicable_parent_baseline"


def test_hst_uses_full_raw_only_for_anatomy_when_carrier_is_partial(monkeypatch, tmp_path):
    from revise.analysis import reconstruction_impact_adapter as adapter

    raw_ids = [f"u{i}" for i in range(5)]
    covered_ids = raw_ids[:3]
    raw = _write_carrier(
        tmp_path / "raw.h5ad",
        raw_ids,
        ["Fibroblast", "Fibroblast", "Fibroblast", "Tumor", "Intestinal Epithelial"],
    )
    spatial = _write_carrier(tmp_path / "spatial.h5ad", covered_ids, ["stale"] * len(covered_ids))
    reference = _reference(tmp_path / "reference.h5ad")
    anatomy_calls = []
    spatial_calls = []
    monkeypatch.setattr(adapter, "run_partition_analysis", _fake_partition([]))
    monkeypatch.setattr(adapter, "map_raw_level2_labels", _fake_mapping([]))

    def fake_anatomy(**kwargs):
        anatomy_calls.append(kwargs)
        return SimpleNamespace(scale_audit={"origin_x_um": 0.0})

    def fake_spatial(**kwargs):
        spatial_calls.append(kwargs)
        return SimpleNamespace()

    monkeypatch.setattr(adapter, "compute_anatomy_regions", fake_anatomy)
    monkeypatch.setattr(adapter, "compute_spatial_impact", fake_spatial)
    monkeypatch.setattr(adapter, "_write_outputs", lambda _output_dir, _bundle: {"audit": {"path": "audit.json"}})

    adapter.run(
        _context(
            tmp_path,
            modality="hST",
            mapping=None,
            raw=raw,
            spatial=spatial,
            parameters=_parameters(
                route_kind="sp_svc",
                parent_labels={"Fibroblast": "Fibroblast"},
            ),
            resources={"raw_level2_reference": reference},
        )
    )

    assert anatomy_calls[0]["full_coordinates"].index.tolist() == raw_ids
    assert anatomy_calls[0]["full_level1_labels"].tolist()[-2:] == ["Tumor", "Intestinal Epithelial"]
    assert spatial_calls[0]["full_coordinates"].index.tolist() == raw_ids
    assert spatial_calls[0]["paired_coordinates"].index.tolist() == covered_ids
    assert spatial_calls[0]["raw_labels"].tolist() == ["r0", "r0", "r0"]


def test_ist_adapter_uses_exact_parent_and_required_level2_resource(monkeypatch, tmp_path):
    from revise.analysis import reconstruction_impact_adapter as adapter

    ids = ["u0", "u1", "u2", "u3"]
    raw = _write_carrier(tmp_path / "raw.h5ad", ids, ["Mono/Macro", "Mono/Macro", "T", "T"])
    spatial = _write_carrier(tmp_path / "spatial.h5ad", ["u0", "u1"], ["stale", "stale"], clusters=["a", "b"])
    reference = _reference(tmp_path / "reference.h5ad")
    partition_calls = []
    mapping_calls = []
    writer_calls = []
    monkeypatch.setattr(adapter, "run_partition_analysis", _fake_partition(partition_calls))
    monkeypatch.setattr(adapter, "map_raw_level2_labels", _fake_mapping(mapping_calls))
    monkeypatch.setattr(
        adapter,
        "_write_outputs",
        lambda output_dir, bundle: writer_calls.append(bundle) or {
            "audit": {"path": "audit.json", "description": "fixture"}
        },
    )

    context = _context(
        tmp_path,
        modality="iST",
        mapping="paired",
        raw=raw,
        spatial=spatial,
        parameters=_parameters(
            route_kind="sc_svc",
            parent_labels={"Mono_Macro": "Mono/Macro"},
            parent_sample_n_units=None,
        ),
        resources={"raw_level2_reference": reference},
    )
    result = adapter.run(context)

    assert result["artifacts"]["audit"]["path"] == "audit.json"
    assert list(writer_calls[0]["scopes"]) == ["Mono_Macro"]
    assert partition_calls[0]["raw"].obs["Level1"].tolist() == ["Mono/Macro", "Mono/Macro"]
    assert partition_calls[0]["reconstructed"].obs["Level1"].tolist() == ["Mono/Macro", "Mono/Macro"]
    assert mapping_calls[0]["kwargs"]["parent_value"] == "Mono/Macro"
    assert mapping_calls[0]["kwargs"]["method"] == "pot"


def test_adapter_requires_raw_level2_reference(tmp_path):
    from revise.analysis import reconstruction_impact_adapter as adapter

    raw = _write_carrier(tmp_path / "raw.h5ad", ["u0", "u1"], ["T", "T"])
    spatial = _write_carrier(tmp_path / "spatial.h5ad", ["u0", "u1"], ["stale", "stale"])
    context = _context(
        tmp_path,
        modality="hST",
        mapping=None,
        raw=raw,
        spatial=spatial,
        parameters=_parameters(route_kind="sp_svc", parent_labels={"T": "T"}),
        resources={},
    )

    with pytest.raises(ValueError, match="raw_level2_reference"):
        adapter.run(context)


@pytest.mark.parametrize("value", [3.9, "3", True])
def test_sampling_limit_requires_exact_positive_integer(value):
    from revise.analysis import reconstruction_impact_adapter as adapter

    with pytest.raises(ValueError, match="positive integer"):
        adapter._positive_limit(value, "sample_n_units")


@pytest.mark.parametrize("value", [3.9, "3", True])
def test_seed_and_spatial_counts_require_exact_integer(value):
    from revise.analysis import reconstruction_impact_adapter as adapter

    with pytest.raises(ValueError, match="explicit integer"):
        adapter._strict_integer(value, "partition_change.random_state")


def test_parent_labels_accept_only_exact_string_values():
    from revise.analysis import reconstruction_impact_adapter as adapter

    with pytest.raises(ValueError, match="exact raw label"):
        adapter._parent_specs(
            {"modality": "hST"},
            {"parent_labels": {"Fibroblast": {"label": "Fibroblast"}}},
        )


def test_writer_bridge_constructs_impact_scope_objects(monkeypatch, tmp_path):
    from revise.analysis import impact_outputs
    from revise.analysis import reconstruction_impact_adapter as adapter

    captured = {}

    def fake_writer(output_dir, *, scopes, audit, calculation):
        captured["output_dir"] = output_dir
        captured["scopes"] = list(scopes)
        captured["audit"] = audit
        captured["calculation"] = calculation
        return {"artifacts": {"audit": {"path": "audit.json"}}, "calculation": calculation}

    monkeypatch.setattr(impact_outputs, "write_impact_outputs", fake_writer)
    anatomy = SimpleNamespace(scale_audit={"origin_x_um": 0.0})
    bundle = {
        "level1_col": "Level1",
        "route_kind": "sp_svc",
        "common_inputs": {
            "record_identity": {"sample_id": "fixture", "cell_type": None},
        },
        "parameters": {"parent_labels": {"T": "T"}},
        "resources": {"raw_level2_reference": str(tmp_path / "reference.h5ad")},
        "anatomy": anatomy,
        "scopes": {
            "global": {
                "audit": {"partition_edge": "raw_to_recon_expression"},
                "partition": object(),
                "spatial": object(),
                "raw_level2": object(),
            }
        },
    }

    result = adapter._write_outputs(tmp_path, bundle)

    assert result["artifacts"]["audit"]["path"] == "audit.json"
    assert len(captured["scopes"]) == 1
    assert isinstance(captured["scopes"][0], impact_outputs.ImpactScope)
    assert json.loads(captured["scopes"][0].comparison_id) == {
        "sample_id": "fixture", "scope": "global", "task_cell_type": None,
    }
    assert captured["scopes"][0].comparison["comparison_edge"] == "raw_to_recon_expression"
