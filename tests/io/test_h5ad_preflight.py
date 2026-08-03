from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, read_h5ad

from revise.backend.policies import ModeValidationPolicy
from revise.config import load_raw_config, merge_unified_config
from revise.config.runner_conf import InputSpec
from revise.config.runner_conf import resolve_input_specs
from revise.config.runner_conf import resolved_input_path
from revise.framework import REVISEPipeline
from revise.io.input_service import REVISEInputService
from revise.recon.context import PipelineContext
from revise.utils import input_identities
from revise.utils.spot_sr_input import ensure_all_cells_in_spot


COLUMNS = {
    "cell_type_col": "Level1",
    "sub_cell_type_col": "Level2",
    "confidence_col": "Confidence",
    "unknown_key": "Unknown",
}

DECLARED_PROFILES = {
    "application_sp",
    "application_sc",
    "application_sc_hyper",
    "application_sc_sr",
    "benchmark_seg",
    "benchmark_bin2cell",
    "benchmark_sr_batch",
    "benchmark_sr_spot_size",
    "benchmark_impute_panel",
    "benchmark_impute_dropout",
}


def _adata(role: str, *, genes=("g1", "g2")) -> AnnData:
    obs = pd.DataFrame(index=pd.Index([f"{role}-1", f"{role}-2"]))
    if role == "sc_ref":
        obs["Level1"] = ["A", "B"]
        obs["Level2"] = ["A1", "B1"]
    if role == "gt":
        obs["cell_id"] = obs.index.astype(str)
        obs["Level1"] = ["A", "B"]
        obs["x"] = [0.0, 1.0]
        obs["y"] = [0.0, 1.0]
    adata = AnnData(
        X=np.ones((len(obs), len(genes)), dtype=np.float64),
        obs=obs,
        var=pd.DataFrame(index=pd.Index(genes)),
    )
    if role in {"st", "gt"}:
        adata.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 1.0]])
    return adata


def _write(path: Path, role: str, *, genes=("g1", "g2")) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    _adata(role, genes=genes).write_h5ad(path)
    return path


def _application_specs(tmp_path: Path):
    runtime = {"mode": "application", "task": "sp_svc"}
    io = {
        "data_root": str(tmp_path),
        "sample_name": "sample",
        "st_file": "st.h5ad",
        "sc_ref_file": "sc.h5ad",
    }
    specs = resolve_input_specs(runtime, io)
    return runtime, io, specs


def _write_application_inputs(tmp_path: Path):
    _write(tmp_path / "sample_st.h5ad", "st")
    _write(tmp_path / "sc.h5ad", "sc_ref")


@pytest.mark.parametrize(
    ("profile", "expected"),
    [
        (
            "application_sp",
            [("st", "sample_Xenium.h5ad"), ("sc_ref", "adata_sc_all_reanno.h5ad")],
        ),
        (
            "application_sc",
            [("st", "sample_Xenium.h5ad"), ("sc_ref", "adata_sc_all_reanno.h5ad")],
        ),
        (
            "application_sc_hyper",
            [("st", "sample_Xenium.h5ad"), ("sc_ref", "adata_sc_all_reanno.h5ad")],
        ),
        (
            "application_sc_sr",
            [("st", "sample_Xenium.h5ad"), ("sc_ref", "adata_sc_all_reanno.h5ad")],
        ),
        (
            "benchmark_seg",
            [
                ("st", "sample/seg_1/Xenium.h5ad"),
                ("sc_ref", "sample/adata_sc_all_reanno.h5ad"),
                ("gt", "sample/selected_xenium.h5ad"),
            ],
        ),
        (
            "benchmark_bin2cell",
            [
                ("st", "sample/bin2cell/Xenium.h5ad"),
                ("sc_ref", "sample/adata_sc_all_reanno.h5ad"),
                ("gt", "sample/selected_xenium.h5ad"),
            ],
        ),
        (
            "benchmark_sr_batch",
            [
                ("st", "sample/spot_50/xenium_spot.h5ad"),
                ("sc_ref", "sample/real_sc_ref_all.h5ad"),
                ("gt", "sample/selected_xenium.h5ad"),
            ],
        ),
        (
            "benchmark_sr_spot_size",
            [
                ("st", "sample/spot_50/xenium_spot.h5ad"),
                ("sc_ref", "sample/real_sc_ref_all.h5ad"),
                ("gt", "sample/selected_xenium.h5ad"),
            ],
        ),
        (
            "benchmark_impute_panel",
            [
                ("st", "sample/selected_xenium.h5ad"),
                ("sc_ref", "sample/real_sc_ref.h5ad"),
                ("gt", "sample/selected_xenium.h5ad"),
            ],
        ),
        (
            "benchmark_impute_dropout",
            [
                ("st", "sample/selected_xenium.h5ad"),
                ("sc_ref", "sample/real_sc_ref.h5ad"),
                ("gt", "sample/selected_xenium.h5ad"),
            ],
        ),
    ],
)
def test_all_declared_profiles_resolve_one_explicit_input_matrix(
    tmp_path,
    profile,
    expected,
):
    raw = load_raw_config("revise/revise.yaml")
    merged = merge_unified_config(
        raw_config=raw,
        profile=profile,
        runtime_overrides={},
        io_overrides={"data_root": str(tmp_path), "sample_name": "sample"},
        algorithm_overrides={},
    )
    runtime = merged["runtime"]

    specs = resolve_input_specs(runtime, merged["io"])

    resolved = [
        (spec.role, Path(spec.path).relative_to(tmp_path).as_posix())
        for spec in specs
    ]
    assert resolved == expected


def test_input_matrix_covers_every_declared_profile():
    assert set(load_raw_config("revise/revise.yaml")["profiles"]) == DECLARED_PROFILES


def test_full_loader_prefers_the_preflight_resolved_path():
    specs = (InputSpec("st", "/shared/st.h5ad"),)

    assert resolved_input_path(specs, "st", "/derived-again/st.h5ad") == "/shared/st.h5ad"


def test_spatialdata_source_path_replaces_the_h5ad_fallback(tmp_path):
    runtime, io, _ = _application_specs(tmp_path)
    io.update(
        input_format="auto",
        spatialdata_path=str(tmp_path / "sample.zarr"),
    )

    specs = resolve_input_specs(runtime, io)

    assert specs[0] == InputSpec("st", str(tmp_path / "sample.zarr"))


def _preflight(tmp_path: Path):
    runtime, io, specs = _application_specs(tmp_path)
    report = REVISEInputService(io).preflight(
        specs,
        runtime=runtime,
        columns=COLUMNS,
    )
    return specs, report


def _write_benchmark_inputs(
    tmp_path: Path,
    task: str,
    *,
    sample_name: str = "sample",
):
    runtime = {"mode": "benchmark", "task": task}
    io = {
        "data_root": str(tmp_path),
        "sample_name": sample_name,
        "st_file": "st.h5ad",
        "sc_ref_file": "sc.h5ad",
        "gt_svc_file": "gt.h5ad",
        "seg_method": "seg_1",
        "spot_size": 50,
    }
    specs = resolve_input_specs(runtime, io)
    paths = {spec.role: Path(spec.path) for spec in specs}
    st = _adata("st")
    if task == "sp_svc":
        st.obs["seg_error"] = [0, 0]
    if task == "sc_svc_impute":
        st.obs["transcript_counts"] = [2.0, 2.0]
    sc_ref = _adata("sc_ref")
    gt = _adata("gt")
    if task in {"sp_svc", "sc_svc_impute"}:
        gt.obs_names = st.obs_names.copy()
    for path, adata in (
        (paths["st"], st),
        (paths["sc_ref"], sc_ref),
        (paths["gt"], gt),
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(path)
    return runtime, io, specs, paths


def test_valid_preflight_reports_resolved_roles_and_proof_boundary(tmp_path):
    _write_application_inputs(tmp_path)

    specs, report = _preflight(tmp_path)

    assert [(spec.role, spec.path) for spec in specs] == [
        ("st", str(tmp_path / "sample_st.h5ad")),
        ("sc_ref", str(tmp_path / "sc.h5ad")),
    ]
    assert [item["role"] for item in report["inputs"]] == ["st", "sc_ref"]
    assert all(item["backed"] is True for item in report["inputs"])
    assert all(item["format"] == "h5ad" for item in report["inputs"])
    assert all(item["shape"] == [2, 2] for item in report["inputs"])
    assert report["gene_overlap"] == 2
    assert report["proof_boundary"] == (
        "metadata_and_required_arrays_only; expression_values_not_fully_scanned"
    )


@pytest.mark.parametrize("task", ["sp_svc", "sc_svc_sr", "sc_svc_impute"])
def test_valid_benchmark_preflight_requires_ground_truth(tmp_path, task):
    runtime, io, specs, _ = _write_benchmark_inputs(tmp_path, task)

    report = REVISEInputService(io).preflight(
        specs,
        runtime=runtime,
        columns=COLUMNS,
    )

    assert [item["role"] for item in report["inputs"]] == ["st", "sc_ref", "gt"]


def test_benchmark_preflight_matches_patient_to_sample_path_prefix(tmp_path):
    runtime, io, specs, paths = _write_benchmark_inputs(
        tmp_path,
        "sp_svc",
        sample_name="P2CRC/cut_part1",
    )
    io["patient_key"] = "Patient"
    sc_ref = read_h5ad(paths["sc_ref"])
    sc_ref.obs["Patient"] = ["P2CRC"] * sc_ref.n_obs
    sc_ref.write_h5ad(paths["sc_ref"])

    report = REVISEInputService(io).preflight(
        specs,
        runtime=runtime,
        columns=COLUMNS,
    )

    assert report["status"] == "ready"
    assert io["sample_name"] == "P2CRC/cut_part1"


def test_benchmark_preflight_rejects_wrong_patient_for_sample_path(tmp_path):
    runtime, io, specs, paths = _write_benchmark_inputs(
        tmp_path,
        "sp_svc",
        sample_name="P2CRC/cut_part1",
    )
    io["patient_key"] = "Patient"
    sc_ref = read_h5ad(paths["sc_ref"])
    sc_ref.obs["Patient"] = ["P3CRC"] * sc_ref.n_obs
    sc_ref.write_h5ad(paths["sc_ref"])

    with pytest.raises(ValueError, match=r"sample 'P2CRC'"):
        REVISEInputService(io).preflight(
            specs,
            runtime=runtime,
            columns=COLUMNS,
        )


def test_application_preflight_keeps_exact_patient_match(tmp_path):
    runtime, io, _ = _application_specs(tmp_path)
    io["sample_name"] = "P2CRC/cut_part1"
    io["patient_key"] = "Patient"
    specs = resolve_input_specs(runtime, io)
    _write(Path(specs[0].path), "st")
    sc_ref_path = tmp_path / "sc.h5ad"
    sc_ref = _adata("sc_ref")
    sc_ref.obs["Patient"] = ["P2CRC"] * sc_ref.n_obs
    sc_ref.write_h5ad(sc_ref_path)

    with pytest.raises(ValueError, match=r"sample 'P2CRC/cut_part1'"):
        REVISEInputService(io).preflight(
            specs,
            runtime=runtime,
            columns=COLUMNS,
        )


def test_batch_effect_preflight_allows_reference_without_target_patient(tmp_path):
    runtime, io, specs, paths = _write_benchmark_inputs(
        tmp_path,
        "sc_svc_sr",
        sample_name="P2CRC/cut_part1",
    )
    runtime["confounding"] = "batch_effect"
    io["patient_key"] = "Patient"
    sc_ref = read_h5ad(paths["sc_ref"])
    sc_ref.obs["Patient"] = ["P3CRC"] * sc_ref.n_obs
    sc_ref.write_h5ad(paths["sc_ref"])

    report = REVISEInputService(io).preflight(
        specs,
        runtime=runtime,
        columns=COLUMNS,
    )

    assert report["status"] == "ready"


def test_impute_preflight_does_not_require_unused_spatial_coordinates(tmp_path):
    runtime, io, specs, paths = _write_benchmark_inputs(
        tmp_path,
        "sc_svc_impute",
    )
    st = read_h5ad(paths["st"])
    del st.obsm["spatial"]
    st.write_h5ad(paths["st"])

    report = REVISEInputService(io).preflight(
        specs,
        runtime=runtime,
        columns=COLUMNS,
    )
    assert report["status"] == "ready"


def test_benchmark_preflight_reports_missing_ground_truth_role(tmp_path):
    runtime, io, specs, paths = _write_benchmark_inputs(tmp_path, "sp_svc")
    paths["gt"].unlink()

    with pytest.raises(FileNotFoundError, match=r"role=gt.*gt\.h5ad"):
        REVISEInputService(io).preflight(
            specs,
            runtime=runtime,
            columns=COLUMNS,
        )


@pytest.mark.parametrize("task", ["sp_svc", "sc_svc_impute"])
def test_benchmark_preflight_reports_ground_truth_alignment_context(tmp_path, task):
    runtime, io, specs, paths = _write_benchmark_inputs(tmp_path, task)
    gt = read_h5ad(paths["gt"])
    gt.obs_names = ["other-1", "other-2"]
    gt.write_h5ad(paths["gt"])

    with pytest.raises(
        ValueError,
        match=r"role=gt.*gt\.h5ad.*obs_names_overlap.*actual=0.*st_path=",
    ):
        REVISEInputService(io).preflight(
            specs,
            runtime=runtime,
            columns=COLUMNS,
        )


@pytest.mark.parametrize(
    ("task", "role", "field"),
    [
        ("sp_svc", "st", "seg_error"),
        ("sc_svc_impute", "st", "transcript_counts"),
        ("sc_svc_sr", "gt", "cell_id"),
        ("sc_svc_sr", "gt", "x"),
        ("sc_svc_sr", "gt", "y"),
    ],
)
def test_benchmark_preflight_rejects_route_specific_missing_fields(
    tmp_path,
    task,
    role,
    field,
):
    runtime, io, specs, paths = _write_benchmark_inputs(tmp_path, task)
    adata = read_h5ad(paths[role])
    del adata.obs[field]
    adata.write_h5ad(paths[role])

    with pytest.raises(KeyError, match=rf"role={role}.*{field}"):
        REVISEInputService(io).preflight(
            specs,
            runtime=runtime,
            columns=COLUMNS,
        )


def test_sr_ground_truth_requires_configured_broad_column(tmp_path):
    runtime, io, specs, paths = _write_benchmark_inputs(tmp_path, "sc_svc_sr")
    gt = read_h5ad(paths["gt"])
    del gt.obs["Level1"]
    gt.write_h5ad(paths["gt"])

    with pytest.raises(KeyError, match=r"role=gt.*Level1"):
        REVISEInputService(io).preflight(
            specs,
            runtime=runtime,
            columns=COLUMNS,
        )


def test_sr_ground_truth_default_accepts_historical_clusters_and_reports_source(
    tmp_path,
):
    runtime, io, specs, paths = _write_benchmark_inputs(tmp_path, "sc_svc_sr")
    gt = read_h5ad(paths["gt"])
    gt.obs["clusters"] = gt.obs.pop("Level1")
    gt.write_h5ad(paths["gt"])

    report = REVISEInputService(io).preflight(
        specs,
        runtime=runtime,
        columns=COLUMNS,
    )

    gt_report = next(item for item in report["inputs"] if item["role"] == "gt")
    assert gt_report["ground_truth_label_source"] == "clusters"


def test_sr_ground_truth_custom_broad_key_does_not_fallback_to_clusters(tmp_path):
    runtime, io, specs, paths = _write_benchmark_inputs(tmp_path, "sc_svc_sr")
    sc_ref = read_h5ad(paths["sc_ref"])
    sc_ref.obs["major_type"] = sc_ref.obs.pop("Level1")
    sc_ref.write_h5ad(paths["sc_ref"])
    gt = read_h5ad(paths["gt"])
    gt.obs["clusters"] = gt.obs.pop("Level1")
    gt.write_h5ad(paths["gt"])

    with pytest.raises(KeyError, match=r"role=gt.*major_type"):
        REVISEInputService(io).preflight(
            specs,
            runtime=runtime,
            columns={**COLUMNS, "cell_type_col": "major_type"},
        )


@pytest.mark.parametrize("problem", ["duplicate_cell_id", "nonfinite_x"])
def test_sr_ground_truth_rejects_invalid_cell_identity_or_coordinates(
    tmp_path,
    problem,
):
    runtime, io, specs, paths = _write_benchmark_inputs(tmp_path, "sc_svc_sr")
    gt = read_h5ad(paths["gt"])
    if problem == "duplicate_cell_id":
        gt.obs["cell_id"] = ["duplicate", "duplicate"]
        expected = "cell_id.*unique"
    else:
        gt.obs["x"] = [np.nan, 1.0]
        expected = "x.*y.*finite"
    gt.write_h5ad(paths["gt"])

    with pytest.raises(ValueError, match=expected):
        REVISEInputService(io).preflight(
            specs,
            runtime=runtime,
            columns=COLUMNS,
        )


def test_sr_pre_allocation_requires_unique_virtual_cell_ids(tmp_path):
    runtime, io, specs, paths = _write_benchmark_inputs(tmp_path, "sc_svc_sr")
    st = read_h5ad(paths["st"])
    st.uns["all_cells_in_spot"] = {
        str(st.obs_names[0]): ["unknown-cell"],
        str(st.obs_names[1]): ["unknown-cell"],
    }
    st.write_h5ad(paths["st"])

    report = REVISEInputService(io).preflight(
        specs,
        runtime=runtime,
        columns=COLUMNS,
    )

    st_report = next(item for item in report["inputs"] if item["role"] == "st")
    assert st_report["sr_mapping"] == {
        "source": "embedded",
        "validation": "pre_allocation",
    }
    with pytest.raises(
        ValueError,
        match=r"all_cells_in_spot.*unique cell ids",
    ):
        ensure_all_cells_in_spot(st, real_adata=read_h5ad(paths["gt"]))


def test_sr_pre_allocation_rejects_mapping_ids_absent_from_ground_truth(tmp_path):
    runtime, io, specs, paths = _write_benchmark_inputs(tmp_path, "sc_svc_sr")
    st = read_h5ad(paths["st"])
    st.uns["all_cells_in_spot"] = {
        str(st.obs_names[0]): ["unknown-1"],
        str(st.obs_names[1]): ["unknown-2"],
    }
    st.write_h5ad(paths["st"])

    report = REVISEInputService(io).preflight(
        specs,
        runtime=runtime,
        columns=COLUMNS,
    )

    assert report["status"] == "ready"
    with pytest.raises(ValueError, match=r"actual_unknown.*unknown-1"):
        ensure_all_cells_in_spot(st, real_adata=read_h5ad(paths["gt"]))


def test_sr_fallback_mapping_uses_ground_truth_cell_id_column():
    st = _adata("st")
    gt = _adata("gt")
    gt.obs_names = ["row-1", "row-2"]
    gt.obs["cell_id"] = ["cell-1", "cell-2"]

    ensure_all_cells_in_spot(st, real_adata=gt)

    mapped = [cell for cells in st.uns["all_cells_in_spot"].values() for cell in cells]
    assert mapped == ["cell-1", "cell-2"]


@pytest.mark.parametrize(
    ("mapping", "message"),
    [
        ({"st-1": ["cell-1"]}, "missing mappings.*active spots"),
        (
            {"st-1": ["cell-1"], "st-2": []},
            "empty cell lists.*active spots",
        ),
    ],
)
def test_sr_pre_allocation_requires_mapping_for_each_nonempty_active_spot(
    mapping,
    message,
):
    st = _adata("st")
    st.uns["all_cells_in_spot"] = mapping

    with pytest.raises((KeyError, ValueError), match=message):
        ensure_all_cells_in_spot(st)


def test_preflight_missing_file_has_role_and_path_context(tmp_path):
    _write(tmp_path / "sc.h5ad", "sc_ref")
    _, _, specs = _application_specs(tmp_path)

    with pytest.raises(FileNotFoundError, match=r"role=st.*sample_st\.h5ad"):
        REVISEInputService().preflight(
            specs,
            runtime={"mode": "application", "task": "sp_svc"},
            columns=COLUMNS,
        )


def test_preflight_corrupt_file_has_role_and_path_context(tmp_path):
    (tmp_path / "sample_st.h5ad").write_text("not hdf5", encoding="utf-8")
    _write(tmp_path / "sc.h5ad", "sc_ref")

    with pytest.raises(ValueError, match=r"role=st.*sample_st\.h5ad.*open"):
        _preflight(tmp_path)


@pytest.mark.parametrize("axis", ["obs", "var"])
def test_preflight_rejects_duplicate_axes(tmp_path, axis):
    _write_application_inputs(tmp_path)
    path = tmp_path / "sample_st.h5ad"
    adata = _adata("st")
    if axis == "obs":
        adata.obs_names = ["duplicate", "duplicate"]
    else:
        adata.var_names = ["duplicate", "duplicate"]
    adata.write_h5ad(path)

    with pytest.raises(ValueError, match=rf"role=st.*{axis}_names.*unique"):
        _preflight(tmp_path)


def test_preflight_rejects_empty_input(tmp_path):
    _write_application_inputs(tmp_path)
    empty = AnnData(
        X=np.empty((0, 2)),
        obs=pd.DataFrame(index=pd.Index([], dtype=str)),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    empty.obsm["spatial"] = np.empty((0, 2))
    empty.write_h5ad(tmp_path / "sample_st.h5ad")

    with pytest.raises(ValueError, match=r"role=st.*shape.*nonempty"):
        _preflight(tmp_path)


@pytest.mark.parametrize("problem", ["missing", "shape", "nonfinite"])
def test_preflight_rejects_invalid_spatial_coordinates(tmp_path, problem):
    _write_application_inputs(tmp_path)
    adata = _adata("st")
    if problem == "missing":
        del adata.obsm["spatial"]
    elif problem == "shape":
        adata.obsm["spatial"] = np.ones((adata.n_obs, 1))
    else:
        adata.obsm["spatial"][0, 0] = np.nan
    adata.write_h5ad(tmp_path / "sample_st.h5ad")

    with pytest.raises((KeyError, ValueError), match=r"role=st.*spatial"):
        _preflight(tmp_path)


def test_preflight_rejects_missing_reference_label(tmp_path):
    _write_application_inputs(tmp_path)
    adata = _adata("sc_ref")
    del adata.obs["Level1"]
    adata.write_h5ad(tmp_path / "sc.h5ad")

    with pytest.raises(KeyError, match=r"role=sc_ref.*Level1"):
        _preflight(tmp_path)


def test_preflight_rejects_reference_labels_that_collide_after_normalization(
    tmp_path,
):
    _write_application_inputs(tmp_path)
    adata = _adata("sc_ref")
    adata.obs["Level1"] = ["A/B", "A_B"]
    adata.write_h5ad(tmp_path / "sc.h5ad")

    with pytest.raises(ValueError, match=r"role=sc_ref.*Level1.*collide"):
        _preflight(tmp_path)


@pytest.mark.parametrize("task", ["sc_svc"])
def test_application_sc_routes_require_the_labels_used_by_full_runners(
    tmp_path,
    task,
):
    _write_application_inputs(tmp_path)
    sc_ref = _adata("sc_ref")
    del sc_ref.obs["Level2"]
    sc_ref.write_h5ad(tmp_path / "sc.h5ad")
    _, io, specs = _application_specs(tmp_path)

    with pytest.raises(KeyError, match=r"role=sc_ref.*Level2"):
        REVISEInputService(io).preflight(
            specs,
            runtime={
                "mode": "application",
                "task": task,
                "platform": task,
            },
            columns=COLUMNS,
        )


def test_sc_svc_custom_columns_are_the_full_runner_contract(tmp_path):
    _write_application_inputs(tmp_path)
    sc_ref = _adata("sc_ref")
    sc_ref.obs["custom_level1"] = sc_ref.obs.pop("Level1")
    sc_ref.obs["custom_level2"] = sc_ref.obs.pop("Level2")
    sc_ref.write_h5ad(tmp_path / "sc.h5ad")
    _, io, specs = _application_specs(tmp_path)

    report = REVISEInputService(io).preflight(
        specs,
        runtime={"mode": "application", "task": "sc_svc", "platform": "sc_svc"},
        columns={
            **COLUMNS,
            "cell_type_col": "custom_level1",
            "sub_cell_type_col": "custom_level2",
        },
    )

    assert report["status"] == "ready"


@pytest.mark.parametrize("mode", ["application", "benchmark"])
def test_sr_preflight_accepts_configured_broad_without_subtype_or_literal_level1(
    tmp_path,
    mode,
):
    columns = {
        **COLUMNS,
        "cell_type_col": "major_type",
        "sub_cell_type_col": "minor_type",
    }
    if mode == "application":
        _write_application_inputs(tmp_path)
        _, io, specs = _application_specs(tmp_path)
        sc_path = tmp_path / "sc.h5ad"
    else:
        runtime, io, specs, paths = _write_benchmark_inputs(
            tmp_path,
            "sc_svc_sr",
        )
        sc_path = paths["sc_ref"]
        gt = read_h5ad(paths["gt"])
        gt.obs["major_type"] = gt.obs.pop("Level1")
        gt.write_h5ad(paths["gt"])
    sc_ref = _adata("sc_ref")
    sc_ref.obs["major_type"] = sc_ref.obs.pop("Level1")
    del sc_ref.obs["Level2"]
    sc_ref.write_h5ad(sc_path)

    report = REVISEInputService(io).preflight(
        specs,
        runtime={
            "mode": mode,
            "task": "sc_svc_sr",
            "platform": "sc_svc_sr" if mode == "application" else "sim2real",
        },
        columns=columns,
    )

    assert report["status"] == "ready"
    if mode == "benchmark":
        gt_report = next(
            item for item in report["inputs"] if item["role"] == "gt"
        )
        assert gt_report["ground_truth_label_source"] == "major_type"


@pytest.mark.parametrize("mode", ["application", "benchmark"])
def test_sr_preflight_rejects_missing_configured_broad_column(tmp_path, mode):
    columns = {
        **COLUMNS,
        "cell_type_col": "major_type",
        "sub_cell_type_col": "minor_type",
    }
    if mode == "application":
        _write_application_inputs(tmp_path)
        _, io, specs = _application_specs(tmp_path)
    else:
        _, io, specs, _ = _write_benchmark_inputs(tmp_path, "sc_svc_sr")

    with pytest.raises(KeyError, match=r"role=sc_ref.*major_type"):
        REVISEInputService(io).preflight(
            specs,
            runtime={
                "mode": mode,
                "task": "sc_svc_sr",
                "platform": (
                    "sc_svc_sr" if mode == "application" else "sim2real"
                ),
            },
            columns=columns,
        )


def test_application_sr_pre_allocation_rejects_mismatched_cell_locations(tmp_path):
    _write_application_inputs(tmp_path)
    st = read_h5ad(tmp_path / "sample_st.h5ad")
    st.uns["all_cells_in_spot"] = {
        str(st.obs_names[0]): ["cell-1"],
        str(st.obs_names[1]): ["cell-2"],
    }
    st.uns["revise_cell_locations"] = pd.DataFrame(
        {"spot_name": [str(st.obs_names[1])], "x": [0.1], "y": [0.2]},
        index=pd.Index(["cell-1"], name="cell_id"),
    )
    st.write_h5ad(tmp_path / "sample_st.h5ad")
    _, io, specs = _application_specs(tmp_path)

    report = REVISEInputService(io).preflight(
        specs,
        runtime={"mode": "application", "task": "sc_svc_sr"},
        columns=COLUMNS,
    )

    assert report["status"] == "ready"
    with pytest.raises(ValueError, match="spot_name disagrees"):
        ensure_all_cells_in_spot(st)


@pytest.mark.parametrize(
    ("mode", "expected_source"),
    [
        ("application", "generated_from_transcript_counts"),
        ("benchmark", "generated_from_ground_truth_coordinates"),
    ],
)
def test_sr_preflight_records_documented_mapping_generation_path(
    tmp_path,
    mode,
    expected_source,
):
    if mode == "application":
        _write_application_inputs(tmp_path)
        _, io, specs = _application_specs(tmp_path)
    else:
        _, io, specs, _ = _write_benchmark_inputs(tmp_path, "sc_svc_sr")

    report = REVISEInputService(io).preflight(
        specs,
        runtime={"mode": mode, "task": "sc_svc_sr"},
        columns=COLUMNS,
    )

    st_report = next(item for item in report["inputs"] if item["role"] == "st")
    assert st_report["sr_mapping"] == {
        "source": expected_source,
        "validation": "pre_allocation",
    }


def test_preflight_rejects_zero_gene_overlap(tmp_path):
    _write(tmp_path / "sample_st.h5ad", "st", genes=("st-gene",))
    _write(tmp_path / "sc.h5ad", "sc_ref", genes=("sc-gene",))

    with pytest.raises(ValueError, match=r"gene overlap.*sample_st.*sc\.h5ad"):
        _preflight(tmp_path)


def test_backed_handles_close_after_success_and_failure(monkeypatch, tmp_path):
    import revise.io.input_service as input_service

    _write_application_inputs(tmp_path)
    opened = []
    original = input_service.read_h5ad

    def tracked_read(*args, **kwargs):
        adata = original(*args, **kwargs)
        opened.append(adata)
        return adata

    monkeypatch.setattr(input_service, "read_h5ad", tracked_read)
    _preflight(tmp_path)
    assert opened and all(not adata.file.is_open for adata in opened)

    opened.clear()
    bad_sc = _adata("sc_ref")
    bad_sc.var_names = ["other-1", "other-2"]
    bad_sc.write_h5ad(tmp_path / "sc.h5ad")
    with pytest.raises(ValueError, match="gene overlap"):
        _preflight(tmp_path)
    assert opened and all(not adata.file.is_open for adata in opened)


def _policy_context(
    tmp_path: Path,
    *,
    ga_solver="pot",
    lr_solver="pot",
    mode="application",
    task="sp_svc",
):
    runtime, io, _ = _application_specs(tmp_path)
    runtime = {**runtime, "mode": mode, "task": task}
    merged = {
        "runtime": runtime,
        "io": {**io, "output_root": str(tmp_path / "output")},
        "columns": COLUMNS,
        "ot": {
            "ga": {"solver": ga_solver},
            "lr": {"solver": lr_solver},
        },
    }
    if task in {"sp_svc", "sc_svc_sr"}:
        merged["local_refinement"] = {
            "strength": 0.2 if task == "sp_svc" else 0.0,
        }
    return PipelineContext(
        merged_config=merged,
        raw_config={},
        config_path="revise/revise.yaml",
        profile="test",
        runtime=runtime,
        route_key="sp_svc:bin2cell",
        run_dir=tmp_path / "run",
        logger=logging.getLogger("test-h5ad-preflight"),
        dry_run=True,
    )


@pytest.mark.parametrize("phase", ["ga", "lr"])
def test_preflight_requires_tacco_before_reconstruction(monkeypatch, tmp_path, phase):
    import revise.backend.policies as policies

    _write_application_inputs(tmp_path)
    solvers = {"ga_solver": "pot", "lr_solver": "pot"}
    solvers[f"{phase}_solver"] = "tacco"
    ctx = _policy_context(tmp_path, **solvers)

    def missing_tacco():
        raise ModuleNotFoundError("TACCO is missing", name="tacco")

    monkeypatch.setattr(policies, "require_tacco", missing_tacco)
    with pytest.raises(ModuleNotFoundError, match="TACCO is missing"):
        ModeValidationPolicy().validate(ctx)


def test_preflight_accepts_tacco_local_refinement_solver(monkeypatch, tmp_path):
    import revise.backend.policies as policies

    _write_application_inputs(tmp_path)
    ctx = _policy_context(
        tmp_path,
        lr_solver="tacco",
    )
    monkeypatch.setattr(policies, "require_tacco", lambda: None)

    ModeValidationPolicy().validate(ctx)

    assert (ctx.run_dir / "preflight.json").exists()

def test_preflight_report_is_persisted_without_scientific_outputs(tmp_path):
    from revise.framework import REVISEPipeline

    _write_application_inputs(tmp_path)
    output_root = tmp_path / "output"

    svc = REVISEPipeline().run(
        profile="application_sp",
        io_overrides={
            "data_root": str(tmp_path),
            "output_root": str(output_root),
            "sample_name": "sample",
            "st_file": "st.h5ad",
            "sc_ref_file": "sc.h5ad",
        },
        dry_run=True,
    )

    run_dir = Path(svc.provenance["run_dir"])
    report = json.loads((run_dir / "preflight.json").read_text())
    assert report["status"] == "ready"
    assert {path.name for path in run_dir.iterdir()} == {
        "merged_config.json",
        "preflight.json",
        "provenance.json",
        "run.log",
    }
    assert not list(run_dir.rglob("*.h5ad"))


def test_pipeline_preflight_failure_persists_terminal_run_truth(tmp_path):
    output_root = tmp_path / "failed-output"

    with pytest.raises(FileNotFoundError, match=r"role=st.*missing-case"):
        REVISEPipeline().run(
            profile="application_sp",
            io_overrides={
                "data_root": str(tmp_path),
                "output_root": str(output_root),
                "sample_name": "missing-case",
            },
            dry_run=True,
        )

    manifest_path = next(output_root.rglob("provenance.json"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["run"]["status"] == "failed"
    assert manifest["stages"][0]["status"] == "failed"
    assert all(
        stage["status"] == "skipped" and stage["reason"] == "upstream_failure"
        for stage in manifest["stages"][1:]
    )
    assert manifest["artifacts"] == []
    expected_inputs = [
        InputSpec("st", str(tmp_path / "missing-case_Xenium.h5ad")),
        InputSpec("sc_ref", str(tmp_path / "adata_sc_all_reanno.h5ad")),
    ]
    assert manifest["input_identities"] == input_identities(expected_inputs)
    assert "data_fingerprint" not in manifest
