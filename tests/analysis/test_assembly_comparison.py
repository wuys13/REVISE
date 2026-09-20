from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import nbformat
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from nbclient import NotebookClient

from revise.analysis.assembly_comparison import (
    DEFAULT_METHODS,
    assert_input_hashes_unchanged,
    compare_assembly_methods,
    input_summary_frame,
    load_assembly_inputs,
    normalize_broad_labels,
    spatial_plot_frame,
    validate_generated_expression,
)

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
NOTEBOOK = REPOSITORY_ROOT / "reproduce/case/assembly_comparison.ipynb"


def _fixture_objects(*, n_obs: int = 24) -> tuple[dict[str, AnnData], AnnData]:
    ids = [f"real-{index:02d}" for index in range(n_obs)]
    broad = ["T"] * n_obs
    subtype = ["T-a"] * (n_obs // 2) + ["T-b"] * (n_obs - n_obs // 2)
    genes_by_method = {
        "mean": ["g5", "g1", "g2", "g3", "g4", "mean-only"],
        "random": ["g4", "g3", "g2", "g1", "g5", "random-only"],
        "within_cluster": ["g1", "g2", "g3", "g4", "g5"],
        "outside_cluster": ["g2", "g3", "g4", "g5", "g1", "outside-only"],
    }
    spatial = np.column_stack([np.arange(n_obs), np.arange(n_obs) % 5]).astype(float)
    methods = {}
    for offset, (method, genes) in enumerate(genes_by_method.items()):
        rng = np.random.default_rng(20260919 + offset)
        x = rng.poisson(2, size=(n_obs, len(genes))).astype(float) + 1
        first = [genes.index(gene) for gene in ("g1", "g2")]
        second = [genes.index(gene) for gene in ("g4", "g5")]
        x[: n_obs // 2, first] += 9
        x[n_obs // 2 :, second] += 9
        obs = pd.DataFrame({"Level1": pd.Categorical(broad)}, index=ids)
        methods[method] = AnnData(
            x,
            obs=obs,
            var=pd.DataFrame(index=genes),
            obsm={"spatial": spatial.copy()},
        )

    baseline_obs = pd.DataFrame(
        {
            "Level1": pd.Categorical(broad),
            "SVC_cluster": pd.Categorical(subtype),
        },
        index=ids,
    )
    baseline_obs.loc[ids[-1], "SVC_cluster"] = np.nan
    baseline = AnnData(
        np.zeros((n_obs, 0), dtype=np.float32),
        obs=baseline_obs,
        obsm={"spatial": spatial.copy()},
    )
    return methods, baseline


def _write_fixture(tmp_path: Path) -> tuple[dict[str, Path], Path]:
    methods, baseline = _fixture_objects()
    method_paths = {}
    for method, adata in methods.items():
        path = tmp_path / f"{method}.h5ad"
        adata.write_h5ad(path)
        method_paths[method] = path
    baseline_path = tmp_path / "baseline.h5ad"
    baseline.write_h5ad(baseline_path)
    return method_paths, baseline_path


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_real_leiden_uses_shared_ids_and_sorted_gene_intersection_without_mutation():
    methods, baseline = _fixture_objects()
    original_matrices = {name: adata.X.copy() for name, adata in methods.items()}

    comparison = compare_assembly_methods(
        methods,
        baseline,
        cell_types=["T"],
        broad_column="Level1",
        resolutions=[0.6, 0.7],
        seed=17,
    )

    result = comparison.by_type["T"]
    assert result.status == "ok"
    assert result.shared_ids.tolist() == sorted(baseline.obs_names.tolist())
    assert result.shared_genes.tolist() == ["g1", "g2", "g3", "g4", "g5"]
    assert set(result.clustered) == set(DEFAULT_METHODS)
    assert len(comparison.metrics) == len(DEFAULT_METHODS) * 2
    assert comparison.metrics["ARI"].notna().all()
    assert comparison.metrics["NMI"].notna().all()
    assert (comparison.metrics["missing_baseline_labels"] == 1).all()
    assert all({"X_pca"} <= set(work.obsm) for work in result.clustered.values())
    assert all(
        {"connectivities", "distances"} <= set(work.obsp)
        for work in result.clustered.values()
    )
    assert all(
        work.uns["assembly_comparison"]["n_pcs"] <= min(work.n_obs - 1, work.n_vars - 1)
        for work in result.clustered.values()
    )
    assert all(
        {"leiden_0.6", "leiden_0.7"} <= set(work.obs)
        for work in result.clustered.values()
    )
    assert all(
        result.contingencies[(method, resolution)].to_numpy().sum() == len(baseline) - 1
        for method in DEFAULT_METHODS
        for resolution in (0.6, 0.7)
    )
    for method, native in methods.items():
        np.testing.assert_array_equal(native.X, original_matrices[method])
        assert not any(column.startswith("leiden_") for column in native.obs)
        assert native.obsm["spatial"].shape == (len(native), 2)

    frame = spatial_plot_frame(result, baseline, method="mean", resolution=0.6)
    assert frame.groupby("panel", observed=True).size().tolist() == [
        len(baseline),
        len(baseline),
    ]
    assert frame.loc[frame["panel"] == "baseline", "label"].isna().sum() == 1


def test_normalizes_slash_labels_preserves_na_and_does_not_change_ids_or_clusters():
    values = pd.Series(["Mono/Macro", pd.NA, "T"], dtype="object")
    normalized = normalize_broad_labels(values)
    assert normalized.tolist()[0] == "Mono_Macro"
    assert pd.isna(normalized.iloc[1])
    assert normalized.tolist()[2] == "T"

    methods, baseline = _fixture_objects()
    original_ids = baseline.obs_names.tolist()
    original_labels = baseline.obs["SVC_cluster"].copy()
    comparison = compare_assembly_methods(
        methods,
        baseline,
        cell_types=["Mono_Macro"],
        broad_column="Level1",
        resolutions=[0.6],
    )
    assert comparison.by_type["Mono_Macro"].status == "unavailable"
    assert baseline.obs_names.tolist() == original_ids
    pd.testing.assert_series_equal(baseline.obs["SVC_cluster"], original_labels)


def test_missing_type_ids_and_genes_are_explicit_and_not_zero_filled():
    methods, baseline = _fixture_objects()
    methods["random"] = methods["random"][:-2, :].copy()

    comparison = compare_assembly_methods(
        methods,
        baseline,
        cell_types=["T", "Fibroblast"],
        broad_column="Level1",
        resolutions=[0.6],
    )

    t_coverage = comparison.coverage.query("broad_type == 'T'").set_index("object")
    assert t_coverage.loc["random", "shared_ids"] == len(baseline) - 2
    assert t_coverage.loc["mean", "excluded_ids"] == ("real-22", "real-23")
    assert t_coverage.loc["mean", "excluded_genes"] == ("mean-only",)
    assert pd.isna(t_coverage.loc["baseline", "original_genes"])

    fibroblast = comparison.by_type["Fibroblast"]
    assert fibroblast.status == "unavailable"
    assert any("broad type absent" in issue for issue in fibroblast.issues)
    assert comparison.coverage.query("broad_type == 'Fibroblast'")["missing_type"].all()
    assert not (comparison.metrics["broad_type"] == "Fibroblast").any()


def test_coordinate_difference_makes_only_type_unavailable_and_is_reported():
    methods, baseline = _fixture_objects()
    methods["mean"].obsm["spatial"][0, 0] += 0.25

    comparison = compare_assembly_methods(
        methods,
        baseline,
        cell_types=["T"],
        broad_column="Level1",
        resolutions=[0.6],
    )
    result = comparison.by_type["T"]
    assert result.status == "unavailable"
    assert any("spatial coordinates differ" in issue for issue in result.issues)
    mean_row = comparison.coverage.query(
        "broad_type == 'T' and object == 'mean'"
    ).iloc[0]
    assert mean_row["spatial_status"] == "unavailable"
    assert int(mean_row["spatial_difference_count"]) == 1
    assert comparison.metrics.empty
    with pytest.raises(ValueError, match="comparison is unavailable"):
        spatial_plot_frame(result, baseline, method="mean", resolution=0.6)


def test_nonshared_coordinate_nan_is_ignored_but_shared_scope_stays_checked():
    methods, baseline = _fixture_objects()

    native = methods["mean"]
    extra_obs = native.obs.iloc[[0]].copy()
    extra_obs.index = ["method-only"]
    extra_obs["Level1"] = "T"
    methods["mean"] = AnnData(
        np.vstack([native.X, native.X[0]]),
        obs=pd.concat([native.obs, extra_obs]),
        var=native.var.copy(),
        obsm={
            "spatial": np.vstack(
                [native.obsm["spatial"], np.array([[np.nan, np.nan]])]
            )
        },
    )

    other_obs = pd.DataFrame(
        {"Level1": ["Fibroblast"], "SVC_cluster": ["other"]},
        index=["baseline-other"],
    )
    baseline = AnnData(
        np.empty((len(baseline) + 1, 0), dtype=np.float32),
        obs=pd.concat([baseline.obs, other_obs]),
        obsm={
            "spatial": np.vstack(
                [baseline.obsm["spatial"], np.array([[np.nan, np.nan]])]
            )
        },
    )

    comparison = compare_assembly_methods(
        methods,
        baseline,
        cell_types=["T"],
        broad_column="Level1",
        resolutions=[0.6],
    )
    result = comparison.by_type["T"]
    assert result.status == "ok"
    assert "method-only" not in result.shared_ids
    assert len(result.shared_ids) == 24
    assert comparison.coverage.query(
        "broad_type == 'T' and object == 'mean'"
    )["spatial_status"].eq("ok").all()


def test_shared_id_nonfinite_coordinate_blocks_type():
    methods, baseline = _fixture_objects()
    baseline.obsm["spatial"][0, 0] = np.nan

    comparison = compare_assembly_methods(
        methods,
        baseline,
        cell_types=["T"],
        broad_column="Level1",
        resolutions=[0.6],
    )
    result = comparison.by_type["T"]
    assert result.status == "unavailable"
    assert any("finite for compared IDs" in issue for issue in result.issues)
    assert comparison.coverage.query(
        "broad_type == 'T' and object == 'mean'"
    )["spatial_status"].eq("unavailable").all()


def test_old_alias_and_fallback_arguments_are_removed():
    methods, baseline = _fixture_objects()
    with pytest.raises(TypeError):
        compare_assembly_methods(
            methods,
            baseline,
            type_aliases={"T": ["T"]},
        )
    with pytest.raises(TypeError):
        compare_assembly_methods(
            methods,
            baseline,
            cell_types=["T"],
            broad_fallback="Level1",
        )


def test_invalid_expression_is_rejected():
    methods, _ = _fixture_objects()
    invalid = methods["mean"].copy()
    invalid.X[0, 0] = -1
    with pytest.raises(ValueError, match="nonnegative unlogged linear expression"):
        validate_generated_expression(invalid, method="mean")

    fractional = methods["mean"].copy()
    fractional.X = fractional.X / 2.5
    validate_generated_expression(fractional, method="fractional-mixture")


def test_h5ad_hashes_remain_unchanged_after_comparison(tmp_path: Path):
    method_paths, baseline_path = _write_fixture(tmp_path)
    loaded = load_assembly_inputs(method_paths, baseline_path)
    before = {name: _digest(path) for name, path in loaded.paths.items()}

    comparison = compare_assembly_methods(
        loaded.methods,
        loaded.baseline,
        cell_types=["T"],
        broad_column="Level1",
        resolutions=[0.6],
    )

    assert_input_hashes_unchanged(loaded)
    assert {name: _digest(path) for name, path in loaded.paths.items()} == before
    assert input_summary_frame(loaded)["sha256"].notna().all()
    assert comparison.assumptions["spatial_coordinates"].startswith("each method")


def test_notebook_executes_fixture_copy_without_overwriting_source(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    method_paths, baseline_path = _write_fixture(tmp_path)
    config = {
        "method_paths": {name: str(path) for name, path in method_paths.items()},
        "baseline_path": str(baseline_path),
        "cell_types": ["T"],
        "broad_column": "Level1",
        "baseline_subtype_column": "SVC_cluster",
        "spatial_key": "spatial",
        "coordinate_unit": "pixel",
        "microns_per_coordinate": 0.2125,
        "seed": 23,
        "resolutions": [0.6],
        "plot_resolution": 0.6,
        "output_dir": str(tmp_path / "report"),
    }
    config_path = tmp_path / "assembly-config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")
    copied_notebook = tmp_path / NOTEBOOK.name
    shutil.copy2(NOTEBOOK, copied_notebook)
    source_hash = _digest(NOTEBOOK)
    monkeypatch.setenv("REVISE_ASSEMBLY_COMPARISON_CONFIG", str(config_path))
    monkeypatch.setenv("MPLCONFIGDIR", str(tmp_path / "mpl"))
    monkeypatch.setenv("NUMBA_CACHE_DIR", str(tmp_path / "numba"))

    notebook = nbformat.read(copied_notebook, as_version=4)
    executed = NotebookClient(
        notebook,
        timeout=180,
        kernel_name="python3",
        allow_errors=False,
        resources={"metadata": {"path": str(REPOSITORY_ROOT)}},
    ).execute()

    assert _digest(NOTEBOOK) == source_hash
    assert all(
        cell.get("execution_count") is not None
        for cell in executed.cells
        if cell.cell_type == "code"
    )
    assert all(
        output.get("output_type") != "error"
        for cell in executed.cells
        if cell.cell_type == "code"
        for output in cell.get("outputs", [])
    )
    rendered_text = "\n".join(
        str(output.get("text", "")) + str(output.get("data", {}).get("text/plain", ""))
        for cell in executed.cells
        for output in cell.get("outputs", [])
    )
    assert "native input hashes unchanged" in rendered_text
    assert (tmp_path / "report" / "effective_parameters.json").exists()
    assert (tmp_path / "report" / "coverage.csv").exists()
    assert (tmp_path / "report" / "metrics.csv").exists()
    assert (tmp_path / "report" / "contingencies" / "T-mean-r0.6.csv").exists()
