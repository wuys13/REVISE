"""Analysis input views are read-only, explicit and lazy."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse


def _write_adata(path, matrix, *, obs_names, var_names, obs=None, spatial=None):
    data = ad.AnnData(matrix, obs=pd.DataFrame(obs or {}, index=obs_names),
                      var=pd.DataFrame(index=var_names))
    if spatial is not None:
        data.obsm["spatial"] = np.asarray(spatial, dtype=float)
    data.write_h5ad(path)
    return data


def _record(tmp_path, *, modality="iST", mapping="paired", output_roles=None):
    raw_path = tmp_path / "raw.h5ad"
    _write_adata(
        raw_path,
        sparse.csr_matrix([[10, 1, 0], [2, 4, 6], [5, 0, 3]], dtype=np.int32),
        obs_names=["spot-a", "spot-b", "spot-c"],
        var_names=["raw-1", "raw-2", "raw-3"],
        spatial=[[1, 1], [2, 2], [3, 3]],
    )
    outputs = {}
    for role, data in (output_roles or {}).items():
        path = tmp_path / f"{role}.h5ad"
        data.write_h5ad(path)
        outputs[role] = {"path": str(path)}
    return {
        "status": "succeeded",
        "modality": modality,
        "ist_mapping": mapping,
        "inputs": {"sources": {"spatial": {"path": str(raw_path)}}},
        "outputs": outputs,
        "pairing": {"status": "available" if mapping == "paired" else "unavailable"},
        "coordinates": {"unit": "um", "microns_per_coordinate": 1.0},
    }


def _paired_record(tmp_path):
    spatial = _write_adata(
        tmp_path / "unused-spatial.h5ad",
        sparse.csr_matrix([[0, 0], [0, 0]], dtype=np.float32),
        obs_names=["spot-b", "spot-a"],
        var_names=["g1", "g2"],
        obs={"SVC_cluster": ["b", "a"]},
        spatial=[[2, 2], [1, 1]],
    )
    donor = _write_adata(
        tmp_path / "unused-expression.h5ad",
        sparse.csr_matrix([[2, 0], [4, 6], [8, 10]], dtype=np.float32),
        obs_names=["donor-2", "donor-1", "donor-3"],
        var_names=["g1", "g2"],
        obs={"SVC_cluster": ["a", "a", "b"]},
    )
    return _record(tmp_path, output_roles={"spatial": spatial, "expression": donor})


def _sst_record(tmp_path):
    output = _write_adata(
        tmp_path / "unused-sst.h5ad",
        sparse.csr_matrix([[99, 98, 97], [96, 95, 94], [93, 92, 91]], dtype=np.int32),
        obs_names=["cell-2", "cell-0", "cell-1"],
        var_names=["generated-1", "generated-2", "generated-3"],
        obs={"spot_name": ["spot-b", "spot-a", "spot-b"]},
        spatial=[[22, 22], [11, 11], [21, 21]],
    )
    return _record(tmp_path, modality="sST", mapping=None, output_roles={"svc": output})


def test_analysis_inputs_are_lazy_until_a_view_is_requested(tmp_path, monkeypatch):
    from revise.batch import inputs as module

    record = _paired_record(tmp_path)
    calls = []
    original = module.ad.read_h5ad

    def read(path, *args, **kwargs):
        calls.append(str(path))
        return original(path, *args, **kwargs)

    monkeypatch.setattr(module.ad, "read_h5ad", read)
    views = module.AnalysisInputs.from_reconstruction(record)
    assert calls == []
    views.raw()
    assert calls == [record["inputs"]["sources"]["spatial"]["path"]]


def test_new_delivery_reads_published_raw_and_legacy_delivery_reads_source(tmp_path):
    from revise.batch.inputs import AnalysisInputs

    record = _paired_record(tmp_path)
    source_path = record["inputs"]["sources"]["spatial"]["path"]
    assert AnalysisInputs.from_reconstruction(record).raw().source == source_path

    published = _write_adata(
        tmp_path / "published-raw.h5ad",
        sparse.csr_matrix([[7, 8]], dtype=np.int32),
        obs_names=["published-spot"],
        var_names=["g1", "g2"],
        spatial=[[4, 5]],
    )
    record["delivery_protocol_version"] = 2
    record["outputs"]["raw"] = {"path": str(tmp_path / "published-raw.h5ad")}
    view = AnalysisInputs.from_reconstruction(record).raw()
    assert view.source == record["outputs"]["raw"]["path"]
    assert list(view.observation_ids) == ["published-spot"]
    np.testing.assert_array_equal(view.adata.X.toarray(), published.X.toarray())


def test_new_delivery_never_falls_back_to_mutable_source_when_raw_is_missing(tmp_path):
    from revise.batch.inputs import AnalysisInputs

    record = _paired_record(tmp_path)
    record["delivery_protocol_version"] = 2

    with pytest.raises(ValueError, match="missing published Raw"):
        AnalysisInputs.from_reconstruction(record)


def test_raw_alignment_uses_ids_and_coordinates_and_records_uncovered_rows(tmp_path):
    from revise.batch.inputs import AnalysisInputs

    record = _paired_record(tmp_path)
    aligned = AnalysisInputs.from_reconstruction(record).aligned_raw()
    assert list(aligned.observation_ids) == ["spot-b", "spot-a"]
    assert list(aligned.gene_ids) == ["raw-1", "raw-2", "raw-3"]
    np.testing.assert_array_equal(aligned.coordinates, [[2, 2], [1, 1]])
    assert list(aligned.uncovered_observations) == ["spot-c"]
    assert aligned.mapping == "native_id"


def test_raw_alignment_rejects_missing_ids_or_changed_coordinates(tmp_path):
    from revise.batch.inputs import AnalysisInputs

    record = _paired_record(tmp_path)
    output_path = tmp_path / "spatial.h5ad"
    output = ad.read_h5ad(record["outputs"]["spatial"]["path"])
    output.obs_names = ["spot-b", "missing"]
    output.write_h5ad(output_path)
    with pytest.raises(ValueError, match="raw observation IDs"):
        AnalysisInputs.from_reconstruction(record).aligned_raw()

    output.obs_names = ["spot-b", "spot-a"]
    output.obsm["spatial"][0, 0] = 999
    output.write_h5ad(output_path)
    with pytest.raises(ValueError, match="coordinates"):
        AnalysisInputs.from_reconstruction(record).aligned_raw()


def test_paired_view_keeps_donor_axis_and_maps_sparse_means_to_spatial_order(tmp_path):
    from revise.batch.inputs import AnalysisInputs

    inputs = AnalysisInputs.from_reconstruction(_paired_record(tmp_path))
    with pytest.raises(ValueError, match="paired iST"):
        inputs.reconstructed()
    paired = inputs.paired()
    assert list(paired.spatial.observation_ids) == ["spot-b", "spot-a"]
    assert list(paired.donor.observation_ids) == ["donor-2", "donor-1", "donor-3"]
    means, clusters = paired.cluster_means()
    assert sparse.issparse(means)
    assert list(clusters) == ["a", "b"]
    np.testing.assert_allclose(means.toarray(), [[3, 3], [8, 10]])
    mapped = paired.map_to_spatial(means)
    assert sparse.issparse(mapped)
    np.testing.assert_allclose(mapped.toarray(), [[8, 10], [3, 3]])
    np.testing.assert_allclose(paired.map_to_spatial(np.array([0.8, 0.2])), [0.2, 0.8])


def test_paired_view_rejects_mismatched_or_null_cluster_labels(tmp_path):
    from revise.batch.inputs import AnalysisInputs

    record = _paired_record(tmp_path)
    donor = ad.read_h5ad(record["outputs"]["expression"]["path"])
    donor.obs["SVC_cluster"] = ["a", "a", "other"]
    donor.write_h5ad(record["outputs"]["expression"]["path"])
    with pytest.raises(ValueError, match="cluster sets"):
        AnalysisInputs.from_reconstruction(record).paired().cluster_means()

    donor.obs["SVC_cluster"] = ["a", None, "b"]
    donor.write_h5ad(record["outputs"]["expression"]["path"])
    with pytest.raises(ValueError, match="null"):
        AnalysisInputs.from_reconstruction(record).paired().cluster_means()


def test_mean_and_random_modes_use_svc_and_preserve_mapping_provenance(tmp_path):
    from revise.batch.inputs import AnalysisInputs

    output = _write_adata(
        tmp_path / "unused-svc.h5ad",
        sparse.csr_matrix([[1, 2], [3, 4]], dtype=np.float32),
        obs_names=["spot-a", "spot-b"],
        var_names=["g1", "g2"],
        obs={"revise_ist_donor_id": ["donor-1", "donor-2"]},
        spatial=[[1, 1], [2, 2]],
    )
    record = _record(tmp_path, mapping="random", output_roles={"svc": output})
    view = AnalysisInputs.from_reconstruction(record).reconstructed()
    assert view.mapping == "random"
    assert view.source == record["outputs"]["svc"]["path"]
    assert list(view.donor_ids) == ["donor-1", "donor-2"]
    assert view.provenance["donor_column"] == "revise_ist_donor_id"


def test_sst_baseline_equal_splits_raw_spots_in_output_order_and_preserves_coordinates(tmp_path):
    from revise.batch.inputs import AnalysisInputs

    record = _sst_record(tmp_path)
    inputs = AnalysisInputs.from_reconstruction(record)
    before = inputs.reconstructed().adata.X.copy()
    baseline = inputs.sst_baseline()
    assert baseline.mapping == "parent_spot_equal_split"
    assert sparse.issparse(baseline.adata.X)
    assert baseline.adata.X.dtype.kind == "f"
    assert list(baseline.observation_ids) == ["cell-2", "cell-0", "cell-1"]
    assert list(baseline.gene_ids) == ["raw-1", "raw-2", "raw-3"]
    np.testing.assert_array_equal(baseline.coordinates, [[22, 22], [11, 11], [21, 21]])
    np.testing.assert_allclose(
        baseline.adata.X.toarray(),
        [[1, 2, 3], [10, 1, 0], [1, 2, 3]],
    )
    assert list(baseline.covered_raw_spots) == ["spot-b", "spot-a"]
    assert list(baseline.uncovered_raw_spots) == ["spot-c"]
    np.testing.assert_array_equal(inputs.reconstructed().adata.X.toarray(), before.toarray())


def test_sst_baseline_accepts_explicit_output_obs_coordinates(tmp_path):
    from revise.batch.inputs import AnalysisInputs

    record = _sst_record(tmp_path)
    output_path = Path(record["outputs"]["svc"]["path"])
    output = ad.read_h5ad(output_path)
    coordinates = np.asarray(output.obsm["spatial"], dtype=float)
    del output.obsm["spatial"]
    output.obs["x"] = coordinates[:, 0]
    output.obs["y"] = coordinates[:, 1]
    output.write_h5ad(output_path)

    inputs = AnalysisInputs.from_reconstruction(record)
    reconstructed = inputs.reconstructed()
    np.testing.assert_array_equal(reconstructed.coordinates, coordinates)
    assert reconstructed.provenance["coordinate_source"] == "output.obs[x,y]"
    baseline = inputs.sst_baseline()

    assert list(baseline.observation_ids) == ["cell-2", "cell-0", "cell-1"]
    assert list(baseline.gene_ids) == ["raw-1", "raw-2", "raw-3"]
    assert sparse.issparse(baseline.adata.X)
    np.testing.assert_allclose(
        baseline.adata.X.toarray(),
        [[1, 2, 3], [10, 1, 0], [1, 2, 3]],
    )
    np.testing.assert_array_equal(baseline.coordinates, coordinates)
    assert baseline.provenance["coordinate_source"] == "output.obs[x,y]"
    assert "parent-spot centers" in baseline.provenance["coordinate_semantics"]


def test_sst_baseline_rejects_output_without_coordinate_fields(tmp_path):
    from revise.batch.inputs import AnalysisInputs

    record = _sst_record(tmp_path)
    output_path = Path(record["outputs"]["svc"]["path"])
    output = ad.read_h5ad(output_path)
    del output.obsm["spatial"]
    output.write_h5ad(output_path)

    with pytest.raises(ValueError, match="spatial coordinates"):
        AnalysisInputs.from_reconstruction(record).sst_baseline()


@pytest.mark.parametrize("problem", ["missing_parent", "unknown_parent", "duplicate_cell_ids"])
def test_sst_baseline_rejects_invalid_actual_output_mapping(tmp_path, problem):
    from revise.batch.inputs import AnalysisInputs

    record = _sst_record(tmp_path)
    output_path = Path(record["outputs"]["svc"]["path"])
    output = ad.read_h5ad(output_path)
    if problem == "missing_parent":
        del output.obs["spot_name"]
    elif problem == "unknown_parent":
        output.obs["spot_name"] = output.obs["spot_name"].astype(object)
        output.obs.iloc[0, output.obs.columns.get_loc("spot_name")] = "not-a-raw-spot"
    else:
        output.obs_names = ["cell-0", "cell-0", "cell-1"]
    output.write_h5ad(output_path)
    with pytest.raises(ValueError):
        AnalysisInputs.from_reconstruction(record).sst_baseline()


def test_sst_baseline_rejects_native_pairing_request(tmp_path):
    from revise.batch.inputs import AnalysisInputs

    record = _sst_record(tmp_path)
    record["pairing"] = {"status": "unavailable", "reason": "generated_spatial_units"}
    with pytest.raises(ValueError, match="spot_name"):
        AnalysisInputs.from_reconstruction(record).paired()
