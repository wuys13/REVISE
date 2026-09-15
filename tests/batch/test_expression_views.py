"""Full-gene expression views preserve spatial and reference axes."""

import numpy as np
import pytest
from scipy import sparse

from test_inputs import _paired_record, _record, _sst_record, _write_adata
from revise.batch.inputs import AnalysisInputs


def test_paired_expression_projects_full_donor_genes_on_spatial_ids(tmp_path):
    record = _paired_record(tmp_path)
    inputs = AnalysisInputs(record)
    before = inputs.paired().donor.X.copy()

    view = inputs.reconstruction_expression()

    assert view.mapping == "cluster_mean_projection"
    assert view.observation_ids.tolist() == ["spot-b", "spot-a"]
    assert view.gene_ids.tolist() == ["g1", "g2"]
    assert not set(view.gene_ids) & set(inputs.raw().gene_ids)
    assert sparse.issparse(view.X)
    np.testing.assert_allclose(view.X.toarray(), [[8, 10], [3, 3]])
    np.testing.assert_array_equal(view.coordinates, [[2, 2], [1, 1]])
    np.testing.assert_array_equal(inputs.paired().donor.X.toarray(), before.toarray())
    assert view.provenance["expression_source"] == "expression_carrier.X_as_is"
    assert view.provenance["donor_id_source"] == "reference_cell"
    assert view.provenance["donor_cells_per_cluster"] == [2, 1]
    assert view.provenance["cluster_labels"] == ["a", "b"]
    assert view.provenance["source_paths"]["expression"] == record["outputs"]["expression"]["path"]
    assert view.provenance["source_identity"]["raw"]["path"] == record["inputs"]["sources"]["spatial"]["path"]
    assert view.provenance["record_identity"]["modality"] == "iST"


def test_expression_strategy_cannot_relabel_paired_projection_as_native(tmp_path):
    inputs = AnalysisInputs(_paired_record(tmp_path))
    with pytest.raises(ValueError, match="strategy"):
        inputs.reconstruction_expression(strategy="native")


def test_paired_expression_rejects_wrong_spatial_coordinates(tmp_path):
    inputs = AnalysisInputs(_paired_record(tmp_path))
    inputs.spatial().adata.obsm["spatial"][0, 0] = 999
    with pytest.raises(ValueError, match="coordinates"):
        inputs.reconstruction_expression()


def test_native_expression_keeps_all_published_genes(tmp_path):
    output = _write_adata(
        tmp_path / "native.h5ad", sparse.csr_matrix([[1, 2, 3, 4]]),
        obs_names=["spot-a"], var_names=["a", "b", "c", "d"],
        spatial=[[1, 1]],
    )
    inputs = AnalysisInputs(_record(tmp_path, modality="hST", output_roles={"svc": output}))
    view = inputs.reconstruction_expression()
    assert view.gene_ids.tolist() == ["a", "b", "c", "d"]
    assert view.observation_ids.tolist() == ["spot-a"]
    assert view.mapping == "native_id"


def test_expression_view_rejects_empty_observation_axis(tmp_path):
    output = _write_adata(
        tmp_path / "empty.h5ad", sparse.csr_matrix((0, 2)),
        obs_names=[], var_names=["a", "b"], spatial=np.empty((0, 2)),
    )
    inputs = AnalysisInputs(_record(tmp_path, modality="hST", output_roles={"svc": output}))
    with pytest.raises(ValueError, match="nonempty"):
        inputs.reconstruction_expression()


@pytest.mark.parametrize("mapping", ["mean", "random"])
def test_native_ist_mapping_is_not_renamed_projection(tmp_path, mapping):
    obs = {"revise_ist_donor_id": ["reference-a"]} if mapping == "random" else {}
    output = _write_adata(
        tmp_path / "native.h5ad", sparse.csr_matrix([[1, 2]]),
        obs_names=["spot-a"], var_names=["a", "b"], obs=obs, spatial=[[1, 1]],
    )
    inputs = AnalysisInputs(_record(tmp_path, mapping=mapping, output_roles={"svc": output}))
    view = inputs.reconstruction_expression()
    assert view.mapping == mapping
    if mapping == "random":
        assert view.donor_ids.tolist() == ["reference-a"]
        assert view.provenance["donor_id_source"] == "reference_cell"


def test_paired_expression_rejects_typed_cluster_mismatch(tmp_path):
    record = _paired_record(tmp_path)
    inputs = AnalysisInputs(record)
    inputs.paired().donor.adata.obs["SVC_cluster"] = np.asarray([1, 1, 2], dtype=object)

    with pytest.raises(ValueError, match="cluster"):
        inputs.reconstruction_expression()


def test_paired_expression_does_not_mutate_spatial_or_donor_axes(tmp_path):
    inputs = AnalysisInputs(_paired_record(tmp_path))
    spatial = inputs.paired().spatial.adata
    donor = inputs.paired().donor.adata
    spatial_obs = spatial.obs.copy(deep=True)
    donor_obs = donor.obs.copy(deep=True)
    spatial_coords = spatial.obsm["spatial"].copy()
    donor_var = donor.var.copy(deep=True)

    inputs.reconstruction_expression()

    assert spatial.obs.equals(spatial_obs)
    assert donor.obs.equals(donor_obs)
    assert donor.var.equals(donor_var)
    np.testing.assert_array_equal(spatial.obsm["spatial"], spatial_coords)


def test_sst_expression_view_is_out_of_scope(tmp_path):
    inputs = AnalysisInputs(_sst_record(tmp_path))
    with pytest.raises(ValueError, match="sST"):
        inputs.reconstruction_expression()
