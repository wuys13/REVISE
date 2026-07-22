from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, read_h5ad

from revise.backend.ops.meta import get_sc_obs
from revise.preprocess.histology_priors import build_histology_prior_from_tables


CELL_LOCATIONS_KEY = "revise_cell_locations"


def test_histology_prior_persists_segmented_centers_and_falls_back_to_spot_centers(
    tmp_path,
):
    st = AnnData(
        X=np.array([[10.0], [20.0]]),
        obs=pd.DataFrame(index=["spot-a", "spot-b"]),
        var=pd.DataFrame(index=["gene"]),
    )
    st.obsm["spatial"] = np.array([[0, 0], [10, 0]])
    cells = pd.DataFrame(
        {"cell_id": ["segmented-1"], "x": [0.25], "y": [0.5]}
    )
    spots = pd.DataFrame(
        {
            "spot_id": ["spot-a", "spot-b"],
            "x": [0.0, 10.0],
            "y": [0.0, 0.0],
        }
    )

    build_histology_prior_from_tables(st, cells, spots, spot_radius=1.0)

    output = tmp_path / "with-centers.h5ad"
    st.write_h5ad(output)
    restored = read_h5ad(output)
    locations = restored.uns[CELL_LOCATIONS_KEY]
    assert locations.index.tolist() == ["segmented-1"]
    assert locations.loc["segmented-1"].to_dict() == {
        "spot_name": "spot-a",
        "x": 0.25,
        "y": 0.5,
    }

    svc_obs = get_sc_obs(
        restored.obs_names,
        restored.uns["all_cells_in_spot"],
        restored.obsm["spatial"],
        cell_locations=locations,
    ).set_index("cell_id")
    fallback_id = restored.uns["all_cells_in_spot"]["spot-b"][0]
    assert svc_obs.loc["segmented-1", ["x", "y"]].tolist() == [0.25, 0.5]
    assert svc_obs.loc[fallback_id, ["x", "y"]].tolist() == [10.0, 0.0]
    assert svc_obs[["x", "y"]].dtypes.tolist() == [np.float64, np.float64]


def test_histology_prior_rejects_spot_coordinates_from_another_frame():
    st = AnnData(
        X=np.array([[10.0]]),
        obs=pd.DataFrame(index=["spot-a"]),
        var=pd.DataFrame(index=["gene"]),
    )
    st.obsm["spatial"] = np.array([[0.0, 0.0]])
    cells = pd.DataFrame(
        {"cell_id": ["segmented-1"], "x": [100.25], "y": [100.5]}
    )
    spots = pd.DataFrame(
        {"spot_id": ["spot-a"], "x": [100.0], "y": [100.0]}
    )

    with pytest.raises(ValueError, match="one coordinate frame"):
        build_histology_prior_from_tables(st, cells, spots, spot_radius=1.0)


def test_histology_prior_coordinate_check_supports_obs_xy_fallback():
    st = AnnData(
        X=np.array([[10.0]]),
        obs=pd.DataFrame({"x": [2.0], "y": [3.0]}, index=["spot-a"]),
        var=pd.DataFrame(index=["gene"]),
    )
    cells = pd.DataFrame(
        {"cell_id": ["segmented-1"], "x": [2.25], "y": [3.25]}
    )
    spots = pd.DataFrame(
        {"spot_id": ["spot-a"], "x": [2.0], "y": [3.0]}
    )

    build_histology_prior_from_tables(st, cells, spots, spot_radius=1.0)

    assert st.uns[CELL_LOCATIONS_KEY].index.tolist() == ["segmented-1"]
