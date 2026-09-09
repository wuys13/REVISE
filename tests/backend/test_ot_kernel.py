from types import SimpleNamespace

import numpy as np
import pandas as pd
from anndata import AnnData

from revise.backend import kernels
from revise.backend.kernels import ot as ot_kernel_module


def test_ot_kernel_exposes_matrix_coupling_contract():
    ot_kernel = getattr(kernels, "OTKernel", None)

    assert ot_kernel is not None
    assert callable(getattr(ot_kernel, "couple", None))


def test_ot_kernel_tacco_annotation_uses_fresh_full_inputs_and_only_publishes_assignment(
    monkeypatch,
):
    target = AnnData(
        X=np.array([[2.0, 1.0], [1.0, 2.0]]),
        obs=pd.DataFrame({"kept": [1, 2]}, index=["spot1", "spot2"]),
        var=pd.DataFrame({"feature": ["x", "y"]}, index=["g1", "g2"]),
        layers={"counts": np.array([[2.0, 1.0], [1.0, 2.0]])},
        obsm={"spatial": np.array([[0.0, 1.0], [2.0, 3.0]])},
    )
    reference = AnnData(
        X=np.array([[2.0, 0.0], [0.0, 2.0]]),
        obs=pd.DataFrame(
            {
                "Level1": pd.Categorical(
                    ["B", "A"], categories=["A", "B"], ordered=True
                )
            },
            index=["cell1", "cell2"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    captured = {}
    expected = pd.DataFrame(
        [[7.0, 3.0], [4.0, 6.0]],
        index=target.obs_names,
        columns=["B", "A"],
    )

    def annotate(
        adata,
        ref,
        annotation_key,
        *,
        result_key,
        return_reference,
        multi_center,
        lamb,
    ):
        captured.update(
            target_obs=adata.obs_names.copy(),
            target_var=adata.var_names.copy(),
            reference_obs=ref.obs_names.copy(),
            reference_var=ref.var_names.copy(),
            reference_categories=list(ref.obs[annotation_key].cat.categories),
            annotation_key=annotation_key,
            result_key=result_key,
            return_reference=return_reference,
            multi_center=multi_center,
            lamb=lamb,
        )
        ref.uns["pollution"] = True
        adata.X[:] = 99
        adata.var["pollution"] = True
        adata.layers["counts"][:] = 99
        adata.obsm["internal"] = np.ones((adata.n_obs, 1))
        adata.obsm[result_key] = expected.copy()
        return adata, ref

    monkeypatch.setattr(
        ot_kernel_module,
        "require_tacco",
        lambda: SimpleNamespace(tl=SimpleNamespace(annotate=annotate)),
    )

    result = kernels.OTKernel.annotate(
        target,
        reference,
        method="tacco",
        annotation_key="Level1",
        confidence_key="Confidence",
        multi_center=1,
        lamb=0.001,
    )

    assert captured["target_obs"].equals(target.obs_names)
    assert captured["target_var"].equals(target.var_names)
    assert captured["reference_obs"].equals(reference.obs_names)
    assert captured["reference_var"].equals(reference.var_names)
    assert captured["reference_categories"] == ["B", "A"]
    assert captured["annotation_key"] == "Level1"
    assert captured["result_key"] != "Level1"
    assert captured["return_reference"] is True
    assert captured["multi_center"] == 1
    assert captured["lamb"] == 0.001
    np.testing.assert_array_equal(result.X, target.X)
    pd.testing.assert_frame_equal(result.var, target.var)
    np.testing.assert_array_equal(result.layers["counts"], target.layers["counts"])
    np.testing.assert_array_equal(result.obsm["spatial"], target.obsm["spatial"])
    assert "internal" not in result.obsm
    assert "pollution" not in reference.uns
    assert list(reference.obs["Level1"].cat.categories) == ["A", "B"]
    pd.testing.assert_frame_equal(result.obsm["Level1"], expected)
    assert result.obs["Level1"].tolist() == ["B", "A"]
    assert result.obs["Confidence"].tolist() == [7.0, 6.0]
