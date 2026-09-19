from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml
from anndata import AnnData


def _write_inputs(tmp_path: Path) -> tuple[Path, Path]:
    st = AnnData(
        X=np.array(
            [
                [8.0, 1.0, 1.0],
                [7.0, 2.0, 1.0],
                [1.0, 8.0, 1.0],
                [2.0, 7.0, 1.0],
            ]
        ),
        obs=pd.DataFrame(index=["spot1", "spot2", "spot3", "spot4"]),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
        obsm={"spatial": np.array([[0, 0], [1, 0], [0, 1], [1, 1]])},
    )
    reference = AnnData(
        X=np.array(
            [
                [9.0, 1.0, 1.0],
                [8.0, 2.0, 1.0],
                [1.0, 9.0, 1.0],
                [2.0, 8.0, 1.0],
                [7.0, 3.0, 1.0],
                [3.0, 7.0, 1.0],
            ]
        ),
        obs=pd.DataFrame(
            {
                "Level1": ["B", "B", "A", "A", "B", "A"],
                "Level2": ["B1", "B2", "A1", "A2", "B1", "A1"],
            },
            index=["cell1", "cell2", "cell3", "cell4", "cell5", "cell6"],
        ),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
    )
    st_path = tmp_path / "st.h5ad"
    reference_path = tmp_path / "candidate.h5ad"
    st.write_h5ad(st_path)
    reference.write_h5ad(reference_path)
    return st_path, reference_path


def _write_config(
    tmp_path: Path,
    *,
    svc_type: str,
    mode: str | None,
    solver: str = "pot",
) -> Path:
    application = {"svc_type": svc_type}
    if mode is not None:
        application["mode"] = mode
    if mode == "cluster":
        local_refinement = {
            "subtype_column": "Level2",
            "alpha": 0.2,
            "resolutions": [0.6],
        }
    elif mode == "sr":
        local_refinement = {
            "strength": 0.0,
            "graph": {
                "method": "pca",
                "alpha": 0.2,
                "n_neighbors": 2,
                "exp_neighbors": 2,
                "spatial_neighbors": 2,
            },
        }
    else:
        local_refinement = {"strength": 0.2}
    document = {
        "schema_version": 1,
        "application": application,
        "paths": {"root_dir": str(tmp_path)},
        "algorithm": {"ot_method": solver},
        "inputs": {
            "st": {"path": "st.h5ad", "format": "h5ad"},
            # Screening must replace this whole-file input and clear this old
            # filter before any reference validation or loading occurs.
            "reference": {
                "path": "configured-reference-does-not-exist.h5ad",
                "format": "h5ad",
                "filter_column": "Patient",
                "filter_value": "configured-only",
            },
        },
        "preprocessing": {
            "spatial": {
                "min_transcript_counts": None,
                "min_counts": None,
                "min_cell_counts": 1,
            },
            "reference": {
                "min_transcript_counts": None,
                "min_genes": None,
                "min_cell_counts": 1,
            },
        },
        "global_anchoring": {"broad_column": "Level1"},
        "local_refinement": local_refinement,
        "output": {"dir": "unused-output", "name": "unused"},
        "execution": {"seed": 17},
    }
    if mode == "cluster":
        document["output"]["ist_mapping"] = "random"
    path = tmp_path / f"{svc_type}-{mode or 'default'}-{solver}.yaml"
    path.write_text(yaml.safe_dump(document, sort_keys=False), encoding="utf-8")
    return path


def _mock_posterior(
    target,
    reference,
    *,
    annotation_key,
    confidence_key,
    method,
    scoring_genes_callback,
    **_,
):
    assert method == "pot"
    scoring_genes_callback(
        [
            str(value)
            for value in target.var_names.intersection(
                reference.var_names,
                sort=False,
            ).tolist()
        ]
    )
    categories = list(pd.unique(reference.obs[annotation_key]))
    assert categories == ["B", "A"]
    values = np.array(
        [[0.8, 0.2], [0.7, 0.3], [0.2, 0.8], [0.3, 0.7]],
        dtype=float,
    )
    result = target.copy()
    result.obsm[annotation_key] = pd.DataFrame(
        values,
        index=target.obs_names.copy(),
        columns=categories,
    )
    result.obs[annotation_key] = result.obsm[annotation_key].idxmax(axis=1)
    result.obs[confidence_key] = values.max(axis=1)
    return result


@pytest.mark.parametrize(
    ("svc_type", "mode", "route"),
    [
        ("sp-SVC", None, "application:sp-SVC"),
        ("sc-SVC", "cluster", "application:sc-SVC:cluster"),
        ("sc-SVC", "sr", "application:sc-SVC:sr"),
    ],
)
def test_mock_solver_runs_all_application_routes_through_ga_only(
    monkeypatch,
    tmp_path,
    svc_type,
    mode,
    route,
):
    from revise.backend.kernels.ot import OTKernel
    from revise.recon.pipeline import UnifiedReconstructionPipeline
    from revise.reference_preparation.host import run_global_anchoring

    st_path, reference_path = _write_inputs(tmp_path)
    config_path = _write_config(tmp_path, svc_type=svc_type, mode=mode)
    monkeypatch.setattr(OTKernel, "annotate", staticmethod(_mock_posterior))

    def forbidden_local_refinement(*_args, **_kwargs):
        raise AssertionError("GA-only screening entered Local Refinement")

    monkeypatch.setattr(
        UnifiedReconstructionPipeline,
        "local_refinement",
        forbidden_local_refinement,
    )
    response = run_global_anchoring(reference_path, config_path)

    assert response.expected_st_unit_ids == ["spot1", "spot2", "spot3", "spot4"]
    assert response.result.st_unit_ids == response.expected_st_unit_ids
    assert response.result.cell_type_labels == ["B", "A"]
    np.testing.assert_allclose(response.result.distribution.sum(axis=1), 1.0)
    assert response.result.metadata["route"] == route
    assert response.result.metadata["effective_seed"] == 17
    assert response.result.metadata["effective_solver"] == "pot"
    assert response.result.metadata["reference_n_obs"] == 6
    assert response.result.metadata["reference_n_vars"] == 3
    assert response.result.metadata["scoring_n_vars"] == 3
    assert response.result.metadata["cell_type_labels"] == ["B", "A"]
    assert not (tmp_path / "unused-output").exists()
    assert st_path.exists()


def test_mock_solver_nonfinite_posterior_fails_core_contract(monkeypatch, tmp_path):
    from revise.backend.kernels.ot import OTKernel
    from revise.reference_preparation.host import run_global_anchoring

    _, reference_path = _write_inputs(tmp_path)
    config_path = _write_config(tmp_path, svc_type="sp-SVC", mode=None)

    def nonfinite(
        target,
        reference,
        *,
        annotation_key,
        scoring_genes_callback,
        **_kwargs,
    ):
        scoring_genes_callback(list(target.var_names))
        result = target.copy()
        result.obsm[annotation_key] = pd.DataFrame(
            [[np.nan, np.nan]] * target.n_obs,
            index=target.obs_names.copy(),
            columns=list(pd.unique(reference.obs[annotation_key])),
        )
        return result

    monkeypatch.setattr(OTKernel, "annotate", staticmethod(nonfinite))
    with pytest.raises(ValueError, match="finite and non-negative"):
        run_global_anchoring(reference_path, config_path)


def test_real_pot_host_adapter_computes_a_small_numerical_ga(tmp_path):
    pytest.importorskip("ot", reason="real POT solver is not installed")
    from revise.reference_preparation.host import run_global_anchoring

    _, reference_path = _write_inputs(tmp_path)
    config_path = _write_config(tmp_path, svc_type="sp-SVC", mode=None)

    response = run_global_anchoring(reference_path, config_path)

    assert response.result.distribution.shape == (4, 2)
    assert np.isfinite(response.result.distribution).all()
    np.testing.assert_allclose(response.result.distribution.sum(axis=1), 1.0)
    assert response.result.cell_type_labels == ["B", "A"]


def test_real_tacco_host_adapter_computes_a_small_numerical_ga(tmp_path):
    pytest.importorskip("tacco", reason="real TACCO solver is not installed")
    from revise.reference_preparation.host import run_global_anchoring

    _, reference_path = _write_inputs(tmp_path)
    config_path = _write_config(
        tmp_path,
        svc_type="sc-SVC",
        mode="cluster",
        solver="tacco",
    )

    response = run_global_anchoring(reference_path, config_path)

    assert response.result.distribution.shape == (4, 2)
    assert np.isfinite(response.result.distribution).all()
    np.testing.assert_allclose(response.result.distribution.sum(axis=1), 1.0)
    assert response.result.metadata["effective_solver"] == "tacco"
    assert response.result.metadata["scoring_n_vars"] > 0
    assert response.result.metadata["scoring_genes_sha256"]
