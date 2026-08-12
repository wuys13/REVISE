from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData


def test_global_anchor_returns_the_tacco_assignment_without_reordering(monkeypatch):
    import batch_ga

    spatial = AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame(index=["spot-1", "spot-2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    reference = AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame({"Level1": ["A", "B"]}, index=["cell-1", "cell-2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    posterior = pd.DataFrame(
        [[0.25, 0.75], [0.9, 0.1]],
        index=spatial.obs_names,
        columns=["B", "A"],
    )

    def annotate(target, _reference, **kwargs):
        assert kwargs["method"] == "tacco"
        assert kwargs["multi_center"] == 1
        assert kwargs["lamb"] == 0.001
        result = target.copy()
        result.obsm["Level1"] = posterior
        result.obs["Level1"] = posterior.idxmax(axis=1).to_numpy()
        return result

    monkeypatch.setattr(batch_ga.OTKernel, "annotate", annotate)

    assignment = batch_ga.global_anchor(spatial, reference, "Level1", seed=42)

    assert assignment.posterior.equals(posterior)
    assert assignment.labels.tolist() == ["A", "B"]


def test_run_batch_ga_runs_all_reference_ids_and_keeps_going_after_failure(
    monkeypatch,
    tmp_path,
):
    import batch_ga
    import yaml

    spatial = AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame({"transcript_counts": [5, 5]}, index=["spot-1", "spot-2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    reference = AnnData(
        X=np.ones((4, 2)),
        obs=pd.DataFrame(
            {"reference_id": ["A", "B", "A", "B"], "Level1": ["T", "T", "T", "T"]},
            index=["a1", "b1", "a2", "b2"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    config_path = tmp_path / "batch.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "schema_version": 1,
                "paths": {"root_dir": str(tmp_path)},
                "inputs": {
                    "st": {"path": "spatial.h5ad", "format": "h5ad"},
                    "reference": {
                        "path": "reference.h5ad",
                        "format": "h5ad",
                        "filter_column": "reference_id",
                        "filter_value": "all",
                    },
                },
                "preprocessing": {
                    "spatial": {"min_transcript_counts": None, "min_counts": None, "min_cell_counts": 0},
                    "reference": {"min_transcript_counts": None, "min_genes": None, "min_cell_counts": 0},
                },
                "global_anchoring": {"broad_column": "Level1"},
                "output": {"dir": "ga_output"},
                "execution": {"seed": 42},
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(batch_ga, "read_h5ad", lambda path: spatial.copy() if "spatial" in str(path) else reference.copy())

    def anchor(target, ref, _column, seed):
        if ref.obs_names[0] == "b1":
            raise ValueError("bad reference")
        posterior = pd.DataFrame([[1.0], [1.0]], index=target.obs_names, columns=["T"])
        return batch_ga.GlobalAssignment(labels=posterior.idxmax(axis=1), posterior=posterior)

    monkeypatch.setattr(batch_ga, "global_anchor", anchor)

    result = batch_ga.run_batch_ga(config_path)

    assert list(result.assignments) == ["A"]
    assert result.failures == {"B": "ValueError: bad reference"}
    assert (tmp_path / "ga_output" / "A.csv").is_file()
    status = __import__("json").loads(result.status_path.read_text())
    assert [task["reference_id"] for task in status["tasks"]] == ["A", "B"]
    assert status["summary"] == {"total": 2, "succeeded": 1, "failed": 1}
