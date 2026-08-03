from types import SimpleNamespace

import pandas as pd


def _context(tmp_path, labels):
    artifacts = {}
    posterior = pd.DataFrame(
        {
            "B": [0.1] * len(labels),
            "T": [0.2] * len(labels),
        },
        index=[f"spot-{i}" for i in range(len(labels))],
    )

    def record_artifact(artifact):
        artifacts.setdefault("records", []).append(artifact)

    return SimpleNamespace(
        runner=SimpleNamespace(
            st_adata=SimpleNamespace(
                obs=pd.DataFrame({"Level1": labels}),
                obsm={"Level1": posterior},
            )
        ),
        columns={"cell_type_col": "Level1"},
        run_dir=tmp_path,
        artifacts=artifacts,
        record_artifact=record_artifact,
    )


def test_ist_selection_assessment_excludes_tumor_and_epi_and_warns_above_threshold(
    tmp_path,
):
    from revise.backend.adapters import _assess_ist_selection

    labels = ["Fibroblast"] * 3 + ["Tumor-like"] * 20_001 + ["EPI progenitor"] * 2
    ctx = _context(tmp_path, labels)

    assessment = _assess_ist_selection(ctx)

    assert assessment["default_candidates"] == ["Fibroblast"]
    assert assessment["excluded_cell_types"] == ["Tumor-like", "EPI progenitor"]
    assert {item["cell_type"] for item in assessment["warnings"] if "cell_type" in item} == {
        "Tumor-like"
    }
    assert assessment["counts"]["Tumor-like"] == 20_001
    assert (tmp_path / "selection_assessment.json").is_file()


def test_ist_selection_assessment_warns_when_no_default_candidates(tmp_path):
    from revise.backend.adapters import _assess_ist_selection

    ctx = _context(tmp_path, ["Tumor", "Epithelial"])

    assessment = _assess_ist_selection(ctx)

    assert assessment["default_candidates"] == []
    assert any(item.get("code") == "no_default_candidates" for item in assessment["warnings"])


def test_ist_selection_threshold_is_strictly_greater_than_20000(tmp_path):
    from revise.backend.adapters import _assess_ist_selection

    ctx = _context(tmp_path, ["Fibroblast"] * 20_000)

    assessment = _assess_ist_selection(ctx)

    assert not any(item.get("code") == "over_20000" for item in assessment["warnings"])


def test_ist_selection_exports_ga_posterior_csv_with_axes_and_values(tmp_path):
    from revise.backend.adapters import _assess_ist_selection

    labels = ["B", "T"]
    ctx = _context(tmp_path, labels)
    ctx.runner.st_adata.obsm["Level1"] = pd.DataFrame(
        [[0.7, 0.3], [0.2, 0.8]],
        index=["spot-2", "spot-1"],
        columns=["B", "T"],
    )

    assessment = _assess_ist_selection(ctx)

    posterior = pd.read_csv(tmp_path / "GA_posterior.csv")
    assert list(posterior.columns) == ["spot_id", "B", "T"]
    assert posterior["spot_id"].tolist() == ["spot-2", "spot-1"]
    assert posterior[["B", "T"]].to_numpy().tolist() == [[0.7, 0.3], [0.2, 0.8]]
    assert assessment["ga_posterior"]["path"] == str(tmp_path / "GA_posterior.csv")
    assert any(record["role"] == "ga_posterior" for record in ctx.artifacts["records"])


def test_ist_explicit_selection_deduplicates_in_first_seen_order():
    from revise.backend.adapters import _normalize_selected_cell_types

    assert _normalize_selected_cell_types(["T", "B", "T", ""]) == ["T", "B"]
