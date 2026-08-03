from types import SimpleNamespace

import pandas as pd


def _context(tmp_path, labels):
    artifacts = {}

    def record_artifact(artifact):
        artifacts.setdefault("records", []).append(artifact)

    return SimpleNamespace(
        runner=SimpleNamespace(
            st_adata=SimpleNamespace(obs=pd.DataFrame({"Level1": labels}))
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


def test_ist_explicit_selection_deduplicates_in_first_seen_order():
    from revise.backend.adapters import _normalize_selected_cell_types

    assert _normalize_selected_cell_types(["T", "B", "T", ""]) == ["T", "B"]
