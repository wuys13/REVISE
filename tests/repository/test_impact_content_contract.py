"""Static contract for the reconstruction-impact content tree."""

from reproduce.case.reconstruction_impact.content_contract import (
    LAYERS,
    NODE_SPECS,
    QUESTION_LABELS,
    QUESTION_LAYERS,
    SCOPES,
)


def test_content_contract_exposes_the_four_layer_tree_and_shared_questions():
    assert SCOPES == ("Fibroblast", "Mono_Macro", "T")
    assert LAYERS == ("foundation", "overall", "localization", "region")

    expected_ids = [
        *(f"{layer}.{index}" for layer, count in ((1, 6), (2, 6), (3, 9), (4, 5)) for index in range(1, count + 1)),
        "summary",
    ]
    assert [spec["node_id"] for spec in NODE_SPECS] == expected_ids
    assert len(NODE_SPECS) == 27
    required = {"node_id", "layer", "title_zh", "title_en", "purpose_zh", "purpose_en", "method", "questions", "reading_placement"}
    assert all(set(spec) == required for spec in NODE_SPECS)
    assert all(isinstance(spec["method"], tuple) and len(spec["method"]) == 2 for spec in NODE_SPECS)
    assert all(isinstance(spec["questions"], tuple) for spec in NODE_SPECS)
    assert set(QUESTION_LAYERS) == set(QUESTION_LABELS)

    by_id = {spec["node_id"]: spec for spec in NODE_SPECS}
    assert by_id["1.1"]["questions"] == ("foundation",)
    assert by_id["2.4"]["questions"] == ("moran_all_valid", "moran_shared_valid")
    assert by_id["3.5"]["questions"] == ("local_vs_raw_leiden", "local_vs_raw_level2")
    assert by_id["summary"]["questions"] == tuple(QUESTION_LAYERS)


def test_reading_placement_keeps_technical_foundation_out_of_main_story():
    by_id = {spec["node_id"]: spec for spec in NODE_SPECS}
    assert all(by_id[f"1.{i}"]["reading_placement"] == "audit" for i in range(1, 7))
    assert by_id["2.3"]["reading_placement"] == "reference"
    assert by_id["2.1"]["reading_placement"] == "main"
