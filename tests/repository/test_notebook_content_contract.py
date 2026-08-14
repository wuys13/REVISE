from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CASE_DIR = ROOT / "reproduce" / "case"
XENIUM_NOTEBOOKS = (
    "Xenium_sc_SVC_T.ipynb",
    "Xenium_sc_SVC_Fibroblast.ipynb",
    "Xenium_sc_SVC_Monocyte.ipynb",
)
CODE_SOURCE_DIGESTS = {
    "Xenium_sc_SVC_T.ipynb": "45938fda09840700e66f2df505759af2c05db10a0f62edcadac93dc91de6aaea",
    "Xenium_sc_SVC_Fibroblast.ipynb": "6b7299ab022f3ce34f4e3ca4d405f5fed0cdf92f4cafa68f2d35aa55f55a53b7",
    "Xenium_sc_SVC_Monocyte.ipynb": "28a2aee83380a77cc1da238e55fc8abda0a5c4bc2a1fd3058c412e826a496818",
    "VisiumHD_sp_SVC.ipynb": "d38bf4308eba8c09a152a77583674a6b35a485eb23d6a0f5271981f177c3c763",
}


def _notebook(name: str) -> dict:
    return json.loads((CASE_DIR / name).read_text(encoding="utf-8"))


def _markdown(name: str) -> str:
    return "\n".join(
        "".join(cell.get("source", []))
        for cell in _notebook(name)["cells"]
        if cell["cell_type"] == "markdown"
    )


def _code_digest(name: str) -> str:
    sources = [
        "".join(cell.get("source", []))
        for cell in _notebook(name)["cells"]
        if cell["cell_type"] == "code"
    ]
    payload = json.dumps(
        sources, ensure_ascii=False, separators=(",", ":")
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def test_curated_notebooks_match_approved_ordered_code_cell_sources():
    for name, expected_digest in CODE_SOURCE_DIGESTS.items():
        assert _code_digest(name) == expected_digest


def test_notebooks_publish_only_the_platform_matched_download_route():
    application = "https://zenodo.org/records/17705737"
    sc_svc = "https://zenodo.org/records/18389211"
    sp_svc = "https://zenodo.org/records/18389835"

    for name in XENIUM_NOTEBOOKS:
        markdown = _markdown(name)
        assert "Two ways to obtain the analysis inputs" in markdown
        assert application in markdown
        assert sc_svc in markdown
        assert sp_svc not in markdown
        assert "spatial-side" in markdown
        assert "expression-side" in markdown
        assert "not observation-paired" in markdown

    visium = _markdown("VisiumHD_sp_SVC.ipynb")
    assert application in visium
    assert sp_svc in visium
    assert sc_svc not in visium
    for stale in ("zenodo.org/uploads/18396759", "[raw_data_path]", "[svc_data_path]"):
        assert stale not in visium


def test_xenium_narratives_separate_direct_observation_and_paper_context():
    for name in XENIUM_NOTEBOOKS:
        markdown = _markdown(name)
        assert "Direct observation" in markdown
        assert "Paper context" in markdown
        assert "bioinfo analysis" not in markdown
        assert "### sc_SVC" not in markdown

    t_markdown = _markdown("Xenium_sc_SVC_T.ipynb")
    for value in ("2287", "123", "90", "64", "padj < 0.01", "padj < 1e-6"):
        assert value in t_markdown

    fib_markdown = _markdown("Xenium_sc_SVC_Fibroblast.ipynb")
    assert "Hypoxia-related pathway score" in fib_markdown
    assert "Glycolysis" not in fib_markdown
    assert "does not establish lineage" in fib_markdown

    mono_markdown = _markdown("Xenium_sc_SVC_Monocyte.ipynb")
    assert "preselected ligand–receptor pairs" in mono_markdown
    assert "does not use spatial distance" in mono_markdown
    assert "does not establish causal communication" in mono_markdown


def test_notebook_labels_explain_their_paper_facing_equivalents():
    t_markdown = _markdown("Xenium_sc_SVC_T.ipynb")
    assert "raw Xenium versus sc-SVC" in t_markdown
    assert "spatial-side panel-space" in t_markdown

    fib_markdown = _markdown("Xenium_sc_SVC_Fibroblast.ipynb")
    assert "CAF_6 (cluster 6 here)" in fib_markdown
    assert "normal-infiltrated comparator" in fib_markdown

    mono_markdown = _markdown("Xenium_sc_SVC_Monocyte.ipynb")
    assert "TAM_1 T-cell-recruiting" in mono_markdown
    for label in (
        "T_5 (paper T_3)",
        "TAM_0 (paper TAM_3)",
        "T_1 (paper T_0)",
        "TAM_1 (paper TAM_1)",
    ):
        assert label in mono_markdown


def test_xenium_retains_executed_snapshots_and_visium_remains_clean():
    for name in XENIUM_NOTEBOOKS:
        code_cells = [
            cell for cell in _notebook(name)["cells"] if cell["cell_type"] == "code"
        ]
        assert all(cell["execution_count"] is not None for cell in code_cells)
        assert all(
            output.get("output_type") != "error"
            for cell in code_cells
            for output in cell["outputs"]
        )

    visium_code = [
        cell
        for cell in _notebook("VisiumHD_sp_SVC.ipynb")["cells"]
        if cell["cell_type"] == "code"
    ]
    assert all(cell["execution_count"] is None for cell in visium_code)
    assert all(cell["outputs"] == [] for cell in visium_code)
