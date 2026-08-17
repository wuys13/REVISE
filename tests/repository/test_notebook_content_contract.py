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
BENCHMARK_NOTEBOOKS = (
    "segmentation.ipynb",
    "bin2cell.ipynb",
    "batch.ipynb",
    "spot.ipynb",
    "imputation_and_dropout.ipynb",
    "plot_imputation_case.ipynb",
)
CODE_SOURCE_DIGESTS = {
    "Xenium_sc_SVC_T.ipynb": "45938fda09840700e66f2df505759af2c05db10a0f62edcadac93dc91de6aaea",
    "Xenium_sc_SVC_Fibroblast.ipynb": "7e260ff47fba229a9742bfd9adcd255b719121a54dd0d49fbc8431821f81c4d1",
    "Xenium_sc_SVC_Monocyte.ipynb": "a47cc825541e70decc7d7a655bbd78929bf3565651c93aca61811274c7255649",
    "VisiumHD_sp_SVC.ipynb": "d38bf4308eba8c09a152a77583674a6b35a485eb23d6a0f5271981f177c3c763",
    "segmentation.ipynb": "0b72170a0fedce72ddeb01f698799f7e9f9e412188fcf84c50d588dc72348a16",
    "bin2cell.ipynb": "2d9082ca947f18036fbf251165decbd4dfd1aca5c054a00b6918cf62294a6274",
    "batch.ipynb": "c56737e6b4c4b042c2c5f78108a341ab1a559b3ec99ecc703bd6ad8a52caf872",
    "spot.ipynb": "37db402cfffe9cac631c9c98d3bdf574a499aa62d638531c7dd9be7d74b67d92",
    "imputation_and_dropout.ipynb": "1c00f690839c1719ba014f98464d05abd813a3e75da0cb8719c37b104a64936c",
    "plot_imputation_case.ipynb": "df15d5935f310583e5299f4e1acc9cd512ec3d52109ef3267b3730cfe408878e",
}
CODE_SNAPSHOT_DIGESTS = {
    "VisiumHD_sp_SVC.ipynb": "03e156ca648a4414dc2a7058b817afa9c27b6dfd526c27512a30a1387fd09c64",
    "Xenium_sc_SVC_T.ipynb": "92e536179f648c92468a005e11e9547b99bbbc525508af39244b1f0f13f703dc",
    "Xenium_sc_SVC_Fibroblast.ipynb": "0a69132b8c167006f3e07abc811fcf4474405a3d7175b32bc778af175f7d5c8d",
    "Xenium_sc_SVC_Monocyte.ipynb": "ada46cc32b632beb08c0af2c1f3baebbe0ee8d0e2691f7c9ccb07649774dbf02",
    "segmentation.ipynb": "d6e35abb508b3b728b8253119a452f77309ca6b3953fb6108ab939ca0cc71ab5",
    "bin2cell.ipynb": "4a6506161921ca09e6fe13bd9184f624010aa24c34f0d1713b898c162fa9395d",
    "batch.ipynb": "e0d1f4a222e9674a5f54e771a179365787feef13332b8e5b66f6887dc11cb6a3",
    "spot.ipynb": "1e05cfa2f85e1d9c2b1548197362466b1eeea4177c53fa687f1c7228355ec67a",
    "imputation_and_dropout.ipynb": "5d5a9940c9f17d4595cd126b11ace7877525bdb8494c2e31a2718dd4aa77f367",
    "plot_imputation_case.ipynb": "56cccedad1a06a5fcfbb688ec9a4046f899be7c0f9e0cf5af2ae427fd5a84251",
}
VISIUM_SC_SVC_NOTEBOOK = "Visium_sc_SVC_mouse_brain.ipynb"
VISIUM_UNCHANGED_CODE_DIGEST = (
    "8ea767f4955eb5ecc507cbbfa16214e72ea5d7c0132db30f44d783a8a329d54c"
)


def _notebook(name: str) -> dict:
    directory = ROOT / "reproduce" / "benchmark" if name in BENCHMARK_NOTEBOOKS else CASE_DIR
    return json.loads((directory / name).read_text(encoding="utf-8"))


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


def _code_snapshot_digest(name: str) -> str:
    payload = [
        {
            "execution_count": cell.get("execution_count"),
            "outputs": cell.get("outputs", []),
        }
        for cell in _notebook(name)["cells"]
        if cell["cell_type"] == "code"
    ]
    return hashlib.sha256(
        json.dumps(payload, ensure_ascii=False, separators=(",", ":")).encode()
    ).hexdigest()


def _code_digest_except(name: str, excluded_indexes: set[int]) -> str:
    sources = [
        "".join(cell.get("source", []))
        for index, cell in enumerate(_notebook(name)["cells"])
        if cell["cell_type"] == "code" and index not in excluded_indexes
    ]
    payload = json.dumps(
        sources, ensure_ascii=False, separators=(",", ":")
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def test_curated_notebooks_match_approved_ordered_code_cell_sources():
    for name, expected_digest in CODE_SOURCE_DIGESTS.items():
        assert _code_digest(name) == expected_digest


def test_active_notebooks_preserve_execution_counts_and_outputs():
    for name, expected_digest in CODE_SNAPSHOT_DIGESTS.items():
        assert _code_snapshot_digest(name) == expected_digest


def test_xenium_notebooks_read_the_current_derived_cluster_directories():
    fib_source = "\n".join(
        "".join(cell.get("source", []))
        for cell in _notebook("Xenium_sc_SVC_Fibroblast.ipynb")["cells"]
    )
    mono_source = "\n".join(
        "".join(cell.get("source", []))
        for cell in _notebook("Xenium_sc_SVC_Monocyte.ipynb")["cells"]
    )

    assert "P2CRC_Xenium/Fibroblast" in fib_source
    assert "P2CRC_Xenium/Mono_Macro" in mono_source
    assert "P2CRC_Xenium/Fib/" not in fib_source
    assert "P2CRC_Xenium/Mono/" not in mono_source


def test_active_notebooks_use_the_current_download_names_and_record():
    material = "https://zenodo.org/records/21921802"
    sc_svc = "https://zenodo.org/records/18389211"
    sp_svc = "https://zenodo.org/records/18389835"

    for name in XENIUM_NOTEBOOKS:
        markdown = _markdown(name)
        assert "Two ways to obtain the analysis inputs" in markdown
        assert material in markdown
        assert "Real-world ST datasets" in markdown
        assert sc_svc in markdown
        assert sp_svc not in markdown
        assert "spatial-side" in markdown
        assert "expression-side" in markdown
        assert "not observation-paired" in markdown

    visium = _markdown("VisiumHD_sp_SVC.ipynb")
    assert material in visium
    assert "Real-world ST datasets" in visium
    assert sp_svc in visium
    assert sc_svc not in visium
    for stale in ("zenodo.org/uploads/18396759", "[raw_data_path]", "[svc_data_path]"):
        assert stale not in visium

    for name in BENCHMARK_NOTEBOOKS:
        notebook = json.loads((ROOT / "reproduce" / "benchmark" / name).read_text(encoding="utf-8"))
        markdown = "\n".join(
            "".join(cell.get("source", []))
            for cell in notebook["cells"]
            if cell["cell_type"] == "markdown"
        )
        assert material in markdown
        assert "Sim2Real-ST benchmark" in markdown
        assert "Reproduced benchmark results" in markdown

    for path in (
        *(CASE_DIR / name for name in CODE_SOURCE_DIGESTS if name in {
            "VisiumHD_sp_SVC.ipynb", *XENIUM_NOTEBOOKS
        }),
        *(ROOT / "reproduce" / "benchmark" / name for name in BENCHMARK_NOTEBOOKS),
    ):
        assert "17705737" not in path.read_text(encoding="utf-8")


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


def test_visium_sc_svc_hard_cut_changes_only_allowed_code_cells_and_names():
    notebook = _notebook(VISIUM_SC_SVC_NOTEBOOK)
    source = "\n".join(
        line
        for cell in notebook["cells"]
        for line in cell.get("source", [])
    )

    assert not (CASE_DIR / "sc_SVC_sr_case_Visium_mouse_brain.ipynb").exists()
    assert _code_digest_except(VISIUM_SC_SVC_NOTEBOOK, {6, 20}) == (
        VISIUM_UNCHANGED_CODE_DIGEST
    )
    assert 'SAMPLE_NAME = "REVISEVisiumMouseBrain_sc-SVC"' in source
    assert 'RECONSTRUCTION_COMMAND = "python reconstruct.py --config configs/application/Visium.yaml"' in source
    assert '("route", "sc-SVC")' in source
    assert '("mode", "sr")' in source
    assert "sc-SVC-sr" not in source

    code_cells = [cell for cell in notebook["cells"] if cell["cell_type"] == "code"]
    assert all(cell["execution_count"] is None for cell in code_cells)
    assert all(cell["outputs"] == [] for cell in code_cells)
