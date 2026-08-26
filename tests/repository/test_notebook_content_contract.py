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
    "Xenium_sc_SVC_T.ipynb": "cfeff26eb68901473b1f59f18e46c09722674e9ab71c8eb75bc0942e1470f6e3",
    "Xenium_sc_SVC_Fibroblast.ipynb": "5076c49f1cbb37f161914c5e1f1267dcd9b5c4ab330b6979d2b8f3a810d3e6d9",
    "Xenium_sc_SVC_Monocyte.ipynb": "2c2c2efd721d2eb07d23908ee862fec58be6c3cb691eb5d64b5c58b29962b95c",
    "VisiumHD_sp_SVC.ipynb": "9833157b60f4cbde93888cf0a61381c2f2403a76b2d42ed51849ca96e9b8adb6",
    "segmentation.ipynb": "0b72170a0fedce72ddeb01f698799f7e9f9e412188fcf84c50d588dc72348a16",
    "bin2cell.ipynb": "2d9082ca947f18036fbf251165decbd4dfd1aca5c054a00b6918cf62294a6274",
    "batch.ipynb": "c56737e6b4c4b042c2c5f78108a341ab1a559b3ec99ecc703bd6ad8a52caf872",
    "spot.ipynb": "37db402cfffe9cac631c9c98d3bdf574a499aa62d638531c7dd9be7d74b67d92",
    "imputation_and_dropout.ipynb": "93f30ef9ef0d142b617d668d6fafd534bcf7150d71603d7a86d2ff9baaa98e93",
    "plot_imputation_case.ipynb": "df15d5935f310583e5299f4e1acc9cd512ec3d52109ef3267b3730cfe408878e",
}
CODE_SNAPSHOT_DIGESTS = {
    "VisiumHD_sp_SVC.ipynb": "76d4c142b43adcb33e5232491ab95a29a470ce599c58c25380eedf98d3b42bef",
    "Xenium_sc_SVC_T.ipynb": "4c1b7e8d9828f37acf686c815e1b65282395647fbc4697606e09d067b9cd6652",
    "Xenium_sc_SVC_Fibroblast.ipynb": "75f4d6eeaed97fa120387fa827c52114ba0aa82ba6b24e525db48d1efc482dfd",
    "Xenium_sc_SVC_Monocyte.ipynb": "94f2e97696698c1d6c7beab12ce00aa8dd9329390283bab2ce95a900f63f7a3f",
    "segmentation.ipynb": "d6e35abb508b3b728b8253119a452f77309ca6b3953fb6108ab939ca0cc71ab5",
    "bin2cell.ipynb": "4a6506161921ca09e6fe13bd9184f624010aa24c34f0d1713b898c162fa9395d",
    "batch.ipynb": "542f8c2deef00b98e2ec5ef72bf80c085ff21e09d46ccdede1f810d67b8a25cb",
    "spot.ipynb": "eabf5ce8a9e4f2e20c6fee182f9a6b6b5e29cd06a8f90f455fc036c7a241aedc",
    "imputation_and_dropout.ipynb": "1dd00c67e23ec7353c462ada7347f881a6298e5430c06dde9975f79f1fffd42b",
    "plot_imputation_case.ipynb": "56cccedad1a06a5fcfbb688ec9a4046f899be7c0f9e0cf5af2ae427fd5a84251",
}
NOTEBOOK_METADATA_DIGESTS = {
    "VisiumHD_sp_SVC.ipynb": "a7762c677ddb239293c8ef52b3a4f641e7c8d67998e7bc3a7902b3fa7dc70b2c",
}
VISIUM_SC_SVC_NOTEBOOK = "Visium_sc_SVC_mouse_brain.ipynb"


def _notebook(name: str) -> dict:
    directory = ROOT / "reproduce" / "benchmark" if name in BENCHMARK_NOTEBOOKS else CASE_DIR
    return json.loads((directory / name).read_text(encoding="utf-8"))


def _markdown(name: str) -> str:
    return "\n".join(
        "".join(cell.get("source", []))
        for cell in _notebook(name)["cells"]
        if cell["cell_type"] == "markdown"
    )


def _code_source(name: str) -> str:
    return "\n".join(
        "".join(cell.get("source", []))
        for cell in _notebook(name)["cells"]
        if cell["cell_type"] == "code"
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


def _metadata_digest(name: str) -> str:
    payload = _notebook(name).get("metadata", {})
    return hashlib.sha256(
        json.dumps(payload, ensure_ascii=False, separators=(",", ":")).encode()
    ).hexdigest()


def test_curated_notebooks_match_approved_ordered_code_cell_sources():
    for name, expected_digest in CODE_SOURCE_DIGESTS.items():
        assert _code_digest(name) == expected_digest


def test_active_notebooks_preserve_execution_counts_and_outputs():
    for name, expected_digest in CODE_SNAPSHOT_DIGESTS.items():
        assert _code_snapshot_digest(name) == expected_digest
    for name, expected_digest in NOTEBOOK_METADATA_DIGESTS.items():
        assert _metadata_digest(name) == expected_digest


def test_xenium_notebooks_read_the_current_derived_cluster_directories():
    t_source = _code_source("Xenium_sc_SVC_T.ipynb")
    fib_source = _code_source("Xenium_sc_SVC_Fibroblast.ipynb")
    mono_source = _code_source("Xenium_sc_SVC_Monocyte.ipynb")

    for source in (t_source, fib_source, mono_source):
        assert 'data_root = "../../results/sc_SVC_case/P2CRC_Xenium"' in source
    assert 'case_dir = f"{data_root}/Fibroblast"' in fib_source
    assert 'case_dir = f"{data_root}/Mono_Macro"' in mono_source
    assert (
        'analysis_output_dir = os.environ.get("REVISE_ANALYSIS_OUTPUT_ROOT", '
        '"../../output/P2CRC_Xenium/Fib")'
    ) in fib_source
    assert (
        'analysis_output_dir = os.environ.get("REVISE_ANALYSIS_OUTPUT_ROOT", '
        '"../../output/P2CRC_Xenium/Mono")'
    ) in mono_source
    assert "results/P2CRC_Xenium/Fibroblast" not in fib_source
    assert "results/P2CRC_Xenium/Mono_Macro" not in mono_source
    assert "../../results/sc_SVC_case/P2CRC_Xenium/T/expr.h5ad" in mono_source
    assert "../../results/sc_SVC_case/P2CRC_Xenium/Mono_Macro/expr.h5ad" in mono_source


def test_active_notebooks_use_the_current_download_names_and_record():
    material = "https://zenodo.org/records/21921802"
    sc_svc = "https://zenodo.org/records/22046001"
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

    imputation_markdown = _markdown("imputation_and_dropout.ipynb")
    imputation_code = _code_source("imputation_and_dropout.ipynb")
    assert "SpaGE_impute.csv" in imputation_markdown
    assert "identity control" in imputation_markdown
    assert 'methods = ["Raw", "Tangram", "gimVI", "SpaGE", "REVISE"]' in imputation_code

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
    assert "TAM_1 recruiting-inflammatory" in mono_markdown
    for label in (
        "T_5 (paper T_3)",
        "TAM_0 (paper TAM_3)",
        "T_1 (paper T_0)",
        "TAM_1 (paper TAM_1)",
    ):
        assert label in mono_markdown


def test_active_notebooks_retain_approved_execution_snapshots():
    for name in XENIUM_NOTEBOOKS:
        code_cells = [
            cell for cell in _notebook(name)["cells"] if cell["cell_type"] == "code"
        ]
        assert any(cell["execution_count"] is not None for cell in code_cells)
        assert any(cell["outputs"] for cell in code_cells)
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
    assert any(cell["execution_count"] is not None for cell in visium_code)
    assert any(cell["outputs"] for cell in visium_code)
    assert all(
        output.get("output_type") != "error"
        for cell in visium_code
        for output in cell["outputs"]
    )
    visium_markdown = _markdown("VisiumHD_sp_SVC.ipynb")
    visium_source = _code_source("VisiumHD_sp_SVC.ipynb")
    assert "sampled independently" in visium_markdown
    assert "observation-paired inference" in visium_markdown
    assert "sample_independently" in visium_source
    assert "sample_paired_by_identity" not in visium_source


def test_visium_sc_svc_case_uses_the_canonical_direct_application_route():
    notebook = _notebook(VISIUM_SC_SVC_NOTEBOOK)
    source = "\n".join(
        "".join(cell.get("source", []))
        for cell in notebook["cells"]
    )

    assert not (CASE_DIR / "sc_SVC_sr_case_Visium_mouse_brain.ipynb").exists()
    assert 'SAMPLE_NAME = "REVISEVisiumMouseBrain_sc-SVC"' in source
    assert 'RECONSTRUCTION_COMMAND = "python reconstruct.py --config configs/application/Visium.yaml"' in source
    assert "sc-SVC-sr" not in source
    assert "urlretrieve" not in source
    assert "run_manifest" not in source

    code_cells = [cell for cell in notebook["cells"] if cell["cell_type"] == "code"]
    assert all(cell["execution_count"] is None for cell in code_cells)
    outputs = [output for cell in code_cells for output in cell.get("outputs", [])]
    assert len(outputs) == 7
    assert all(output.get("output_type") == "display_data" for output in outputs)
    assert all(set(output.get("data", {})) == {"image/png"} for output in outputs)
