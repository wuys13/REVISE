"""Behavioral style contracts for the streamlined root case-study notebooks.

These checks intentionally validate notebook behavior and structure rather than
pinning a whole source digest.  The builders are called only to obtain their
in-memory cell definitions; the published notebooks are never rewritten.
"""

from __future__ import annotations

import ast
import base64
import importlib.util
import io
import json
import re
from pathlib import Path
from types import ModuleType

import nbformat
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
CASE_DIR = ROOT / "reproduce" / "case"

NOTEBOOKS = (
    "CosMx_SMI_267T_not_sp_SVC.ipynb",
    "MERFISH_Allen_VISp_sc_SVC_cluster.ipynb",
    "SlideSeq_mouse_colon_sp_SVC.ipynb",
    "SlideSeq_mouse_olfactory_bulb_sp_SVC.ipynb",
    "StereoSeq_zebrafish_5hpf_sp_SVC.ipynb",
    "Visium_sc_SVC_mouse_brain.ipynb",
    "osmFISH_sc_SVC_cluster.ipynb",
)

PNG_CELL_COUNTS = {
    "CosMx_SMI_267T_not_sp_SVC.ipynb": {5: 1},
    "MERFISH_Allen_VISp_sc_SVC_cluster.ipynb": {8: 2, 10: 1},
    "SlideSeq_mouse_colon_sp_SVC.ipynb": {9: 1},
    "SlideSeq_mouse_olfactory_bulb_sp_SVC.ipynb": {8: 1},
    "StereoSeq_zebrafish_5hpf_sp_SVC.ipynb": {8: 1},
    "Visium_sc_SVC_mouse_brain.ipynb": {
        11: 1,
        16: 1,
        18: 1,
        24: 1,
        26: 1,
    },
    "osmFISH_sc_SVC_cluster.ipynb": {10: 3},
}

BUILDER_SPECS = {
    "CosMx_SMI_267T_not_sp_SVC.ipynb": (
        CASE_DIR / "cosmx" / "build_notebook.py",
        "build_notebook",
    ),
    "MERFISH_Allen_VISp_sc_SVC_cluster.ipynb": (
        CASE_DIR / "tacco" / "build_notebooks.py",
        "build_allen_visp_merfish",
    ),
    "SlideSeq_mouse_colon_sp_SVC.ipynb": (
        CASE_DIR / "tacco" / "build_notebooks.py",
        "build_colon",
    ),
    "SlideSeq_mouse_olfactory_bulb_sp_SVC.ipynb": (
        CASE_DIR / "tacco" / "build_notebooks.py",
        "build_olfactory",
    ),
    "StereoSeq_zebrafish_5hpf_sp_SVC.ipynb": (
        CASE_DIR / "stereo_seq" / "build_notebook.py",
        "build_zesta_zf5",
    ),
    "Visium_sc_SVC_mouse_brain.ipynb": (
        CASE_DIR / "visium_mouse_brain" / "build_notebook.py",
        "build_notebook",
    ),
    "osmFISH_sc_SVC_cluster.ipynb": (
        CASE_DIR / "tacco" / "build_notebooks.py",
        "build_osmfish",
    ),
}

MARKDOWN_REQUIREMENTS = (
    "Case study",
    "Notebook Guide",
    "Question",
    "Method",
    "Direct observation",
    "Interpretation boundary",
    "Phase 1 snapshot",
    "has not been rerun",
)

MACHINE_PATTERNS = (
    ("absolute user path", re.compile(r"/Users/", re.IGNORECASE)),
    ("temporary path", re.compile(r"/tmp/")),
    (
        "IPv4 address",
        re.compile(
            r"(?<![\d.])"
            r"(?:25[0-5]|2[0-4]\d|1?\d?\d)"
            r"(?:\.(?:25[0-5]|2[0-4]\d|1?\d?\d)){3}"
            r"(?![\d.])"
        ),
    ),
    ("ssh command", re.compile(r"\bssh\b", re.IGNORECASE)),
)

FORBIDDEN_CONTROL_PATTERNS = (
    ("FORCE_ control", re.compile(r"FORCE_")),
    ("cache control", re.compile(r"\bcache\b", re.IGNORECASE)),
    ("reuse control", re.compile(r"\breus(?:e|ing)\b", re.IGNORECASE)),
    ("fallback control", re.compile(r"\bfallback\b", re.IGNORECASE)),
)

ALLOWED_NOTEBOOK_IMPORT_ROOTS = {
    # Python standard library
    "os",
    "pathlib",
    "subprocess",
    "sys",
    # revise-svc runtime dependencies used by these cases
    "anndata",
    "matplotlib",
    "numpy",
    "pandas",
    "scanpy",
    "scipy",
    # Repository modules reached through the documented repository-root entrypoint
    "reconstruct",
    "revise",
}

DIRECT_UMAP_NOTEBOOKS = (
    "CosMx_SMI_267T_not_sp_SVC.ipynb",
    "MERFISH_Allen_VISp_sc_SVC_cluster.ipynb",
    "SlideSeq_mouse_colon_sp_SVC.ipynb",
    "SlideSeq_mouse_olfactory_bulb_sp_SVC.ipynb",
    "StereoSeq_zebrafish_5hpf_sp_SVC.ipynb",
    "osmFISH_sc_SVC_cluster.ipynb",
)


def _notebook_path(name: str) -> Path:
    return CASE_DIR / name


def _read_notebook(name: str) -> dict:
    return json.loads(_notebook_path(name).read_text(encoding="utf-8"))


def _source(cell: dict) -> str:
    source = cell.get("source", "")
    return "".join(source) if isinstance(source, list) else source


def _ordered_cell_sources(notebook: dict) -> list[tuple[str, str]]:
    return [
        (cell["cell_type"], _source(cell))
        for cell in notebook["cells"]
    ]


def _code_source(notebook: dict) -> str:
    return "\n".join(
        _source(cell)
        for cell in notebook["cells"]
        if cell["cell_type"] == "code"
    )


def _load_builder(path: Path, module_name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(module_name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _output_strings(value):
    """Yield textual output values while excluding binary image payloads."""

    if isinstance(value, str):
        yield value
    elif isinstance(value, list):
        for item in value:
            yield from _output_strings(item)
    elif isinstance(value, dict):
        for item in value.values():
            yield from _output_strings(item)


def _cell_texts(cell: dict):
    yield _source(cell)
    for output in cell.get("outputs", []):
        output_type = output.get("output_type")
        if output_type == "error":
            yield from _output_strings(output.get("traceback", []))
        elif output_type == "stream":
            yield from _output_strings(output.get("text", ""))
        for mime, value in output.get("data", {}).items():
            if mime != "image/png":
                yield from _output_strings(value)


def _png_outputs(cell: dict):
    for output in cell.get("outputs", []):
        if "image/png" in output.get("data", {}):
            yield output


def _decode_png(encoded) -> bytes:
    if isinstance(encoded, list):
        encoded = "".join(encoded)
    encoded = "".join(str(encoded).split())
    return base64.b64decode(encoded, validate=True)


def test_root_case_notebooks_are_nbformat_valid():
    for name in NOTEBOOKS:
        notebook = _read_notebook(name)
        nbformat.validate(notebook)


def test_root_case_code_cells_parse_without_asserts_or_display_calls():
    for name in NOTEBOOKS:
        notebook = _read_notebook(name)
        for index, cell in enumerate(notebook["cells"]):
            if cell["cell_type"] != "code":
                continue
            source = _source(cell)
            tree = ast.parse(source, filename=f"{name}:cell-{index}")
            assert not any(
                isinstance(node, ast.Assert) for node in ast.walk(tree)
            ), f"{name} cell {index} contains an assert statement"
            for node in ast.walk(tree):
                if not isinstance(node, ast.Call):
                    continue
                function = node.func
                is_display = (
                    isinstance(function, ast.Name) and function.id == "display"
                ) or (
                    isinstance(function, ast.Attribute) and function.attr == "display"
                )
                assert not is_display, (
                    f"{name} cell {index} explicitly calls display()"
                )


def test_root_case_notebooks_use_only_direct_packaged_imports():
    for name in NOTEBOOKS:
        source = _code_source(_read_notebook(name))
        assert "notebook_utils" not in source
        assert "sys.path.insert" not in source

        imported_roots = set()
        for node in ast.walk(ast.parse(source, filename=name)):
            if isinstance(node, ast.Import):
                imported_roots.update(
                    alias.name.split(".", 1)[0] for alias in node.names
                )
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported_roots.add(node.module.split(".", 1)[0])
        assert imported_roots <= ALLOWED_NOTEBOOK_IMPORT_ROOTS, (
            f"{name} imports undeclared packages: "
            f"{sorted(imported_roots - ALLOWED_NOTEBOOK_IMPORT_ROOTS)}"
        )


def test_root_case_code_keeps_runtime_controls_without_audit_scaffolding():
    forbidden_fragments = (
        "importlib.metadata",
        "platform.",
        "sys.version",
        "__version__",
        "hashlib",
        "sha256",
        "disk_usage",
        "provenance",
        "audit_table",
        "load_json",
        "ipython.display",
        "locate_repo_root",
        "run_log",
        "savefig(",
        "print(",
    )
    for name in NOTEBOOKS:
        source = _code_source(_read_notebook(name))
        assert "SEED = 42" in source, f"{name} lost seed 42"
        assert '"OMP_NUM_THREADS"' in source, f"{name} lost single-thread setup"
        assert '"NUMBA_NUM_THREADS"' in source, f"{name} lost single-thread setup"
        assert "plt.show()" in source, f"{name} no longer shows its figure directly"
        compact_source = source.casefold()
        for fragment in forbidden_fragments:
            assert fragment.casefold() not in compact_source, (
                f"{name} retains audit or logging scaffold {fragment!r}"
            )


def test_root_case_code_cells_are_unexecuted():
    for name in NOTEBOOKS:
        notebook = _read_notebook(name)
        for index, cell in enumerate(notebook["cells"]):
            if cell["cell_type"] == "code":
                assert cell.get("execution_count") is None, (
                    f"{name} cell {index} retains execution_count"
                )


def test_root_case_outputs_are_only_the_approved_retained_pngs():
    for name in NOTEBOOKS:
        notebook = _read_notebook(name)
        observed = {}
        for index, cell in enumerate(notebook["cells"]):
            outputs = cell.get("outputs", [])
            for output in outputs:
                assert output.get("output_type") == "display_data", (
                    f"{name} cell {index} has a non-display output"
                )
                data = output.get("data", {})
                assert set(data) == {"image/png"}, (
                    f"{name} cell {index} has non-PNG output data"
                )
                encoded = data["image/png"]
                png = _decode_png(encoded)
                with Image.open(io.BytesIO(png)) as image:
                    image.verify()
                observed[index] = observed.get(index, 0) + 1
        assert observed == PNG_CELL_COUNTS[name], (
            f"{name} retained PNG cells changed"
        )


def test_root_case_cells_contain_no_machine_paths_or_execution_controls():
    for name in NOTEBOOKS:
        notebook = _read_notebook(name)
        for index, cell in enumerate(notebook["cells"]):
            text = "\n".join(_cell_texts(cell))
            for description, pattern in (*MACHINE_PATTERNS, *FORBIDDEN_CONTROL_PATTERNS):
                assert not pattern.search(text), (
                    f"{name} cell {index} contains forbidden {description}"
                )


def test_root_case_builders_match_published_ordered_cell_sources():
    for ordinal, name in enumerate(NOTEBOOKS):
        builder_path, function_name = BUILDER_SPECS[name]
        builder = _load_builder(builder_path, f"revise_case_builder_{ordinal}")
        expected = getattr(builder, function_name)()
        observed = _read_notebook(name)
        assert _ordered_cell_sources(observed) == _ordered_cell_sources(expected), (
            f"{name} differs from {function_name}() cell sources"
        )


def test_case_markdown_exposes_the_same_reading_contract():
    for name in NOTEBOOKS:
        notebook = _read_notebook(name)
        markdown = "\n".join(
            _source(cell)
            for cell in notebook["cells"]
            if cell["cell_type"] == "markdown"
        )
        for requirement in MARKDOWN_REQUIREMENTS:
            assert requirement in markdown, (
                f"{name} markdown is missing {requirement!r}"
            )


def test_global_anchoring_diagnostics_are_scoped_by_case():
    no_complete_ga_rerun = (
        "SlideSeq_mouse_colon_sp_SVC.ipynb",
        "SlideSeq_mouse_olfactory_bulb_sp_SVC.ipynb",
        "StereoSeq_zebrafish_5hpf_sp_SVC.ipynb",
        "Visium_sc_SVC_mouse_brain.ipynb",
    )
    complete_ga_diagnostics = (
        "MERFISH_Allen_VISp_sc_SVC_cluster.ipynb",
        "osmFISH_sc_SVC_cluster.ipynb",
    )
    for name in no_complete_ga_rerun:
        source = _code_source(_read_notebook(name))
        assert "preprocess_data(" not in source, f"{name} reruns complete GA"
        assert "OTKernel.annotate(" not in source, f"{name} reruns complete GA"
    for name in complete_ga_diagnostics:
        source = _code_source(_read_notebook(name))
        assert "preprocess_data(" in source, f"{name} lost complete GA setup"
        assert "OTKernel.annotate(" in source, f"{name} lost complete GA annotation"


def test_case_specific_input_and_figure_contracts_are_direct():
    for name in DIRECT_UMAP_NOTEBOOKS:
        source = _code_source(_read_notebook(name))
        for call in (
            "sc.pp.normalize_total(",
            "sc.pp.log1p(",
            "sc.pp.pca(",
            "sc.pp.neighbors(",
            "sc.tl.umap(",
            "sc.pl.umap(",
        ):
            assert call in source, f"{name} lost direct {call}"
        assert "target_sum=1e4" in source
        assert "n_comps=30" in source
        assert "n_neighbors=15" in source
        assert "n_pcs=30" in source
        assert "random_state=SEED" in source
        for removed_helper in (
            "assigned_labels",
            "highly_variable_genes",
            "independent_umap",
            "MAX_UMAP_OBSERVATIONS",
            "max_observations",
            "plot_independent_panels",
        ):
            assert removed_helper not in source, (
                f"{name} retains {removed_helper}"
            )

    cosmx = _code_source(_read_notebook("CosMx_SMI_267T_not_sp_SVC.ipynb"))
    assert 'os.environ["REVISE_COSMX_CASE_ROOT"]' in cosmx
    assert 'os.environ.get("REVISE_COSMX_CASE_ROOT"' not in cosmx
    assert "run_reconstruction" not in cosmx

    for name in (
        "SlideSeq_mouse_colon_sp_SVC.ipynb",
        "SlideSeq_mouse_olfactory_bulb_sp_SVC.ipynb",
    ):
        source = _code_source(_read_notebook(name))
        assert source.count("subprocess.run(") == 2
        assert 'svc_level1 = svc.obs["Level1"].astype(str)' in source
        assert "svc_level1.loc[raw_for_umap.obs_names]" in source
        assert "formal_ga_level1" not in source

    stereo = _code_source(
        _read_notebook("StereoSeq_zebrafish_5hpf_sp_SVC.ipynb")
    )
    assert stereo.count("subprocess.run(") == 2
    assert 'svc.obs.loc[' in stereo
    assert 'raw_for_umap.obs_names, "celltype_new"' in stereo
    assert stereo.count("sc.pl.umap(") == 3

    visium = _code_source(
        _read_notebook("Visium_sc_SVC_mouse_brain.ipynb")
    )
    assert visium.count("subprocess.run(") == 1
    assert "urlretrieve" not in visium
    assert "run_manifest" not in visium
    assert 'raw_reference = ad.read_h5ad(RAW_REFERENCE_PATH)' in visium
    assert 'spatial.uns["all_cells_in_spot"] = all_cells_in_spot' in visium
    assert 'pm_on_cell.to_csv(PM_ON_CELL_PATH)' in visium
    assert (
        '"python reconstruct.py --config configs/application/Visium.yaml"'
        in visium
    )
    assert "sc.tl.umap(" not in visium
    assert "sc.pl.umap(" not in visium
    assert "independent_umap" not in visium

    allen = _code_source(
        _read_notebook("MERFISH_Allen_VISp_sc_SVC_cluster.ipynb")
    )
    assert "= ga_spatial.copy()" in allen
    assert 'for column, key in enumerate(("Level2", "SVC_cluster"))' in allen

    osmfish_notebook = _read_notebook("osmFISH_sc_SVC_cluster.ipynb")
    osmfish_markdown = "\n".join(
        _source(cell)
        for cell in osmfish_notebook["cells"]
        if cell["cell_type"] == "markdown"
    ).casefold()
    assert "sidecar" in osmfish_markdown
    assert "not passed" in osmfish_markdown
    assert "not reconstruction truth" in osmfish_markdown
    osmfish_source = _code_source(osmfish_notebook)
    assert "= ga_spatial.copy()" in osmfish_source
