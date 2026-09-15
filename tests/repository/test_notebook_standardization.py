"""Presentation refactors must retain visible steps without long helper cells."""
import ast
from pathlib import Path
from reproduce.case.reconstruction_impact import build_notebooks as builder


def test_notebook_helpers_are_imported_not_defined_in_startup():
    for route in ('visium', 'xenium'):
        code = [s for kind, s in builder.cells_for(route) if kind == 'code']
        assert not any(isinstance(n, (ast.FunctionDef, ast.ClassDef)) for n in ast.walk(ast.parse(code[0])))
        assert 'notebook_helpers' in code[0]


def test_region_reliability_is_presented_before_masks():
    for route in ('visium', 'xenium'):
        code = [s for kind, s in builder.cells_for(route) if kind == 'code']
        reliability = next((i for i,s in enumerate(code) if 'plot_threshold_reliability(' in s), None)
        masks = next(i for i,s in enumerate(code) if 'plot_high_diversity(' in s)
        assert reliability is not None and reliability < masks


def test_long_analysis_stages_are_short_visible_delegations():
    """Long audits and plots live in helpers while notebook stages stay readable."""

    assert len(builder.MORAN_SOURCE_CELL.splitlines()) <= 35
    assert len(builder.MORAN_AUDIT_CELL.splitlines()) <= 35
    assert len(builder.MORAN_DISTRIBUTION_PLOT_CELL.splitlines()) <= 35
    assert len(builder.MORAN_SCATTER_PLOT_CELL.splitlines()) <= 35
    assert len(builder.AUCELL_OVERALL_CODE.splitlines()) <= 35

    assert "prepare_moran" in builder.MORAN_SOURCE_CELL
    assert "compute_moran" in builder.MORAN_COMPUTE_CELL
    assert "summarize_moran" in builder.MORAN_AUDIT_CELL
    assert "plot_moran" in builder.MORAN_Q75_PLOT_CELL
    assert "plot_moran" in builder.MORAN_DISTRIBUTION_PLOT_CELL
    assert "plot_moran" in builder.MORAN_SCATTER_PLOT_CELL
    assert "plot_aucell" in builder.AUCELL_OVERALL_CODE


def test_setup_cell_keeps_only_notebook_surface_imports():
    """Heavy analysis and plotting imports belong to the existing helper modules."""

    setup = builder.HELPERS
    for import_name in ("scanpy", "squidpy", "seaborn", "scipy", "metadata", "gc"):
        assert f"import {import_name}" not in setup
    assert "from matplotlib.colors" not in setup
    assert len(setup.splitlines()) <= 40


def test_notebook_block_template_is_consumed_with_source_relative_method_links():
    template = Path(builder.NOTEBOOK_BLOCK_TEMPLATE_PATH)
    assert template.is_file()
    assert "$heading" in template.read_text(encoding="utf-8")
    assert "$purpose" in template.read_text(encoding="utf-8")
    assert "$method_links" in template.read_text(encoding="utf-8")
    assert "../../../docs/design/reconstruction-impact/gene-and-function.md#moran" in builder.MORAN_MARKDOWN
    assert "../../../docs/design/reconstruction-impact/gene-and-function.md#emt-score" in builder.AU_CELL_MARKDOWN


def test_each_route_keeps_the_full_38_figure_inventory():
    for route in ("visium", "xenium"):
        figure_calls = [
            source
            for kind, cell in builder.cells_for(route)
            if kind == "code"
            for source in cell.splitlines()
            if "save_figure(" in source
        ]
        assert len(figure_calls) == 38


def test_foundation_long_blocks_delegate_input_projection_and_pairing_audit():
    """Input preparation stays visible as short calls with a durable audit."""

    helper_source = Path(builder.__file__).with_name("notebook_helpers.py").read_text(
        encoding="utf-8"
    )
    assert "def prepare_route_inputs(" in helper_source
    assert "def prepare_expression_views(" in helper_source
    assert "pairing_audit.csv" in helper_source
    for route in ("visium", "xenium"):
        code = [source for kind, source in builder.cells_for(route) if kind == "code"]
        input_cells = [source for source in code if "prepare_route_inputs(" in source]
        expression_cells = [source for source in code if "prepare_expression_views(" in source]
        assert input_cells and expression_cells
        assert all(len(source.splitlines()) <= 15 for source in input_cells + expression_cells)


def test_notebook_nodes_use_post_run_placeholders_and_number_region_evidence():
    for route in ("visium", "xenium"):
        cells = builder.cells_for(route)
        source = "\n".join(text for _kind, text in cells)
        assert "impact-records:" not in source
        assert "Overall observation" not in source
        assert "impact-node:1.6" in source
        assert "impact-node:2.4" in source
        assert "impact-node:3.9" in source
        assert "impact-node:4.5" in source
        assert "display(Markdown(f\"" not in source
        region = "\n".join(text for _kind, text in cells if "Regions:" in text or "4." in text)
        assert "### 4.1 Threshold reliability" in region
        assert "### 4.2 Continuous reconstructed-Neff state and mask" in region
        assert "### 4.3 Region extent by anatomy" in region
        assert "### 4.4 Scale sensitivity" in region
        assert "### 4.5 Region conclusion" in region


def test_generated_notebooks_keep_science_and_shorten_heading_prose():
    """Reader cleanup changes Markdown only, with no alternate calculation path."""
    import json
    for route, filename in (("visium", "VisiumHD_sp_SVC_Reconstruction_Impact.ipynb"),
                            ("xenium", "Xenium_sc_SVC_Fibroblast_Reconstruction_Impact.ipynb")):
        notebook = json.loads((Path(builder.__file__).parent / filename).read_text())
        expected = [text for kind, text in builder.cells_for(route) if kind == "code"]
        actual = ["".join(cell["source"]) for cell in notebook["cells"] if cell["cell_type"] == "code"]
        assert expected == actual
        headings = ["".join(cell["source"]) for cell in notebook["cells"]
                    if cell["cell_type"] == "markdown" and "".join(cell["source"]).startswith("### ")]
        assert headings
        assert all("Batch output:" not in text for text in headings)
        assert all(len(text.splitlines()) <= 5 for text in headings)
