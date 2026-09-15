"""Structural contracts for the notebook-first reconstruction-impact workflow.

These checks deliberately inspect the ordered notebook cells instead of only
checking that a few words occur somewhere in the combined source.  The
notebooks are the scientific record for this iteration: the reader must be
able to follow the four-layer argument, and the calculations which feed a
later layer must occur before that layer is presented.
"""

from __future__ import annotations

import ast
import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
NOTEBOOK_DIR = ROOT / "reproduce" / "case" / "reconstruction_impact"
NOTEBOOKS = (
    "VisiumHD_sp_SVC_Reconstruction_Impact.ipynb",
    "Xenium_sc_SVC_Fibroblast_Reconstruction_Impact.ipynb",
)


def _notebook(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _source(cell: dict) -> str:
    source = cell.get("source", "")
    return "".join(source) if isinstance(source, list) else source


def _cells(path: Path) -> list[dict]:
    return _notebook(path)["cells"]


def _all_source(cells: list[dict]) -> str:
    return "\n".join(_source(cell) for cell in cells)


def _heading_index(cells: list[dict], pattern: str, *, level: int = 2) -> int:
    expression = re.compile(
        rf"^{'#' * level}\s+.*(?:{pattern})",
        flags=re.IGNORECASE | re.MULTILINE,
    )
    for index, cell in enumerate(cells):
        if cell.get("cell_type") == "markdown" and expression.search(_source(cell)):
            return index
    raise AssertionError(f"missing level-{level} heading matching {pattern!r}")


def _section(cells: list[dict], heading: str) -> str:
    """Return the source between one top-level section and the next."""

    start = _heading_index(cells, heading)
    next_top_level = re.compile(r"^##\s+", re.MULTILINE)
    for index in range(start + 1, len(cells)):
        if cells[index].get("cell_type") == "markdown" and next_top_level.search(
            _source(cells[index])
        ):
            return "\n".join(_source(cell) for cell in cells[start:index])
    return "\n".join(_source(cell) for cell in cells[start:])


def _first_code_cell(cells: list[dict], pattern: str) -> int:
    expression = re.compile(pattern, flags=re.IGNORECASE | re.MULTILINE)
    for index, cell in enumerate(cells):
        if cell.get("cell_type") == "code" and expression.search(_source(cell)):
            return index
    raise AssertionError(f"missing executable code matching {pattern!r}")


def _heading_positions(
    cells: list[dict], pattern: str, *, levels: tuple[int, ...] = (3,)
) -> list[tuple[int, int]]:
    """Return positions of matching heading lines, not prose occurrences."""

    expression = re.compile(
        rf"^(?:{'|'.join('#' * level for level in levels)})\s+.*(?:{pattern})",
        flags=re.IGNORECASE,
    )
    positions: list[tuple[int, int]] = []
    for cell_index, cell in enumerate(cells):
        if cell.get("cell_type") != "markdown":
            continue
        for line_index, line in enumerate(_source(cell).splitlines()):
            if expression.search(line):
                positions.append((cell_index, line_index))
    return positions


def _assert_ordered_headings(
    cells: list[dict], tokens: tuple[tuple[str, ...], ...], *, levels: tuple[int, ...] = (3,)
) -> None:
    positions: list[tuple[int, int]] = []
    for aliases in tokens:
        matches = []
        for alias in aliases:
            matches.extend(_heading_positions(cells, alias, levels=levels))
        assert matches, f"missing heading matching one of {aliases!r}"
        positions.append(min(matches))
    assert positions == sorted(positions), f"out-of-order headings: {positions}"


def _call_cells(cells: list[dict], names: set[str]) -> list[int]:
    """Find calls executed by notebook cells, ignoring helper definitions.

    Notebook loops are deliberately traversed: a call inside a ``for`` is
    still an executed scientific block.  Function and lambda bodies are
    skipped because helper definitions must not satisfy execution-order
    checks.
    """

    class Calls(ast.NodeVisitor):
        def __init__(self) -> None:
            self.found = False

        def visit_FunctionDef(self, node: ast.FunctionDef) -> None:  # noqa: N802
            return

        def visit_AsyncFunctionDef(self, node: ast.AsyncFunctionDef) -> None:  # noqa: N802
            return

        def visit_Lambda(self, node: ast.Lambda) -> None:  # noqa: N802
            return

        def visit_Call(self, node: ast.Call) -> None:  # noqa: N802
            function = node.func
            name = function.id if isinstance(function, ast.Name) else function.attr if isinstance(function, ast.Attribute) else None
            if name in names:
                self.found = True
            self.generic_visit(node)

    found_cells: list[int] = []
    for index, cell in enumerate(cells):
        if cell.get("cell_type") != "code":
            continue
        visitor = Calls()
        visitor.visit(ast.parse(_source(cell), filename=f"cell-{index}"))
        if visitor.found:
            found_cells.append(index)
    return found_cells


def test_reconstruction_impact_notebooks_have_four_layers_in_reading_order():
    """The notebook argument is foundation -> change -> localization -> Region."""

    for name in NOTEBOOKS:
        cells = _cells(NOTEBOOK_DIR / name)
        top_level = (
            "comparison foundation|comparison basis|what are we comparing",
            "overall changes|what changed",
            "spatial localization|where .*change",
            "regions?",
        )
        positions = [_heading_index(cells, pattern) for pattern in top_level]
        assert positions == sorted(positions), (name, positions)

        # These are semantic section contracts.  They make it possible to
        # rename a chart without silently moving Moran/AUCell behind Region.
        foundation = _section(cells, top_level[0])
        _assert_ordered_headings(
            cells,
            (
                ("route semantics", "carrier audit"),
                ("cohort and coordinate", "coordinate audit"),
                ("expression views", "gene availability", "full-gene"),
                ("analysis-specific preprocessing", "partition preprocessing"),
            ),
        )

        overall = _section(cells, top_level[1])
        _assert_ordered_headings(
            cells,
            (
                ("partition complexity",),
                ("matched-complexity", "matched-k", "matched K"),
                ("change by cell type", "parent-internal", "level1"),
                ("moran",),
                ("aucell",),
                ("summary of observed changes", "overall change summary"),
            ),
        )

        localization = _section(cells, top_level[2])
        _assert_ordered_headings(
            cells,
            (
                ("anatomy context", "level1 anatomy"),
                ("changed-unit spatial", "assignment.*localization"),
                ("emt.*spatial", "spatial fields"),
                ("internal-state baselines", "window decision"),
                ("kobs",),
                ("neff",),
                ("evenness",),
                ("anatomy-stratified", "anatomy summary"),
            ),
        )

        region = _section(cells, top_level[3])
        _assert_ordered_headings(
            cells,
            (
                ("threshold reliability",),
                ("continuous reconstructed-neff state and mask",),
                ("region extent", "extent by anatomy"),
                ("scale sensitivity", "sensitivity"),
                ("region conclusion",),
            ),
            levels=(2, 3),
        )

        summary_positions = _heading_positions(
            cells, "cross-parent summary|evidence boundary", levels=(2, 3)
        )
        assert summary_positions
        assert min(summary_positions)[0] >= positions[-1]


def test_scientific_calculation_order_is_explicit_and_not_hidden_in_setup():
    """Each layer's executable work is separate and feeds the next layer."""

    for name in NOTEBOOKS:
        cells = _cells(NOTEBOOK_DIR / name)
        partition = min(_call_cells(cells, {"run_partition_analysis"}))
        moran = min(_call_cells(cells, {"compute_moran"}))
        aucell = min(_call_cells(cells, {"compute_emt"}))
        anatomy = min(_call_cells(cells, {"compute_anatomy_regions"}))
        window = min(_call_cells(cells, {"compute_rarefied_window_diversity"}))
        threshold = min(_call_cells(cells, {"select_region_threshold"}))

        assert partition < moran < anatomy < window < threshold
        assert partition < aucell < anatomy

        # A single setup cell must not hide all scientific stages.  The
        # ordered calls should be inspectable as separate notebook blocks.
        assert len({partition, moran, aucell, anatomy, window, threshold}) >= 5


def test_moran_and_aucell_are_overall_results_before_spatial_location():
    for name in NOTEBOOKS:
        cells = _cells(NOTEBOOK_DIR / name)
        overall = _heading_index(cells, "overall changes|what changed")
        localization = _heading_index(cells, "spatial localization|where .*change")
        region = _heading_index(cells, "regions?")
        for pattern in (
            {"compute_moran"},
            {"compute_emt"},
        ):
            calculation = min(_call_cells(cells, pattern))
            assert overall < calculation < localization < region


def test_existing_evidence_formats_are_retained_and_expanded():
    for name in NOTEBOOKS:
        cells = _cells(NOTEBOOK_DIR / name)
        source = _all_source(cells)

        # Existing 2 x 3 state/delta matrices and 1 x 2 Neff/State panels are
        # part of the notebook evidence format, not optional report styling.
        helpers = (NOTEBOOK_DIR / "notebook_helpers.py").read_text()
        assert re.search(r"subplots\(\s*2\s*,\s*3", helpers)
        assert re.search(r"subplots\(\s*1\s*,\s*2", helpers)
        for label in (
            "Raw Leiden",
            "Reconstructed",
            "Recon − Raw Leiden",
            "Raw Level2",
            "Recon − Raw Level2",
        ):
            assert label in source or label in helpers
        assert "Reconstructed Neff" in source or "State Region" in source
        assert "High-diversity Region" in source or "State Region" in source

        # There must still be one explicit comparison for each metric, even if
        # the implementation factors the plotting helper out.
        for metric in ("k_obs", "neff", "evenness"):
            assert metric in source


def test_notebooks_use_full_gene_views_without_an_adapter_or_cache_bypass():
    for name in NOTEBOOKS:
        cells = _cells(NOTEBOOK_DIR / name)
        source = _all_source(cells)
        lowered = source.lower()
        assert "full gene" in lowered or "full_gene" in lowered
        assert "gene availability" in lowered or "gene_status" in lowered
        assert "raw" in lowered and "recon" in lowered

        # Science stays visible in the notebook.  It must not be replaced by
        # the batch adapter or a serialized checkpoint branch.
        for forbidden in (
            "reconstruction_impact_adapter",
            "pathway_activity_adapter",
            "batch_analyze",
            "analysisinputs",
            "reconstruction_impact_cache_path",
            "cache_path",
            "pickle",
        ):
            assert forbidden not in lowered

        assert "common_genes" not in lowered or "delta" in lowered
        assert "notebook_development" not in lowered
        assert "GENE_PANEL" not in source
        for application_gene in ("SFRP4", "SULF1", '"COMP"'):
            assert application_gene not in source
        assert "revise_analysis_output_root" in lowered


def test_notebook_cells_have_identity_and_each_plot_cell_has_one_save():
    for name in NOTEBOOKS:
        cells = _cells(NOTEBOOK_DIR / name)
        assert all(cell.get("id") for cell in cells)
        assert len({cell["id"] for cell in cells}) == len(cells)
        for cell in cells:
            if cell.get("cell_type") == "code":
                assert _source(cell).count("save_figure(") <= 1


def test_route_specific_baselines_and_sampling_parameters_remain_visible():
    visium_source = _all_source(_cells(NOTEBOOK_DIR / NOTEBOOKS[0]))
    xenium_source = _all_source(_cells(NOTEBOOK_DIR / NOTEBOOKS[1]))

    for marker in (
        "deterministic_same_id_sample",
        "30_000",
        "global_ids",
        "Raw-QC retained units",
        "raw_qc_min_genes",
        "raw_qc_min_cells",
        "Raw-derived",
        "parent-internal",
        "candidate_window_sides_um",
        "rarefaction_draws",
        "threshold_bootstraps",
        "min_parent_units",
    ):
        assert marker in visium_source

    assert "representation_audit" in xenium_source
    assert all(parent in xenium_source for parent in ("Fibroblast", "Mono_Macro", "T"))
    assert "expression_h5ad" in xenium_source


def test_reconstruction_impact_configs_keep_the_existing_route_inputs():
    design_dir = ROOT / "docs" / "design" / "reconstruction-impact"
    assert (design_dir / "README.md").is_file()
    assert (design_dir / "outputs-and-test-plan.md").is_file()
    for name in (
        "reconstruction_impact_visiumhd_p1crc.yaml",
        "reconstruction_impact_xenium_p2crc_fibroblast.yaml",
    ):
        path = ROOT / "configs" / "analysis" / name
        assert path.is_file()
        text = path.read_text(encoding="utf-8")
        assert "schema_version: 1" in text
        assert "final_notebook:" in text
        assert "candidate_window_sides_um" in text
        assert "min_parent_units: 4" in text
        assert "raw_level2_mapping:" in text
        assert "adata_sc_all_reanno.h5ad" in text

    xenium = (
        ROOT / "configs" / "analysis" /
        "reconstruction_impact_xenium_p2crc_fibroblast.yaml"
    ).read_text(encoding="utf-8")
    assert all(parent in xenium for parent in ("Fibroblast", "Mono_Macro", "T"))
    assert "expr.h5ad" in xenium
    assert "method: tacco" in xenium
    assert "value: P2CRC" in xenium

    visium = (
        ROOT / "configs" / "analysis" /
        "reconstruction_impact_visiumhd_p1crc.yaml"
    ).read_text(encoding="utf-8")
    assert "method: pot" in visium
    assert "raw_qc_min_genes: 50" in visium
    assert "raw_qc_min_cells: 3" in visium
    assert "feature_selection: raw_canonical_seurat_v3_shared" in visium


def test_xenium_complexity_names_fixed_final_clusters():
    source = _all_source(_cells(NOTEBOOK_DIR / NOTEBOOKS[1]))
    assert 'Fixed final-cluster K' in source
    assert 'raw_reference_vs_final' in source
    assert 'Recon K at Raw resolution' not in source
