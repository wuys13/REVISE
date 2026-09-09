# Reconstruction Impact output and test contract (internal v4)

Each route writes to `output/reconstruction_impact/<sample_id>/` or the
`REVISE_ANALYSIS_OUTPUT_ROOT` override. The root has the resolved config,
manifest (input SHA, sampling, carrier role, selected scale and threshold
status), the VisiumHD Raw-QC/HVG and Leiden contract, and a compact input audit.

`output/reconstruction_impact/<sample_id>/notebook/<route>.ipynb` is the sole
final notebook for a route: it is the executed result and is the artifact to
open, review, or hand off. `reproduce/case/reconstruction_impact/` contains
only the version-controlled source notebooks and their renderer; it is never a
second final-result location. `run_route_notebook.py` executes a source
notebook, writes it to the configured final path, and verifies that its code
cells match the source.

For each global or parent scope, `st_unit/` contains the matched-K partition
summary, mapping, normalized and absolute contingency, unit assignments,
`matched_k_resolution_sweep.csv`, `complexity_raw_resolution_sweep.csv`, and
`change_by_level1.csv`. Same-resolution complexity diagnostic tables are
separate from the matched-K headline tables.

Each parent `spatial/` directory contains:

- occupancy scale-knee audit and selected physical side;
- paired unit/window assignments, full-Level1 anatomy assignments and map;
- rarefied Raw-Leiden, Raw-Level2, Recon and paired-delta Kobs/Neff/evenness
  window metrics and anatomy summaries;
- State and Gain threshold bootstrap tables and separate anatomy-stratified
  extent tables;
- all-scale continuous diversity sensitivity.

`region_extent_by_anatomy.csv` refers to State and
`gain_region_extent_by_anatomy.csv` refers to Gain. If a threshold is not
stable, their Region numerators/fractions remain `NaN` rather than being forced
to zero. Area uses exact um2 in artifacts and mm2 in notebooks.

Each parent also contains `raw_level2/assignments.csv.gz`,
`raw_level2/posterior.csv.gz`, and `raw_level2/audit.json`. These files prove
that Level2 was mapped from original Raw expression using the configured
single-cell reference and route-native method. The Gain files remain audit
artifacts; only the reconstructed-Neff High-diversity Region appears in the
notebook headline.

Tests cover strict observation pairing, matching count selection and its tie/
failure status, complexity-vs-change separation, Hungarian metric complements,
Raw-defined paired sp-SVC QC, Raw-derived shared HVGs, Wilson Level1 intervals,
coordinate conversion, occupancy-only scale selection,
paired three-assignment rarefaction, uniform-Level1 metrics, identical Delta,
unstable thresholds, anatomy isolation, artifact naming, notebook hierarchy,
and unchanged canonical case notebooks. Real execution verifies both source
notebooks without errors before their executed copies replace prior outputs.
