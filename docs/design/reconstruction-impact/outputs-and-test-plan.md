# Reconstruction Impact output and test contract (internal v2)

Each route writes to `output/reconstruction_impact/<sample_id>/` or the
`REVISE_ANALYSIS_OUTPUT_ROOT` override. The root has the resolved config,
manifest (input SHA, sampling, carrier role, selected scale and threshold
status), and a compact input audit.

For each global or parent scope, `st_unit/` contains the matched-K partition
summary, mapping, normalized and absolute contingency, unit assignments,
`matched_k_resolution_sweep.csv`, `complexity_raw_resolution_sweep.csv`, and
`change_by_level1.csv`. Same-resolution complexity diagnostic tables are
separate from the matched-K headline tables.

Each parent `spatial/` directory contains:

- occupancy scale-knee audit and selected physical side;
- paired unit/window assignments, full-Level1 anatomy assignments and map;
- rarefied Raw/Recon/Delta window metrics and anatomy summaries;
- State and Gain threshold bootstrap tables and separate anatomy-stratified
  extent tables;
- all-scale continuous diversity sensitivity.

`region_extent_by_anatomy.csv` refers to State and
`gain_region_extent_by_anatomy.csv` refers to Gain. If a threshold is not
stable, their Region numerators/fractions remain `NaN` rather than being forced
to zero. Area uses exact um2 in artifacts and mm2 in notebooks.

Tests cover strict observation pairing, matching count selection and its tie/
failure status, complexity-vs-change separation, Hungarian metric complements,
Wilson Level1 intervals, coordinate conversion, occupancy-only scale selection,
paired deterministic rarefaction, uniform-Level1 Neff, identical Delta,
unstable thresholds, anatomy isolation, artifact naming, notebook hierarchy,
and unchanged canonical case notebooks. Real execution verifies both source
notebooks without errors before their executed copies replace prior outputs.
