# Reconstruction Impact output and test contract

Each route writes to `output/reconstruction_impact/<sample_id>/` (or the
`REVISE_ANALYSIS_OUTPUT_ROOT` override). The root contains
`resolved_config.yaml`, `manifest.json`, and `input_audit.csv`.

`st_unit/` contains a partition summary, one Hungarian mapping and contingency
table per scientifically applicable comparison edge, and per-unit assignments. Every partition output
contains `comparison_edge`; the summary also records the chosen resolution,
resolution source, pairing audit, and shared-gene count.

`spatial/` contains scale and support-selection audits, a matched-cohort
unit-to-window map, a separate full-Level1 anatomy unit-to-window map and
anatomy candidate map, and per-window Raw / Recon / delta `Neff`. Its visible
metric layers have separate summaries:

- `anatomy_context_summary.csv` for full-Level1 units and tissue-window area;
- `cluster_change_by_anatomy.csv` for globally matched unit change and its
  window distribution;
- `diversity_by_anatomy.csv` for Raw, reconstructed, and delta `Neff` median/IQR;
- `region_extent_by_anatomy.csv` for exact Region windows, area, and unit
  coverage.

Full threshold and scale sensitivity tables remain available. Main Region area
uses non-overlapping 40 um windows, so each selected window contributes 1600
um2; notebooks additionally display area in mm2.

Automated tests cover strict pairing/reordering, label permutation, Hungarian
metric complements, rare and unmatched clusters, Level1-ARI resolution choice,
fixed within-Level1 resolution, deterministic feature selection, coordinate to
micron conversion, boundary window assignment, support-selection knees, `Neff`,
invalid windows, anatomy candidates/context denominators, median/IQR summaries,
globally mapped anatomy-stratified changes, Region area/fractions, config
carrier-role validation, and notebook content boundaries. Existing
canonical case notebooks are protected separately and are not modified by this
work.

The default real-data acceptance run is P2CRC Xenium Fibroblast full cohort plus
P1CRC VisiumHD's deterministic same-ID 30k partition and spatial-impact cohort,
with full-tissue Level1 anatomy context. The VisiumHD notebook switch can run
the full paired cohort when required. Successful execution demonstrates reproducible
representation/partition and Region calculations only; it is not biological or
mechanistic validation.
