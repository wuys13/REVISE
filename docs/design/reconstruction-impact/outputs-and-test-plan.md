# Reconstruction Impact output and test contract

Each route writes to `output/reconstruction_impact/<sample_id>/` (or the
`REVISE_ANALYSIS_OUTPUT_ROOT` override). The root contains
`resolved_config.yaml`, `manifest.json`, and `input_audit.csv`.

`st_unit/` contains a partition summary, one Hungarian mapping and contingency
table per comparison edge, and per-unit assignments. Every partition output
contains `comparison_edge`; the summary also records the chosen resolution,
resolution source, pairing audit, and shared-gene count.

`spatial/` contains a matched-cohort unit-to-window map, a separate
full-Level1 anatomy unit-to-window map and anatomy-region map, per-window Raw /
Recon / delta `Neff`, Region summaries, threshold and scale sensitivity tables.
The maps are Raw `Neff`, Recon `Neff`, centered delta `Neff`, the
high-diversity Region mask, and the separate anatomy map. Main Region area uses
non-overlapping 40 um windows, so each selected window contributes 1600 um2.

Automated tests cover strict pairing/reordering, label permutation, Hungarian
metric complements, rare and unmatched clusters, Level1-ARI resolution choice,
fixed within-Level1 resolution, deterministic feature selection, boundary window
assignment, `Neff`, invalid windows, Region area/fractions, config carrier-role
validation, and notebook content boundaries. Existing canonical case notebooks
are protected separately and are not modified by this work.

The real-data acceptance run is P2CRC Xenium Fibroblast full cohort plus P1CRC
VisiumHD's deterministic 30k partition sample and full matched-cohort spatial
window calculation. Successful execution demonstrates reproducible
representation/partition and Region calculations only; it is not biological or
mechanistic validation.
