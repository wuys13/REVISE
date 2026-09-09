# Reconstruction Impact analysis (internal v4)

This post-reconstruction analysis describes paired representation and local
spatial-pattern changes. It does not alter the reconstruction engine and does
not establish biological mechanism, ground truth, or clinical relevance.

Two notebooks cover P1CRC VisiumHD `sp-SVC` and P2CRC Xenium `sc-SVC`. Both
evaluate the internal diversity of Fibroblast, Mono_Macro and T. Xenium's
`expr.h5ad` is recorded as an expression-side carrier only: it is excluded from
observation pairing, Leiden, spatial coordinates, and diversity calculations.

## Partition evidence

The analysis deliberately distinguishes two questions.

For VisiumHD sp-SVC, QC is defined only on Raw counts (`min_genes=50`,
`min_cells=3`, mitochondrial genes excluded) and the retained observation IDs
and genes are applied to both paired carriers. The canonical VisiumHD
normalize/log/HVG sequence is run on Raw only, and that Raw-derived feature set
is shared with reconstructed expression. Reconstruction therefore cannot
define the Raw feature space. The manifest records the feature rule, explicit
igraph Leiden backend, two iterations, and seed 42. Xenium keeps its existing
spatial-carrier partition contract.

- The **complexity diagnostic** holds the Raw resolution definition fixed. It
  shows whether reconstructed expression produces more clusters and altered
  cluster-size/split patterns. VisiumHD global Raw chooses from `0.3/0.5/0.8`
  by its Level1 ARI; parent Raw uses `0.5`.
- The **matched-K comparison** independently scans the other graph (or Raw for
  sc-SVC) to match the target cluster count. It then performs one Hungarian
  alignment and reports ST-unit change, balanced change, and ARI. The global
  VisiumHD comparison additionally reports change fractions for every Raw
  Level1 label with Wilson intervals. If the closest count differs by more than
  one, the result is retained as `unmatched_cluster_complexity`, not a headline
  change estimate.

Thus “partition is finer” and “assignment changed” cannot be conflated. Global
Level1 and spatial results reuse one global matched-K mapping. The three
parent blocks instead report **parent-internal reassignment** between internal
clusters; those percentages are not Level1 identity changes. No result is
re-matched inside an anatomy region.

## Parent diversity and Regions

Each parent has three Raw descriptions with different roles.

- **Raw Level1** is one uniform parent label, so Kobs, Neff and evenness are one
  by construction. It is retained only as an unexpanded descriptive baseline.
- **Raw Leiden** is expression-derived and matched as closely as possible to
  reconstructed/final K. It is the controlled partition baseline.
- **Raw Level2** is mapped from the original Raw expression after Level1 parent
  selection. VisiumHD uses its route-native POT assignment; Xenium uses TACCO
  with the P2CRC reference subset. Reconstructed-carrier Level2 labels are not
  used to construct this baseline.

Kobs counts labels observed locally, Neff is the abundance-aware effective
cluster count, and evenness (`Neff / Kobs`) describes balance conditional on
richness. Each metric is shown as a 2x3 matrix: Raw Leiden, reconstructed, and
their delta on the first row; Raw Level2, the same reconstructed state, and
their delta on the second row.

Candidate square sides are 16, 24, 32, 40, 56, and 80 um. A parent chooses one
only from its coordinates and occupancy: it retains windows with at least four
parent units and takes the chord-distance knee of retained units against log
window side. The grid is anchored at the full Raw Level1 coordinate minimum.
Within every valid window, Raw Leiden, Raw Level2 and Recon use the same four
sampled parent units for each of 200 deterministic draws. Kobs, entropy, Neff
and evenness are calculated per draw before summarisation. This controls
density-driven diversity inflation and keeps both delta baselines paired.

The notebook headline Region derives from rarefied reconstructed Neff. Its
threshold is estimated from a survival-curve segmented breakpoint and 500
window bootstraps. Insufficient, edge, or wide-interval estimates are labelled
`no_stable_threshold` and produce no binary mask. Delta-threshold artifacts are
retained for audit compatibility but are not shown as a second headline Region.
Region fractions are therefore not cross-platform effect sizes.

Anatomy remains a separate fixed Level1 context: Tumor, Normal (current source
`Intestinal Epithelial`), Interface, and Other. It never tunes K, scale, or a
threshold. VisiumHD coordinates use 0.273808 um per coordinate and Xenium uses
0.2125 um per morphology pixel ([Xenium file format](https://cf.10xgenomics.com/supp/xenium/xenium_documentation.html),
[Space Ranger spatial outputs](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/spatial-outputs)).

The VisiumHD notebook exposes `USE_FULL_VISIUMHD_COHORT`. Its default first
samples 30,000 global paired IDs and at most 30,000 IDs per parent, then reports
the Raw-QC retained count explicitly; setting it to `True` removes the sampling
limits but not Raw QC. Mono_Macro and T are already below that sampling limit.

See [outputs-and-test-plan.md](outputs-and-test-plan.md) for the artifact and
test contract.
