# Reconstruction Impact analysis (internal v2)

This post-reconstruction analysis describes paired representation and local
spatial-pattern changes. It does not alter the reconstruction engine and does
not establish biological mechanism, ground truth, or clinical relevance.

Two notebooks cover P1CRC VisiumHD `sp-SVC` and P2CRC Xenium `sc-SVC`. Both
evaluate the internal diversity of Fibroblast, Mono_Macro and T. Xenium's
`expr.h5ad` is recorded as an expression-side carrier only: it is excluded from
observation pairing, Leiden, spatial coordinates, and diversity calculations.

## Partition evidence

The analysis deliberately distinguishes two questions.

- The **complexity diagnostic** holds the Raw resolution definition fixed. It
  shows whether reconstructed expression produces more clusters and altered
  cluster-size/split patterns. VisiumHD global Raw chooses from `0.3/0.5/0.8`
  by its Level1 ARI; parent Raw uses `0.5`.
- The **matched-K comparison** independently scans the other graph (or Raw for
  sc-SVC) to match the target cluster count. It then performs one Hungarian
  alignment and reports ST-unit change, balanced change, ARI, and Level1 change
  fractions with Wilson intervals. If the closest count differs by more than
  one, the result is retained as `unmatched_cluster_complexity`, not a headline
  change estimate.

Thus “partition is finer” and “assignment changed” cannot be conflated. The
Level1 and spatial results reuse that one global matched-K mapping; they never
re-match inside a region.

## Parent diversity and Regions

Each parent has two separate diversity baselines.

- **Revealed subtype diversity:** Raw is one uniform Level1 label, so its Neff
  is one by construction. `Raw-Level1 / Recon subtype / Recon−1` describes
  internal structure exposed after reconstruction, not a fair effect size.
- **Reconstruction-associated diversity change:** Raw and reconstructed/final
  labels are the matched-K partitions. `Raw-Leiden / Recon / Delta` is the
  comparable local change view and supplies the Gain Region input.

Candidate square sides are 16, 24, 32, 40, 56, and 80 um. A parent chooses one
only from its coordinates and occupancy: it retains windows with at least four
parent units and takes the chord-distance knee of retained units against log
window side. The grid is anchored at the full Raw Level1 coordinate minimum.
Within every valid window, Raw and Recon use the same four sampled parent units
for each of 200 deterministic draws; Neff is their draw mean. This controls
density-driven diversity inflation.

State Region derives from rarefied reconstructed Neff. Gain Region derives from
positive matched-K Delta Neff. Each threshold is estimated independently from
a survival-curve segmented breakpoint and 500 window bootstraps. Insufficient,
edge, or wide-interval estimates are labelled `no_stable_threshold` and produce
no binary mask. Region fractions are therefore not cross-platform effect sizes.

Anatomy remains a separate fixed Level1 context: Tumor, Normal (current source
`Intestinal Epithelial`), Interface, and Other. It never tunes K, scale, or a
threshold. VisiumHD coordinates use 0.273808 um per coordinate and Xenium uses
0.2125 um per morphology pixel ([Xenium file format](https://cf.10xgenomics.com/supp/xenium/xenium_documentation.html),
[Space Ranger spatial outputs](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/spatial-outputs)).

The VisiumHD notebook exposes `USE_FULL_VISIUMHD_COHORT`. Its default samples
30,000 global paired IDs and at most 30,000 IDs per parent; setting it to `True`
removes both limits. Mono_Macro and T are already below that parent limit.

See [outputs-and-test-plan.md](outputs-and-test-plan.md) for the artifact and
test contract.
