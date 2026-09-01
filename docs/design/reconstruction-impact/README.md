# Reconstruction Impact analysis

This post-reconstruction analysis describes two endpoints for observation-paired
spatial units: **ST-unit Expression-Space Change** and **Region**. It is not a
reconstruction route, does not change the reconstruction engine, and does not
establish biological mechanism or ground-truth biological validity.

The v0.1 branch covers P1CRC VisiumHD `sp-SVC` and P2CRC Xenium Fibroblast
`sc-SVC`. Raw and reconstructed spatial-side inputs must have exactly the same
observation IDs. The sc-SVC expression-side H5AD is not observation-paired and
is excluded from ST-unit matching.

## Partition definition

For a multi-Level1 scope, candidate Leiden resolutions `0.3`, `0.5`, and `0.8`
are evaluated on the Raw representation. The selected value has the highest ARI
against Raw `Level1`; equal ARIs choose the lower resolution. That one value is
then applied to Raw and reconstructed expression graphs. Within a single Level1
parent, the fixed resolution is `0.5`.

sc-SVC reports Raw-to-Recon-expression, Recon-expression-to-Final-SVC, and
Raw-to-Final-SVC edges. sp-SVC reports only Raw-to-Recon-expression because it
has no separate final SVC cluster carrier. Hungarian matching makes cluster
labels comparable; matched accuracy/F1 are primary and ARI/AMI/NMI/VI are
agreement diagnostics.

## Spatial definition

Spatial windows are non-overlapping squares anchored at the full Raw-context
coordinate minimum; the paired diversity cohort and the anatomy map use this
same grid. Both routes use an explicit 8 um base scale, with a 40 um main
window and 8/16/24/40/56/80 um sensitivity. The Region is the set of valid
windows with reconstructed `Neff >= 2`; threshold sensitivity spans 1.5--3.5.
Its area is the number of windows times the window area, with valid windows as
the area-fraction denominator. Invalid windows remain `NaN`.

Anatomy is a separate Level1-derived window assignment on the same spatial
grid: Tumor, Normal, Interface, or Other. A window with both Tumor and
`Intestinal Epithelial` is Interface; either label alone is Tumor or Normal.
It is written and plotted independently of reconstruction-driven diversity,
and never changes partitions, `Neff`, or high-diversity Region membership.
`Intestinal Epithelial` is the current source label for Normal; this is a
spatial convention, not a new single-cell annotation claim.

See [outputs-and-test-plan.md](outputs-and-test-plan.md) for the output
contract and verification boundary.
