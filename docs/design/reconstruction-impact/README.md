# Reconstruction Impact analysis

This post-reconstruction analysis describes two endpoints for observation-paired
spatial units: **ST-unit Expression-Space Change** and **Region**. It is not a
reconstruction route, does not change the reconstruction engine, and does not
establish biological mechanism or ground-truth biological validity.

The branch has exactly two route notebooks: P1CRC VisiumHD `sp-SVC` and P2CRC
Xenium Fibroblast `sc-SVC`. Raw and reconstructed spatial-side inputs must have
exactly the same observation IDs. The sc-SVC expression-side H5AD is not
observation-paired and is excluded from ST-unit matching.

## Partition definition

For a multi-Level1 scope, candidate Leiden resolutions `0.3`, `0.5`, and `0.8`
are evaluated on the Raw representation. The selected value has the highest ARI
against Raw `Level1`; equal ARIs choose the lower resolution. That one value is
then applied to Raw and reconstructed expression graphs. Within a single Level1
parent, the fixed resolution is `0.5`.

sp-SVC reports only Raw-to-Reconstructed-expression. sc-SVC audits the Raw and
spatial-carrier expression matrices; when they are identical, it reports only
Raw Leiden-to-Final-SVC rather than labelling the identity edge as an impact.
Hungarian matching makes cluster labels comparable; matched accuracy/F1 are
primary and ARI/AMI/NMI/VI are agreement diagnostics.

## Spatial definition

Spatial windows are non-overlapping squares anchored at the full Raw-context
coordinate minimum; the paired diversity cohort and the anatomy map use this
same grid. Coordinates are first converted to microns. VisiumHD uses the
verified 0.273808 um per coordinate and Xenium uses 0.2125 um per morphology
pixel ([Xenium file format](https://cf.10xgenomics.com/supp/xenium/xenium_documentation.html);
[Space Ranger spatial outputs](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/spatial-outputs)).
Both routes use an explicit 8 um cell-equivalent scale, with a 5x = 40 um main
window and 1/2/3/5/7/10x sensitivity. Parent support is selected at the main
scale from the retained-unit curve's chord-distance knee, then held fixed for
scale sensitivity. The Region is valid windows with reconstructed `Neff >= 2`;
threshold sensitivity spans 1.5--3.5. Its area is selected-window count times
window area, with valid windows as the area-fraction denominator. Invalid
windows remain `NaN`.

Anatomy is a separate Level1-derived window assignment on the same spatial
grid. It records binary Tumor and Normal (`Intestinal Epithelial`) candidates,
their 0/1/2 union score, and the categorical Tumor, Normal, Interface, or Other
context. It is written and plotted independently of reconstruction-driven
diversity, and never changes partitions, `Neff`, or high-diversity Region
membership. `Intestinal Epithelial` is the current source label for Normal;
this is a spatial convention, not a new single-cell annotation claim.

See [outputs-and-test-plan.md](outputs-and-test-plan.md) for the output
contract and verification boundary.
