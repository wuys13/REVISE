Concepts
========

REVISE uses one orchestration lifecycle for three reconstruction types. The
route describes the computational contract; it does not by itself prove that a
reconstructed object is a true cell or biologically valid result.

Routes and result types
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 1 2 2

   * - Public type
     - Input rows
     - Public 1.x result type
   * - ``sp-SVC``
     - high-definition bins or pseudo-cells
     - ``sp-SVC``
   * - ``sc-SVC``
     - segmented imaging-ST cells
     - ``sc-SVC``
   * - ``sc-SVC-sr``
     - spot-level observations
     - ``sc-SVC-sr``

sp-SVC and sc-SVC-sr publish ``<output-root>/<sample-name>/SVC.h5ad``.
Standard sc-SVC publishes separate spatial and reference-expression H5AD files
under ``<output-root>/<sample-name>/sc-SVC/<cell-type>/``. Its
``provenance.json`` records both logical roles separately from the internal
route.

Shared lifecycle
----------------

Every full pipeline request follows the same ordered stages:

1. input validation;
2. Global Anchoring (GA);
3. local-unit preparation;
4. graph construction;
5. OT problem construction and Local Refinement (LR);
6. expression update;
7. SVC finalization;
8. optional benchmark evaluation.

``ot.ga.solver`` and ``ot.lr.solver`` select POT or TACCO for their respective
stages. ``--ot-method`` is the convenience control that sets both stages.

Spot super-resolution cell locations and random assignment
-----------------------------------------------------------

A sc-SVC-sr route needs a count of virtual cells per spot. If a curated
``uns["all_cells_in_spot"]`` mapping is absent, the input adapter estimates
counts from spot transcript totals and creates virtual-cell rows.

Segmentation-derived cell centers use the optional standardized
``uns["revise_cell_locations"]`` DataFrame. Its index contains unique
``cell_id`` values, its columns are ``spot_name``, ``x``, and ``y``, and its
cell-to-spot assignments must agree with ``uns["all_cells_in_spot"]``. The
histology-prior preprocessor writes this table from segmented-cell centroids.
Its x/y values must use the same coordinate system and scale as
``obsm["spatial"]``. Rows without a supplied center use the source spot
coordinate; REVISE does not invent a sub-spot coordinate for them.

The LR contribution proportions are converted with ``np.round`` and then
repaired in stable order until each spot has an exact quota equal to its virtual
cell count. For an application route, the only resolved score-matrix location
is the case-sensitive ``<data-root>/PM_on_cell.csv``; there is no CLI path override.
Its rows must cover every current virtual-cell ID and its columns must cover
every requested normalized cell-type label. Extra Patient/case rows and extra
class columns are permitted; REVISE strictly subsets and reorders the matrix to
the current case before requiring the active values to be finite numbers. A
score-maximizing assignment then places the quota slots. If that file is
missing, one seeded random permutation assigns those slots to the existing
virtual-cell rows inside each spot.

The tested invariant is exact per-spot composition, plus repeatability for the
same seed. Which virtual-cell row receives a type can change with the seed.
This random assignment changes which existing virtual-cell row receives a cell
type; it does not generate an x/y coordinate. It is not a nucleus or
cell-localization result and cannot establish which inferred type belongs to a
supplied segmentation center.

Inputs
------

The default input carrier is AnnData/H5AD. ST input provides an expression
matrix and spatial coordinates; the reference provides an expression matrix
and cell-type labels. Optional SpatialData input normalizes a selected ST table
into that same AnnData contract. Reference and benchmark ground-truth inputs
remain AnnData/H5AD.

The pipeline validates input roles, axes, required fields, coordinate shape,
and gene overlap. It aligns benchmark observations by shared identifiers before
metrics are called. Validation does not convert an arbitrary or mismatched
dataset into a scientifically comparable one.

Evidence versus interpretation
------------------------------

Automated tests establish route selection, array/label invariants, solver
events, deterministic identities, failure states, output publication, and
small synthetic execution. They do not establish biological validation,
cross-solver biological parity, or production-scale suitability. See
:doc:`limitations` before interpreting an SVC or benchmark table.
