Concepts
========

REVISE uses one orchestration lifecycle for three kinds of spatial input. The
route describes the computational contract; it does not by itself prove that a
reconstructed object is a true cell or biologically valid result.

Platforms and results
---------------------

.. list-table::
   :header-rows: 1
   :widths: 1 2 2

   * - Platform
     - Input rows
     - Public CLI result
   * - ``hST``
     - high-definition bins or pseudo-cells
     - ``hST-SVC.h5ad``
   * - ``iST``
     - segmented imaging-ST cells
     - ``iST-SVC.h5ad``
   * - ``sST``
     - spot-level observations
     - ``sST-SVC.h5ad``

The common result location is
``<output-root>/<sample-name>/<platform>-SVC.h5ad``. hST selects the sp-SVC
route. iST and sST select sc-SVC routes; the CLI merges route-internal objects
into one platform result when necessary.

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

sST virtual cells and random placement
--------------------------------------

An sST route needs a count of virtual cells per spot. If a curated
``uns["all_cells_in_spot"]`` mapping is absent, the input adapter estimates
counts from spot transcript totals and creates virtual-cell rows. Those rows
share the source spot coordinate; the fallback does not estimate sub-spot
coordinates.

The LR contribution proportions are converted with ``np.round`` and then
repaired in stable order until each spot has an exact quota equal to its virtual
cell count. For an application route, the only resolved score-matrix location
is the case-sensitive ``<data-root>/PM_on_cell.csv``; there is no CLI path override.
Its row labels must exactly match the virtual-cell IDs; columns must exactly
match the normalized cell-type labels; values must be
finite numbers. A score-maximizing assignment then places the quota slots. If
that file is missing, one seeded random permutation assigns those slots to the
existing virtual-cell rows inside each spot.

The tested invariant is exact per-spot composition, plus repeatability for the
same seed. Which virtual-cell row receives a type can change with the seed.
This random assignment is not a nucleus or cell-localization result and cannot
establish true within-spot cell identity or morphology-aware position.

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

Candidate tests establish route selection, array/label invariants, solver
events, deterministic identities, failure states, output publication, and
small synthetic execution. They do not establish biological validation,
cross-solver biological parity, or production-scale suitability. See
:doc:`limitations` before interpreting an SVC or benchmark table.
