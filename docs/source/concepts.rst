Concepts
========

REVISE reconstructs Spatially-inferred Virtual Cells (SVCs) with one shared
reconstruction lifecycle. Application users choose a reconstruction type from
the rows represented by their spatial dataset. Benchmark users choose a
confounding-factor family. Those frontends prepare data differently, but both
reach the same pipeline and shared OT implementation.

Choose from the data shape
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 2 3 3 3

   * - Public type
     - Spatial observations
     - Typical platform or task
     - Public result
   * - ``sp-SVC``
     - high-definition bins or pseudo-cells
     - Visium HD; segmentation or bin-to-cell artifacts
     - one spatially refined ``AnnData``
   * - ``sc-SVC``
     - segmented imaging-ST cells
     - Xenium, CosMx, or MERFISH; one selected broad cell type
     - spatial and expression ``AnnData`` objects
   * - ``sc-SVC-sr``
     - multi-cell spots
     - Visium; spot super-resolution
     - one reconstructed ``AnnData``

This is a data-shape decision, not an automatic platform detector. Users are
responsible for supplying data that match the selected type and for choosing
the correct annotations and mode.

Inputs and outputs
------------------

The normal carrier is AnnData/H5AD. The ST object contains expression and
finite two-dimensional coordinates in ``obsm["spatial"]``. The matched
single-cell reference contains expression and the broad annotation named by
``global_anchoring.broad_column``. Standard ``sc-SVC`` also requires its
configured subtype annotation and one concrete selected broad type. Inputs
must share at least one gene.

``sp-SVC`` and ``sc-SVC-sr`` publish ``<output-dir>/<output-name>.h5ad`` and
return one ``AnnData``. Standard ``sc-SVC`` publishes
``<output-dir>/<output-name>_spatial.h5ad`` and
``<output-dir>/<output-name>_expr.h5ad`` and returns the fixed
``(spatial, expression)`` pair. Each public H5AD links to the run's
``provenance.json``.

Unified lifecycle and OT
------------------------

Application and Benchmark requests converge on this fixed stage order:

1. validate inputs;
2. Global Anchoring;
3. Local Refinement, including route-owned local-unit, graph, and OT work;
4. finalize the SVC;
5. evaluate only for an enabled Benchmark request.

POT and TACCO share the same bottom-level OT surface. In an Application YAML,
``algorithm.ot_method`` sets both the GA solver and the LR solver. REVISE does
not switch solvers automatically when a selected implementation is missing or
fails.

Spot super-resolution assignment
--------------------------------

``sc-SVC-sr`` needs a virtual-cell count for every spot. A prepared object may
provide ``uns["all_cells_in_spot"]``. Otherwise the input service estimates
counts from spot transcript totals and creates virtual-cell rows.

Optional segmentation-derived centers use the standardized
``uns["revise_cell_locations"]`` table. Its index contains unique cell IDs;
the columns are ``spot_name``, ``x``, and ``y``; and its assignments must agree
with ``uns["all_cells_in_spot"]``. Coordinates must already use the same scale
and coordinate system as ``obsm["spatial"]``. Rows without a center retain the
source spot coordinate; REVISE does not invent a sub-spot location.

The assignment converts spot-level Global Anchoring posterior proportions with
``np.round`` and repairs them in stable order to an exact quota. Application
can declare one exact
``inputs.pm_on_cell.path``. This sample-local score matrix must have rows
exactly equal to the active virtual-cell IDs, columns exactly equal to the
active normalized broad labels, and numeric finite values in ``[0, 1]``.
REVISE only reorders matching
axes; it never clips or normalizes the values. PM is not a case table, cohort
registry, or generic assignment posterior.

If that field is omitted, no sidecar is probed. One seeded random permutation
assigns the exact quota to existing virtual-cell rows within each spot. The
tested contract is exact composition and same-seed repeatability. It is not a
nucleus or cell-localization result.

Evidence boundary
-----------------

Tests establish routing, input and axis contracts, deterministic identities,
failure states, output publication, and small synthetic execution. A rendered
notebook or passing software test does not establish biological validation,
clinical validity, cross-solver biological parity, or production-scale
suitability.
