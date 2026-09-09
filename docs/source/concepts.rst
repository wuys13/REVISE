Concepts
========

REVISE reconstructs Spatially-inferred Virtual Cells (SVCs) with one shared
reconstruction lifecycle. Application users choose from the shape of their
spatial observations; Sim2Real-ST Benchmark users choose a confounding-factor
family. The frontends prepare different data, then meet at the same engine and
optimal-transport implementation.

Choose from the data shape
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 2 3 3 3

   * - Public type and mode
     - Spatial rows
     - Typical platform
     - Public result
   * - ``sp-SVC``
     - high-definition bins or pseudo-cells
     - Visium HD
     - one spatially refined ``AnnData``
   * - ``sc-SVC``, cluster mode
     - segmented cells
     - Xenium, CosMx, MERFISH
     - spatial and expression ``AnnData`` objects for one selected broad type
   * - ``sc-SVC``, sr mode
     - multi-cell spots
     - Visium
     - one virtual-cell reconstruction ``AnnData``

This is a data-shape decision, not an automatic platform detector. The user is
responsible for supplying the matching input and annotation columns.

sc-SVC modes
------------

Cluster mode starts from segmented-cell observations. Its selected broad cell
type identifies the cohort that receives subtype and expression refinement;
the two output carriers keep spatial-side and expression-side analysis
separate.

SR mode starts from a multi-cell spot. It constructs the virtual-cell rows
needed for each spot and assigns the spot-level broad composition to those
rows. A supplied PM-on-cell score matrix can guide that assignment; otherwise
the seeded quota allocation is used. SR mode reconstructs cell-type composition
and expression within a spot. It does **not** by itself prove physical
sub-spot cell locations or a nucleus/localization result. Segmentation-derived
centers, when supplied, are retained; missing centers remain at the source
spot coordinate.

For exact input axes, PM semantics, field constraints, and output paths, see
:doc:`application-reference`.

Unified lifecycle and OT
------------------------

Application and Benchmark requests converge on this fixed stage order:

1. validate inputs;
2. Global Anchoring;
3. Local Refinement, including route-owned local-unit, graph, and OT work;
4. finalize the SVC;
5. evaluate only for an enabled Benchmark request.

``algorithm.ot_method`` selects both the Global Anchoring and Local Refinement
solver in an Application YAML. POT and TACCO share the same lower-level OT
surface. A missing or failed selected solver is an error; REVISE never switches
to the other solver automatically.

Evidence boundary
-----------------

Tests establish routing, input and axis contracts, deterministic identities,
failure states, output publication, and small synthetic execution. A notebook
snapshot or a passing software test does not establish biological validation,
clinical validity, cross-solver biological parity, or production-scale
suitability.
