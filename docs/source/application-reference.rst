Application Reference
=====================

This is the canonical reference for an Application YAML. It describes the
public request accepted by ``reconstruct.py`` and ``revise-reconstruct``;
Benchmark YAML uses a separate ``cf`` interface.

Route and mode
--------------

Every request has ``schema_version: 1`` and an ``application`` mapping. Its
public fields are ``application.svc_type`` and ``application.mode``:

.. code-block:: yaml

   application:
     svc_type: sc-SVC
     mode: cluster

Use exactly one of these combinations:

.. list-table::
   :header-rows: 1
   :widths: 2 2 4

   * - ``svc_type``
     - ``mode``
     - Use
   * - ``sp-SVC``
     - omitted
     - high-resolution bins or pseudo-cells
   * - ``sc-SVC``
     - ``cluster``
     - segmented imaging-ST cells; whole-sample delivery by default, or one concrete broad cell type for a compatibility call
   * - ``sc-SVC``
     - ``sr``
     - virtual cells reconstructed from multi-cell spots

``mode`` is required for ``sc-SVC`` and forbidden for ``sp-SVC``.

Paths and input formats
-----------------------

``paths.root_dir`` is either literal ``.`` (the command launch directory) or
an existing absolute directory. All ``inputs`` and ``output`` paths are
non-empty relative children of that root; absolute paths and ``..`` traversal
are rejected.

``inputs.st`` requires ``path`` and ``format``. ``format`` is ``h5ad``,
``spatialdata``, or ``auto``. A SpatialData/auto request can set
``inputs.st.spatialdata.table`` and ``inputs.st.spatialdata.element``. The
reference requires ``inputs.reference.path`` and ``format: h5ad``.

Both AnnData inputs must have non-empty expression, unique observation and
variable names, and at least one shared gene. The ST input also needs finite
two-dimensional ``obsm["spatial"]`` coordinates. The reference must contain
the configured annotation columns after any requested filter.

Reference filter and annotations
--------------------------------

Set ``inputs.reference.filter_column`` and ``filter_value`` together, or omit
both. Filtering happens before preprocessing; REVISE does not infer a patient,
cohort, or sample filter.

``global_anchoring.broad_column`` names the reference ``obs`` column used for
broad cell types. Cluster mode additionally requires
``local_refinement.subtype_column``. SR mode uses the broad assignment and has
no subtype-column field.

Preprocessing and OT
--------------------

``preprocessing.spatial`` requires ``min_transcript_counts`` and
``min_cell_counts``. ``preprocessing.spatial.min_counts`` is optional:

.. code-block:: yaml

   min_transcript_counts: null
   min_counts: 20
   min_cell_counts: 30

``preprocessing.reference`` requires ``min_transcript_counts`` and
``min_cell_counts``. ``preprocessing.reference.min_genes`` is optional. The
optional ``min_counts`` and ``min_genes`` fields accept non-negative integers;
``null`` disables either optional threshold. Required count thresholds also
accept ``null`` only where their route permits it and otherwise are
non-negative integers.

``algorithm.ot_method`` is optional and, when supplied, is ``pot`` or
``tacco``. It selects both Global Anchoring and Local Refinement solvers. A
missing or failing solver stops the request; REVISE does not substitute the
other solver automatically. sc-SVC SR defaults to ``tacco``; its Visium
template records that choice explicitly.

Cluster-mode fields
-------------------

Cluster mode accepts only this ``local_refinement`` shape:

.. code-block:: yaml

   local_refinement:
     subtype_column: Level2
     alpha: 0.2
     resolutions: [0.6, 0.7, 0.8]

When ``select_cell_type`` is omitted, cluster mode runs whole-sample delivery:
Global Anchoring runs once, then actual broad labels are considered for Local
Refinement. A type is eligible only when its reference subset has more than one
valid, non-empty subtype label. Ineligible types are recorded as skipped; a
failure in an eligible type fails the sample publication. Whole-sample cluster
mode accepts ``output.ist_mapping: random | mean`` and publishes one assembled
``SVC.h5ad`` together with ``raw.h5ad`` and ``sample.yaml``.

``select_cell_type`` is one concrete non-empty broad label when a compatibility
call needs a single type. The CLI may override it with ``--select-ct VALUE``.
The override wins over the YAML but is rejected for ``sp-SVC`` and SR mode, as
are empty values, ``all``, and wildcards. A selected single type keeps the
legacy single-type scope and is not a whole-sample delivery. If
``output.ist_mapping`` is omitted, its legacy result is paired; an explicit
``mean`` or ``random`` requests one assembled ``SVC.h5ad`` for that single
type. Its existing compatibility eligibility is preserved: one valid subtype
is allowed, while no valid subtype still fails.

Whole-sample delivery can optionally declare the sample identity and coordinate
metadata:

.. code-block:: yaml

   delivery:
     sample_id: CRC/P2CRC_Xenium
     coordinates:
       key: spatial
       unit: pixel
       microns_per_coordinate: 0.2125

``delivery.sample_id`` is a normalized relative hierarchical ID; it cannot be
absolute, contain ``..`` or backslashes, or be rewritten by taking only the
last path component. ``delivery.coordinates.unit`` accepts ``um``,
``micron``, ``pixel`` or ``unknown``. ``um`` is recorded as ``micron``;
``micron`` requires a scale of 1, and an unknown unit cannot declare a physical
scale. If ``delivery`` is omitted, the output directory name supplies the
sample identity and coordinate metadata remains unknown unless supplied by the
route.

The whole-sample delivery snapshots the original spatial input before
preprocessing and verifies that the source does not change before publication.
``raw.h5ad`` keeps the original expression matrix, observation/variable axes,
coordinates, and original annotations, including units excluded by QC. The
delivery adds ``revise_Level1`` and ``revise_Level2`` only for inference
obtained for real observation IDs; missing inference remains missing. Original
annotation columns are never overwritten. If an input already uses one of the
reserved inference column names, publication fails with an explicit error.
The generated ``sample.yaml`` points its broad and subtype declarations to
these inference columns and records the original sample ID and source
identity. It does not claim that a normalized reconstruction matrix is raw
counts.

SR-mode fields
--------------

SR mode accepts ``strength`` and optional graph controls:

.. code-block:: yaml

   local_refinement:
     strength: 0.0
     graph:
       method: pca
       alpha: 0.2
       n_neighbors: 10
       exp_neighbors: 10
       spatial_neighbors: 10

SR always applies parent-spot per-gene correction against internally normalized
spot expression. It does not finally scale each generated cell to 10,000 or
restore raw counts. Internal spot/reference normalization and graph-copy
normalization/log1p remain unchanged. The removed
``local_refinement.match_spot_sum`` and engine ``sc.match_spot_sum`` keys raise
a migration error for either boolean value; delete the key from old configs.

SR mode also accepts one cell-count method under ``algorithm``:

.. code-block:: yaml

   algorithm:
     ot_method: tacco
     sr_cell_count_method: cyto_linear_v1

``cyto_linear_v1`` is the default. Before spatial or gene filtering, it uses
the full-gene spot matrix ``X`` to compute
``R_s = sum_g(log2(1 + 1e6 * X_sg / sum_g(X_sg)))`` and then
``N_s = max(1, rint(-13.2645109457 + 0.000681111936 * R_s))``. Sparse input
stays sparse and the estimate has no upper cap. Set
``sr_cell_count_method: transcript_heuristic`` only to reproduce the previous
median-four, 1--12-cell behavior. If ``uns["all_cells_in_spot"]`` is already
present, that supplied mapping remains authoritative and neither estimator is
run.

``strength: 0.0`` does not skip Local Refinement; it turns off posterior
blending in its cost. SR inputs can optionally declare one exact
``inputs.pm_on_cell.path``. When present, its rows must equal the active
virtual-cell IDs, its columns must equal normalized broad labels, and values
must be finite in ``[0, 1]``. REVISE may reorder matching axes but never clips
or normalizes the scores. When it is absent, no sidecar is searched and the
seeded quota allocation is used.

sp-SVC fields
-------------

``sp-SVC`` accepts only optional ``local_refinement.strength``. It has no
mode, subtype, selected-cell-type, graph, or PM field.

Output and publication
----------------------

``output.dir`` is required. ``output.name`` is an optional filename stem; do
not include ``.h5ad`` or a path separator. For a selected cluster type,
``output.dir`` is the base directory and the final directory is
``<output.dir>/<normalized selected cell-type label>``. For whole-sample
cluster mode, the output directory itself is the sample delivery directory.
Labels are trimmed and ``/`` becomes ``_`` (for example, ``Mono/Macro`` becomes
``Mono_Macro``). A selected label must be safe for an output directory: empty
values, traversal, backslashes, wildcards, and control characters are
rejected.

The maintained templates use ``raw_data/`` for reconstruction inputs and set
``output.dir`` under ``results/`` for reconstruction publication. The
``output`` mapping is a schema name, not a requirement that reconstruction
files live under an ``output/`` directory; case-notebook analysis artifacts
use ``output/`` separately.

``sp-SVC`` and ``sc-SVC`` SR mode each publish one H5AD. A selected single-type
``sc-SVC`` cluster call publishes a fixed pair by default
(``output.ist_mapping: paired``). Whole-sample cluster mode defaults to
``output.ist_mapping: random`` and publishes one assembled ``SVC.h5ad`` (or
``<name>.h5ad``), plus ``raw.h5ad`` and ``sample.yaml``. It also accepts
``mean``. ``mean`` assigns cluster means, while ``random`` samples
within-cluster donors using the execution seed and records donor IDs. GA/LR are
unchanged by the final assembly choice. Other routes reject this key. See
`batch protocol <https://github.com/wuys13/REVISE/blob/revise-2.0/docs/design/batch/protocol.md#output-roles-and-current-result-semantics>`_ for carrier semantics
and batch directory conventions. The table below distinguishes compatibility
single-type output from whole-sample delivery.

.. list-table::
   :header-rows: 1
   :widths: 3 5

   * - Route
     - Published H5AD(s)
   * - ``sp-SVC``
     - ``<dir>/<name>.h5ad`` or ``<dir>/svc.h5ad``
   * - ``sc-SVC`` cluster
     - ``<dir>/<name>_spatial.h5ad`` and ``<dir>/<name>_expr.h5ad``; without a name, ``spatial.h5ad`` and ``expr.h5ad``
   * - ``sc-SVC`` cluster, whole-sample
     - ``<dir>/<name>.h5ad`` or ``<dir>/SVC.h5ad``, plus ``raw.h5ad`` and ``sample.yaml``
   * - ``sc-SVC`` sr
     - ``<dir>/<name>.h5ad`` or ``<dir>/svc.h5ad``

The effective request and published H5AD metadata record the normalized
``svc_type``, route, mode, and selected cell type. ``provenance.json`` records
the resolved request, input identities, output roles, stages, and terminal
state. Treat a succeeded manifest plus the route's promised artifacts as the
success contract.

Execution
---------

``execution.seed`` is optional and defaults to 42. It must be an integer from
0 through ``2**32 - 1``. The effective seed is recorded in the request and
provenance metadata.


Prepared reference override
---------------------------

Reference preparation runs independently::

   python -m revise.reference_preparation --config preparation.yaml
   python reconstruct.py --config application.yaml --reference-config prepared/reference.yaml

The Python equivalent is ``run_application(config_path, reference_config=path)``.
The optional reference file replaces only reference path/format and clears old
reference filters before validation. Other preprocessing, solver, annotation,
output and seed settings remain controlled by the Application YAML. With no
external reference config, existing behavior is unchanged. A failed override
never falls back to the old reference.

The preparation report records candidate rankings and input identities; consumer
checks reject a selected file that changed after preparation. Screening and
reconstruction each run GA; this version does not cache the screening matrix.
See ``docs/development/reconstruction-analysis/plans/reference-preparation.md``
for paired/screen examples, batch usage and scientific interpretation limits.

OT assembly and linear expression declarations
---------------------------------------------

``output.ist_mapping`` additionally accepts ``within_cluster`` and
``outside_cluster``. Both use TACCO for the new assembly step; upstream GA/LR
continues to obey ``algorithm.ot_method``. The default remains ``random``.
``output.ist_ot`` accepts ``spatial_weight`` (default 0.2),
``max_cost_entries`` (default 2000000), and ``gene_block_size`` (default 256).
Spatial mixing uses same-SVC_cluster neighbours only. Within mode matches
individual donors in that cluster; outside mode matches cluster mean profiles
from the same broad type. The full valid row-normalized weights reconstruct
reference-wide expression. Exceeding the cost-size limit fails without fallback.

Optional ``inputs.st.expression`` and ``inputs.reference.expression`` mappings
contain ``identity`` and ``scale: untransformed_nonnegative``. They explicitly
attest to non-log linear input history; finite/nonnegative checks do not infer
this history. Missing declarations stay unknown. Alternatively, the source H5AD
may explicitly carry these fields in ``uns['revise_expression']``. A reference
configuration override discards the old reference declaration. The actual new
file's explicit source metadata may supply its declaration.

See ``docs/development/reconstruction-analysis/plans/ot-assembly.md`` for the
configuration, confidence boundaries and comparison Notebook.
