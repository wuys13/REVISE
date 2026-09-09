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
     - segmented imaging-ST cells; one concrete broad cell type
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
other solver automatically.

Cluster-mode fields
-------------------

Cluster mode accepts only this ``local_refinement`` shape:

.. code-block:: yaml

   local_refinement:
     subtype_column: Level2
     select_cell_type: T
     alpha: 0.2
     resolutions: [0.6, 0.7, 0.8]

``select_cell_type`` is one concrete non-empty broad label. The CLI may
override it with ``--select-ct VALUE``. The override wins over the YAML but is
rejected for ``sp-SVC`` and SR mode, as are empty values, ``all``, and
wildcards.

SR-mode fields
--------------

SR mode accepts ``strength`` and optional graph/mass controls:

.. code-block:: yaml

   local_refinement:
     strength: 0.0
     graph:
       method: pca
       alpha: 0.2
       n_neighbors: 10
       exp_neighbors: 10
       spatial_neighbors: 10
     match_spot_sum: true

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
not include ``.h5ad`` or a path separator. For cluster mode, ``output.dir`` is
the base directory and the final directory is
``<output.dir>/<normalized selected cell-type label>``. Labels are trimmed and
``/`` becomes ``_`` (for example, ``Mono/Macro`` becomes ``Mono_Macro``). The
selected label must be safe for an output directory: empty values, traversal,
backslashes, wildcards, and control characters are rejected.

The maintained templates use ``raw_data/`` for reconstruction inputs and set
``output.dir`` under ``results/`` for reconstruction publication. The
``output`` mapping is a schema name, not a requirement that reconstruction
files live under an ``output/`` directory; case-notebook analysis artifacts
use ``output/`` separately.

``sp-SVC`` and ``sc-SVC`` SR mode each publish one H5AD. ``sc-SVC`` cluster
mode publishes a fixed pair.

.. list-table::
   :header-rows: 1
   :widths: 3 5

   * - Route
     - Published H5AD(s)
   * - ``sp-SVC``
     - ``<dir>/<name>.h5ad`` or ``<dir>/svc.h5ad``
   * - ``sc-SVC`` cluster
     - ``<dir>/<name>_spatial.h5ad`` and ``<dir>/<name>_expr.h5ad``; without a name, ``spatial.h5ad`` and ``expr.h5ad``
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
