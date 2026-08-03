Quick Start
===========

REVISE has two public entry points: Benchmark and Application. Paper notebooks
require a repository checkout and their corresponding external data.

Benchmark
---------

Run one Sim2Real-ST confounding family from the repository root:

.. code-block:: bash

   python reproduce/benchmark_main.py \
     --confounding segmentation \
     --data-root raw_data/Sim2Real-ST \
     --sample-name P2CRC/cut_part1 \
     --dataset-task segmentation \
     --output-root output/benchmark

This command runs one confounding family, which may contain multiple leaf
cases.

``--confounding`` accepts ``segmentation``, ``bin2cell``, ``batch_effect``,
``spot_size``, ``gene_panel``, or ``gene_dropout``. Run the bounded multi-family
launcher with:

.. code-block:: bash

   bash reproduce/benchmark_main.sh

Benchmark analysis notebooks are tracked under ``reproduce/benchmark/``.

Application inputs
------------------

For ``--sample-name sample``, ``--st-file st.h5ad``, and
``--sc-ref-file sc_ref.h5ad``, use the flat input layout resolved by the CLI:

.. code-block:: text

   data/
   |-- sample_st.h5ad
   `-- sc_ref.h5ad

The resolved paths are ``data/sample_st.h5ad`` and ``data/sc_ref.h5ad``. Both
inputs require non-empty ``X``, unique ``obs_names`` and unique ``var_names``.
The ST input requires two spatial coordinate columns in ``obsm["spatial"]``.
Every route requires the configured broad annotation in reference ``obs``
(``Level1`` by default). Only standard sc-SVC requires the configured subtype
annotation (``Level2`` by default). sc-SVC-sr composition and expression
allocation use the broad assignment and do not require a subtype column. The
inputs must share at least one gene.

If the reference has the default ``Patient`` column, its values are matched to
``--sample-name``. Select another column with ``--patient-key``.

Application command
-------------------

.. code-block:: bash

   python reconstruct.py \
     --svc-type sp-SVC \
     --sample-name sample \
     --data-root data \
     --st-file st.h5ad \
     --sc-ref-file sc_ref.h5ad \
     --output-root output

Use ``--svc-type sc-SVC`` for molecular completion and ``--svc-type
sc-SVC-sr`` for spot super-resolution. ``--cell-type-col`` selects the broad
reference annotation on every route. ``--sub-cell-type-col`` selects the
refined annotation required only by standard sc-SVC. sc-SVC-sr uses the broad
assignment for composition and expression allocation. sc-SVC also accepts
``--select-ct``. Its application profile defaults to TACCO; install it with
``python -m pip install ".[tacco]"`` from source or ``python -m pip install
"revise-svc[tacco]"`` from a published package. If TACCO is unavailable and a
different algorithm is acceptable, explicitly add ``--ot-method pot``. REVISE
never switches algorithms automatically.

For sc-SVC-sr, optional segmentation-derived centers use a DataFrame in
``st_adata.uns["revise_cell_locations"]`` with a unique ``cell_id`` index and
``spot_name/x/y`` columns. Its assignments agree with
``uns["all_cells_in_spot"]`` and its coordinates use the same coordinate system
and scale as ``obsm["spatial"]``; rows without centers remain at the spot center.
The optional sample-local probability prior is resolved from the prepared ST
path as ``<st-parent>/<st-stem>_PM_on_cell.csv``. Its rows must exactly equal
the active virtual-cell IDs and its columns must exactly equal the active
normalized cell-type labels. Values must be numeric and finite within
``[0, 1]``; every row must sum to one with zero relative tolerance and an
absolute tolerance of ``1e-6``. REVISE only reorders exact axes and never clips
or normalizes PM. It is not a case table, cohort registry, or generic assignment
posterior. Without that file, those coordinates are retained while inferred
cell types are assigned to the existing rows by a seeded random permutation.

After installation, the equivalent package command is:

.. code-block:: bash

   revise-reconstruct \
     --svc-type sp-SVC \
     --sample-name sample \
     --data-root data \
     --st-file st.h5ad \
     --sc-ref-file sc_ref.h5ad \
     --output-root output

Append ``--dry-run`` to validate the resolved route, inputs, and dependencies
without running reconstruction.

Application output
------------------

``sp-SVC`` and ``sc-SVC-sr`` publish:

.. code-block:: text

   <output-root>/<sample-name>/SVC.h5ad

Standard ``sc-SVC`` publishes:

.. code-block:: text

   <output-root>/<sample-name>/sc-SVC/<cell-type>/sc_SVC_spatial.h5ad
   <output-root>/<sample-name>/sc-SVC/<cell-type>/sc_SVC_expr.h5ad

Each file links to the canonical run's ``provenance.json``. The manifest
records the result role, route, stages, configuration, per-role input
identities, software identity, and artifacts.

Paper reproduction notebooks
----------------------------

Curated application notebooks are tracked under ``reproduce/case/``. Standard
sc-SVC now preserves the notebook spatial and reference-expression carriers as
separate public outputs.

Application utilities
---------------------

Build optional morphology-derived priors with the installed package command:

.. code-block:: bash

   revise-build-histology-priors \
     --st-h5ad st.h5ad \
     --mask segmented_cells.tif \
     --out-h5ad st_with_histology_priors.h5ad

Compute biology-facing post-reconstruction metrics through the package-owned
analysis layer:

.. code-block:: bash

   revise-compute-biological-metrics \
     --input-h5ad output/sample/SVC.h5ad \
     --output-dir output/sample/biological_metrics
