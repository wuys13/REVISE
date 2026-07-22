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
sp-SVC requires ``Level1`` in reference ``obs``; sc-SVC and sc-SVC-sr require
both ``Level1`` and ``Level2``. The inputs must share at least one gene.

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
sc-SVC-sr`` for spot super-resolution. sc-SVC also accepts ``--select-ct``,
``--cell-type-col``, ``--sub-cell-type-col``, and ``--sc-mapping mean|random``.

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

All application routes publish:

.. code-block:: text

   <output-root>/<sample-name>/SVC.h5ad

The file links to the canonical run's ``provenance.json``. The manifest's
``result.type`` identifies the result as ``sp-SVC``, ``sc-SVC``, or
``sc-SVC-sr`` and records the route, stages, configuration, inputs, and
artifacts.

Paper reproduction notebooks
----------------------------

Curated application notebooks are tracked under ``reproduce/case/``. They may
refer to historical output layouts, but the current application command has
only the canonical ``SVC.h5ad`` output contract.

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
