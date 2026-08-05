REVISE Documentation
====================

REVISE reconstructs Spatially-inferred Virtual Cells (SVCs) from spatial
transcriptomics data and a matched single-cell reference. It provides two
public workflows: Sim2Real-ST benchmark reproduction and SVC reconstruction
for real-data applications.

The paper data and reproduced results are published at
``https://zenodo.org/records/17705737``. Curated notebooks are tracked under
``reproduce/``; installation does not download research data.

Start here
----------

- :doc:`source/installation` describes the base package and optional analysis
  capabilities.
- :doc:`source/quickstart` shows the Benchmark and Application entry points.
- :doc:`source/concepts` defines sp-SVC, sc-SVC, inputs, and evidence limits.
- :doc:`source/configuration` documents detailed runtime and OT configuration.
- :doc:`source/benchmark` records the implemented metric formulas and their
  proof limits.
- :doc:`source/limitations` states the current data, scale, and scientific
  evidence boundaries.

Public result contract
----------------------

sp-SVC and sc-SVC-sr publish:

.. code-block:: text

   <output-root>/<sample-name>/SVC.h5ad

Standard sc-SVC publishes its paired carriers:

.. code-block:: text

   <output-root>/<sample-name>/sc-SVC/<cell-type>/sc_SVC_spatial.h5ad
   <output-root>/<sample-name>/sc-SVC/<cell-type>/sc_SVC_expr.h5ad

Each result links to the run's ``provenance.json``. Single-output runs use
``result.type``; sc-SVC records both roles under ``results``. The manifest also
records the application request, internal route, engine configuration, stages,
inputs, and artifacts.

.. toctree::
   :caption: Start Here
   :maxdepth: 2
   :hidden:

   source/concepts
   source/quickstart
   source/installation
   source/limitations

.. toctree::
   :caption: Run REVISE
   :maxdepth: 2
   :hidden:

   Reconstruction routes <source/case>
   Benchmark reproduction <source/benchmark>
   source/configuration

.. toctree::
   :caption: Reproduction
   :maxdepth: 2
   :hidden:

   Curated notebooks <source/gallery>

.. toctree::
   :caption: Reference
   :maxdepth: 2
   :hidden:

   source/architecture
   source/api/index
