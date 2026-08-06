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

Application output contract
---------------------------

Every Application invocation uses exactly one run directory:

.. code-block:: text

   <output.dir>/<output.name>/application__<svc-type>/<timestamp_uuid>/

The run envelope contains ``provenance.json``, ``merged_config.json``, and
``run.log``; completed input validation adds ``preflight.json``. Successful
formal sp-SVC and sc-SVC-sr runs add
``<output.name>.h5ad``; standard sc-SVC adds
``<output.name>_spatial.h5ad`` and ``<output.name>_expr.h5ad`` as siblings in
the same directory. Dry runs add no H5AD. Repeating a run creates a new
timestamp/UUID leaf and never overwrites an earlier run.

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
