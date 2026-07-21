REVISE Documentation
====================

REVISE reconstructs Spatially-inferred Virtual Cells (SVCs) from spatial
transcriptomics data and a matched single-cell reference. Version |release| is
a release candidate whose public reconstruction paths use one installed
command, one pipeline lifecycle, and one platform result file.

Candidate status
----------------

The candidate matrix targets Python 3.10/3.11 packaging, the installed command,
small synthetic POT runs, small TACCO solver smokes, failure provenance, and
strict documentation checks. Real-data end-to-end scientific validation is
deferred. The documentation therefore makes no biological-validation or
production-scale claim.

External data identities, publication identifiers, project roles, and release
contacts are pending owner confirmation. Historical notebooks and research
assets remain in ``REVISE-legacy`` and are located through the exact source
commit recorded in ``legacy-assets.json``.

Start here
----------

- :doc:`source/installation` distinguishes candidate Wheel/source installs from
  optional scientific domains.
- :doc:`source/quickstart` runs hST, iST, or sST through
  ``revise-reconstruct`` or the source compatibility wrapper.
- :doc:`source/configuration` explains GA/LR solver selection and why TACCO
  failures never fall back to POT.
- :doc:`source/benchmark` records the exact metric formulas and their proof
  limits.
- :doc:`source/limitations` states the candidate's data, scale, solver, and
  scientific evidence boundaries.

Public result contract
----------------------

The installed/source CLI publishes one file:

.. code-block:: text

   <output-root>/<sample-name>/<platform>-SVC.h5ad

The concrete names are ``hST-SVC.h5ad``, ``iST-SVC.h5ad``, and
``sST-SVC.h5ad``. Each result points to a canonical run directory containing
``provenance.json`` and route-specific evidence. See :doc:`source/architecture`
for the lifecycle and :doc:`source/case` for route selection.

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
   :caption: Historical Assets
   :maxdepth: 2
   :hidden:

   Legacy research index <source/gallery>

.. toctree::
   :caption: Reference
   :maxdepth: 2
   :hidden:

   source/architecture
   source/api/index
