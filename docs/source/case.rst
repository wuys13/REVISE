Reconstruction Types
====================

REVISE 2.0 exposes exactly three public selectors through ``--svc-type``.
Every route publishes one route-qualified public H5AD:
``<output-root>/<sample-name>/<svc-type>/SVC.h5ad``.

hST-SVC
-------

.. code-block:: bash

   revise-reconstruct \
     --svc-type hST-SVC \
     --sample-name sample \
     --data-root data \
     --st-file st.h5ad \
     --sc-ref-file sc_ref.h5ad \
     --output-root output

The public result is ``output/sample/hST-SVC/SVC.h5ad`` and the manifest
records ``result.type`` as ``hST-SVC``.

iST-SVC
-------

.. code-block:: bash

   revise-reconstruct \
     --svc-type iST-SVC \
     --sample-name sample \
     --data-root data \
     --st-file st.h5ad \
     --sc-ref-file sc_ref.h5ad \
     --output-root output \
     --ist-mapping mean

``mean`` is the default. It preserves the spatial carrier's observations and
coordinates, preserves the expression carrier's genes, and assigns each
spatial row its cluster's mean expression. ``--ist-mapping random`` instead
selects a donor row from the same cluster with the effective seed and records
the donor IDs and hash. Both modes publish only
``output/sample/iST-SVC/SVC.h5ad``. This default route requires TACCO; pass
``--ot-method pot`` only when deliberately selecting a different algorithm.

sST-SVC
-------

.. code-block:: bash

   revise-reconstruct \
     --svc-type sST-SVC \
     --sample-name sample \
     --data-root data \
     --st-file st.h5ad \
     --sc-ref-file sc_ref.h5ad \
     --output-root output

The public result is ``output/sample/sST-SVC/SVC.h5ad``. Missing curated
spot-to-cell and PM-on-cell inputs use the quota and seeded random allocation
described in :doc:`concepts`.

Shared provenance
-----------------

The canonical run's ``provenance.json`` records the public AnnData. Its
``result`` contains exactly ``filename`` and ``type``. Only iST-SVC adds the
top-level ``assembly`` evidence. Inspect the manifest rather than inferring
success from directory existence.

Historical 1.x notebook compatibility
--------------------------------------

The files under ``reproduce/case/`` are 1.x historical reproduction material,
not current 2.0 output. Their carrier names, including ``sp_SVC.h5ad``,
``sc_SVC_expr.h5ad``, and ``sc_SVC_spatial.h5ad``, remain unchanged so the
historical notebooks can still describe their original workflows. New 2.0
runs must not treat those names or paired carriers as public output contracts.
