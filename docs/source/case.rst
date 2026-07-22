Reconstruction Platforms
========================

The installed/source interface selects one public 2.0 input class with
``--platform``. Every route publishes the stable filename ``SVC.h5ad``.

hST
---

.. code-block:: bash

   revise-reconstruct \
     --platform hST \
     --sample-name sample \
     --data-root data \
     --st-file st.h5ad \
     --sc-ref-file sc_ref.h5ad \
     --output-root output

The manifest records ``result.type`` as ``hST-SVC``.

iST
---

The internal strategy returns spatial and expression carriers. The application
service validates their cluster sets and publishes one result:

.. code-block:: bash

   revise-reconstruct \
     --platform iST \
     --sample-name sample \
     --data-root data \
     --st-file st.h5ad \
     --sc-ref-file sc_ref.h5ad \
     --output-root output \
     --select-ct T \
     --sc-mapping mean

The manifest records ``result.type`` as ``iST-SVC``. ``--sc-mapping random``
uses the command seed to choose a same-cluster reference expression row instead
of the cluster mean.

sST
---

.. code-block:: bash

   revise-reconstruct \
     --platform sST \
     --sample-name sample \
     --data-root data \
     --st-file st.h5ad \
     --sc-ref-file sc_ref.h5ad \
     --output-root output

The manifest records ``result.type`` as ``sST-SVC``. Missing curated
spot-to-cell and PM-on-cell inputs use the quota and seeded random allocation
described in :doc:`concepts`.

Shared provenance
-----------------

The public AnnData links to the canonical run's ``provenance.json``. The
manifest records requested, attempted, and completed stages. Inspect it rather
than inferring success from directory existence.

Paper notebook compatibility
----------------------------

Historical notebooks may consume copies such as ``sp_SVC.h5ad``,
``sc_SVC_expr.h5ad``, and ``sc_SVC_spatial.h5ad``. These are reproduction
material and do not change the current ``SVC.h5ad`` contract.
