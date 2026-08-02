Reconstruction Types
====================

The installed/source interface selects one public 1.x reconstruction type
with ``--svc-type``. sp-SVC and sc-SVC-sr publish ``SVC.h5ad``; standard
sc-SVC publishes separate spatial and reference-expression carriers.

sp-SVC
------

.. code-block:: bash

   revise-reconstruct \
     --svc-type sp-SVC \
     --sample-name sample \
     --data-root data \
     --st-file st.h5ad \
     --sc-ref-file sc_ref.h5ad \
     --output-root output

The manifest records ``result.type`` as ``sp-SVC``.

sc-SVC
------

The internal strategy returns spatial and expression carriers. The application
service validates and publishes both without merging their expression spaces:

.. code-block:: bash

   revise-reconstruct \
     --svc-type sc-SVC \
     --sample-name sample \
     --data-root data \
     --st-file st.h5ad \
     --sc-ref-file sc_ref.h5ad \
     --output-root output \
     --select-ct T

The outputs are ``output/sample/sc-SVC/T/sc_SVC_spatial.h5ad`` and
``output/sample/sc-SVC/T/sc_SVC_expr.h5ad``. The manifest records both logical
roles. This default route requires TACCO. Install the ``tacco`` extra from a
published package with ``python -m pip install "revise-svc[tacco]"`` or from a
source checkout with ``python -m pip install ".[tacco]"``. If a different
algorithm is acceptable, append ``--ot-method pot`` explicitly; REVISE does
not fall back automatically.

sc-SVC-sr
---------

.. code-block:: bash

   revise-reconstruct \
     --svc-type sc-SVC-sr \
     --sample-name sample \
     --data-root data \
     --st-file st.h5ad \
     --sc-ref-file sc_ref.h5ad \
     --output-root output

The manifest records ``result.type`` as ``sc-SVC-sr``. Missing curated
spot-to-cell and PM-on-cell inputs use the quota and seeded random allocation
described in :doc:`concepts`.

Shared provenance
-----------------

Each public AnnData artifact links to the canonical run's ``provenance.json``.
The manifest records requested, attempted, and completed stages. Inspect it
rather than inferring success from directory existence.

Paper notebook compatibility
----------------------------

Historical notebooks may consume copies such as ``sp_SVC.h5ad``. Standard
sc-SVC intentionally preserves the notebook-compatible
``sc_SVC_expr.h5ad``/``sc_SVC_spatial.h5ad`` pair as its current public
contract.
