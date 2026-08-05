Reconstruction Types
====================

The installed/source interface selects one public 1.x reconstruction type
through ``application.svc_type`` in a strict application YAML. sp-SVC and
sc-SVC-sr publish ``SVC.h5ad``; standard sc-SVC publishes separate spatial and
reference-expression carriers.

sp-SVC
------

.. code-block:: bash

   revise-reconstruct --config configs/application/VisiumHD.yaml

The manifest records ``result.type`` as ``sp-SVC``.

sc-SVC
------

The internal strategy returns spatial and expression carriers. The application
service validates and publishes both without merging their expression spaces:

.. code-block:: bash

   revise-reconstruct --config configs/application/Xenium_T.yaml

The outputs are
``<output-root>/<sample-name>/sc-SVC/<cell-type>/sc_SVC_spatial.h5ad`` and
``<output-root>/<sample-name>/sc-SVC/<cell-type>/sc_SVC_expr.h5ad``. The manifest records both logical
roles. This default route requires TACCO. Install the ``tacco`` extra from a
published package with ``python -m pip install "revise-svc[tacco]"`` or from a
source checkout with ``python -m pip install ".[tacco]"``. Choose the Xenium
template whose fixed ``local_refinement.select_cell_type`` matches the broad
label in the reference. If a different algorithm is acceptable, set
``algorithm.ot_method: pot`` explicitly; REVISE does not fall back automatically.

sc-SVC-sr
---------

.. code-block:: bash

   revise-reconstruct --config configs/application/Visium.yaml

The manifest records ``result.type`` as ``sc-SVC-sr``. Missing curated
spot-to-cell and PM-on-cell inputs use the quota and seeded random allocation
described in :doc:`concepts`.

Shared provenance
-----------------

Each public AnnData artifact links to the canonical run's ``provenance.json``.
The manifest records stage status and errors, not solver-event telemetry.
Inspect it rather than inferring success from directory existence.

Paper notebook compatibility
----------------------------

Historical notebooks may consume copies such as ``sp_SVC.h5ad``. Standard
sc-SVC intentionally preserves the notebook-compatible
``sc_SVC_expr.h5ad``/``sc_SVC_spatial.h5ad`` pair as its current public
contract.
