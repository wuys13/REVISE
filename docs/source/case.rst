Reconstruction Types
====================

The installed/source interface selects one public 1.x reconstruction type
through ``application.svc_type`` in a strict application YAML. sp-SVC and
sc-SVC-sr publishes ``<output-name>.h5ad``; standard sc-SVC publishes separate
spatial and reference-expression carriers.

sp-SVC
------

.. code-block:: bash

   revise-reconstruct --config configs/application/VisiumHD.yaml

The manifest records the logical output name and ``sp-SVC`` route.

sc-SVC
------

The internal strategy returns spatial and expression carriers. The application
service validates and publishes both without merging their expression spaces:

.. code-block:: bash

   revise-reconstruct --config configs/application/Xenium_T.yaml

The outputs are
``<output-dir>/<output-name>_spatial.h5ad`` and
``<output-dir>/<output-name>_expr.h5ad``. The manifest records both logical
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

The manifest records the logical output name and ``sc-SVC-sr`` route. Omitting
curated spot-to-cell and PM-on-cell inputs uses the quota and seeded random
allocation described in :doc:`concepts`; a configured input path must exist.

Shared provenance
-----------------

Each public AnnData artifact links to the canonical run's ``provenance.json``.
The manifest records stage status and errors, not solver-event telemetry.
Inspect it rather than inferring success from directory existence.

Paper notebook compatibility
----------------------------

Historical notebooks may consume copies with older names. Standard sc-SVC
publishes the configured ``<output-name>_expr.h5ad``/``<output-name>_spatial.h5ad``
pair as its current public contract.

The canonical downstream analysis notebooks are
``Xenium_sc_SVC_T.ipynb``, ``Xenium_sc_SVC_Fibroblast.ipynb``,
``Xenium_sc_SVC_Monocyte.ipynb``, and ``VisiumHD_sp_SVC.ipynb`` under
``reproduce/case/``. They document reconstruction-to-analysis handoff but do
not replace the public reconstruction CLI or its provenance contract.
