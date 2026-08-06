Reconstruction Types
====================

The installed/source interface selects one public 1.x reconstruction type
through ``application.svc_type`` in a strict application YAML. Every
Application invocation writes to one run directory:
``<output.dir>/<output.name>/application__<svc-type>/<timestamp_uuid>/``.
sp-SVC and sc-SVC-sr write ``<output.name>.h5ad`` there; standard sc-SVC writes
its spatial and reference-expression carriers as sibling files in that same
directory.

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

The outputs are ``<output.name>_spatial.h5ad`` and
``<output.name>_expr.h5ad`` in the run directory. The manifest records both
logical roles. This default route requires TACCO. Install the ``tacco`` extra from a
published package with ``python -m pip install "revise-svc[tacco]"`` or from a
source checkout with ``python -m pip install ".[tacco]"``. Choose the Xenium
template whose fixed ``local_refinement.select_cell_type`` matches the broad
label in the reference. If a different algorithm is acceptable, set
``algorithm.ot_method: pot`` explicitly; REVISE does not fall back automatically.

sc-SVC-sr
---------

.. code-block:: bash

   revise-reconstruct --config configs/application/Visium.yaml

The manifest records the logical output name and ``sc-SVC-sr`` route. Missing curated
spot-to-cell and PM-on-cell inputs use the quota and seeded random allocation
described in :doc:`concepts`.

Shared provenance
-----------------

Each run envelope contains ``provenance.json``, ``merged_config.json``, and
``run.log``; completed input validation adds ``preflight.json``. A dry run gets
its own run directory and contains no H5AD; repeating a formal run creates a
new ``<timestamp_uuid>`` leaf. The manifest records stage status and errors,
not solver-event telemetry. Inspect it rather than inferring success from
directory existence.

Paper notebook compatibility
----------------------------

Historical notebooks may consume copies with older names. The current standard
sc-SVC contract is the configured ``<output.name>_expr.h5ad``/
``<output.name>_spatial.h5ad`` pair inside the unique Application run directory;
there is no shared fixed-name publication copy.
