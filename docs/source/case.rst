Reconstruction Routes
=====================

The installed/source interface selects an internal route with ``--platform``.
The public 1.x result remains an sp-SVC or sc-SVC and always uses the stable
filename ``SVC.h5ad``.

hST route
---------

hST selects the sp-SVC application profile:

.. code-block:: bash

   revise-reconstruct \
     --platform hST \
     --sample-name sample \
     --data-root data \
     --st-file st.h5ad \
     --sc-ref-file sc_ref.h5ad \
     --output-root output

The public result is ``output/sample/SVC.h5ad`` and its manifest records
``result.type`` as ``sp-SVC``.

iST route
---------

iST selects the sc-SVC imaging route. The internal strategy returns spatial
and expression carriers; the CLI validates their cluster sets and publishes
one result:

.. code-block:: bash

   revise-reconstruct \
     --platform iST \
     --sample-name sample \
     --data-root data \
     --st-file st.h5ad \
     --sc-ref-file sc_ref.h5ad \
     --output-root output \
     --select-ct T \
     --ist-mapping mean

The public result is ``output/sample/SVC.h5ad`` and its manifest records
``result.type`` as ``sc-SVC``. ``--ist-mapping random`` uses the CLI seed to
choose a same-cluster reference expression row instead of the cluster mean.

sST route
---------

sST selects the spot super-resolution sc-SVC route:

.. code-block:: bash

   revise-reconstruct \
     --platform sST \
     --sample-name sample \
     --data-root data \
     --st-file st.h5ad \
     --sc-ref-file sc_ref.h5ad \
     --output-root output

The public result is ``output/sample/SVC.h5ad`` and its manifest records
``result.type`` as ``sc-SVC``. Missing curated spot-to-cell and PM-on-cell
inputs use the quota and seeded random assignment described in
:doc:`concepts`.

Shared provenance
-----------------

The public AnnData links to the canonical run's ``provenance.json``. The
manifest distinguishes public result type from internal route and records
requested, attempted, and completed stages. Inspect it rather than inferring
success from directory existence.

Paper notebook compatibility
----------------------------

Historical wrappers retain notebook-specific copies such as ``sp_SVC.h5ad``,
``sc_SVC_expr.h5ad``, and ``sc_SVC_spatial.h5ad``. These are checkout-only
compatibility paths and do not change the installed ``SVC.h5ad`` contract.
