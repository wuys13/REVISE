Reconstruction Routes
=====================

The canonical installed/source interface selects a route with ``--platform``.
The route names describe input and output contracts, not a validated biological
application claim.

hST
---

hST selects the sp-SVC application profile and the bin-to-cell route:

.. code-block:: bash

   revise-reconstruct \
     --platform hST \
     --sample-name sample \
     --data-root data \
     --st-file st.h5ad \
     --sc-ref-file sc_ref.h5ad \
     --output-root output \
     --ot-method pot

The public result is ``output/sample/hST-SVC.h5ad``.

iST
---

iST selects the sc-SVC imaging route. The internal strategy returns separate
spatial and expression carriers; the CLI validates their cluster sets and
publishes one merged result:

.. code-block:: bash

   revise-reconstruct \
     --platform iST \
     --sample-name sample \
     --data-root data \
     --st-file st.h5ad \
     --sc-ref-file sc_ref.h5ad \
     --output-root output \
     --select-ct T \
     --ist-mapping mean \
     --ot-method pot

The public result is ``output/sample/iST-SVC.h5ad``. ``--ist-mapping random``
uses the CLI seed to choose a same-cluster reference expression row instead of
the same-cluster mean. This is a reproducible mapping rule, not evidence that a
particular reference cell is the true biological match.

sST
---

sST selects the spot super-resolution route:

.. code-block:: bash

   revise-reconstruct \
     --platform sST \
     --sample-name sample \
     --data-root data \
     --st-file st.h5ad \
     --sc-ref-file sc_ref.h5ad \
     --output-root output \
     --ot-method pot

The public result is ``output/sample/sST-SVC.h5ad``. If the input has no curated
spot-to-cell mapping and no ``pm_on_cell`` score matrix, the route uses the
count/quota and seeded random assignment described in :doc:`concepts`.

Shared output and provenance
----------------------------

All three commands follow:

.. code-block:: text

   <output-root>/<sample-name>/<platform>-SVC.h5ad

The public AnnData records a relative link to the canonical run's
``provenance.json``. The manifest distinguishes requested, attempted, and
completed OT stages and records failures or interruption. Inspect it together
with the public file rather than inferring success from an output directory.

Paper notebook compatibility
----------------------------

The source repository retains legacy application reconstruction scripts for
the historical notebooks. They are checkout-only compatibility paths and
publish notebook-specific copies; they are not installed commands and do not
change the single public result contract above.

Their current compatibility copies are:

.. code-block:: text

   output/sp_SVC_case/<sample-name>/sp_SVC.h5ad
   output/sc_SVC_case/<sample-name>_<st-file-stem>/<select-ct>/sc_SVC_expr.h5ad
   output/sc_SVC_case/<sample-name>_<st-file-stem>/<select-ct>/sc_SVC_spatial.h5ad

The notebooks and their data/results identities are not current RC scientific
evidence. They remain in ``REVISE-legacy`` and can be recovered through the
exact commit and paths described in :doc:`gallery`.
