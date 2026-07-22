Configuration
=============

``revise/revise.yaml`` contains defaults, route profiles, routing rules, and
locked low-level values. ``REVISEPipeline.run()`` resolves one configuration in
this order:

1. package defaults;
2. the selected profile;
3. runtime and IO overrides;
4. repeated dotted ``--set KEY=VALUE`` overrides.

Unknown keys and incomplete OT sections fail validation rather than being
silently ignored. The merged configuration is written into the canonical run
evidence.

GA and LR OT selection
----------------------

REVISE exposes two OT stages:

.. list-table::
   :header-rows: 1
   :widths: 1 2 2

   * - Stage
     - Meaning
     - Solver key
   * - GA
     - Global Anchoring between ST and the reference
     - ``ot.ga.solver``
   * - LR
     - Local Refinement inside the platform strategy
     - ``ot.lr.solver``

Each solver key accepts ``pot`` or ``tacco``. The installed/source CLI offers
one high-level switch:

.. code-block:: bash

   revise-reconstruct ... --ot-method tacco

This appends both ``ot.ga.solver=tacco`` and ``ot.lr.solver=tacco`` to the
resolved request. ``--ot-method pot`` does the same for POT. When the option is
omitted, the configuration retains the two stage values independently. A mixed
diagnostic request can therefore use:

.. code-block:: bash

   revise-reconstruct ... \
     --set ot.ga.solver=tacco \
     --set ot.lr.solver=pot

TACCO behavior
--------------

TACCO 0.5.0 is an optional dependency. Choosing it performs a dependency and
version check during preflight. A missing package, incompatible version,
unsupported conditioning combination, invalid result, or solver exception
fails the run. REVISE does not fall back to POT and records requested,
attempted, and completed solver events separately.

Reference-mode posterior conditioning is not compatible with TACCO LR. Select
a supported conditioning mode or POT LR; the configuration/preflight layer
rejects that combination before reconstruction.

POT parameter example
---------------------

This is an illustrative stage-scoped request, not the default for every route.
The resolved route profile remains the source of truth:

.. code-block:: yaml

   ot:
     ga:
       solver: pot
       pot:
         reg: 0.1
         reg_m: 0.0
         reg_type: entropy
     lr:
       solver: pot
       pot:
         reg: 0.1
         reg_m: 0.0
         reg_type: entropy

Locked parameters cannot be overridden by ``--set``. Change the governed
profile and its tests when a paper-facing low-level value needs to change.

Runtime and IO
--------------

The canonical CLI manages ``runtime.seed`` and ``io.save_outputs`` itself so a
conflicting ``--set`` cannot make the command summary disagree with the actual
run. Common IO values are:

.. code-block:: yaml

   io:
     data_root: data
     output_root: output
     sample_name: sample
     st_file: st.h5ad
     sc_ref_file: sc_ref.h5ad

The input resolver applies the same route-specific path rules in preflight and
full execution. The public result remains
``<output-root>/<sample-name>/SVC.h5ad``; the canonical run directory contains
``provenance.json`` and internal evidence.

Python API
----------

Programmatic callers use the same merge and validation path:

.. code-block:: python

   from revise.framework import REVISEPipeline

   pipeline = REVISEPipeline()
   result = pipeline.run(
       profile="application_sp",
       runtime_overrides={"platform": "sp_svc", "confounding": "bin2cell"},
       io_overrides={
           "data_root": "data",
           "output_root": "output",
           "sample_name": "sample",
           "st_file": "st.h5ad",
           "sc_ref_file": "sc_ref.h5ad",
       },
       set_overrides=["ot.ga.solver=pot", "ot.lr.solver=pot"],
   )

Direct API use returns an ``SVC`` carrier. The single public result file is a
contract of ``reconstruct.py``/``revise-reconstruct``, not of every
low-level pipeline call.
