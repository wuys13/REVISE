Configuration
=============

``revise/revise.yaml`` contains defaults, route profiles, routing rules, and
locked low-level values. ``REVISEPipeline.run()`` resolves one configuration in
this order:

1. package defaults;
2. the selected profile;
3. runtime and IO overrides;
4. dedicated high-level options such as ``--ot-method``.

Unknown keys and incomplete OT sections fail validation rather than being
silently ignored. The merged configuration is written into the canonical run
evidence. For advanced settings, copy ``revise/revise.yaml``, edit the relevant
profile, and select that file with ``--config``. There is no generic CLI
key/value override surface.

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

This selects TACCO for both ``ot.ga.solver`` and ``ot.lr.solver`` in the
resolved request. ``--ot-method pot`` does the same for POT. When the option is
omitted, the configuration retains the two stage values independently. A mixed
diagnostic request can therefore be written in a custom configuration:

.. code-block:: yaml

   ot:
     ga:
       solver: tacco
     lr:
       solver: pot

The ``application_sc`` profile defaults both stages to TACCO and reads the
high-level annotation arguments from the profile itself:

.. code-block:: yaml

   sc:
     resolutions: [0.6, 0.7, 0.8]
     tacco_annotate:
       multi_center: 1
       lamb: 0.001

These two values are forwarded to the Level1, Level2, and ``SVC_cluster``
``tacco.tl.annotate`` calls. They are validated and locked against runtime
algorithm overrides.

Run the normal command with ``--config path/to/custom-revise.yaml``.

TACCO behavior
--------------

TACCO 0.5.0 is an optional dependency. Choosing it performs a dependency and
version check during preflight. A missing package, incompatible version,
unsupported conditioning combination, invalid result, or solver exception
fails the run. REVISE does not fall back to POT and records requested,
attempted, and completed solver events separately. If TACCO is unavailable and
a different algorithm is acceptable, application users may explicitly rerun
with ``--ot-method pot``; this is a solver change, not an automatic fallback.

Assignment and local refinement
-------------------------------

Global anchoring produces a validated posterior ``Q`` and ``argmax(Q)`` labels.
The only public local-refinement option is:

.. code-block:: yaml

   local_refinement:
     strength: 0.2

sp-SVC defaults to ``0.2``. sc-SVC-sr defaults to ``0.0`` and accepts an
explicit non-negative finite strength. Standard sc-SVC and imputation reject
this option and omit ``local_refinement`` from resolved configuration. sp-SVC
uses ``Q`` in local OT; sc-SVC-sr first projects ``Q`` to virtual cells;
standard sc-SVC uses only ``argmax(Q)`` for cohort routing. There are no
policy, compatibility-mode, reference-mode, or one-hot-fallback settings.
Application and benchmark entrypoints expose the same single override as
``--local-refinement-strength``.

The manifest records exactly ``route``, ``applied``, and ``strength`` under
``local_refinement``. Mandatory sc-SVC-sr expression allocation is recorded
separately under ``sr_allocation``.

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

Locked parameters have no generic CLI bypass. Change the governed profile and
its tests when a paper-facing low-level value needs to change.

Runtime and IO
--------------

The canonical CLI manages task identity, ``runtime.seed``, and its public
inputs and outputs through dedicated parameters. Common IO values are:

.. code-block:: yaml

   io:
     data_root: data
     output_root: output
     sample_name: sample
     st_file: st.h5ad
     sc_ref_file: sc_ref.h5ad

The input resolver applies the same route-specific path rules in preflight and
full execution. sp-SVC and sc-SVC-sr publish
``<output-root>/<sample-name>/SVC.h5ad``. Standard sc-SVC publishes its pair
under ``<output-root>/<sample-name>/sc-SVC/<cell-type>/``. The canonical run
directory contains ``provenance.json`` and internal evidence.

Python API
----------

Programmatic callers use the same configuration-resolution and validation path:

.. code-block:: python

   from revise.framework import REVISEPipeline

   pipeline = REVISEPipeline(config_path="path/to/custom-revise.yaml")
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
   )

Direct API use returns an ``SVC`` carrier. Route-specific public result
publication is a contract of ``reconstruct.py``/``revise-reconstruct``, not of
every low-level pipeline call.
