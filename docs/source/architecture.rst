Architecture
============

REVISE has one public orchestration API and one fixed reconstruction lifecycle.
Profiles select route behavior; registries select platform adapters and local
strategies; the configuration selects OT solvers for GA and LR.

System map
----------

.. code-block:: text

   revise-reconstruct / reconstruct.py
   `-- revise.cli
       `-- REVISEPipeline.run()
           |-- load and validate merged configuration
           |-- resolve and preflight route inputs
           |-- create canonical run envelope
           `-- UnifiedReconstructionPipeline
               |-- validate_inputs
               |-- global_anchoring
               |-- prepare_local_units
               |-- build_graph
               |-- build_ot_problem
               |-- solve_ot
               |-- update_expression
               |-- finalize_svc
               `-- evaluate_if_needed

``REVISEPipeline`` owns configuration, route resolution, input preflight,
deterministic setup, run/provenance lifecycle, and strategy dispatch. The
unified pipeline owns stage order. Strategies change stage internals without
creating a second lifecycle.

Configuration and routing
-------------------------

``revise/revise.yaml`` contains:

- ``defaults`` for runtime, IO, columns, preprocessing, graph, OT, and route
  behavior;
- ``profiles`` for declared application and benchmark requests;
- ``router`` mappings from platform/confounding to strategy and adapter;
- ``locked_params`` for governed low-level values.

GA uses ``ot.ga.solver`` and LR uses ``ot.lr.solver``. The plugin registry does
not choose an OT solver. Adapter translations preserve the same two-stage
selection even when a legacy runner has a different internal call shape.

Run evidence
------------

A full run allocates a unique canonical directory beneath the route-specific
output tree. It contains at least the merged configuration, logs, and
``provenance.json``; successful stages may add hashed artifacts and benchmark
metric tables. The exact directory leaf is unique and should be discovered
through the public result link or manifest rather than reconstructed from a
hard-coded timestamp pattern.

The canonical CLI separately publishes one stable-facing result:

.. code-block:: text

   <output-root>/<sample-name>/SVC.h5ad

The manifest records ``result.type`` as ``hST``, ``iST``, or ``sST`` and records
the internal platform route separately. Strategy artifacts remain in the
canonical run but are not additional public output contracts.

Failure model
-------------

Every lifecycle stage is recorded as pending, running, succeeded, failed,
skipped, or interrupted. A stage exception marks later stages skipped because
of upstream failure. Dry-run marks non-validation stages skipped. Abrupt process
death can leave incomplete states, which is evidence that the run did not
complete—not permission to infer success.

Solver telemetry records requested, attempted, and completed events separately.
A requested TACCO run cannot be reported as completed POT because TACCO does not
fall back to POT.

Extension boundary
------------------

Add a new route by extending the existing profile/router/registry surfaces and
their focused tests. Do not create another orchestration entrypoint or output
alias layer. Candidate evidence is intentionally limited to tested routes and
scales; see :doc:`limitations`.
