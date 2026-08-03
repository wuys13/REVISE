Architecture
============

REVISE has one public orchestration API and one fixed reconstruction lifecycle.
Profiles and the router select one reconstruction strategy; the configuration
selects OT solvers for GA and LR.

System map
----------

.. code-block:: text

   reconstruct.py / revise-reconstruct
   `-- revise.application.cli
       `-- revise.application.service
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
- ``router`` mappings from internal route/confounding to strategy;
- ``locked_params`` for governed low-level values.

GA uses ``ot.ga.solver`` and LR uses ``ot.lr.solver``. Runner translations
preserve the same two-stage selection even when a legacy runner has a different
internal call shape.

Assignment boundary
-------------------

Global anchoring produces one validated posterior ``Q`` and its ``argmax(Q)``
labels. Downstream ownership is explicit per route:

- sp-SVC conditions each local OT cost with ``Q``;
- sc-SVC-sr projects spot-level ``Q`` to virtual cells, then conditions local
  OT cost;
- standard sc-SVC uses only ``argmax(Q)`` to split broad cohorts and does not
  reweight GraphCluster with ``Q``;
- imputation does not expose assignment-based local refinement.

There is no optional policy or fallback state machine. ``Q`` must already have
the expected observation/category axes, finite non-negative values, and
positive row mass. sc-SVC-sr composition and closed-form expression allocation
remain mandatory algorithm steps independent of local OT conditioning.

Run evidence
------------

A full run allocates a unique canonical directory beneath the route-specific
output tree. It contains at least the merged configuration, logs, and
``provenance.json``; successful stages may add hashed artifacts and benchmark
metric tables. The exact directory leaf is unique and should be discovered
through the public result link or manifest rather than reconstructed from a
hard-coded timestamp pattern.

The canonical CLI publishes stable-facing application results. sp-SVC and
sc-SVC-sr use:

.. code-block:: text

   <output-root>/<sample-name>/SVC.h5ad

Standard sc-SVC publishes ``sc_SVC_spatial.h5ad`` and ``sc_SVC_expr.h5ad``
under a cell-type subdirectory. The manifest records each public result role
and the internal route separately. Strategy artifacts remain in the canonical
run but are not additional public output contracts.

For both the single-file and paired 1.x outputs, the publisher writes a
same-directory temporary H5AD, reloads it before replacement, and provides
best-effort caught-exception rollback. The pair is not reader-atomic or
crash-atomic. The caller must guarantee one writer per stable public target;
violating that precondition is undefined.

``provenance.json.local_refinement`` is the minimal route-level evidence:
``route``, ``applied``, and ``strength``. For standard sc-SVC and imputation,
``strength`` is ``null`` because those routes do not accept the option.
``sr_allocation`` remains adjacent durable evidence for mandatory SR
allocation.

Failure model
-------------

Run status is limited to ``running``, ``succeeded``, and ``failed``. Captured
SIGTERM and KeyboardInterrupt become failed with their error evidence, as do
captured stage exceptions; later stages are skipped because of upstream
failure. Dry-run marks non-validation stages skipped. An uncatchable
termination leaves the last manifest running, which is evidence that the run
did not complete—not permission to infer success.

``input_identities`` records one content identity per external role; there is
no aggregate data fingerprint. Software identity is collected once per run.
The manifest retains the resolved ``ot_config`` and minimal
``local_refinement`` evidence, but has no OT or Assignment event state machine.

The ``local_refinement.applied`` flag changes to true only after a local OT
conditioning call succeeds. Failure and interruption continue through the
normal stage error and publication rollback; stage/run errors remain the
authoritative failure explanation.

Extension boundary
------------------

Add a new route by extending the existing profile, router, and strategy registry and
their focused tests. Do not create another orchestration entrypoint or output
alias layer. Candidate evidence is intentionally limited to tested routes and
scales; see :doc:`limitations`.
