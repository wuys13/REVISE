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

- hST-SVC conditions each local OT cost with ``Q``;
- sST-SVC projects spot-level ``Q`` to virtual cells, then conditions local
  OT cost;
- iST-SVC uses only ``argmax(Q)`` to split broad cohorts and does not
  reweight GraphCluster with ``Q``;
- imputation does not expose assignment-based local refinement.

There is no optional policy or fallback state machine. ``Q`` must already have
the expected observation/category axes, finite non-negative values, and
positive row mass. sST-SVC composition and closed-form expression allocation
remain mandatory algorithm steps independent of local OT conditioning.

Run evidence
------------

A full run allocates a unique canonical directory beneath the route-specific
output tree. It contains at least the merged configuration, logs, and
``provenance.json``; successful stages may add hashed artifacts and benchmark
metric tables. The exact directory leaf is unique; locate its
``provenance.json`` beneath the route-specific output tree rather than
reconstructing a hard-coded timestamp pattern.

The canonical CLI publishes one stable-facing application result for every
public selector:

.. code-block:: text

   <output-root>/<sample-name>/<svc-type>/SVC.h5ad

The manifest ``result`` contains exactly ``filename`` and ``type``. Only
iST-SVC adds the top-level ``assembly`` record for mean/random construction.
Strategy carriers and other artifacts remain in the canonical run but are not
additional public output contracts.

The single-file publisher writes a same-directory temporary H5AD, reloads it
before replacement, and provides best-effort caught-exception rollback. It is
not reader-atomic or crash-atomic. The caller must guarantee one writer per
stable public target; violating that precondition is undefined.

For the 2.0 iST-SVC route, ``select_ct: null`` runs validation and GA, then
writes ``selection_assessment.json`` and returns ``needs_review``. The report
excludes labels containing ``tumor`` or ``epi`` (case-insensitive), warns for
any label with more than 20,000 GA spots, and never silently selects a largest
class. Human confirmation is supplied by repeating ``--select-ct`` for each
requested cell type; only then does iST refinement and ``SVC.h5ad`` publication
run. The assessment records the resolved input identities for the current
reference; this version does not rank multiple references automatically.

Inapplicable reconstruction metadata values that are ``None`` are omitted by
H5AD serialization; it does not invent sentinel values. Applicable keys remain
exact.

``provenance.json.local_refinement`` is the minimal route-level evidence:
``route``, ``applied``, and ``strength``. For iST-SVC and imputation,
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

The ``local_refinement.applied`` flag changes to true only after at least one
route-owned local refinement unit completes successfully. It is independent of
posterior conditioning and its strength: a completed local OT refinement with
strength zero is still applied. Failure and interruption continue through the
normal stage error and publication rollback; stage/run errors remain the
authoritative failure explanation.

Extension boundary
------------------

Add a new route by extending the existing profile, router, and strategy registry and
their focused tests. Do not create another orchestration entrypoint or output
alias layer. Candidate evidence is intentionally limited to tested routes and
scales; see :doc:`limitations`.
