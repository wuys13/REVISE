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

Assignment-guidance boundary
----------------------------

The assignment subsystem separates inference, mandatory use, and optional
local-refinement guidance. Routes share one Assignment State carrier, one
label-aligned compatibility definition, and one ``off|prefer|require`` policy.
They do not share a single local-problem implementation:

- sp-SVC injects compatibility into neighbor or replacement OT cost;
- standard sc-SVC reweights existing Graph edges;
- sc-SVC-sr injects projected spot assignment into virtual-cell OT cost;
- imputation injects spatial-to-reference-subcluster compatibility into OT.

Argmax labels are represented as one-hot Assignment States rather than a
parallel mechanism. Standard sc-SVC uses Level1 for cohort routing and the
newly inferred Level2 soft state for Graph guidance. It selects resolution on
the unguided Graph and only then performs fixed-resolution guided clustering.
sc-SVC-sr's composition and closed-form expression allocation are
algorithm-defining steps and remain active when optional guidance is off.

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

The manifest records ``result.type`` as ``sp-SVC``, ``sc-SVC``, or
``sc-SVC-sr`` and records the internal route separately. Strategy artifacts remain in the
canonical run but are not additional public output contracts.

``provenance.json.assignment_guidance`` records schema version, configured
request, resolved policy, resolution source, and every local invocation. Each
event records route, operator, solver, compatibility numerics, bilateral
assignment source/level/value semantics/lineage, attempted state, terminal
outcome, and stable reason code. ``sr_allocation`` is adjacent durable evidence
because mandatory allocation is not a guidance outcome.

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

Guidance outcomes are likewise explicit: ``not_started``, ``not_applicable``,
``off``, ``applied``, ``fallback``, ``mixed``, ``failed``, or ``interrupted``.
Failure/interruption manifests retain earlier completed events. A required
guidance failure follows the normal publication rollback and cannot publish a
successful result.

Extension boundary
------------------

Add a new route by extending the existing profile, router, and strategy registry and
their focused tests. Do not create another orchestration entrypoint or output
alias layer. Candidate evidence is intentionally limited to tested routes and
scales; see :doc:`limitations`.
