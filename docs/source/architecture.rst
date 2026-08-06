Architecture
============

REVISE has one public orchestration API and one fixed reconstruction lifecycle.
Profiles and the router select one reconstruction strategy; the configuration
selects OT solvers for GA and LR.

System map
----------

.. code-block:: text

   reconstruct.py / revise-reconstruct
       `-- reconstruct.run_application
           `-- REVISEPipeline.run(svc_type=..., cf=None)
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

``reconstruct.py`` owns the public argument parser, YAML overrides, engine
mapping, and H5AD output writing. The package ``revise.application.config``
module only compiles YAML into a validated configuration. Application and
Benchmark meet at ``REVISEPipeline.run``; the engine router resolves the
profile/task/strategy for the selector supplied by each frontend.

Configuration and routing
-------------------------

``revise/revise.yaml`` contains:

- ``defaults`` for runtime, IO, columns, preprocessing, graph, OT, and route
  behavior;
- ``profiles`` for declared application and benchmark requests;
- ``router.application`` mappings from SVC data type to profile and strategy;
- ``router.benchmark`` mappings from confounding family to profile and strategy;
- ``locked_params`` for governed low-level values.

GA uses ``ot.ga.solver`` and LR uses ``ot.lr.solver``. Runner translations
preserve the same two-stage selection even when a legacy runner has a different
internal call shape.

Application and Benchmark keep separate request/case preparation. They meet at
one engine execution Seam: Application supplies ``svc_type`` and no ``cf``;
Benchmark supplies ``cf`` and no ``svc_type``. If both are supplied,
``svc_type`` wins with a warning. Application provenance records
``application_route`` and never records Benchmark ``confounding`` semantics.

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

Every Application invocation, including ``--dry-run``, allocates exactly one
run directory at

.. code-block:: text

   <output.dir>/<output.name>/application__<svc-type>/<timestamp_uuid>/

The run envelope contains ``provenance.json``, ``merged_config.json``, and
``run.log``; completed input validation adds ``preflight.json``. A successful
formal run writes its final H5AD artifact(s) directly into that same directory:
sp-SVC and sc-SVC-sr write
``<output.name>.h5ad``; standard sc-SVC writes
``<output.name>_spatial.h5ad`` and ``<output.name>_expr.h5ad`` as sibling
files. The sc-SVC pair therefore shares one run directory and provenance
record. A dry run receives an independent run directory with no H5AD. Every
repeated formal run receives a new ``<timestamp_uuid>`` leaf, leaving previous
runs untouched. Application creates no temporary H5AD files or shared
fixed-name copies and does not use ``os.replace`` for H5AD output.

The manifest records each result role and the internal route separately.
Strategy artifacts, when present, remain in the same canonical run directory
but are not additional public output contracts.

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
termination leaves the last manifest running; a failure while persisting
terminal provenance can do the same. That is evidence that the run did not
complete—not permission to infer success.

``input_identities`` records one content identity per external role; there is
no aggregate data fingerprint. Software identity is collected once per run.
The manifest retains the resolved ``ot_config`` and minimal
``local_refinement`` evidence, but has no OT or Assignment event state machine.

The ``local_refinement.applied`` flag changes to true only after at least one
route-owned local refinement unit completes successfully. It is independent of
posterior conditioning and its strength: a completed local OT refinement with
strength zero is still applied. Failure and interruption continue through the
normal stage error path; stage/run errors remain the authoritative failure
explanation.

Extension boundary
------------------

Add a new route by extending the existing profile, router, built-in route
contract, and strategy registry together with their focused tests. The loader
rejects route fields that disagree with that contract before execution. Do not
create another orchestration entrypoint or output alias layer. Candidate
evidence is intentionally limited to tested routes and scales; see
:doc:`limitations`.
