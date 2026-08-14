Architecture
============

REVISE is unified at the reconstruction engine and bottom-level optimal
transport implementation. Application and Benchmark keep separate,
task-appropriate request and preprocessing layers; both converge on the same
``REVISEPipeline`` and fixed reconstruction lifecycle.

User-facing flow
----------------

.. code-block:: text

   Application YAML
       -> reconstruct.run_application
       -> load -> filter -> preprocess -> reconstruct -> publish
                              |
                              v
                         REVISEPipeline
                              ^
                              |
   Benchmark YAML + CLI
       -> revise.benchmark.run_benchmark
       -> expand route cases -> route-specific preparation

                         REVISEPipeline
                              |
                              v
                UnifiedReconstructionPipeline
       validate -> Global Anchoring -> Local Refinement -> finalize
                                                    -> evaluate (Benchmark only)
                              |
                              v
                 route strategy -> shared OTKernel

For new Application data, ``reconstruct.py`` deliberately keeps preprocessing
visible. It compiles the YAML, loads spatial and reference ``AnnData``, applies
an explicit reference filter when requested, preprocesses both inputs, prepares
the standard ``sc-SVC`` pair or normalizes reference labels, and only then
calls the shared engine. Benchmark preprocessing remains route-specific
because its six confounding families have different case layouts and
ground-truth roles.

Routing and lifecycle
---------------------

``revise.config.authority`` owns engine defaults and the typed route mapping.
Application selects one ``application.svc_type``; Benchmark selects one YAML
``route``. The frontends translate only their validated request into runtime,
IO, and algorithm values. Users do not pass the internal engine authority as a
run configuration.

``REVISEPipeline.run`` resolves the route, merges package defaults and the
validated frontend values, creates the run directory and manifest, sets the
run seed, and dispatches one strategy. ``UnifiedReconstructionPipeline`` then
keeps the stage order fixed:

1. validate inputs;
2. Global Anchoring;
3. Local Refinement: prepare local units, build the graph and OT problem,
   solve OT, and update expression;
4. finalize the ``SVC``;
5. evaluate when an enabled Benchmark request has aligned ground truth.

Strategies may implement the work inside a stage differently for
``sp-SVC``, ``sc-SVC``, ``sc-SVC-sr``, and the six Benchmark families. They do
not create a second lifecycle.

One OT implementation
---------------------

POT and TACCO share ``revise.backend.kernels.ot.OTKernel`` as the bottom-level
facade:

.. code-block:: text

   Global Anchoring
       GlobalAnchoringKernel -> OTKernel.annotate -> POT or TACCO

   Local Refinement
       standard sc-SVC: LocalAnchoringKernel -> OTKernel.annotate
       other local routes: route runner -> OTKernel.couple
                                      -> POT or TACCO

The Local Refinement stage does not construct or call a
``GlobalAnchoringKernel``. Standard ``sc-SVC`` may use its distinct
``LocalAnchoringKernel`` for local subtype annotation; the other local routes
use the shared coupling operation directly.

An Application ``algorithm.ot_method`` sets both ``ot.ga.solver`` and
``ot.lr.solver``. The selected implementation must be installed and satisfy
its contract. TACCO and POT are never substituted for one another on failure.

Assignment ownership
--------------------

Global Anchoring produces a validated posterior ``Q`` and ``argmax(Q)`` broad
labels. Each route owns how that result enters Local Refinement:

- sp-SVC conditions each local OT cost with ``Q``;
- sc-SVC-sr first projects spot-level ``Q`` to virtual cells, then conditions
  the local OT cost;
- standard sc-SVC uses only ``argmax(Q)`` to choose broad cohorts and does not
  reweight GraphCluster with ``Q``;
- gene-panel and gene-dropout imputation do not expose assignment-conditioned
  Local Refinement.

There is no user-selectable fallback or policy state machine around this
boundary. TACCO and POT share one Global Assignment contract: a wholly-NaN
posterior row is published as ``Unknown`` with NaN confidence, while every
assigned row must satisfy the required axes, finite numeric values and argmax
label. Consumers that numerically condition on ``Q`` require a complete finite
posterior; routes that use only labels may retain the unassigned observation.

Results and publication
-----------------------

The engine returns a canonical ``SVC`` carrier containing expression/spatial
objects, assignment information, artifacts, metrics, and provenance. The
Application publication callback exposes only the route's promised objects:

- ``sp-SVC`` and ``sc-SVC-sr`` publish one H5AD and return one ``AnnData``;
- standard ``sc-SVC`` publishes and returns the fixed ``(spatial,
  expression)`` pair.

The entry point writes same-directory temporary H5AD files before replacing
the public targets. Paired outputs are not reader-atomic or crash-atomic, and
replacement does not provide rollback. The caller must guarantee one writer
per stable public target; violating that precondition is undefined.

Run evidence and failure
------------------------

Every canonical run writes ``merged_config.json`` and ``provenance.json``.
The manifest status is one of ``running``, ``succeeded``, and ``failed``.
Captured SIGTERM and KeyboardInterrupt become failed with error evidence;
ordinary stage exceptions are recorded and re-raised. An uncatchable
termination leaves the last manifest running, which is incomplete evidence,
not success.

``input_identities`` records one content identity per external role; there is
no aggregate data fingerprint. Software identity is collected once per run.
The manifest keeps resolved OT configuration and minimal Local Refinement
evidence: ``route``, ``applied``, and ``strength``. It has no OT or Assignment
event state machine. A stage error in the manifest is the authoritative
failure explanation.

For public fields and signatures, see :doc:`api/index`.
