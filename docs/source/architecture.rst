Architecture
============

REVISE has one reconstruction engine and a shared optimal-transport layer.
Application and Sim2Real-ST Benchmark deliberately keep separate request,
preprocessing, and publication boundaries.

Frontend boundary
-----------------

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
       -> expand experimental cases -> route-specific preparation

                         REVISEPipeline
                              |
                              v
                UnifiedReconstructionPipeline
       validate -> Global Anchoring -> Local Refinement -> finalize
                                                    -> evaluate (Benchmark only)
                              |
                              v
                 route strategy -> shared OTKernel

``reconstruct.py`` keeps Application preparation visible: it compiles one
YAML, filters and preprocesses the two inputs, prepares a cluster-mode pair or
normalizes SR/reference labels, and then invokes the engine. Benchmark owns
its six experimental case layouts and ground-truth roles; Application mode is
not added to Benchmark YAML.

Application routing
-------------------

``revise.config.authority`` owns package defaults and typed route resolution.
The public Application request resolves to one of these internal records:

.. list-table::
   :header-rows: 1
   :widths: 2 2 3 3 3

   * - Public request
     - Selector
     - Profile / task
     - Strategy
     - Runner configuration
   * - ``sp-SVC``
     - ``sp-SVC``
     - ``application_sp`` / ``sp_svc``
     - ``SpSvcApplicationStrategy``
     - ``ApplicationSpConf``
   * - ``sc-SVC`` cluster
     - ``sc-SVC:cluster``
     - ``application_sc`` / ``sc_svc``
     - ``ScSvcApplicationStrategy``
     - ``ApplicationScConf``
   * - ``sc-SVC`` sr
     - ``sc-SVC:sr``
     - ``application_sc_super_resolution`` / ``sc_svc_super_resolution``
     - ``ScSvcSuperResolutionApplicationStrategy``
     - ``ApplicationScSuperResolutionConf``

Runtime metadata records ``application_route: sc-SVC`` plus either
``application_mode: cluster`` or ``application_mode: sr`` for the two sc-SVC
cases. The mode is
also part of the route key, so cluster and SR runs do not share a run directory.
The App SR runner is
``sc_svc_super_resolution_application.ScSVCSuperResolution``.

Benchmark retains its own ``sc_svc_sr`` task, profile, runner, and compatible
result naming. It is a Benchmark implementation detail, not an Application
SVC category.

Assignment and Local Refinement
-------------------------------

Global Anchoring produces a validated posterior ``Q`` and broad labels. Each
route owns how it enters Local Refinement:

- sp-SVC conditions each local OT cost with ``Q``;
- sc-SVC sr mode projects spot-level ``Q`` to virtual cells, then conditions
  the local OT cost;
- sc-SVC cluster mode uses broad labels to choose cohorts and does not reweight
  GraphCluster with ``Q``;
- Benchmark gene-panel and gene-dropout imputation do not expose
  assignment-conditioned Local Refinement.

POT and TACCO meet at ``revise.backend.kernels.ot.OTKernel``. Application
``algorithm.ot_method`` sets both ``ot.ga.solver`` and ``ot.lr.solver``. There
is no user-selectable solver fallback.

Publication and failure contract
--------------------------------

The engine returns a canonical ``SVC`` carrier. Application publication
exposes only the promised route artifacts: one H5AD for sp-SVC/SR and the
fixed spatial/expression pair for cluster mode. It writes same-directory
temporary H5AD files before replacing public targets. Paired outputs are not
reader-atomic or crash-atomic, but catchable replacement failures attempt
rollback; one writer per stable public target is a caller precondition.

Each canonical run writes ``merged_config.json`` and ``provenance.json``. The
manifest state is ``running``, ``succeeded``, or ``failed``. Catchable errors,
SIGTERM, and KeyboardInterrupt record failure evidence and re-raise. An
uncatchable termination can leave the last manifest ``running``; that is
incomplete evidence, not success. The stage error in the manifest is the
authoritative failure explanation.

For public YAML fields use :doc:`application-reference`; for Python function
signatures use :doc:`api/index`.
