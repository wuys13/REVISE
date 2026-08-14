Python API
==========

Use the workflow-level functions for new code. They preserve the same YAML
contracts as the command-line entry points and keep route-specific preparation
inside the package.

Application
-----------

.. code-block:: python

   from reconstruct import run_application

   result = run_application("configs/application/VisiumHD.yaml")

``run_application(config_path)`` compiles one Application YAML, loads its two
inputs, runs the explicit preprocessing sequence, reconstructs, publishes, and
returns the same in-memory ``AnnData`` object(s) that it writes to H5AD; it does
not reload them from disk. ``sp-SVC`` and ``sc-SVC-sr`` return one ``AnnData``.
Standard ``sc-SVC`` returns ``(spatial_adata, expression_adata)``.

.. autosummary::
   :toctree: generated
   :nosignatures:

   reconstruct.run_application

Sim2Real Benchmark
------------------

.. code-block:: python

   from revise.benchmark import run_benchmark

   report = run_benchmark(
       "configs/benchmark/segmentation.yaml",
       "raw_data/Sim2Real-ST",
       "P2CRC/cut_part1",
       "output/benchmark",
       dataset_task="segmentation",
       evaluate=True,
   )

``run_benchmark`` runs every leaf declared by one of the six Benchmark YAML
routes and returns a report mapping. Inspect ``report["ok"]``,
``total_runs``, ``passed_runs``, and each member of ``results``; do not infer
suite success from one output directory. A route/YAML error fails before the
report is created, while a leaf execution error is represented in its result.

.. autosummary::
   :toctree: generated
   :nosignatures:

   revise.benchmark.run_benchmark

Advanced engine objects
-----------------------

``REVISEPipeline`` is the shared engine interface used by both workflow
functions. It resolves one Application ``svc_type`` or Benchmark ``cf`` and
returns the canonical ``SVC`` carrier. Direct callers are responsible for the
correct route-specific IO and algorithm overrides; this is not a replacement
for the simpler Application YAML interface.

``SVC`` holds the expression/spatial carriers, assignment information,
provenance, quality metrics, and named artifacts produced inside the unified
engine.

.. autosummary::
   :toctree: generated
   :nosignatures:

   revise.framework.REVISEPipeline
   revise.svc.SVC

Metrics
-------

``compute_metric`` produces the Benchmark-compatible per-gene ``PCC``,
``SSIM``, ``MSE``, and ``NRMSE`` table. The caller must align observations;
the function aligns genes but does not join cells by identifier.
``compute_clustering_metrics`` returns ``(ARI, NMI)`` from two columns in
``adata.obs``.

.. autosummary::
   :toctree: generated
   :nosignatures:

   revise.analysis.compute_metric
   revise.analysis.compute_clustering_metrics

Internal runner configuration dataclasses, backend runner classes, pipeline
contexts, and compatibility output switches are implementation surfaces. They
remain importable where required by repository notebooks and tests, but they
are not recommended entry points for new users.
