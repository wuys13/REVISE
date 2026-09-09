Python API
==========

Use the workflow-level functions for new code. They preserve the same YAML
contracts as the command-line entry points and keep route-specific preparation
inside the package.

Application
-----------

.. code-block:: python

   from reconstruct import run_application

   spatial, expression = run_application(
       "configs/application/Xenium.yaml",
       select_ct="T",
   )

``run_application(config_path, *, select_ct=None)`` compiles one Application
YAML, loads and preprocesses its inputs, reconstructs, publishes, and returns
the same in-memory objects that it writes. ``select_ct`` is an optional
cluster-mode-only override; it has priority over the YAML value. The current
Application output rules are defined in :doc:`../application-reference`.

``sp-SVC`` and ``sc-SVC`` SR mode return one ``AnnData``. ``sc-SVC`` cluster
mode returns ``(spatial_adata, expression_adata)``.

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

``run_benchmark`` expands the selected Benchmark YAML's cases and returns a
report mapping. Inspect ``report["ok"]``, ``total_runs``, ``passed_runs``, and
each member of ``results``; one output directory is not proof that a suite
succeeded.

.. autosummary::
   :toctree: generated
   :nosignatures:

   revise.benchmark.run_benchmark

Advanced engine boundary
------------------------

``REVISEPipeline`` is the shared engine interface. Direct Application callers
must pass ``svc_type="sc-SVC"`` with
``application_mode="cluster"`` or ``application_mode="sr"``. They are
responsible for correct route-specific IO and algorithm overrides, so this is
not a replacement for the YAML entry point.

``SVC`` holds the canonical expression/spatial carriers, assignment
information, artifacts, metrics, and provenance produced by the engine.

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
``compute_clustering_metrics`` returns ``(ARI, NMI)`` from two ``adata.obs``
columns.

.. autosummary::
   :toctree: generated
   :nosignatures:

   revise.analysis.compute_metric
   revise.analysis.compute_clustering_metrics
