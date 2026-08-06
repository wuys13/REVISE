API
===

The stable public Application Python entry point is ``run_application``. It
accepts the same YAML contract as the CLI. ``REVISEPipeline`` remains the
engine-level interface; Benchmark callers select one confounding family with
``cf``. The package-owned router then resolves the profile and strategy.

Start with :doc:`../quickstart` for runnable examples and :doc:`../architecture`
for the full lifecycle. This page is the reference surface for classes and
extension points.

.. figure:: classes_revise.svg
   :alt: REVISE API architecture
   :align: center

   Current API architecture.

.. code-block:: python

    from reconstruct import run_application

    execution = run_application(
        "configs/application/Xenium_T.yaml",
        dry_run=True,
    )

Pipeline
~~~~~~~~

Application entry point
~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: generated
    :nosignatures:

    reconstruct.run_application
    reconstruct.ApplicationExecution

The engine pipeline remains available for Benchmark and advanced integrations:

.. autosummary::
    :toctree: generated
    :nosignatures:

    revise.framework.REVISEPipeline
    revise.recon.context.PipelineContext
    revise.recon.pipeline.UnifiedReconstructionPipeline
    revise.svc.SVC

Configuration
~~~~~~~~~~~~~

``revise/revise.yaml`` is the package-owned engine configuration and routing
authority; application YAML is the external request surface.
``revise.config.runner_conf`` contains internal runner contracts used by
backend adapters and compatibility notebooks.

.. autosummary::
    :toctree: generated
    :nosignatures:

    revise.config.runner_conf.ApplicationSpConf
    revise.config.runner_conf.ApplicationScConf
    revise.config.runner_conf.ApplicationScSrConf
    revise.config.runner_conf.BenchmarkSrConf
    revise.config.runner_conf.BenchmarkSegConf
    revise.config.runner_conf.BenchmarkImputeConf

Analysis Services
~~~~~~~~~~~~~~~~~

Analysis services consume the unified ``SVC`` result carrier and provide
notebook-compatible downstream helpers.

.. autosummary::
    :toctree: generated
    :nosignatures:

    revise.analysis.ScSVCAnalysisService
    revise.analysis.SpSVCAnalysisService
    revise.analysis.compute_metric
    revise.analysis.compute_clustering_metrics

Strategy contract and registry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These classes are the extension points used by the unified backend.

.. autosummary::
    :toctree: generated
    :nosignatures:

    revise.backend.registry.StrategyRegistry
    revise.backend.contracts.LocalRefinementStrategy

Backend Compatibility Runners
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These classes are kept for notebooks, parity checks, and low-level debugging.
New Application code should import ``run_application`` from ``reconstruct``.
Benchmark and low-level integration code may use ``REVISEPipeline`` or the
benchmark entrypoint.
Direct sp-SVC and sc-SVC-sr runner callers receive the resolved
``local_refinement_strength`` value. Standard sc-SVC and imputation callers do
not provide it. Assignment ``Q`` is validated at the global-assignment boundary
and is not synthesized or repaired inside a runner.

.. autosummary::
    :toctree: generated
    :nosignatures:

    revise.backend.runners.sp_svc_application.SpSVC
    revise.backend.runners.sc_svc_application.ScSVC
    revise.backend.runners.sc_svc_sr_application.ScSVCSr
    revise.backend.runners.sp_svc_benchmark.SpSVC
    revise.backend.runners.sc_svc_sr_benchmark.ScSVCSr
    revise.backend.runners.sc_svc_impute_benchmark.ScSVCImpute
