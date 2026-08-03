API
===

The stable public entry point is ``REVISEPipeline``. It resolves
``revise/revise.yaml`` profiles, applies runtime and IO overrides, and routes
each run to the appropriate backend strategy.

Start with :doc:`../quickstart` for runnable examples and :doc:`../architecture`
for the full lifecycle. This page is the reference surface for classes and
extension points.

.. figure:: classes_revise.svg
   :alt: REVISE API architecture
   :align: center

   Current API architecture.

.. code-block:: python

    from revise.framework import REVISEPipeline

    pipeline = REVISEPipeline()
    svc = pipeline.run(
        profile="application_sc",
        runtime_overrides={"platform": "sc_svc", "confounding": "segmentation"},
        io_overrides={
            "data_root": "raw_data/Real_application",
            "output_root": "output/sc_SVC_case",
            "sample_name": "P2CRC",
            "st_file": "Xenium.h5ad",
            "sc_ref_file": "adata_sc_all_reanno.h5ad",
        },
    )

Pipeline
~~~~~~~~

.. autosummary::
    :toctree: generated
    :nosignatures:

    revise.framework.REVISEPipeline
    revise.recon.context.PipelineContext
    revise.recon.pipeline.UnifiedReconstructionPipeline
    revise.svc.SVC

Facade Helpers
~~~~~~~~~~~~~~

The facade helpers are thin wrappers around ``REVISEPipeline.run`` for callers
that want named sp-SVC or sc-SVC functions while still using the unified config
surface.

.. autosummary::
    :toctree: generated
    :nosignatures:

    revise.recon.facade.sp_svc
    revise.recon.facade.sc_svc

Configuration
~~~~~~~~~~~~~

``revise/revise.yaml`` is the external configuration and routing surface.
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
New application and benchmark code should prefer ``REVISEPipeline`` or the
application or benchmark entrypoints.
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
