.. _batch-framework:

Batch code framework
====================

The public directory and configuration contract is defined in
:ref:`batch-input-protocol`. This page explains its implementation and the
extension point for analysis. Start with :ref:`batch-reconstruction` to run it.

Execution and module boundaries
-------------------------------

.. code-block:: text

   batch_reconstruct.py / revise-batch-reconstruct
     → runner.run_batch
       → config.discover_samples: find spatial.h5ad sample directories
       → config.resolve_sample: inherit settings and record their sources
       → sample.read_sample: check standard H5AD without writing inputs
       → runner.application_document: translate to the single-run contract
       → reconstruct subprocess: existing QC, GA/LR and publication
       → reconstruction.json and batch_status.json

   batch_analyze.py / revise-batch-analyze
     → analysis.run_analysis_batch
       → current batch inventory and shared configuration resolver
       → verify reconstruction and its inputs/configuration
       → run requested aspect adapters and publish their artifacts
       → analysis/analysis.json and analysis_status.json

The modules are in ``revise/batch/``. ``config.py`` owns discovery, inheritance
and declaration-relative paths. ``sample.py`` owns read-only data checks and file
identity. ``runner.py`` owns task expansion, existing single-run invocation,
state and handoff. ``analysis.py`` owns the adapter interface and lifecycle.
``cli.py`` provides both installed commands. The batch layer does not duplicate
reconstruction algorithms or implement a second preprocessing pipeline.

A task is one hST/sST sample or one iST sample/type pair. The runner executes tasks
sequentially in child processes; failures are recorded independently. There is
no parallel/distributed scheduler or cross-type GA cache. Shared references are
read by path rather than materialized for every task.

.. _batch-result-identity:

Identity and recovery
---------------------

Successful reuse requires matching standard input identities, effective
configuration, configuration source chain, runtime code/dependency identity and
published artifact hashes. The source chain includes possible override locations,
so adding or deleting a YAML also invalidates affected results. Reconstruction
and analysis use the same resolver; analysis cannot accept a stale handoff merely
because the files that existed during the old run still match. The fingerprint
covers the complete effective configuration, including analysis settings. After
changing an adapter specification or its configuration source, rerun batch
reconstruction before batch analysis; an analysis-only configuration edit is
not exempt from reconstruction invalidation.

Each task saves its actual single-run YAML, log and state under ``.revise/``.
The handoff references standard files directly; it does not reference a
preparation manifest. Data and configuration are rechecked before results are
accepted. Disabled/removed tasks retain files but are marked inactive. A running
marker from an interrupted process is not a successful cache entry.

Output publication and framework-owned obsolete-file cleanup retain the
existing rollback behavior. A failed rerun may preserve older scientific files
while invalidating their status. Read current state rather than inferring
success from file existence. A shared output-root lock prevents concurrent
reconstruction and analysis from publishing into the same tree.

.. _batch-analysis-interface:

Connect an analysis module
--------------------------

Configure adapters at project, group or sample level, using the same inheritance
rules as reconstruction. The supported aspects are ``partition``,
``spatial_diversity`` and ``spatial_regions``. Without a selection, the default
is all three with ``not_implemented`` status. A nonempty effective mapping selects
its listed aspects; an empty override mapping does not clear inherited adapters.

.. code-block:: yaml

   analysis:
     partition:
       entrypoint: project_impact.partition:run
       version: "1"
       requires_pairing: true
       parameters: {}

``entrypoint`` is an importable ``module:function`` in the same Python environment.
``version`` identifies the adapter and dependencies; bump it when dependencies
outside the adapter's own source change. ``requires_pairing`` defaults to false.
Set it true when the algorithm needs verified raw/reconstructed observation
pairing. Coordinate units, carrier semantics and other scientific prerequisites
must also be checked by the adapter.

The existing ``AnalysisContext`` contains:

* ``reconstruction``: verified handoff dictionary with inputs, output roles,
  coordinate metadata and pairing semantics.
* ``output_dir``: temporary directory for this aspect's artifacts.
* ``parameters``: configured parameter mapping.

The function writes only beneath ``output_dir`` and returns a nonempty mapping
from artifact role to relative file path. This runnable interface example creates
an inventory, not a scientific partition analysis:

.. code-block:: python

   import csv
   from revise.batch.analysis import AnalysisContext

   def run(context: AnalysisContext) -> dict[str, str]:
       with (context.output_dir / 'carriers.csv').open('w', newline='') as stream:
           writer = csv.writer(stream)
           writer.writerow(['role', 'observations', 'genes', 'source_path'])
           for role, carrier in context.reconstruction['outputs'].items():
               writer.writerow([role, *carrier['shape'], carrier['path']])
       return {'carrier_inventory': 'carriers.csv'}

Adapters execute in-process as trusted project code. Exceptions fail that aspect
and allow others to continue; process termination requires restarting the batch.
No matrix is copied by the adapter runner. It validates returned files, publishes
them together and records hashes. Failed computation preserves prior artifacts
but marks the aspect failed. Nonempty directories not owned by the framework are
protected. Reuse checks the reconstruction, parameters, adapter source, version,
runtime identity and artifact hashes.

Logs and states are in ``.revise/analysis/<aspect>.log`` and ``<aspect>.json``.
``analysis/analysis.json`` records the current reconstruction fingerprint and
aspect status. The handoff's initial ``analysis.status: not_run`` is not a mutable
completion flag. Analysis status distinguishes succeeded, reused, failed,
blocked reconstruction, unavailable pairing and not_implemented. Reconstruction
completion does not establish analysis completion or biological validity.

.. _batch-verification-map:

Verification and development boundaries
---------------------------------------

.. list-table:: Code-to-test map
   :header-rows: 1
   :widths: 40 60

   * - Tests
     - Contract covered
   * - ``tests/batch/test_config.py``
     - Hierarchy, reference replacement, paths, discovery and invalid settings.
   * - ``tests/batch/test_sample.py``
     - Standard H5AD checks and unchanged/shared inputs.
   * - ``tests/batch/test_runner.py``
     - Task expansion, failure isolation, configuration invalidation and reuse.
   * - ``tests/batch/test_analysis.py``
     - Adapter prerequisites, configuration freshness, publication and recovery.
   * - ``tests/application/test_ist_publication.py``
     - iST paired/mean/random semantics and publication rollback.
   * - ``tests/integration/batch/test_real_sample_parity.py``
     - Bounded hST/iST/sST batch versus same-configuration single-run parity.

The framework implements input protocol checks, orchestration and analysis
handoff. It does not implement data conversion, label harmonization or the
reconstruction-impact algorithms developed on another branch. Real-sample call
parity is a technical regression check, not a full biological validation.

The one-time implementation record is
``docs/plans/2026-09-09-hierarchical-batch-inputs.md``. It is separate from this
maintained protocol and framework documentation.
