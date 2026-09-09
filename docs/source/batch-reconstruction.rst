.. _batch-reconstruction:

Batch reconstruction and analysis
=================================

Place each standard ST sample in a directory containing ``spatial.h5ad``.
Configure shared settings once in the project's ``batch.yaml`` and add another
``batch.yaml`` only where a group or sample needs overrides. References are
shared by path. Reconstruction and analysis write to a separate output tree.

.. toctree::
   :maxdepth: 2

   batch-protocol
   batch-framework

.. _batch-getting-started:

Set up and run
--------------

Copy the example tree under ``configs/batch/`` to your project. Supply real
``ST/CRC/S01/spatial.h5ad``, ``ST/CRC/S02/spatial.h5ad`` and reference files in
``SC/``. The template directories alone are not samples: discovery requires
``spatial.h5ad``. Remove the S02 override if it should share the CRC reference.
Review modality, coordinate units, reference label columns and protocol settings;
template QC thresholds are examples, not universal defaults.

.. code-block:: text

   project/
   ├── batch.yaml
   ├── ST/
   │   └── CRC/
   │       ├── batch.yaml
   │       ├── S01/spatial.h5ad
   │       └── S02/
   │           ├── spatial.h5ad
   │           └── batch.yaml       # optional override
   ├── SC/
   │   ├── crc_reference.h5ad
   │   └── special_reference.h5ad
   └── output/

Install this checkout in the reconstruction environment with
``python -m pip install -e .``, then run:

.. code-block:: bash

   revise-batch-reconstruct --config /data/project/batch.yaml
   revise-batch-analyze --config /data/project/batch.yaml

   # Equivalent source-checkout entrypoints:
   python /path/to/REVISE/batch_reconstruct.py --config /data/project/batch.yaml
   python /path/to/REVISE/batch_analyze.py --config /data/project/batch.yaml

The launch directory does not affect configured paths. See
:ref:`batch-input-protocol` for fields and inheritance, and
:ref:`batch-framework` for execution and adapter integration.

.. _batch-recovery:

Resume and inspect
------------------

Rerun the same command to resume. Successful tasks are reused only when their
inputs, effective configuration, configuration sources, runtime identity and
artifacts still match. Failed or interrupted tasks rerun. A changed, newly added
or removed configuration file invalidates affected results, including changes
limited to analysis settings. Rerun reconstruction after such an edit before
invoking analysis. Disabled samples and
removed iST types become inactive without deleting their scientific files.

Consult ``output/batch_status.json`` and each task's ``.revise/reconstruction.log``.
A failed rerun can leave older scientific files in place; require a current
successful ``reconstruction.json`` instead of checking filenames alone.
``output/analysis_status.json`` and ``analysis/analysis.json`` report analysis
separately. An empty analysis directory never establishes analysis completion.

Reconstruction exits 0 on success, 1 for task failures and 2 for invalid invocation.
Analysis exits 0 only when requested tasks succeeded or were reused; missing
adapters and unmet prerequisites are incomplete outcomes. Execution is sequential;
one task failure does not stop subsequent tasks. Concurrent runs against the same
output root are rejected.

.. _batch-migration:

Migrate the experimental input protocol
---------------------------------------

Schema version 2 replaces schema version 1, per-sample ``sample.yaml`` discovery
and ``revise-prepare-sample``. Old batch configurations produce a migration error.
Single-run reconstruction configuration remains supported.

1. Choose separate input and output roots and set ``schema_version: 2`` in the
   project ``batch.yaml``.
2. Put each already-standard ST file at ``<input_root>/<sample>/spatial.h5ad``.
   Select existing standardized files deliberately; this release does not move
   or transform originals.
3. Move common reconstruction settings into the project or group ``batch.yaml``.
   Move ``sample.modality`` to ``modality`` and coordinate metadata to
   ``coordinates``. Reference paths can point to one shared file.
4. Keep sample overrides only where needed. Remove ``inputs.st`` and
   ``preparation`` settings from the new configuration. Convert nonstandard data
   before batch execution; the new entrypoint does not perform those conversions.
5. Run the new configuration. Old cached successes are not accepted directly as
   schema-2 successes. Existing data, preparation files and results are never
   automatically relocated or deleted.

The implemented scope is reconstruction orchestration and analysis adapters.
A separate data processing stage is planned but not implemented here.
Reconstruction-impact algorithms remain in their development branch; an aspect
without an adapter reports ``not_implemented``.
