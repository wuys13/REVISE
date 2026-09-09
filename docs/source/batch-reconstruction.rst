Batch reconstruction and analysis handoff
========================================

A batch contains spatial samples. Each sample has one ``sample.yaml``; optional
parent directories (cancer, study, or another classification) do not change its
meaning. Samples explicitly select hST, iST, or sST and a reference. REVISE does
not select references, infer platform routes, or exclude cell types by name.

Prepare and run
---------------

Copy a template from ``configs/batch/sample_iST.yaml``, ``sample_hST.yaml``, or
``sample_sST.yaml`` into each sample directory as ``sample.yaml``. Edit its
paths, expression source, coordinate metadata, reference labels, and protocol
parameters. Template QC values are examples, not universal batch defaults.

Copy ``configs/batch/batch.yaml`` and point ``data_root`` to the collection.
All sample paths resolve relative to their ``sample.yaml``; ``data_root``
resolves relative to the batch YAML. Absolute input paths are also accepted.
The command's launch directory does not affect these paths.

.. code-block:: bash

   # Optional: inspect input preparation separately.
   revise-prepare-sample --config /data/CRC/sample_001/sample.yaml

   # Preparation is also performed automatically by the batch runner.
   revise-batch-reconstruct --config /data/batch.yaml

   # Equivalent source-checkout batch entrypoint.
   python /path/to/REVISE/batch_reconstruct.py --config /data/batch.yaml

The installed entrypoints require installing this checkout, for example with
``python -m pip install -e .`` in the project's Python environment. The existing
single-run ``reconstruct.py`` and ``revise-reconstruct`` remain available.

Input contract
--------------

``inputs.st`` accepts explicit ``h5ad`` or ``spatialdata``; references use H5AD.
For SpatialData, also supply ``inputs.st.spatialdata.table`` and ``element``.
These use the existing REVISE SpatialData reader. This interface does not parse
platform-native matrix/coordinate downloads.

Preparation requires explicit ``matrix: X`` or ``matrix: layers/<name>`` for
both inputs, and ``spatial_key`` plus ``coordinate_unit: um | pixel`` for ST.
It checks nonempty unique observation/gene IDs, finite nonnegative expression,
finite coordinates, shared genes, and the required reference label columns.
Choosing a matrix does not establish that it contains counts; supply the
matrix appropriate to the chosen reconstruction protocol. No normalization,
species inference, or gene-name conversion occurs here.

Pixel coordinates can declare ``microns_per_coordinate``. Without calibration,
physical-distance availability remains false in the handoff. Coordinates are
copied into ``obsm['spatial']`` without numerical scaling. Input sources remain
unchanged; standardized spatial data are written to ``inputs/spatial.h5ad``.
An unchanged reference is shared by path. A changed reference matrix or explicit
label mapping produces ``inputs/reference.h5ad`` for this sample.

For example, only when the source labels have this meaning:

.. code-block:: yaml

   preparation:
     spatial:
       matrix: layers/counts
       spatial_key: spatial
       coordinate_unit: pixel
       microns_per_coordinate: 0.2125
     reference:
       matrix: X
       label_mapping:
         Level1:
           Mono_Macro: Macro

The original reference labels remain in ``Level1_original``. Without a mapping,
``Macro`` and ``Mono_Macro`` are distinct. Reference filtering is specified by
``inputs.reference.filter_column`` and ``filter_value`` and performed by the
single-run application. The preparation manifest records sources, transforms,
calibration, and output hashes. Keep the source references valid for subsequent
validation; do not point input paths at files the framework will overwrite.

Output modes and layout
-----------------------

hST maps to ``sp-SVC``; iST maps to ``sc-SVC`` cluster mode; sST maps to
``sc-SVC`` sr mode. Batch labels do not rename these existing single-run APIs.

For iST, ``local_refinement.cell_types`` defaults to
``[T, Macro, Fibroblast]``. Each selected type is an independent task. Missing
types fail explicitly while the other tasks continue. Directory labels must
be unique safe path components; use explicit label mapping where needed.

``output.ist_mapping`` is available in both batch and single-run iST YAML:

* ``paired`` (default): publish ``spatial.h5ad`` and ``expr.h5ad``.
* ``mean``: assign cluster-mean reference expression to spatial observations
  and publish ``SVC.h5ad``.
* ``random``: sample a reference donor within each cluster, using the configured
  seed, and publish ``SVC.h5ad`` with donor IDs and assembly provenance.

The key controls assembly/publication only. GA and LR remain unchanged.
The batch output names are fixed; single-run custom output names remain
supported. A successful mode change removes only identified framework-owned
alternatives. Failed publication restores previous files.

.. code-block:: text

   data/CRC/sample_001/
     sample.yaml
     inputs/
       spatial.h5ad
       preparation.json
     T/
       spatial.h5ad
       expr.h5ad
       reconstruction.json
       .revise/
         application.yaml
         reconstruction.log
         task.json
       analysis/
         partition/
         spatial_diversity/
         spatial_regions/
     Macro/
     Fibroblast/

For hST and sST, ``SVC.h5ad``, ``reconstruction.json``, ``.revise/``, and
``analysis/`` are at sample level. Engine provenance directories also remain
available and are referenced by the handoff.

Resume and inspect failures
---------------------------

Run the same command again to resume. Only successful tasks whose input,
preparation, effective configuration, runtime code/dependency identity and
recorded artifacts still match are reused. Missing or changed results, failed
tasks and interrupted tasks rerun. Execution is sequential in separate child
processes, so a solver process failure is recorded independently. There is no
cross-type GA cache or distributed scheduler. Concurrent batches using the same
data root are rejected.

``data_root/batch_status.json`` records task outcomes and succeeded/failed/reused
counts. Exit status is 0 for complete success, 1 for task failures, and 2 for
an invalid batch invocation. Each task records its resolved single-run YAML,
log and state in ``.revise/``. Preparation failures are recorded at sample level.
An interrupted batch leaves a running marker; the next invocation retries it.
Removed cell types are marked inactive without deleting their scientific data.

Before reading a result, require ``reconstruction.json.status == 'succeeded'``.
Old files can remain after a failed rerun and must not be treated as current
merely because a filename exists. Consult the current batch report for sample
preparation errors and the task log for solver failures.

Analysis handoff
----------------

This release prepares the handoff; it does **not** run reconstruction-impact
analysis. Empty aspect directories do not mean analysis has completed.
``reconstruction.json`` records ``analysis.status: not_run``, actual output
roles, input identities, calibration, observation pairing checks, and engine
provenance. The reconstruction-impact branch continues to own its metrics,
thresholds, scientific interpretation and notebook execution.

Organize outputs by analysis aspect, not by three repeated raw/reconstruct/
comparison directory trees. Store values measured on the same observations or
windows in one paired table, with explicitly named baselines and deltas. Store
independent axes (such as reference-expression donors) in separate files. Keep
large H5AD carriers once and reference them. A notebook may read several aspects
and present raw, reconstructed and comparison panels together.

Suggested future artifacts include ``partition/unit_assignments.csv.gz``,
``spatial_diversity/window_metrics.csv.gz`` and
``spatial_regions/region_summary.csv``. They are conventions for handoff, not
new files generated by this release. Raw Leiden and Raw Level2 are distinct
baselines and should remain explicitly named in table columns.

.. code-block:: python

   import json
   from pathlib import Path
   import anndata as ad

   folder = Path('/data/CRC/sample_001/T')
   record = json.loads((folder / 'reconstruction.json').read_text())
   assert record['status'] == 'succeeded'
   assert record['ist_mapping'] == 'paired'
   assert record['pairing']['status'] == 'available'
   spatial = ad.read_h5ad(record['outputs']['spatial']['path'])
   # Expression-side observations are not spatial observations.
   expression_path = record['outputs']['expression']['path']

The iST expression carrier does not participate in spatial observation pairing.
Mean/random assembled expression has different semantics from the existing
paired spatial carrier and requires an appropriate downstream analysis contract.
sST generated units are not declared one-to-one paired with raw spots. A valid
file or successful reconstruction does not establish biological validation.
