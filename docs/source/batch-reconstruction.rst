Batch reconstruction and analysis
=================================

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

Copy ``configs/batch/batch.yaml`` and set ``input_root`` to the input collection
and ``output_root`` to a separate result tree. The roots must not overlap.
All sample paths resolve relative to their ``sample.yaml``; both batch roots
resolve relative to the batch YAML. Absolute input paths are also accepted.
The command's launch directory does not affect these paths.

.. code-block:: bash

   # Optional: inspect input preparation separately.
   revise-prepare-sample --config /data/input/CRC/sample_001/sample.yaml

   # Preparation is also performed automatically by the batch runner.
   revise-batch-reconstruct --config /data/batch.yaml

   # Run connected analysis modules; missing adapters are reported explicitly.
   revise-batch-analyze --config /data/batch.yaml

   # Equivalent source-checkout batch entrypoint.
   python /path/to/REVISE/batch_reconstruct.py --config /data/batch.yaml
   python /path/to/REVISE/batch_analyze.py --config /data/batch.yaml

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

   data/
     batch.yaml
     input/                           # input_root
       references/
         shared_reference.h5ad        # optional shared reference
       CRC/sample_001/
         sample.yaml
         source.h5ad                  # optional local original; external paths also work
         inputs/
           spatial.h5ad               # normalized ST
           reference.h5ad             # only if reference transformation is needed
           preparation.json
     output/                          # output_root
       batch_status.json
       analysis_status.json           # written by the analysis runner
       CRC/sample_001/
         .revise/sample.json
         T/
           spatial.h5ad
           expr.h5ad
           reconstruction.json
           .revise/
             application.yaml
             reconstruction.log
             task.json
             analysis/                # per-aspect states and logs
           analysis/
             analysis.json            # current analysis status and reconstruction identity
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
output root are rejected, including concurrent reconstruction and analysis.

``output_root/batch_status.json`` records task outcomes and succeeded/failed/reused
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

This release supplies an executable analysis adapter runner. Reconstruction-impact
algorithms are **not bundled**; empty aspect directories do not mean analysis
has completed.
``reconstruction.json`` records the initial ``analysis.status: not_run``, actual output
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

   folder = Path('/data/output/CRC/sample_001/T')
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

Separate input and output responsibilities
------------------------------------------

Preparation writes normalized inputs and ``preparation.json`` inside the input
sample. Original files are retained. Once prepared, an unchanged input package
is reused without rewriting it. Reconstruction never places H5AD results, engine
provenance, task logs or task state in the input tree. Result paths mirror the
sample directory relative to ``input_root``; category names have no algorithmic
meaning. Shared reference paths may be absolute or relative to ``sample.yaml``.
For the illustrated tree, ``../../references/shared_reference.h5ad`` selects the
shared reference from ``CRC/sample_001/sample.yaml``.

The earlier experimental ``data_root`` batch key is replaced by explicit
``input_root`` and ``output_root``; old batch YAML is rejected with a configuration
error. Set the new roots and rerun into an empty output tree. Existing combined
folders are not automatically moved or deleted. Single-run YAML is unchanged.

Connect and run analysis modules
-------------------------------

``revise.batch.analysis.run_analysis_batch`` and ``revise-batch-analyze`` use the
same batch YAML as reconstruction. They consume its latest completed
``batch_status.json`` inventory, verify reconstruction artifacts and input hashes,
and run each requested aspect sequentially. They do not discover unrelated old
results by filename. The sample YAML must still match the reconstruction.

.. code-block:: yaml

   schema_version: 1
   input_root: ./input
   output_root: ./output
   analysis:
     partition:
       entrypoint: project_impact.partition:run
       version: "1"
       requires_pairing: true
       parameters: {}
     spatial_diversity: {}
     spatial_regions: {}

An aspect without an entrypoint reports ``not_implemented``. With no analysis
selection (or an empty mapping), all three aspects are reported this way. A
nonempty mapping selects only its listed aspects. ``requires_pairing`` defaults
to false; set it true for algorithms requiring the verified raw/reconstructed
observation pairing. Other scientific prerequisites, including coordinate units,
carrier semantics and applicable baselines, remain the adapter's responsibility.

An adapter is an importable Python function accepting one ``AnalysisContext``:

* ``reconstruction``: verified handoff dictionary, including input references,
  output roles, units and pairing semantics.
* ``output_dir``: temporary directory for this aspect's artifacts.
* ``parameters``: its explicit YAML parameter mapping.

Return a nonempty dictionary mapping artifact roles to relative file paths.
Write all outputs beneath ``output_dir`` and do not modify input or reconstruction
files. Functions run in the current Python process; exceptions are recorded and
other aspects continue, while process termination requires restarting the batch.
Adapters are trusted project code, not sandboxed plugins. Install their package
in the same environment or make it importable with ``PYTHONPATH``.

The following minimal adapter demonstrates file handoff using an inventory of
existing carriers. It performs no partition or impact analysis. Save it as
``my_analysis.py`` and use ``entrypoint: my_analysis:run`` for an interface check:

.. code-block:: python

   import csv
   from revise.batch.analysis import AnalysisContext

   def run(context: AnalysisContext) -> dict[str, str]:
       path = context.output_dir / 'carriers.csv'
       with path.open('w', newline='') as stream:
           writer = csv.writer(stream)
           writer.writerow(['role', 'observations', 'genes', 'source_path'])
           for role, carrier in context.reconstruction['outputs'].items():
               writer.writerow([role, *carrier['shape'], carrier['path']])
       return {'carrier_inventory': 'carriers.csv'}

No large matrix is copied by the runner. Each adapter's logs and state live in
``.revise/analysis/<aspect>.log`` and ``<aspect>.json``. Successful artifacts are
published together; failed computation leaves prior artifacts in place but marks
the aspect failed. Only previously published framework directories can be
replaced; unowned, nonempty aspect directories are protected. Reuse requires matching
reconstruction, parameters, adapter source, declared version, runtime identity,
and artifact hashes. Bump ``version`` when code dependencies outside the adapter
source file change. Interrupted aspects rerun.

Read ``analysis/analysis.json`` for current status, and require its
``reconstruction_fingerprint`` to match the current successful handoff. The
handoff's initial ``analysis.status`` is not a mutable analysis completion flag.
Starting reconstruction again resets the analysis summary to ``not_run``.
The batch-level ``analysis_status.json`` distinguishes ``succeeded``, ``reused``,
``failed``, ``blocked`` (invalid reconstruction), ``unavailable`` (required pairing
absent), and ``not_implemented``. Exit status is 0 only when every requested task
succeeded or was reused, 1 for incomplete tasks, and 2 for invalid invocation.
