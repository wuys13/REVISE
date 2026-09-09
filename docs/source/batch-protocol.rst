.. _batch-input-protocol:

Batch input and output protocol
===============================

This page defines schema version 2. The :ref:`batch-reconstruction` page is the
usage entrypoint; :ref:`batch-framework` maps this protocol to code and tests.

Configuration and inheritance
-----------------------------

The project ``batch.yaml`` declares ``schema_version: 2``, ``input_root`` and
``output_root``. These keys are project-only. Roots must be separate and must
not contain one another. The input tree contains ST samples; shared SC files may
live beside that tree or elsewhere and are referenced by path.

.. code-block:: yaml

   schema_version: 2
   input_root: ST
   output_root: output
   modality: iST
   coordinates:
     unit: pixel
   inputs:
     reference:
       path: SC/shared_reference.h5ad

Settings inherit in order: project, input root, each intermediate directory,
then sample. Each optional override file is named ``batch.yaml``. Missing files
leave inherited settings unchanged. Directory names have no algorithmic meaning.

Mappings merge recursively; scalars and lists replace inherited values. An empty
mapping does not clear inherited fields. ``inputs.reference`` is the exception:
a lower-level reference replaces the entire block, including any filtering.
A replacement must supply its own path and any desired filter pair. This prevents
filters for one reference from leaking into another.

Every configured filesystem path resolves relative to the YAML declaring it,
before inheritance. Absolute paths are allowed. Adapter parameters are passed to
the adapter unchanged; they do not define additional framework filesystem fields.
For example, ``ST/CRC/batch.yaml`` can contain:

.. code-block:: yaml

   inputs:
     reference:
       path: ../../SC/crc_reference.h5ad
       filter_column: Patient
       filter_value: P2CRC

``ST/CRC/S02/batch.yaml`` can replace it with:

.. code-block:: yaml

   inputs:
     reference:
       path: ../../../SC/special_reference.h5ad

S02 then uses the special reference without the Patient filter. Other CRC
samples share the original reference file; no per-sample copy is created.

.. list-table:: Configuration fields
   :header-rows: 1
   :widths: 28 72

   * - Field
     - Meaning
   * - ``modality``
     - Required hST, iST or sST; maps to the existing single-run route.
   * - ``coordinates``
     - Required ``unit: um | pixel``; optional positive ``microns_per_coordinate``.
   * - ``enabled``
     - Boolean, defaults to true; lower levels may override it.
   * - ``inputs.reference``
     - H5AD path and optional ``filter_column`` / ``filter_value`` pair.
   * - ``inputs.pm_on_cell``
     - Optional existing single-run mapping input, with a path relative to its YAML.
   * - ``algorithm``, ``preprocessing``
     - Existing single-run solver and QC settings.
   * - ``global_anchoring``, ``local_refinement``
     - Existing GA/LR settings; batch adds iST ``cell_types``.
   * - ``output``
     - iST ``ist_mapping`` only; locations and filenames are batch-owned.
   * - ``execution``
     - Existing execution settings, including seed.
   * - ``analysis``
     - Per-aspect adapter specifications; see :ref:`batch-analysis-interface`.

Required protocol fields cannot be null. ``inputs.pm_on_cell`` is optional but
cannot be null when declared; it must specify a path, as in the single-run
interface. Other single-run fields retain their own null rules, such as null QC
thresholds. Missing reference, modality or units
is an error; none is inferred from platform names, directory names or labels.
The full example configuration is under ``configs/batch/``. Refer to the existing
application reference for scientific parameter meanings.

.. _batch-standard-h5ad:

Standard ST and reference files
-------------------------------

A sample is a directory containing ``spatial.h5ad`` inside ``input_root``. Its
identity is the relative directory path, such as ``CRC/S01``. Two different
categories may contain S01 without collision. A sample cannot contain another
sample. Discovery does not follow directory symlinks and does not scan SC or the
output root. No sample list, sample ID field or per-sample YAML is required.

Both H5AD files must have nonempty unique observation and gene identifiers and
finite nonnegative expression in ``X``. ST requires finite coordinates in
``obsm['spatial']``. Reference labels use configured GA/LR columns. Reference
filtering must select usable observations, and the inputs must share genes.
Choosing X does not prove that it contains counts; supply the representation
appropriate to the reconstruction protocol.

Checks are read-only. The runner does not write input H5AD, transform labels,
select layers, normalize coordinates or create preparation manifests. Prepare
nonstandard data separately. No platform converter or new preprocessing stage
is implemented in this release. Existing single-run QC, reference filtering and
gene alignment still run as configured.

Coordinate metadata describes existing values without scaling them. Pixel units
without calibration do not support physical-distance claims. Micrometre units
already imply ``microns_per_coordinate: 1``; another scale is invalid. Label aliases are
not inferred: ``Macro`` and ``Mono_Macro`` remain distinct unless already
harmonized before batch execution.

.. _batch-output-protocol:

Output roles and current-result semantics
-----------------------------------------

Output directories mirror sample paths relative to ``input_root``. hST uses
``sp-SVC``; iST uses ``sc-SVC`` cluster mode; sST uses ``sc-SVC`` sr mode.
For iST, ``local_refinement.cell_types`` defaults to ``[T, Macro, Fibroblast]``;
selected types become independent tasks in type subdirectories. Labels must be
safe, unique path components. Missing types fail explicitly; there is no
name-based exclusion.

.. code-block:: text

   output/
   ├── batch_status.json
   ├── analysis_status.json
   └── CRC/S01/
       ├── .revise/sample.json
       └── T/
           ├── spatial.h5ad
           ├── expr.h5ad
           ├── reconstruction.json
           ├── .revise/
           │   ├── application.yaml
           │   ├── reconstruction.log
           │   ├── task.json
           │   └── analysis/          # per-aspect state and logs
           └── analysis/
               ├── analysis.json
               ├── partition/
               ├── spatial_diversity/
               └── spatial_regions/

hST/sST place ``SVC.h5ad``, their handoff, controls and analysis at sample level.
Engine provenance remains available through references in the handoff.
iST ``output.ist_mapping`` selects assembly only:

* ``paired`` (default): spatial and expression H5AD carriers.
* ``mean``: sparse cluster-mean expression on spatial observations in ``SVC.h5ad``.
* ``random``: seeded within-cluster donor assignment in ``SVC.h5ad``, retaining
  donor identities and assembly provenance.

GA/LR are unchanged. Successful mode switches remove identified framework-owned
obsolete alternatives; failed publication restores prior files.

``reconstruction.json`` records task identity, mode, actual carrier roles and
hashes, standard input sources, configuration provenance, coordinate meanings,
pairing checks and engine provenance. Require a successful current handoff before
analysis. iST expression observations are reference donors, not spatial
observations. sST generated units are not paired one-to-one with raw spots.
Unavailable pairing is explicit, not silently approximated.

Analysis is organized by aspect. Paired values on the same observations/windows
belong in one table with named baselines and deltas; independent axes belong in
separate files. Large carriers are referenced, not copied. There is no mandatory
raw/reconstruct/comparison directory layer. Specific scientific table schemas
remain the responsibility of the analysis modules.
