Quick Start
===========

REVISE has separate entry points for benchmark reproduction and application
reconstruction. Paper notebooks require a source checkout and external data.

Benchmark
---------

From the repository root, ``python reproduce/benchmark_main.py`` runs one Benchmark route YAML;
``bash reproduce/benchmark_main.sh`` launches the bounded multi-family workflow.
Benchmark notebooks are under
``reproduce/benchmark/``.

Choose and copy a template
--------------------------

Source-checkout users select one maintained file:

.. list-table::
   :header-rows: 1

   * - ST data / sample
     - Template
     - Declared inputs
   * - Xenium T / ``P2CRC``
     - ``configs/application/Xenium_T.yaml``
     - ``raw_data/Real_application/P2CRC_Xenium.h5ad`` and
       ``raw_data/Real_application/adata_sc_all_reanno.h5ad``
   * - Xenium fibroblast / ``P2CRC``
     - ``configs/application/Xenium_Fib.yaml``
     - ``raw_data/Real_application/P2CRC_Xenium.h5ad`` and
       ``raw_data/Real_application/adata_sc_all_reanno.h5ad``
   * - Xenium mono/macro / ``P2CRC``
     - ``configs/application/Xenium_Mono.yaml``
     - ``raw_data/Real_application/P2CRC_Xenium.h5ad`` and
       ``raw_data/Real_application/adata_sc_all_reanno.h5ad``
   * - Visium HD / ``P1CRC``
     - ``configs/application/VisiumHD.yaml``
     - ``raw_data/Real_application/P1CRC_HD.h5ad`` and
       ``raw_data/Real_application/adata_sc_all_reanno.h5ad``
   * - Visium / ``REVISEVisiumMouseBrain``
     - ``configs/application/Visium.yaml``
     - ``raw_data/visium_mouse_brain/ST_mouse_brain_prepared.h5ad`` and
       ``raw_data/visium_mouse_brain/scRNA_mouse_brain_prepared.h5ad``

Installed-package users can copy the matching file from the packaged
``revise.application/templates`` resource into their project, then edit it.
Alternatively, a bare official name first checks
``configs/application/<name>.yaml`` in the launch directory and then uses the
packaged resource. This is the standard ``importlib.resources`` interface, not
another REVISE API:

.. code-block:: python

   from importlib.resources import as_file, files
   from shutil import copyfile

   resource = files("revise.application").joinpath("templates", "VisiumHD.yaml")
   with as_file(resource) as source:
       copyfile(source, "VisiumHD.yaml")

The real P1CRC/P2CRC and mouse-brain H5AD files are not distributed in the package.
Replace paths or stage those files beneath the chosen run root before use.

Application YAML
----------------

The complete schema contains ``application``, ``paths``, ``algorithm``,
``inputs``, ``preprocessing``, ``global_anchoring``, ``local_refinement``, ``output``, and
``execution`` (plus ``schema_version``). A complete Visium HD request is:

The two exact input keys are ``inputs.st.path`` and
``inputs.reference.path``.

.. code-block:: yaml

   schema_version: 1

   application:
     svc_type: sp-SVC

   paths:
     root_dir: .

   algorithm:
     ot_method: pot

   inputs:
     st:
       path: data/sample_st.h5ad
       format: h5ad
     reference:
       path: data/sc_ref.h5ad
       format: h5ad

   preprocessing:
     spatial:
       min_transcript_counts: 60
       min_cell_counts: 100
     reference:
       min_transcript_counts: null
       min_cell_counts: 100

   global_anchoring:
     broad_column: Level1

   local_refinement:
     strength: 0.2

   output:
     dir: output
     name: sample_sp-SVC

   execution:
     seed: 42

``paths.root_dir: .`` means the launch current working directory, not the
application YAML directory. The other accepted form is an existing absolute
directory. Input and output paths must be relative children of
``paths.root_dir`` and cannot escape it with ``..``.

Application inputs are always exact direct paths. ST and reference resolve as
``<root-dir>/<st-path>`` and ``<root-dir>/<reference-path>``. An optional
sc-SVC-sr prior is also an exact ``inputs.pm_on_cell.path``; when omitted, no
sidecar is searched and seeded random quota-slot assignment is used.

The full engine configuration in ``revise.config.authority`` is package-owned
and authoritative. ``algorithm`` contains only optional ``ot_method``.
``algorithm.ot_method`` controls both GA and LR; omitting it keeps the selected
engine profile authoritative.

Template-specific fields are intentionally small:

- ``Xenium_T.yaml``, ``Xenium_Fib.yaml``, and ``Xenium_Mono.yaml`` select
  ``sc-SVC`` and require ``local_refinement.subtype_column``. They carry the
  fixed selections ``T``, ``Fibroblast``, and ``Mono_Macro``; each selected
  label must exist in the full reference.
- Visium HD selects ``sp-SVC`` with ``local_refinement.strength: 0.2``.
- Visium selects ``sc-SVC-sr`` with ``local_refinement.strength: 0.0``. Zero
  does not skip LR; it only disables assignment-posterior strengthening of the
  LR cost.

Both H5AD inputs require non-empty ``X``, unique ``obs_names`` and unique
``var_names``, and at least one shared gene. ST requires finite two-dimensional
coordinates in ``obsm["spatial"]``. Every route requires the configured broad
annotation in reference ``obs``. Only standard sc-SVC requires the configured
subtype annotation. sc-SVC-sr composition and expression allocation use the
broad assignment and do not require a subtype column. Application does not
filter the reference by patient; prepare the reference before running if it
contains multiple cohorts.

Run
---

From a source checkout:

.. code-block:: bash

   python reconstruct.py --config configs/application/VisiumHD.yaml

The installed command is equivalent and also accepts the official external
template path:

.. code-block:: bash

   revise-reconstruct --config configs/application/VisiumHD.yaml

The Application command runs the YAML request and publishes its result H5AD.

Standard sc-SVC selects TACCO by profile. Install
``python -m pip install "revise-svc[tacco]"`` for a published package or
``python -m pip install ".[tacco]"`` from source. A missing or incompatible
TACCO fails with installation guidance; REVISE never switches algorithms
automatically and does not fall back to POT.

Application output and evidence
-------------------------------

``sp-SVC`` and ``sc-SVC-sr`` publish:

.. code-block:: text

   <output-dir>/<output-name>.h5ad

Standard ``sc-SVC`` publishes two carriers:

.. code-block:: text

   <output-dir>/<output-name>_spatial.h5ad
   <output-dir>/<output-name>_expr.h5ad

Each file links to the canonical run's ``provenance.json``. Application request
identity is namespaced under ``application_config`` as ``source_path``,
``source_sha256``, ``declared_root``, ``resolved_root``, ``cwd``,
``resolved_inputs``, ``output_paths``, ``effective_request``,
``effective_request_hash``. The top-level engine
configuration identity remains separate in ``engine_defaults_hash``,
``authority_hash``, ``algorithm_config_hash``, and ``effective_config_hash``. There is no solver-event
telemetry.

Paper reproduction notebooks
----------------------------

Curated application notebooks are under ``reproduce/case/``. Their presence
does not prove that the current package reran the real datasets.

The canonical application-analysis notebooks are
``Xenium_sc_SVC_T.ipynb``, ``Xenium_sc_SVC_Fibroblast.ipynb``,
``Xenium_sc_SVC_Monocyte.ipynb``, and ``VisiumHD_sp_SVC.ipynb``. The Visium
mouse-brain sc-SVC-sr notebook is preserved separately.

Application utilities
---------------------

Optional morphology priors and biology-facing metrics remain available through
``revise-build-histology-priors`` and ``revise-compute-biological-metrics``.
