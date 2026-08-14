Quick Start
===========

For new data, REVISE has one direct Application workflow: choose a template,
edit its exact input paths and annotations, then run the YAML. Benchmark
reproduction notebooks are listed directly under Sim2Real Benchmark in the
documentation navigation.

1. Choose the reconstruction type
---------------------------------

.. list-table::
   :header-rows: 1
   :widths: 3 2 3

   * - Your spatial observations
     - Type
     - Start from
   * - Visium HD bins or pseudo-cells
     - ``sp-SVC``
     - ``configs/application/VisiumHD.yaml``
   * - segmented Xenium cells, T-cell analysis
     - ``sc-SVC``
     - ``configs/application/Xenium_T.yaml``
   * - segmented Xenium cells, fibroblast/CAF analysis
     - ``sc-SVC``
     - ``configs/application/Xenium_Fib.yaml``
   * - segmented Xenium cells, mono/macro analysis
     - ``sc-SVC``
     - ``configs/application/Xenium_Mono.yaml``
   * - Visium multi-cell spots
     - ``sc-SVC-sr``
     - ``configs/application/Visium.yaml``

REVISE does not infer this choice from a filename or repair a mismatched mode.
Choose the type that matches the rows in the ST object.

2. Install the required solver
------------------------------

The base install includes POT. Standard ``sc-SVC`` selects TACCO 0.5.0 in its
maintained Xenium templates, so install the matching extra when using them:

.. code-block:: bash

   python -m pip install .
   python -m pip install ".[tacco]"

For a published package, use ``revise-svc`` and ``revise-svc[tacco]`` instead.
``algorithm.ot_method`` controls both GA and LR. A missing or failed selected
solver stops the run; REVISE never switches algorithms automatically and does
not fall back to POT.

3. Prepare the two inputs
-------------------------

Every Application request names exact direct paths with
``inputs.st.path`` and ``inputs.reference.path``. Both objects require
non-empty expression, unique ``obs_names`` and unique ``var_names``, and at
least one shared gene. ST requires finite two-dimensional coordinates in
``obsm["spatial"]``. Every route requires the configured broad annotation in
reference ``obs``. Only standard sc-SVC requires the configured subtype
annotation.

``sc-SVC-sr`` composition and expression allocation use the broad assignment
and do not require a subtype column. Its optional probability prior is one
exact ``inputs.pm_on_cell.path``; when omitted, no sidecar is searched and the
seeded quota assignment described in :doc:`concepts` is used.

The package does not install the paper datasets. A convenient project-local
layout is ``raw_data/Real_application/`` for the Xenium and Visium HD cases
and ``raw_data/visium_mouse_brain/`` for the Visium mouse-brain case. You may
use another external data root by changing ``paths.root_dir`` and the two input
paths. Research data are inputs only: do not add ``raw_data/`` to Git.

4. Edit one YAML
----------------

Copy a maintained source template, or copy the matching installed resource
from ``revise.application.templates`` with the standard
``importlib.resources`` interface:

.. code-block:: python

   from importlib.resources import as_file, files
   from shutil import copyfile

   resource = files("revise.application").joinpath("templates", "VisiumHD.yaml")
   with as_file(resource) as source:
       copyfile(source, "VisiumHD.yaml")

A minimal ``sp-SVC`` request is:

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
       min_transcript_counts: null
       min_counts: 20
       min_cell_counts: 30
     reference:
       min_transcript_counts: null
       min_genes: 20
       min_cell_counts: 50
   global_anchoring:
     broad_column: Level1
   local_refinement:
     strength: 0.2
   output:
     dir: output
     name: sample_sp-SVC
   execution:
     seed: 42

``paths.root_dir: .`` means the launch current working directory, not the YAML
directory. The other accepted form is an existing absolute directory. Input
and output paths must be relative children of ``paths.root_dir`` and cannot
escape it with ``..``.

REVISE applies no implicit patient filter. A filter runs only when the YAML
declares both ``inputs.reference.filter_column`` and ``filter_value``; the
maintained Xenium YAML files explicitly select ``Patient == P2CRC`` before
spatial and reference preprocessing.

5. Run
------

From a source checkout:

.. code-block:: bash

   python reconstruct.py --config configs/application/VisiumHD.yaml

The installed command accepts the same YAML:

.. code-block:: bash

   revise-reconstruct --config configs/application/VisiumHD.yaml

The source entry point deliberately keeps the processing sequence visible:

.. code-block:: text

   load YAML and AnnData
       -> filter the reference when the YAML declares a filter
       -> preprocess spatial and reference inputs
       -> prepare the sc-SVC pair or normalize reference labels
       -> run unified reconstruction
       -> publish the returned AnnData object(s)

There is no second public reconstruction command and no required low-level
kernel setup.

6. Inspect the result and manifest
----------------------------------

``sp-SVC`` and ``sc-SVC-sr`` return one ``AnnData`` and publish:

.. code-block:: text

   <output-dir>/<output-name>.h5ad

If ``output.name`` is omitted, the filename is ``svc.h5ad``. Standard
``sc-SVC`` returns ``(spatial_adata, expression_adata)`` and publishes:

.. code-block:: text

   <output-dir>/<output-name>_spatial.h5ad
   <output-dir>/<output-name>_expr.h5ad

Without a name, those filenames are ``spatial.h5ad`` and ``expr.h5ad``. Each
H5AD contains ``uns["revise_reconstruction"]`` with its output role and a link
to the canonical ``provenance.json``. A successful result requires the expected
public artifact(s) and a succeeded manifest; directory existence alone is not
success evidence.

Python callers use the same YAML:

.. code-block:: python

   from reconstruct import run_application
   result = run_application("configs/application/VisiumHD.yaml")

Next steps
----------

- Sim2Real Benchmark lists its six reproduction notebooks directly.
- The Gallery preserves static downstream analysis notebooks. They are not the
  current reconstruction interface.
- Optional morphology priors and biology-facing metrics are exposed through
  ``revise-build-histology-priors`` and
  ``revise-compute-biological-metrics``.
