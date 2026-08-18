Quick Start
===========

For new data, copy one Application template locally, edit its data-specific
fields, and run it. REVISE does not infer a route from a filename or convert
one mode to another.

1. Choose the template
----------------------

.. list-table::
   :header-rows: 1
   :widths: 3 3 2

   * - Spatial observations
     - Template
     - Public type
   * - Visium HD bins or pseudo-cells
     - ``VisiumHD.yaml``
     - ``sp-SVC``
   * - Segmented Xenium cells
     - ``Xenium.yaml``
     - ``sc-SVC``, ``cluster`` mode
   * - Multi-cell Visium spots
     - ``Visium.yaml``
     - ``sc-SVC``, ``sr`` mode

``sc-SVC`` requires ``application.mode``. The maintained Xenium request uses
``cluster`` mode and one selected broad cell type. The Visium request uses
``sr`` mode to reconstruct virtual cells from multi-cell spots. ``sp-SVC``
does not accept a mode.

For an installed package, use the :ref:`Application templates <application-templates>`
copy step to create the local YAML that you will edit. A source checkout keeps
the same three templates under ``configs/application/``.

2. Quick run
------------

Before every run, update ``inputs.st.path`` and ``inputs.reference.path``, the
optional reference ``filter_column``/``filter_value``,
``global_anchoring.broad_column``, cluster-mode
``local_refinement.subtype_column``, preprocessing thresholds, and
``output.dir``. The full schema and output rules are owned by
:doc:`application-reference`.

For Visium HD bins or pseudo-cells, set the ST/reference paths and the broad
annotation column in your local ``VisiumHD.yaml`` copy.

.. code-block:: bash

   revise-reconstruct --config VisiumHD.yaml

For Xenium segmented cells, configure your local ``Xenium.yaml`` copy and
select one concrete broad cell type. The final output directory is determined
automatically; see :doc:`application-reference` for that contract.

.. code-block:: bash

   revise-reconstruct --config Xenium.yaml --select-ct T
   revise-reconstruct --config Xenium.yaml --select-ct Fibroblast
   revise-reconstruct --config Xenium.yaml --select-ct "Mono/Macro"

``--select-ct`` has priority over the YAML's
``local_refinement.select_cell_type``. It is valid only for ``sc-SVC`` cluster
mode and must name one non-empty concrete broad cell type; values such as
``all`` and wildcards are rejected.

For multi-cell Visium spots, configure your local ``Visium.yaml`` copy; an
exact PM prior is optional for this route.

.. code-block:: bash

   revise-reconstruct --config Visium.yaml

From a source checkout, the equivalent commands use the maintained repository
templates:

.. code-block:: bash

   python reconstruct.py --config configs/application/VisiumHD.yaml
   python reconstruct.py --config configs/application/Xenium.yaml --select-ct T
   python reconstruct.py --config configs/application/Xenium.yaml --select-ct Fibroblast
   python reconstruct.py --config configs/application/Xenium.yaml --select-ct "Mono/Macro"
   python reconstruct.py --config configs/application/Visium.yaml

``paths.root_dir: .`` means the command launch directory, not the YAML
directory. Input and output paths are relative children of that root.

3. Inspect the promised artifacts
---------------------------------

``sp-SVC`` and ``sc-SVC`` sr mode publish one
``<output-dir>/<output-name>.h5ad`` (or ``svc.h5ad`` without a name). ``sc-SVC``
cluster mode publishes and returns a fixed pair in the selected cell type's
final directory:

.. code-block:: text

   <output-dir>/<output-name>_spatial.h5ad
   <output-dir>/<output-name>_expr.h5ad

Without an output name, the pair is ``spatial.h5ad`` and ``expr.h5ad``. Each
public H5AD records its route/mode and points to the run's
``provenance.json``. A directory alone is not success evidence: check the
expected artifact(s) and a succeeded manifest.

Next steps
----------

- :doc:`concepts` explains what the three data shapes and two sc-SVC modes
  mean scientifically.
- :doc:`gallery` lists the preserved benchmark and Application notebook
  snapshots; it is not a replacement for this current reconstruction entry.
