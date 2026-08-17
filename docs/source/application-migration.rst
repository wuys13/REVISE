Application migration
=====================

This release is a hard cut for the Application interface. It does not read old
Application YAML, aliases, internal imports, notebook links, or output names.
Benchmark routes and their historical internal SR task names are outside this
migration.

Public request changes
----------------------

.. list-table::
   :header-rows: 1
   :widths: 4 5

   * - Previous form
     - Replace with
   * - ``application.svc_type: sc-SVC`` without a mode
     - Add ``application.mode: cluster`` for segmented cells, or ``sr`` for multi-cell spots.
   * - ``application.svc_type: sc-SVC-sr``
     - ``application.svc_type: sc-SVC`` and ``application.mode: sr``.
   * - ``Xenium_T.yaml``, ``Xenium_Fib.yaml``, or ``Xenium_Mono.yaml``
     - ``Xenium.yaml`` plus ``--select-ct T``, ``Fibroblast``, or ``Mono/Macro`` respectively. Use the cluster output-directory contract in :doc:`application-reference`.
   * - ``REVISEPipeline.run(svc_type="sc-SVC-sr")``
     - ``REVISEPipeline.run(svc_type="sc-SVC", application_mode="sr")``.
   * - ``REVISEPipeline.run(svc_type="sc-SVC")`` for segmented cells
     - ``REVISEPipeline.run(svc_type="sc-SVC", application_mode="cluster")``.

``--select-ct`` only applies to cluster mode. It is not a replacement for
checking the selected cell type and annotation columns in a new reference.

Visium and Gallery names
------------------------

The Visium application notebook moved from
``reproduce/case/sc_SVC_sr_case_Visium_mouse_brain.ipynb`` to
``reproduce/case/Visium_sc_SVC_mouse_brain.ipynb``. Its ReadTheDocs bridge has
the corresponding new name. The public result name changed from
``REVISEVisiumMouseBrain_sc-SVC-sr.h5ad`` to
``REVISEVisiumMouseBrain_sc-SVC.h5ad``. Existing local output files are not
renamed automatically; old and new outputs can coexist.

Application internal imports
----------------------------

Application-only code moved to names that describe SR as an sc-SVC mode:

.. list-table::
   :header-rows: 1
   :widths: 4 4

   * - Removed Application name
     - Current Application name
   * - ``application_sc_sr``
     - ``application_sc_super_resolution``
   * - ``sc_svc_sr_application``
     - ``sc_svc_super_resolution_application``
   * - ``ScSvcSrApplicationStrategy``
     - ``ScSvcSuperResolutionApplicationStrategy``
   * - ``ApplicationScSrConf``
     - ``ApplicationScSuperResolutionConf``
   * - ``ScSVCSr`` from the old Application runner
     - ``ScSVCSuperResolution`` from the new Application runner

These old Application imports intentionally fail. Do not mechanically rename
Benchmark internals: its ``sc_svc_sr`` task, runner, strategy, profile, CF
routing, and compatible result names remain Benchmark-owned.

What to verify after migration
------------------------------

1. Confirm the template order is VisiumHD, Xenium, then Visium, and that
   ``sc-SVC`` specifies its mode.
2. Confirm a Xenium cell type is selected through the YAML or ``--select-ct``
   and verify outputs against :doc:`application-reference`.
3. Confirm output consumers expect one H5AD for sp/SR and two H5ADs for
   cluster mode.
4. Update links to the renamed Visium notebook and output before removing any
   old local artifacts.

For the full current schema, use :doc:`application-reference`; for the
scientific distinction between modes, use :doc:`concepts`.
