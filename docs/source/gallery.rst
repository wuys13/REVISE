:orphan:

Curated Reproduction Notebooks
==============================

The documentation navigation exposes a single layer of curated notebook
pages under ``docs/benchmark/`` and ``docs/case/``. The source notebooks live
under ``reproduce/benchmark/`` and ``reproduce/case/``; the ``.nblink`` stubs
keep those paper-facing workflows and embedded historical outputs visible
without copying or editing the notebooks.

This page is an unlisted explanation page rather than a navigation parent.
For the current input and execution contract, see :doc:`Quick Start
<quickstart>`.

Evidence boundary
-----------------

Notebook pages are static historical snapshots and are not executed during
documentation builds. Their presence proves only that the workflow is
preserved in the repository; it does not establish that the current source
checkout reran a displayed result or biologically validated a downstream
pattern. Use the current CLI and ``provenance.json`` contract for new runs.

Canonical application-analysis set
----------------------------------

The maintained application-analysis notebooks are
``Xenium_sc_SVC_T.ipynb``, ``Xenium_sc_SVC_Fibroblast.ipynb``,
``Xenium_sc_SVC_Monocyte.ipynb``, and ``VisiumHD_sp_SVC.ipynb`` under
``reproduce/case/``. ``sc_SVC_sr_case_Visium_mouse_brain.ipynb`` remains a
separate preserved workflow outside that four-notebook set.
