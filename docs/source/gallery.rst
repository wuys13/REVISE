:orphan:

Curated Reproduction Notebooks
==============================

The documentation navigation exposes one curated layer of notebook pages under
``docs/benchmark/`` and ``docs/case/``. The source notebooks live under
``reproduce/benchmark/`` and ``reproduce/case/``. Their ``.nblink`` stubs keep
the paper-facing workflows and embedded historical outputs visible without
copying the notebooks.

Application gallery
-------------------

The Application cases appear in this order in the website navigation:

1. ``VisiumHD_sp_SVC.ipynb`` — VisiumHD, sp-SVC.
2. ``Xenium_sc_SVC_T.ipynb`` — Xenium, sc-SVC T cells.
3. ``Xenium_sc_SVC_Fibroblast.ipynb`` — Xenium, sc-SVC Fibroblast.
4. ``Xenium_sc_SVC_Monocyte.ipynb`` — Xenium, sc-SVC Mono/Macro.
5. ``Visium_sc_SVC_mouse_brain.ipynb`` — Visium, sc-SVC mouse brain in SR
   mode.

The three Xenium notebooks preserve downstream analysis for different selected
cell types. They share the one current ``Xenium.yaml`` reconstruction template;
the current CLI chooses the type with ``--select-ct``. The Visium case is an
sc-SVC SR-mode example, not a separate public SVC category.

Evidence boundary
-----------------

Notebook pages are static historical snapshots and are not executed during
documentation builds. Their presence proves that a workflow and its displayed
outputs are preserved in the repository; it does not establish that the
current source checkout reran a result or biologically validated a downstream
pattern. Use the current CLI, its H5AD artifacts, and ``provenance.json`` for
new runs.

For the current input and execution contract, use :doc:`quickstart`; for
download and optional-analysis dependencies, see
`reproduce/README.md`_.

.. _reproduce/README.md: https://github.com/wuys13/REVISE/blob/main/reproduce/README.md
