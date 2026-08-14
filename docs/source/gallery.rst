Curated Reproduction Notebooks
==============================

Curated notebooks are included under ``reproduce/benchmark/`` and
``reproduce/case/``. They preserve the paper-facing workflows and embedded
historical outputs but are not part of the installed Python package.

Benchmark notebooks inspect Sim2Real-ST metric outputs. Application notebooks
cover reconstruction and downstream cell-state, pathway, spatial-pattern, and
other analyses. Some require optional package extras and external reference
resources.

Evidence boundary
-----------------

The presence of a notebook proves only that the workflow is preserved in the
repository. It does not establish that the current source checkout reran a
displayed result or biologically validated a downstream pattern. Use the
current CLI and ``provenance.json`` contract for new runs.

Canonical application-analysis set
----------------------------------

The maintained application-analysis notebooks are
``Xenium_sc_SVC_T.ipynb``, ``Xenium_sc_SVC_Fibroblast.ipynb``,
``Xenium_sc_SVC_Monocyte.ipynb``, and ``VisiumHD_sp_SVC.ipynb`` under
``reproduce/case/``. ``sc_SVC_sr_case_Visium_mouse_brain.ipynb`` remains a
separate preserved workflow outside that four-notebook set.
