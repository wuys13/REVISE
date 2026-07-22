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
