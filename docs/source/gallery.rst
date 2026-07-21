Curated Notebooks and Legacy Assets
===================================

Curated notebooks are included in the clean repository under
``reproduce/benchmark/`` and ``reproduce/case/``, together with the root
application analysis notebook. They preserve the paper-facing workflows and
embedded historical outputs, but they are not part of the installed package.

Other historical assets remain in ``REVISE-legacy``. These include stored
outputs, internal implementation plans, duplicate documentation links, release
observations, and supporting analysis files that are not required by the
maintained source tree.

The root ``legacy-assets.json`` records the exact source commit and, for every
excluded tracked path, its SHA-256, role, exclusion reason, and retrieval
command. This keeps the remaining historical material recoverable from the
exact source commit without treating it as part of the maintained product.

Evidence boundary
-----------------

The presence of a notebook or indexed legacy item proves only where its bytes
came from. It does not establish that version ``0.1.0rc1`` reran a displayed
result or biologically validated a downstream pattern. Use the current CLI and
provenance contracts for new runs.
