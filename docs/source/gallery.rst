Legacy Research Assets
======================

Executed notebooks, stored outputs, research figures, and supporting analysis
files remain in ``REVISE-legacy``. They are intentionally not included in the
clean repository and are not part of the installed package.

The root ``legacy-assets.json`` records the exact source commit and, for every
excluded tracked path, its SHA-256, role, exclusion reason, and retrieval
command. This keeps historical material recoverable without coupling the public
documentation build to notebook files or notebook-specific dependencies.

Evidence boundary
-----------------

The presence of an item in the legacy index proves only where the historical
bytes came from. It does not establish that version ``0.1.0rc1`` reproduced a
displayed result, that the inputs match an official release dataset, or that a
downstream pattern is biologically validated. Use the current CLI and
provenance contracts for new runs.
