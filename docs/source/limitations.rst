Current Limitations
===================

These boundaries distinguish implemented and tested software behavior from
scientific conclusions.

Data and biology
----------------

- Real-data end-to-end reconstruction and paper-result reproduction have not
  yet been rerun for the current source checkout.
- No current gate proves biological validation, clinical validity, cell-level
  localization accuracy, or biological parity across OT implementations.
- Historical notebook outputs are preserved material, not results reproduced
  by the current checkout.
- Route-specific assignment tests establish software behavior, not that posterior
  compatibility improves reconstruction, subtype recovery, or imputation.
- Cost conditioning is the fixed product mechanism for sp-SVC and sc-SVC-sr
  because it is portable across supported local solvers; this is not evidence
  that it is scientifically superior to another biological model.
- Paper data and reproduced results are available at
  ``https://zenodo.org/records/17705737``.

Scale
-----

The installed end-to-end synthetic smoke uses 52 observations and 52 genes. A
50,000-observation sparse graph check and a 200,000-cell quota check exercise
individual components only. These tests do not establish a complete
production-scale pipeline limit, memory requirement, or throughput guarantee.

Spot super-resolution localization
-----------------------------------

When segmentation-derived centers are present in
``uns["revise_cell_locations"]``, sc-SVC-sr uses those x/y coordinates. Other
virtual-cell rows share the source spot centroid. With no ``PM_on_cell.csv``
matrix, proportions are converted with ``np.round``, repaired to an exact quota,
and assigned to the existing virtual-cell rows by a seeded random permutation.
Supplied centers must already use the same coordinate system and scale as
``obsm["spatial"]``. The histology-prior preprocessor checks shared spot
coordinates before writing the table, but a manually constructed H5AD must
preserve this alignment itself.

This proves quota composition and same-seed repeatability. The random assignment
is not a nucleus or cell-localization result and does not establish sub-spot
morphology or which inferred cell type belongs to a supplied center.

When a PM prior is present, REVISE requires exact active axes, numeric finite
values in ``[0, 1]``, and row sums within an absolute tolerance of ``1e-6``. It
does not clip, normalize, or interpret this software contract as biological
validation. PM is not a case table, cohort registry, or generic assignment
posterior.

Publication
-----------

The 1.x single-file and paired publishers stage a same-directory temporary
H5AD, reload it before replacement, and use best-effort caught-exception
rollback. Paired outputs are not reader-atomic or crash-atomic. The caller must
guarantee one writer per stable public target; violating that precondition is
undefined.

Metrics
-------

Normalized NRMSE is directional because its denominator is the normalized
ground-truth mean. Constant or zero-mean genes can yield ``NaN`` or positive
infinity. SSIM uses aligned expression vectors in row order and is not a
coordinate-aware spatial-image metric. See :doc:`benchmark`.

Dependency and platform support
-------------------------------

Optional capabilities require their corresponding installation extras.
Selecting a missing or incompatible implementation fails explicitly. The
tested dependency constraints target Linux with Python 3.10 and 3.11; other
systems require their own validation.

The CI TACCO gate downloads the candidate wheel produced by the package job,
installs it with TACCO 0.5.0, verifies the import from outside the source
checkout, copies the solver tests outside that checkout, and runs them with
pytest importlib mode. The test process also asserts that ``revise.__file__``
is under the isolated candidate-wheel environment before running the base
TACCO and local-refinement solver smokes. An ad-hoc local run may still
report the TACCO checks as skipped when the exact optional dependency is
absent; a synthetic stub or source-checkout import is not installed-solver
evidence.
