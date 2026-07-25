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
- Assignment-guidance tests establish software behavior, not that posterior
  compatibility improves reconstruction, subtype recovery, or imputation.
- Cost guidance is the unified product mechanism because it is portable across
  supported local solvers; this is not evidence that it is scientifically
  superior to the POT-only reference ablation.
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
TACCO and assignment-guidance solver smokes. An ad-hoc local run may still
report the TACCO checks as skipped when the exact optional dependency is
absent; a synthetic stub or source-checkout import is not installed-solver
evidence.
