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
Benchmark runs retain the historical ``<data-root>/PM_on_cell.csv`` convention.
Application runs declare an exact ``inputs.pm_on_cell.path`` instead, so the
logical output name and input YAML make the selected prior explicit.

Output persistence
------------------

Application output persistence is scoped to one unique run directory per invocation:
``<output.dir>/<output.name>/application__<svc-type>/<timestamp_uuid>/``. Final
H5AD artifact(s) are written directly into that directory beside the run
envelope: ``provenance.json``, ``merged_config.json``, and ``run.log``.
Completed input validation adds ``preflight.json``. Standard sc-SVC keeps its
spatial and expression files in the same run directory. Dry runs allocate an
independent run directory and write no H5AD; repeated formal runs allocate a
new leaf. Application creates no temporary H5AD files or shared fixed-name
copies and does not use ``os.replace`` for H5AD output.

A caught write or later lifecycle failure removes H5AD files created by that
run while retaining its failed manifest and log. If the filesystem refuses a
cleanup operation, the failed manifest records ``output_cleanup_errors``
instead of hiding the original run failure. If terminal provenance itself
cannot be persisted, the last manifest can remain ``running``. An uncatchable
process death or power loss can leave a partial H5AD inside the isolated
directory and can additionally leave its run lock. Neither state is a
successful result; because every Application leaf must be fresh, it cannot mix
files with an earlier or concurrent run.

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
