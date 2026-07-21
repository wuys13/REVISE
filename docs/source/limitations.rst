Candidate Limitations
=====================

This page defines what the ``0.1.0rc1`` evidence does and does not support.
These boundaries are release claims, not a list of planned features.

Data and biology
----------------

- Real-data end-to-end reconstruction and paper-result reproduction are
  deferred until the non-data candidate gates are complete and the owner
  authorizes the download.
- No current gate proves biological validation, clinical validity, cell-level
  localization accuracy, or biological parity between POT and TACCO.
- Historical executed notebooks are preserved material. They are not treated
  as results reproduced by this candidate.
- Official dataset/archive identifiers, a DOI/BibTeX record, named project
  roles, release contacts, and the final canonical repository are pending owner
  confirmation.

Scale
-----

The installed end-to-end POT smoke uses small synthetic inputs with 52
observations and 52 genes. A 50,000-observation sparse graph check and a
200,000-cell quota check exercise individual components only. The TACCO
coverage is a small synthetic solver smoke.

These tests do not establish a complete production-scale pipeline limit,
whole-section memory requirement, or throughput guarantee. Current routes can
materialize dense intermediate arrays. Users evaluating larger inputs should
start with representative subsets and measure memory in their own environment;
the candidate does not publish a safe maximum.

Spot SR localization
--------------------

When the sST fallback creates virtual-cell rows, they share the source spot
centroid. With no ``pm_on_cell`` matrix, proportions are converted with
``np.round``, repaired to an exact quota, and placed into those virtual-cell
rows by a seeded random permutation.

The evidence proves quota composition and same-seed repeatability. It is not a
nucleus or cell-localization result and does not prove sub-spot morphology,
true cell identity, or biological position.

Metrics
-------

Normalized NRMSE is directional because its denominator is the normalized
ground-truth mean. Constant/zero-mean genes can yield ``NaN`` or positive
infinity. SSIM uses aligned expression vectors in row order and is not a
coordinate-aware spatial-image metric. See :doc:`benchmark` for the formula and
preprocessing sequence.

Solver and platform support
---------------------------

POT is a base dependency. TACCO 0.5.0 is optional, and selecting TACCO with a
missing/incompatible dependency or unsupported configuration fails. REVISE
does not fall back to POT.

The declared installed-artifact CI matrix targets Linux with Python 3.10 and
3.11. Other operating systems and interpreters are not claimed by this
candidate until matching gate evidence is attached.

Packaging and release identity
------------------------------

The release constraints pin release-critical packages but are not complete
transitive environment locks. Hosted CI and final publication must still prove
the exact tagged artifacts. Until publication and repository cutover are
completed, this document describes a candidate built from source rather than a
verified current PyPI release.
