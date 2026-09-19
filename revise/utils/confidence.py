"""Small, shared confidence metrics for posterior-like matrices.

The reference-preparation scorer and OT annotation publication use the same
row-wise definitions.  The functions intentionally keep the matrix library's
existing reduction semantics: in particular, ``max_confidence`` delegates to
``matrix.max(axis=1)`` so a pandas DataFrame keeps pandas' default ``skipna``
behaviour.
"""

from __future__ import annotations

from typing import Any

import numpy as np


CONFIDENCE_METADATA_KEY = "revise_confidence"


def max_confidence(matrix: Any):
    """Return the row-wise maximum candidate confidence.

    This is deliberately a thin reduction.  It preserves the current
    ``pandas.DataFrame.max(axis=1)`` behaviour used when publishing OT
    annotations as well as NumPy's existing behaviour used by reference
    scoring.
    """

    return matrix.max(axis=1)


def entropy_certainty(matrix: Any) -> np.ndarray:
    """Return normalized certainty from each row's entropy.

    The arithmetic is kept identical to the reference-selection scorer,
    including the zero-safe logarithm and the denominator based on the number
    of candidate columns.
    """

    values = np.asarray(matrix)
    log_values = np.zeros_like(values)
    np.log(values, out=log_values, where=values > 0)
    entropy = -(values * log_values).sum(axis=1)
    return 1.0 - entropy / np.log(values.shape[1])


def top2_margin(matrix: Any) -> np.ndarray:
    """Return the difference between the two largest values in each row."""

    ordered = np.sort(np.asarray(matrix), axis=1)
    return ordered[:, -1] - ordered[:, -2]


# Names used by the reference-preparation public metrics remain available as
# neutral aliases for downstream callers that need the same row definitions.
normalized_certainty = entropy_certainty
margin = top2_margin


__all__ = [
    "CONFIDENCE_METADATA_KEY",
    "max_confidence",
    "entropy_certainty",
    "normalized_certainty",
    "top2_margin",
    "margin",
]
