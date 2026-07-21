from __future__ import annotations

import importlib
import sys
import types

import numpy as np


def test_bhattacharyya_distance_accepts_integer_count_matrices(monkeypatch):
    numba = types.ModuleType("numba")
    numba.njit = lambda *args, **kwargs: (
        (lambda function: function) if not args else args[0]
    )
    numba.prange = range
    numba.get_num_threads = lambda: 1
    numba.set_num_threads = lambda count: None
    monkeypatch.setitem(sys.modules, "numba", numba)
    monkeypatch.delitem(sys.modules, "revise.backend.ops.distance", raising=False)
    distance = importlib.import_module("revise.backend.ops.distance")

    left = np.array([[2, 1], [1, 3]], dtype=np.int64)
    right = np.array([[3, 1], [1, 2]], dtype=np.int64)

    observed = distance.bhattacharyya_distance(left, right)
    expected = distance.bhattacharyya_distance(
        left.astype(np.float64),
        right.astype(np.float64),
    )

    np.testing.assert_allclose(observed, expected)
    assert np.isfinite(observed).all()
