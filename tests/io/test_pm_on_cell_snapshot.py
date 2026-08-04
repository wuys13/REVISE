from __future__ import annotations

import hashlib
import logging
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from revise.backend.kernels.spot_sr import SpotSrKernel
from revise.io.input_service import REVISEInputService


def _kernel(pm_on_cell, *, seed=17):
    config = SimpleNamespace(
        pm_on_cell=pm_on_cell,
        svc_completeness=True,
        sr_assignment_seed=seed,
        cell_type_col="Level1",
    )
    return SpotSrKernel(config, logging.getLogger("test-pm-snapshot"))


def _svc_obs():
    return pd.DataFrame(
        {
            "spot_name": ["spot-a", "spot-a"],
            "cell_id": ["c2", "c1"],
        }
    )


def _quota():
    return pd.DataFrame(
        [[1, 1]],
        index=["spot-a"],
        columns=["Mono_Macro", "T"],
    )


def test_snapshot_hashes_and_parses_the_same_single_read_bytes(monkeypatch, tmp_path):
    path = tmp_path / "PM_on_cell.csv"
    original = b",T,Mono/Macro\nc1,0.8,0.2\nc2,0.1,0.9\n"
    replacement = b",T,Mono/Macro\nc1,0.4,0.6\nc2,0.7,0.3\n"
    path.write_bytes(original)
    real_read_bytes = Path.read_bytes
    calls = []

    def replace_after_read(self):
        payload = real_read_bytes(self)
        calls.append(self)
        self.write_bytes(replacement)
        return payload

    monkeypatch.setattr(Path, "read_bytes", replace_after_read)

    frame, identity = REVISEInputService().snapshot_pm_on_cell(path)

    assert calls == [path]
    assert identity == {
        "role": "pm_on_cell",
        "path": str(path),
        "sha256": hashlib.sha256(original).hexdigest(),
    }
    expected = pd.DataFrame(
        [[0.8, 0.2], [0.1, 0.9]],
        index=pd.Index(["c1", "c2"]),
        columns=["T", "Mono_Macro"],
    )
    pd.testing.assert_frame_equal(frame, expected)

    kernel = _kernel(frame)
    kernel.assign_cell_types(_svc_obs(), _quota())
    pd.testing.assert_frame_equal(
        kernel.pm_on_cell,
        expected.loc[["c2", "c1"], ["Mono_Macro", "T"]],
    )


@pytest.mark.parametrize(
    "frame",
    [
        pd.DataFrame([[0.2, 0.2]], index=["c1"], columns=["A", "B"]),
        pd.DataFrame([[-0.1, 1.1]], index=["c1"], columns=["A", "B"]),
        pd.DataFrame([[np.nan, np.nan]], index=["c1"], columns=["A", "B"]),
        pd.DataFrame([[np.inf, -np.inf]], index=["c1"], columns=["A", "B"]),
        pd.DataFrame([["x", "y"]], index=["c1"], columns=["A", "B"]),
    ],
)
def test_snapshot_rejects_invalid_probability_matrices(frame, tmp_path):
    path = tmp_path / "PM_on_cell.csv"
    frame.to_csv(path)

    with pytest.raises(ValueError, match="pm_on_cell"):
        REVISEInputService().snapshot_pm_on_cell(path)


@pytest.mark.parametrize(
    "payload",
    [
        b",A,A\nc1,0.5,0.5\n",
        b",,B\nc1,0.5,0.5\n",
    ],
)
def test_snapshot_rejects_duplicate_or_blank_raw_headers(payload, tmp_path):
    path = tmp_path / "PM_on_cell.csv"
    path.write_bytes(payload)

    with pytest.raises(ValueError, match="header"):
        REVISEInputService().snapshot_pm_on_cell(path)


def test_missing_snapshot_stays_none_after_file_appears(tmp_path):
    path = tmp_path / "PM_on_cell.csv"

    frame, identity = REVISEInputService().snapshot_pm_on_cell(path)
    path.write_text(",A,B\nc1,1,0\nc2,0,1\n", encoding="utf-8")

    assert frame is None
    assert identity is None
    kernel = _kernel(frame, seed=23)
    assigned = kernel.assign_cell_types_random(
        _svc_obs(),
        pd.DataFrame([[1, 1]], index=["spot-a"], columns=["A", "B"]),
    )
    expected = np.random.default_rng(23).permutation(["A", "B"]).tolist()
    assert assigned["cell_type"].tolist() == expected


def test_kernel_does_not_discover_pm_from_a_path(tmp_path):
    path = tmp_path / "PM_on_cell.csv"
    path.write_text(",A,B\nc1,1,0\nc2,0,1\n", encoding="utf-8")
    config = SimpleNamespace(
        pm_on_cell=None,
        pm_on_cell_file=str(path),
        svc_completeness=True,
        sr_assignment_seed=5,
        cell_type_col="Level1",
    )

    kernel = SpotSrKernel(config, logging.getLogger("test-no-path-discovery"))

    assert kernel.pm_on_cell is None


def test_kernel_requires_exact_active_cell_and_type_sets():
    extra_cell = pd.DataFrame(
        [[0.5, 0.5], [0.5, 0.5], [0.5, 0.5]],
        index=["c1", "c2", "c3"],
        columns=["Mono_Macro", "T"],
    )
    with pytest.raises(ValueError, match="extra"):
        _kernel(extra_cell).assign_cell_types(_svc_obs(), _quota())

    extra_type = pd.DataFrame(
        [[0.4, 0.4, 0.2], [0.4, 0.4, 0.2]],
        index=["c1", "c2"],
        columns=["Mono_Macro", "T", "B"],
    )
    with pytest.raises(ValueError, match="extra"):
        _kernel(extra_type).assign_cell_types(_svc_obs(), _quota())
