from __future__ import annotations

import numpy as np
import pytest

from revise.backend.ops.assignment import AssignmentState
from revise.backend.ops.assignment_guidance import (
    assignment_compatibility,
    ot_cost_guidance,
)
from revise.backend.ops.local_ot import solve_local_ot


def _conditioned_cost() -> np.ndarray:
    state = AssignmentState(
        values=np.array([[0.9, 0.1], [0.1, 0.9]]),
        observation_labels=("left", "right"),
        category_labels=("A", "B"),
        source="solver_smoke",
        level="Level1",
        value_semantics="soft",
        lineage=[{"operation": "synthetic_solver_smoke"}],
    )
    affinity = assignment_compatibility(
        state,
        state,
        beta=1.0,
        min_affinity=0.05,
    )
    return ot_cost_guidance(
        np.array([[0.0, 1.0], [1.0, 0.0]]),
        affinity,
        strength=0.2,
    )


def _solve_with_installed_solver(method: str) -> np.ndarray:
    return solve_local_ot(
        [0.5, 0.5],
        [0.5, 0.5],
        _conditioned_cost(),
        method=method,
        pot_reg=0.1,
        pot_reg_m=1.0,
        pot_reg_type="kl",
    )


def test_real_pot_executes_conditioned_cost_candidate():
    pytest.importorskip("ot", reason="installed POT smoke requires POT")

    coupling = _solve_with_installed_solver("pot")

    assert coupling.shape == (2, 2)
    assert np.isfinite(coupling).all()
    assert float(coupling.sum()) > 0.0


def test_real_tacco_executes_conditioned_cost_candidate():
    pytest.importorskip("tacco", reason="optional TACCO smoke requires tacco")

    coupling = _solve_with_installed_solver("tacco")

    assert coupling.shape == (2, 2)
    assert np.isfinite(coupling).all()
    assert float(coupling.sum()) > 0.0


def test_reference_tacco_is_rejected_before_any_solver_fallback(
    tmp_path,
    monkeypatch,
):
    import builtins

    from revise.backend.ops import local_ot
    from revise.framework import REVISEPipeline

    pytest.importorskip("tacco", reason="optional TACCO smoke requires tacco")
    pot_imports = []
    pot_solves = []
    original_import = builtins.__import__

    def import_spy(name, *args, **kwargs):
        if name == "ot" or name.startswith("ot."):
            pot_imports.append(name)
        return original_import(name, *args, **kwargs)

    def solve_spy(*args, **kwargs):
        pot_solves.append((args, kwargs))
        raise AssertionError("no local solver may run after incompatible preflight")

    monkeypatch.setattr(builtins, "__import__", import_spy)
    monkeypatch.setattr(local_ot, "solve_local_ot", solve_spy)
    data_root = tmp_path / "data"
    data_root.mkdir()
    pipeline = REVISEPipeline(config_path="revise/revise.yaml")

    with pytest.raises(ValueError, match="TACCO"):
        pipeline._run_with_algorithm_overrides(
            profile="benchmark_seg",
            runtime_overrides={
                "platform": "sim2real",
                "confounding": "segmentation",
            },
            io_overrides={
                "data_root": str(data_root),
                "output_root": str(tmp_path / "output"),
                "sample_name": "sample",
                "st_file": "st.h5ad",
                "gt_svc_file": "gt.h5ad",
                "sc_ref_file": "ref.h5ad",
                "seg_method": "seg_1",
            },
            algorithm_overrides={
                "local_refinement": {
                    "guidance": "prefer",
                    "compatibility": {"mode": "reference"},
                },
                "ot": {"lr": {"solver": "tacco"}},
            },
            dry_run=True,
        )

    assert pot_imports == []
    assert pot_solves == []
