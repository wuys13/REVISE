from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd


def _require_omicverse():
    try:
        import omicverse
    except ImportError as exc:
        raise ImportError(
            "Palantir trajectory inference requires OmicVerse; install it with "
            '`python -m pip install "revise-svc[trajectory]"`.'
        ) from exc
    return omicverse


def infer_palantir(
    adata,
    *,
    groupby: str,
    basis: str,
    use_rep: str,
    n_comps: int,
    origin_cells: str,
    terminal_cells: Sequence[str],
    num_waypoints: int,
) -> pd.Series:
    """Run the notebook's OmicVerse Palantir inference and return pseudotime."""
    if groupby not in adata.obs:
        raise ValueError(f"groupby '{groupby}' not found in adata.obs")

    labels = set(adata.obs[groupby].dropna().astype(str))
    if str(origin_cells) not in labels:
        raise ValueError(f"origin_cells label '{origin_cells}' not found in '{groupby}'")
    if isinstance(terminal_cells, (str, bytes)) or not terminal_cells:
        raise ValueError("terminal_cells must be a non-empty sequence of labels")
    missing_terminal = [str(label) for label in terminal_cells if str(label) not in labels]
    if missing_terminal:
        raise ValueError(
            f"terminal_cells labels not found in '{groupby}': {missing_terminal}"
        )

    omicverse = _require_omicverse()
    trajectory = omicverse.single.TrajInfer(
        adata,
        basis=basis,
        groupby=groupby,
        use_rep=use_rep,
        n_comps=n_comps,
    )
    trajectory.set_origin_cells(origin_cells)
    trajectory.set_terminal_cells(terminal_cells)
    trajectory.inference(method="palantir", num_waypoints=num_waypoints)

    result_key = "palantir_pseudotime"
    if result_key not in adata.obs:
        raise RuntimeError(
            f"Palantir trajectory inference did not create '{result_key}'"
        )
    pseudotime = adata.obs[result_key]
    if not pd.api.types.is_numeric_dtype(pseudotime):
        raise ValueError("palantir_pseudotime must be numeric")
    if not np.isfinite(pseudotime.to_numpy()).all():
        raise ValueError("palantir_pseudotime must contain only finite values")
    if not pseudotime.index.equals(adata.obs_names):
        raise ValueError("palantir_pseudotime index must match adata.obs_names")
    return pseudotime.copy()


__all__ = ["infer_palantir"]
