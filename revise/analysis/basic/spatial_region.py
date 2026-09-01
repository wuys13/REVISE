"""Non-overlapping spatial-window diversity and Region helpers."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd


def effective_number(labels: Sequence[object]) -> float:
    """Return Shannon's effective number of observed labels (``exp(H)``)."""
    values = pd.Series(labels).dropna().astype(str)
    if values.empty:
        return float("nan")
    probabilities = values.value_counts(normalize=True).to_numpy()
    return float(np.exp(-(probabilities * np.log(probabilities)).sum()))


def assign_square_windows(
    coordinates: pd.DataFrame,
    *,
    window_side_length: float,
    origin: tuple[float, float] | None = None,
) -> pd.DataFrame:
    """Assign every coordinate to one deterministic, non-overlapping square window."""
    if not coordinates.index.is_unique:
        raise ValueError("Coordinates must have unique unit IDs")
    if window_side_length <= 0 or not np.isfinite(window_side_length):
        raise ValueError("window_side_length must be positive and finite")
    if not {"x", "y"} <= set(coordinates.columns):
        raise KeyError("Coordinates must contain x and y columns")
    work = coordinates.loc[:, ["x", "y"]].copy().astype(float)
    if not np.isfinite(work.to_numpy()).all():
        raise ValueError("Coordinates must be finite")
    x0, y0 = origin if origin is not None else (float(work["x"].min()), float(work["y"].min()))
    work["window_x_index"] = np.floor((work["x"] - x0) / window_side_length).astype(int)
    work["window_y_index"] = np.floor((work["y"] - y0) / window_side_length).astype(int)
    work["window_id"] = (
        work["window_x_index"].astype(str) + "_" + work["window_y_index"].astype(str)
    )
    work["window_side_length"] = float(window_side_length)
    return work


def _align_labels(labels: pd.Series, unit_ids: pd.Index, name: str) -> pd.Series:
    if not labels.index.is_unique:
        raise ValueError(f"{name} labels must have unique unit IDs")
    if set(labels.index) != set(unit_ids):
        raise ValueError(f"{name} labels must contain the same unit IDs as window assignments")
    labels = labels.reindex(unit_ids)
    if labels.isna().any():
        raise ValueError(f"{name} labels must not contain missing values")
    return labels.astype(str)


def compute_window_diversity(
    window_assignments: pd.DataFrame,
    raw_labels: pd.Series,
    reconstructed_labels: pd.Series,
    *,
    min_units_per_window: int = 1,
) -> pd.DataFrame:
    """Compute Raw, reconstructed and delta effective diversity per spatial window."""
    if min_units_per_window < 1:
        raise ValueError("min_units_per_window must be at least one")
    if "window_id" not in window_assignments:
        raise KeyError("window_assignments must contain window_id")
    if not window_assignments.index.is_unique:
        raise ValueError("window_assignments must have unique unit IDs")
    raw = _align_labels(raw_labels, window_assignments.index, "Raw")
    recon = _align_labels(reconstructed_labels, window_assignments.index, "Reconstructed")
    rows: list[dict[str, object]] = []
    groups = window_assignments.groupby("window_id", sort=True)
    for window_id, frame in groups:
        labels_raw = raw.reindex(frame.index)
        labels_recon = recon.reindex(frame.index)
        n_units = int(frame.shape[0])
        valid = n_units >= min_units_per_window
        neff_raw = effective_number(labels_raw) if valid else float("nan")
        neff_recon = effective_number(labels_recon) if valid else float("nan")
        rows.append(
            {
                "window_id": window_id,
                "n_units": n_units,
                "valid_window": valid,
                "k_obs_raw": int(labels_raw.nunique()) if valid else np.nan,
                "k_obs_recon": int(labels_recon.nunique()) if valid else np.nan,
                "neff_raw": neff_raw,
                "neff_recon": neff_recon,
                "delta_neff": neff_recon - neff_raw if valid else float("nan"),
            }
        )
    return pd.DataFrame(rows)


def assign_anatomy_regions(
    window_assignments: pd.DataFrame,
    level1_labels: pd.Series,
    *,
    tumor_label: str = "Tumor",
    normal_source_label: str = "Intestinal Epithelial",
) -> pd.DataFrame:
    """Assign Level1-derived anatomy labels to spatial windows."""
    if "window_id" not in window_assignments:
        raise KeyError("window_assignments must contain window_id")
    labels = _align_labels(level1_labels, window_assignments.index, "Level1")
    rows: list[dict[str, object]] = []
    for window_id, frame in window_assignments.groupby("window_id", sort=True):
        values = labels.reindex(frame.index)
        tumor_units = int((values == tumor_label).sum())
        normal_units = int((values == normal_source_label).sum())
        if tumor_units and normal_units:
            region_label = "Interface"
        elif tumor_units:
            region_label = "Tumor"
        elif normal_units:
            region_label = "Normal"
        else:
            region_label = "Other"
        rows.append(
            {
                "window_id": window_id,
                "tumor_units": tumor_units,
                "normal_units": normal_units,
                "tumor_present": bool(tumor_units),
                "normal_present": bool(normal_units),
                "selected_type_count": int(bool(tumor_units)) + int(bool(normal_units)),
                "region_label": region_label,
            }
        )
    return pd.DataFrame(rows)


def summarize_region(
    window_metrics: pd.DataFrame,
    *,
    threshold: float,
    window_side_length: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Flag high-diversity Region windows and report absolute area/fractions."""
    required = {"window_id", "valid_window", "neff_recon", "n_units"}
    missing = required - set(window_metrics.columns)
    if missing:
        raise KeyError(f"Window metrics are missing columns: {sorted(missing)}")
    if threshold <= 0 or window_side_length <= 0:
        raise ValueError("threshold and window_side_length must be positive")
    flagged = window_metrics.copy()
    flagged["in_region"] = flagged["valid_window"] & (flagged["neff_recon"] >= threshold)
    valid = flagged.loc[flagged["valid_window"]]
    region = flagged.loc[flagged["in_region"]]
    n_valid = int(valid.shape[0])
    valid_units = int(valid["n_units"].sum())
    region_units = int(region["n_units"].sum())
    window_area = float(window_side_length**2)
    summary = pd.DataFrame(
        [
            {
                "threshold": float(threshold),
                "window_side_length": float(window_side_length),
                "window_area": window_area,
                "n_valid_windows": n_valid,
                "n_region_windows": int(region.shape[0]),
                "region_area": float(region.shape[0] * window_area),
                "region_area_fraction": float(region.shape[0] / n_valid) if n_valid else np.nan,
                "n_valid_units": valid_units,
                "n_region_units": region_units,
                "region_unit_fraction": float(region_units / valid_units) if valid_units else np.nan,
            }
        ]
    )
    return summary, flagged
