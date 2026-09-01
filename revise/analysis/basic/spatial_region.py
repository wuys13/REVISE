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


def convert_coordinates_to_microns(
    coordinates: pd.DataFrame, *, microns_per_coordinate: float
) -> pd.DataFrame:
    """Return a scaled coordinate copy in microns without mutating the input."""
    if microns_per_coordinate <= 0 or not np.isfinite(microns_per_coordinate):
        raise ValueError("microns_per_coordinate must be positive and finite")
    if not {"x", "y"} <= set(coordinates.columns):
        raise KeyError("Coordinates must contain x and y columns")
    converted = coordinates.loc[:, ["x", "y"]].copy().astype(float)
    if not np.isfinite(converted.to_numpy()).all():
        raise ValueError("Coordinates must be finite")
    converted.loc[:, ["x", "y"]] *= float(microns_per_coordinate)
    return converted


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
    work = window_assignments.loc[:, ["window_id", "x", "y"]].copy()
    work["raw_label"] = raw.reindex(work.index).to_numpy()
    work["recon_label"] = recon.reindex(work.index).to_numpy()
    metrics = work.groupby("window_id", sort=True).agg(
        window_x=("x", "mean"), window_y=("y", "mean"), n_units=("x", "size")
    )

    def label_metrics(label_column: str, prefix: str) -> pd.DataFrame:
        counts = work.groupby(["window_id", label_column], sort=True).size().rename("count")
        proportions = counts / counts.groupby(level=0).sum()
        neff = np.exp(-(proportions * np.log(proportions)).groupby(level=0).sum()).rename(f"neff_{prefix}")
        observed = counts.groupby(level=0).size().rename(f"k_obs_{prefix}")
        return pd.concat([observed, neff], axis=1)

    metrics = metrics.join(label_metrics("raw_label", "raw")).join(label_metrics("recon_label", "recon"))
    metrics["valid_window"] = metrics["n_units"] >= min_units_per_window
    invalid = ~metrics["valid_window"]
    metrics.loc[invalid, ["k_obs_raw", "k_obs_recon", "neff_raw", "neff_recon"]] = np.nan
    metrics["delta_neff"] = metrics["neff_recon"] - metrics["neff_raw"]
    return metrics.reset_index()


def select_min_parent_support(
    window_assignments: pd.DataFrame,
) -> tuple[dict[str, int | str | None], pd.DataFrame]:
    """Select parent-window support from retained-unit curve's chord-distance knee."""
    if "window_id" not in window_assignments:
        raise KeyError("window_assignments must contain window_id")
    occupancies = window_assignments.groupby("window_id", sort=True).size()
    if occupancies.empty:
        return {"status": "insufficient_support", "min_parent_units": None, "p95_occupancy": 0}, pd.DataFrame()
    p95 = int(np.ceil(np.percentile(occupancies.to_numpy(), 95)))
    if p95 < 2:
        return {"status": "insufficient_support", "min_parent_units": None, "p95_occupancy": p95}, pd.DataFrame()

    total_units = int(occupancies.sum())
    rows: list[dict[str, float | int]] = []
    for minimum in range(1, p95 + 1):
        valid = occupancies[occupancies >= minimum]
        rows.append(
            {
                "min_parent_units": minimum,
                "n_valid_windows": int(valid.size),
                "valid_window_fraction": float(valid.size / occupancies.size),
                "n_retained_parent_units": int(valid.sum()),
                "retained_unit_fraction": float(valid.sum() / total_units),
            }
        )
    sensitivity = pd.DataFrame(rows)
    x = (sensitivity["min_parent_units"] - 1) / max(p95 - 1, 1)
    y = sensitivity["retained_unit_fraction"].to_numpy(dtype=float)
    x0, y0, x1, y1 = float(x.iloc[0]), float(y[0]), float(x.iloc[-1]), float(y[-1])
    denominator = float(np.hypot(y1 - y0, x1 - x0))
    if denominator == 0:
        distances = np.zeros(len(sensitivity))
    else:
        distances = np.abs((y1 - y0) * x.to_numpy() - (x1 - x0) * y + x1 * y0 - y1 * x0) / denominator
    sensitivity["knee_distance"] = distances
    eligible = sensitivity.loc[sensitivity["min_parent_units"] >= 2]
    maximum = float(eligible["knee_distance"].max())
    minimum = int(eligible.loc[np.isclose(eligible["knee_distance"], maximum), "min_parent_units"].min())
    return {"status": "ok", "min_parent_units": minimum, "p95_occupancy": p95}, sensitivity


def assign_anatomy_regions(
    window_assignments: pd.DataFrame,
    level1_labels: pd.Series,
    *,
    tumor_label: str = "Tumor",
    normal_source_label: str = "Intestinal Epithelial",
) -> pd.DataFrame:
    """Backward-compatible alias for the explicit anatomy candidate assignment."""
    return assign_anatomy_candidates(
        window_assignments,
        level1_labels,
        tumor_label=tumor_label,
        normal_source_label=normal_source_label,
    ).rename(columns={"level1_region": "region_label"})


def assign_anatomy_candidates(
    window_assignments: pd.DataFrame,
    level1_labels: pd.Series,
    *,
    tumor_label: str = "Tumor",
    normal_source_label: str = "Intestinal Epithelial",
) -> pd.DataFrame:
    """Assign binary Tumor/Normal candidates and categorical anatomy context."""
    if "window_id" not in window_assignments:
        raise KeyError("window_assignments must contain window_id")
    labels = _align_labels(level1_labels, window_assignments.index, "Level1")
    rows: list[dict[str, object]] = []
    for window_id, frame in window_assignments.groupby("window_id", sort=True):
        values = labels.reindex(frame.index)
        tumor_units = int((values == tumor_label).sum())
        normal_units = int((values == normal_source_label).sum())
        tumor_candidate = int(bool(tumor_units))
        normal_candidate = int(bool(normal_units))
        union_candidate = tumor_candidate + normal_candidate
        if union_candidate == 2:
            level1_region = "Interface"
        elif tumor_units:
            level1_region = "Tumor"
        elif normal_units:
            level1_region = "Normal"
        else:
            level1_region = "Other"
        rows.append(
            {
                "window_id": window_id,
                "window_x": float(frame["x"].mean()),
                "window_y": float(frame["y"].mean()),
                "tumor_units": tumor_units,
                "normal_units": normal_units,
                "tumor_candidate": tumor_candidate,
                "normal_candidate": normal_candidate,
                "union_candidate": union_candidate,
                "level1_region": level1_region,
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


def summarize_region_by_anatomy(
    window_metrics: pd.DataFrame,
    anatomy_windows: pd.DataFrame,
    *,
    window_side_length: float,
) -> pd.DataFrame:
    """Summarize diversity and Region coverage within fixed Level1 anatomy context."""
    required_metrics = {"window_id", "valid_window", "in_region", "n_units", "neff_recon", "delta_neff"}
    missing = required_metrics - set(window_metrics.columns)
    if missing:
        raise KeyError(f"Window metrics are missing columns: {sorted(missing)}")
    if "level1_region" not in anatomy_windows:
        raise KeyError("anatomy_windows must contain level1_region")
    if window_side_length <= 0:
        raise ValueError("window_side_length must be positive")
    anatomy_context = anatomy_windows.loc[:, ["window_id", "level1_region"]]
    if "level1_region" in window_metrics:
        merged = window_metrics.copy()
        expected = anatomy_context.set_index("window_id")["level1_region"].reindex(merged["window_id"])
        if not merged["level1_region"].reset_index(drop=True).equals(expected.reset_index(drop=True)):
            raise ValueError("Existing Level1 anatomy context does not match anatomy_windows")
    else:
        merged = window_metrics.merge(anatomy_context, on="window_id", how="left", validate="one_to_one")
    if merged["level1_region"].isna().any():
        raise ValueError("Every diversity window must have Level1 anatomy context")
    rows: list[dict[str, float | int | str]] = []
    for region, frame in merged.groupby("level1_region", sort=True):
        valid = frame.loc[frame["valid_window"]]
        selected = valid.loc[valid["in_region"]]
        valid_units = int(valid["n_units"].sum())
        selected_units = int(selected["n_units"].sum())
        rows.append(
            {
                "level1_region": region,
                "n_windows": int(frame.shape[0]),
                "n_valid_windows": int(valid.shape[0]),
                "n_region_windows": int(selected.shape[0]),
                "mean_neff_recon": float(valid["neff_recon"].mean()) if not valid.empty else np.nan,
                "mean_delta_neff": float(valid["delta_neff"].mean()) if not valid.empty else np.nan,
                "region_area": float(selected.shape[0] * window_side_length**2),
                "region_area_fraction": float(selected.shape[0] / valid.shape[0]) if not valid.empty else np.nan,
                "n_valid_units": valid_units,
                "n_region_units": selected_units,
                "region_unit_fraction": float(selected_units / valid_units) if valid_units else np.nan,
            }
        )
    return pd.DataFrame(rows)
