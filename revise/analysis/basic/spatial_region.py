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


def select_window_scale(
    coordinates: pd.DataFrame,
    *,
    candidate_window_sides: Sequence[float],
    min_parent_units: int = 4,
    origin: tuple[float, float] | None = None,
) -> tuple[dict[str, float | int | str | None], pd.DataFrame]:
    """Select a parent-specific square side from occupancy-only retention curves."""
    if min_parent_units < 1:
        raise ValueError("min_parent_units must be at least one")
    candidates = sorted({float(side) for side in candidate_window_sides})
    if not candidates or any(not np.isfinite(side) or side <= 0 for side in candidates):
        raise ValueError("candidate_window_sides must contain positive finite values")
    if coordinates.empty:
        return {"status": "insufficient_support", "window_side_length": None, "min_parent_units": min_parent_units}, pd.DataFrame()

    rows: list[dict[str, float | int]] = []
    total_units = int(coordinates.shape[0])
    for side in candidates:
        assigned = assign_square_windows(coordinates, window_side_length=side, origin=origin)
        occupancy = assigned.groupby("window_id", sort=True).size()
        retained = occupancy.loc[occupancy >= min_parent_units]
        rows.append(
            {
                "window_side_length": side,
                "n_tissue_windows": int(occupancy.size),
                "n_valid_windows": int(retained.size),
                "valid_window_fraction": float(retained.size / occupancy.size) if occupancy.size else np.nan,
                "retained_parent_units": int(retained.sum()),
                "retained_parent_unit_fraction": float(retained.sum() / total_units),
            }
        )
    sensitivity = pd.DataFrame(rows)
    if not (sensitivity["n_valid_windows"] > 0).any():
        return {
            "status": "insufficient_support",
            "window_side_length": None,
            "min_parent_units": min_parent_units,
        }, sensitivity

    x = np.log(sensitivity["window_side_length"].to_numpy(dtype=float))
    y = sensitivity["retained_parent_unit_fraction"].to_numpy(dtype=float)
    x0, y0, x1, y1 = x[0], y[0], x[-1], y[-1]
    denominator = float(np.hypot(y1 - y0, x1 - x0))
    if denominator == 0:
        distances = np.zeros_like(x)
    else:
        distances = np.abs((y1 - y0) * x - (x1 - x0) * y + x1 * y0 - y1 * x0) / denominator
    sensitivity["knee_distance"] = distances
    maximum = float(sensitivity["knee_distance"].max())
    selected_side = float(
        sensitivity.loc[np.isclose(sensitivity["knee_distance"], maximum), "window_side_length"].min()
    )
    return {
        "status": "ok",
        "window_side_length": selected_side,
        "min_parent_units": min_parent_units,
    }, sensitivity


def compute_rarefied_window_diversity(
    window_assignments: pd.DataFrame,
    raw_labels: pd.Series,
    reconstructed_labels: pd.Series,
    *,
    raw_level2_labels: pd.Series | None = None,
    min_parent_units: int = 4,
    n_draws: int = 200,
    random_state: int = 42,
) -> pd.DataFrame:
    """Estimate local partition diversity from identical within-window draws."""
    if min_parent_units < 1 or n_draws < 1:
        raise ValueError("min_parent_units and n_draws must be at least one")
    if "window_id" not in window_assignments:
        raise KeyError("window_assignments must contain window_id")
    raw = _align_labels(raw_labels, window_assignments.index, "Raw")
    recon = _align_labels(reconstructed_labels, window_assignments.index, "Reconstructed")
    level2 = (
        _align_labels(raw_level2_labels, window_assignments.index, "Raw Level2")
        if raw_level2_labels is not None
        else None
    )
    raw_codes = pd.Series(pd.factorize(raw, sort=True)[0], index=raw.index)
    recon_codes = pd.Series(pd.factorize(recon, sort=True)[0], index=recon.index)
    level2_codes = (
        pd.Series(pd.factorize(level2, sort=True)[0], index=level2.index)
        if level2 is not None
        else None
    )

    def composition_from_codes(codes: np.ndarray) -> tuple[float, float, float, float]:
        counts = np.bincount(codes)
        probabilities = counts[counts > 0] / codes.size
        k_obs = float(probabilities.size)
        entropy = float(-(probabilities * np.log(probabilities)).sum())
        neff = float(np.exp(entropy))
        return k_obs, entropy, neff, float(neff / k_obs)

    generator = np.random.default_rng(random_state)
    rows: list[dict[str, float | int | bool | str]] = []
    for fallback_index, (window_id, frame) in enumerate(
        window_assignments.groupby("window_id", sort=True)
    ):
        raw_values_for_window = raw_codes.reindex(frame.index).to_numpy(dtype=int)
        recon_values_for_window = recon_codes.reindex(frame.index).to_numpy(dtype=int)
        level2_values_for_window = (
            level2_codes.reindex(frame.index).to_numpy(dtype=int)
            if level2_codes is not None
            else None
        )
        n_units = int(frame.shape[0])
        row: dict[str, float | int | bool | str] = {
            "window_id": str(window_id),
            "window_x": float(frame["x"].mean()),
            "window_y": float(frame["y"].mean()),
            "window_x_index": int(frame["window_x_index"].iloc[0])
            if "window_x_index" in frame
            else fallback_index,
            "window_y_index": int(frame["window_y_index"].iloc[0])
            if "window_y_index" in frame
            else 0,
            "n_units": n_units,
            "valid_window": n_units >= min_parent_units,
        }
        if n_units < min_parent_units:
            for metric in ("k_obs", "entropy", "neff", "evenness"):
                for suffix in ("raw", "recon", "sd_raw", "sd_recon"):
                    row[f"{metric}_{suffix}"] = np.nan
                row[f"delta_{metric}"] = np.nan
                row[f"delta_{metric}_sd"] = np.nan
                row[f"{metric}_level2"] = np.nan
                row[f"{metric}_sd_level2"] = np.nan
                row[f"delta_{metric}_vs_raw_leiden"] = np.nan
                row[f"delta_{metric}_vs_raw_leiden_sd"] = np.nan
                row[f"delta_{metric}_vs_raw_level2"] = np.nan
                row[f"delta_{metric}_vs_raw_level2_sd"] = np.nan
            rows.append(row)
            continue
        raw_values: list[tuple[float, float, float, float]] = []
        recon_values: list[tuple[float, float, float, float]] = []
        level2_values: list[tuple[float, float, float, float]] = []
        for _ in range(n_draws):
            sampled = generator.choice(n_units, size=min_parent_units, replace=False)
            raw_values.append(composition_from_codes(raw_values_for_window[sampled]))
            recon_values.append(composition_from_codes(recon_values_for_window[sampled]))
            if level2_values_for_window is not None:
                level2_values.append(composition_from_codes(level2_values_for_window[sampled]))
        raw_array = np.asarray(raw_values, dtype=float)
        recon_array = np.asarray(recon_values, dtype=float)
        level2_array = np.asarray(level2_values, dtype=float) if level2_values else None
        for index, metric in enumerate(("k_obs", "entropy", "neff", "evenness")):
            delta = recon_array[:, index] - raw_array[:, index]
            row.update(
                {
                    f"{metric}_raw": float(raw_array[:, index].mean()),
                    f"{metric}_recon": float(recon_array[:, index].mean()),
                    f"delta_{metric}": float(delta.mean()),
                    f"{metric}_sd_raw": float(raw_array[:, index].std(ddof=0)),
                    f"{metric}_sd_recon": float(recon_array[:, index].std(ddof=0)),
                    f"delta_{metric}_sd": float(delta.std(ddof=0)),
                    f"delta_{metric}_vs_raw_leiden": float(delta.mean()),
                    f"delta_{metric}_vs_raw_leiden_sd": float(delta.std(ddof=0)),
                }
            )
            if level2_array is None:
                row[f"{metric}_level2"] = np.nan
                row[f"{metric}_sd_level2"] = np.nan
                row[f"delta_{metric}_vs_raw_level2"] = np.nan
                row[f"delta_{metric}_vs_raw_level2_sd"] = np.nan
            else:
                level2_delta = recon_array[:, index] - level2_array[:, index]
                row[f"{metric}_level2"] = float(level2_array[:, index].mean())
                row[f"{metric}_sd_level2"] = float(level2_array[:, index].std(ddof=0))
                row[f"delta_{metric}_vs_raw_level2"] = float(level2_delta.mean())
                row[f"delta_{metric}_vs_raw_level2_sd"] = float(level2_delta.std(ddof=0))
        rows.append(row)
    return pd.DataFrame(rows)


def select_region_threshold(
    values: Sequence[float] | np.ndarray,
    *,
    n_bootstrap: int = 500,
    random_state: int = 42,
) -> tuple[dict[str, float | int | str | None], pd.DataFrame]:
    """Find a stable survival-curve breakpoint, or explicitly decline a Region mask."""
    observed = np.asarray(values, dtype=float)
    observed = observed[np.isfinite(observed)]
    if observed.size < 200 or n_bootstrap < 1:
        return {
            "status": "no_stable_threshold",
            "threshold": None,
            "n_windows": int(observed.size),
            "n_valid_bootstrap": 0,
        }, pd.DataFrame()

    def breakpoint(sample: np.ndarray) -> float | None:
        unique = np.unique(np.quantile(sample, np.linspace(0.0, 1.0, 81)))
        if unique.size < 5:
            return None
        survival = np.asarray([(sample >= value).mean() for value in unique])
        low = int(np.ceil(unique.size * 0.10))
        high = int(np.floor(unique.size * 0.90))
        candidates = range(max(1, low), min(unique.size - 1, high) + 1)
        best: tuple[float, float] | None = None
        for index in candidates:
            left_x, left_y = unique[: index + 1], np.log(survival[: index + 1])
            right_x, right_y = unique[index:], np.log(survival[index:])
            if left_x.size < 2 or right_x.size < 2:
                continue
            left_fit = np.polyfit(left_x, left_y, deg=1)
            right_fit = np.polyfit(right_x, right_y, deg=1)
            error = float(((left_y - np.polyval(left_fit, left_x)) ** 2).sum() + ((right_y - np.polyval(right_fit, right_x)) ** 2).sum())
            candidate = (error, float(unique[index]))
            if best is None or candidate < best:
                best = candidate
        return None if best is None else best[1]

    point = breakpoint(observed)
    if point is None:
        return {
            "status": "no_stable_threshold",
            "threshold": None,
            "n_windows": int(observed.size),
            "n_valid_bootstrap": 0,
        }, pd.DataFrame()
    generator = np.random.default_rng(random_state)
    bootstrap_values = [breakpoint(generator.choice(observed, size=observed.size, replace=True)) for _ in range(n_bootstrap)]
    valid = np.asarray([value for value in bootstrap_values if value is not None], dtype=float)
    audit = pd.DataFrame({"bootstrap_threshold": valid})
    if valid.size < int(np.ceil(n_bootstrap * 0.80)):
        return {
            "status": "no_stable_threshold",
            "threshold": None,
            "n_windows": int(observed.size),
            "n_valid_bootstrap": int(valid.size),
        }, audit
    lower, upper = np.quantile(valid, [0.025, 0.975])
    value_range = float(observed.max() - observed.min())
    stable = value_range > 0 and (upper - lower) <= 0.25 * value_range
    return {
        "status": "ok" if stable else "no_stable_threshold",
        "threshold": float(np.median(valid)) if stable else None,
        "point_threshold": float(point),
        "ci_lower": float(lower),
        "ci_upper": float(upper),
        "n_windows": int(observed.size),
        "n_valid_bootstrap": int(valid.size),
    }, audit


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
                "window_x_index": int(frame["window_x_index"].iloc[0])
                if "window_x_index" in frame
                else 0,
                "window_y_index": int(frame["window_y_index"].iloc[0])
                if "window_y_index" in frame
                else 0,
                "tumor_units": tumor_units,
                "normal_units": normal_units,
                "tumor_candidate": tumor_candidate,
                "normal_candidate": normal_candidate,
                "union_candidate": union_candidate,
                "level1_region": level1_region,
            }
        )
    return pd.DataFrame(rows)


def summarize_anatomy_context(
    anatomy_unit_assignments: pd.DataFrame,
    anatomy_windows: pd.DataFrame,
    *,
    window_side_length: float,
) -> pd.DataFrame:
    """Summarize full-Level1 units and tissue area for each anatomy context."""
    if "window_id" not in anatomy_unit_assignments:
        raise KeyError("anatomy_unit_assignments must contain window_id")
    if not {"window_id", "level1_region"} <= set(anatomy_windows.columns):
        raise KeyError("anatomy_windows must contain window_id and level1_region")
    if anatomy_windows["window_id"].duplicated().any():
        raise ValueError("anatomy_windows must contain one row per window")
    if window_side_length <= 0:
        raise ValueError("window_side_length must be positive")

    region_by_window = anatomy_windows.set_index("window_id")["level1_region"]
    unit_regions = anatomy_unit_assignments["window_id"].map(region_by_window)
    if unit_regions.isna().any():
        raise ValueError("Every full-Level1 unit must map to an anatomy window")
    total_windows = int(anatomy_windows.shape[0])
    window_area_um2 = float(window_side_length**2)
    rows = []
    for region, frame in anatomy_windows.groupby("level1_region", sort=True):
        area_um2 = float(frame.shape[0] * window_area_um2)
        rows.append(
            {
                "level1_region": region,
                "full_level1_units": int((unit_regions == region).sum()),
                "tissue_windows": int(frame.shape[0]),
                "area_um2": area_um2,
                "area_mm2": area_um2 / 1_000_000.0,
                "area_fraction": float(frame.shape[0] / total_windows) if total_windows else np.nan,
            }
        )
    return pd.DataFrame(rows)


def _distribution(values: pd.Series, prefix: str) -> dict[str, float]:
    numeric = pd.to_numeric(values, errors="coerce").dropna()
    if numeric.empty:
        return {
            f"median_{prefix}": np.nan,
            f"q1_{prefix}": np.nan,
            f"q3_{prefix}": np.nan,
        }
    return {
        f"median_{prefix}": float(numeric.median()),
        f"q1_{prefix}": float(numeric.quantile(0.25)),
        f"q3_{prefix}": float(numeric.quantile(0.75)),
    }


def _region_frames(frame: pd.DataFrame) -> list[tuple[str, pd.DataFrame]]:
    return [("Overall", frame)] + [
        (str(region), group)
        for region, group in frame.groupby("level1_region", sort=True)
    ]


def summarize_cluster_change_by_anatomy(
    unit_assignments: pd.DataFrame,
    window_metrics: pd.DataFrame,
) -> pd.DataFrame:
    """Summarize globally matched unit changes and their window distribution."""
    required_units = {"level1_region", "unit_changed"}
    required_windows = {
        "level1_region",
        "valid_window",
        "unit_change_fraction",
    }
    if missing := required_units - set(unit_assignments.columns):
        raise KeyError(f"Unit assignments are missing columns: {sorted(missing)}")
    if missing := required_windows - set(window_metrics.columns):
        raise KeyError(f"Window metrics are missing columns: {sorted(missing)}")

    unit_groups = dict(_region_frames(unit_assignments))
    window_groups = dict(_region_frames(window_metrics))
    rows = []
    for region, units in unit_groups.items():
        windows = window_groups.get(region, window_metrics.iloc[0:0])
        valid = windows.loc[windows["valid_window"]]
        changed_units = int(units["unit_changed"].astype(bool).sum())
        paired_units = int(units.shape[0])
        rows.append(
            {
                "level1_region": region,
                "paired_units": paired_units,
                "changed_units": changed_units,
                "change_fraction": float(changed_units / paired_units) if paired_units else np.nan,
                "n_valid_windows": int(valid.shape[0]),
                **_distribution(valid["unit_change_fraction"], "window_change"),
            }
        )
    return pd.DataFrame(rows)


def summarize_diversity_by_anatomy(window_metrics: pd.DataFrame) -> pd.DataFrame:
    """Summarize every available diversity metric and Raw baseline with IQR."""
    required = {"level1_region", "valid_window"}
    if missing := required - set(window_metrics.columns):
        raise KeyError(f"Window metrics are missing columns: {sorted(missing)}")
    metrics = ("k_obs", "entropy", "neff", "evenness")
    source_columns = []
    for metric in metrics:
        source_columns.extend(
            [
                f"{metric}_raw",
                f"{metric}_recon",
                f"delta_{metric}",
                f"{metric}_level2",
                f"delta_{metric}_vs_raw_level2",
            ]
        )
    source_columns = [column for column in source_columns if column in window_metrics]
    if not source_columns:
        raise KeyError("Window metrics contain no diversity metric columns")
    rows = []
    for region, frame in _region_frames(window_metrics):
        valid = frame.loc[frame["valid_window"]]
        row = {"level1_region": region, "n_valid_windows": int(valid.shape[0])}
        for column in source_columns:
            row.update(_distribution(valid[column], column))
        rows.append(row)
    return pd.DataFrame(rows)


def summarize_region_extent_by_anatomy(
    window_metrics: pd.DataFrame,
    *,
    window_side_length: float,
) -> pd.DataFrame:
    """Summarize high-diversity Region area and unit coverage by anatomy."""
    required = {"level1_region", "valid_window", "in_region", "n_units"}
    if missing := required - set(window_metrics.columns):
        raise KeyError(f"Window metrics are missing columns: {sorted(missing)}")
    if window_side_length <= 0:
        raise ValueError("window_side_length must be positive")
    window_area_um2 = float(window_side_length**2)
    rows = []
    for region, frame in _region_frames(window_metrics):
        valid = frame.loc[frame["valid_window"]]
        available = valid["in_region"].notna().any()
        selected = valid.loc[valid["in_region"].fillna(False).astype(bool)]
        valid_units = int(valid["n_units"].sum())
        region_units = int(selected["n_units"].sum()) if available else np.nan
        region_area_um2 = float(selected.shape[0] * window_area_um2) if available else np.nan
        rows.append(
            {
                "level1_region": region,
                "valid_windows": int(valid.shape[0]),
                "region_windows": int(selected.shape[0]) if available else np.nan,
                "region_area_um2": region_area_um2,
                "region_area_mm2": region_area_um2 / 1_000_000.0,
                "area_fraction": float(selected.shape[0] / valid.shape[0]) if available and not valid.empty else np.nan,
                "valid_units": valid_units,
                "region_units": region_units,
                "unit_fraction": float(region_units / valid_units) if available and valid_units else np.nan,
                "region_available": bool(available),
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
