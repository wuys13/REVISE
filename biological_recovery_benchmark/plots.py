from __future__ import annotations

from collections.abc import Mapping, Sequence
from itertools import cycle

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes
from matplotlib.figure import Figure

_COLORS = ["#0077BB", "#33BBEE", "#009988", "#EE7733", "#CC3311"]


def _figure_axis(ax: Axes | None, figsize: tuple[float, float]) -> tuple[Figure, Axes]:
    if ax is None:
        return plt.subplots(figsize=figsize)
    return ax.figure, ax


def _method_palette(
    data: pd.DataFrame,
    method_col: str,
    method_order: Sequence[str] | None,
    palette: Mapping[str, str] | Sequence[str] | None,
) -> Mapping[str, str] | Sequence[str]:
    if palette is not None:
        return palette
    methods = (
        list(method_order)
        if method_order is not None
        else list(pd.unique(data[method_col]))
    )
    colors = cycle(_COLORS)
    return {str(method): next(colors) for method in methods}


def plot_metric_comparison(
    data: pd.DataFrame,
    *,
    method_col: str = "method",
    metric_col: str = "metric",
    value_col: str = "value",
    method_order: Sequence[str] | None = None,
    palette: Mapping[str, str] | Sequence[str] | None = None,
    ax: Axes | None = None,
) -> tuple[Figure, Axes]:
    """Plot already-computed scalar metrics across methods."""
    fig, axis = _figure_axis(ax, (7.0, 4.0))
    sns.stripplot(
        data=data,
        x=metric_col,
        y=value_col,
        hue=method_col,
        hue_order=method_order,
        dodge=True,
        jitter=False,
        size=7,
        palette=_method_palette(data, method_col, method_order, palette),
        ax=axis,
    )
    axis.set_xlabel("")
    axis.set_ylabel("Metric value")
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    return fig, axis


def plot_spatial_metric_comparison(
    data: pd.DataFrame,
    *,
    metric: str,
    method_col: str = "method",
    method_order: Sequence[str] | None = None,
    palette: Mapping[str, str] | Sequence[str] | None = None,
    ax: Axes | None = None,
) -> tuple[Figure, Axes]:
    """Show the descriptive per-gene distribution of MISC, MIDC or Moran's I."""
    if metric not in data:
        raise KeyError(f"Column not found: {metric!r}")
    fig, axis = _figure_axis(ax, (6.0, 4.0))
    sns.boxplot(
        data=data,
        x=method_col,
        y=metric,
        order=method_order,
        hue=method_col,
        hue_order=method_order,
        palette=_method_palette(data, method_col, method_order, palette),
        legend=False,
        showfliers=False,
        ax=axis,
    )
    axis.set_xlabel("")
    axis.set_ylabel(metric)
    axis.tick_params(axis="x", rotation=30)
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    return fig, axis


def plot_moran_heatmap(
    data: pd.DataFrame,
    *,
    method_col: str = "method",
    group_col: str = "group",
    value_col: str = "MoranI",
    statistic: str = "mean",
    ax: Axes | None = None,
) -> tuple[Figure, Axes]:
    """Plot a method-by-group heatmap from precomputed gene-level Moran values."""
    if statistic not in {"mean", "median"}:
        raise ValueError("statistic must be 'mean' or 'median'")
    grouped = data.groupby([group_col, method_col], observed=True)[value_col]
    summary = grouped.mean() if statistic == "mean" else grouped.median()
    matrix = summary.unstack(method_col)
    fig, axis = _figure_axis(
        ax, (max(5.0, 1.1 * matrix.shape[1]), max(3.0, 0.35 * matrix.shape[0]))
    )
    sns.heatmap(
        matrix,
        cmap="viridis",
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": f"{statistic.capitalize()} Moran's I"},
        ax=axis,
    )
    axis.set_xlabel("")
    axis.set_ylabel("")
    fig.tight_layout()
    return fig, axis
