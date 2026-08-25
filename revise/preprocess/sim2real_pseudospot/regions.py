"""Candidate-region selection for Sim2Real-ST."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING, Mapping

import numpy as np
import pandas as pd
from anndata import AnnData
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from scipy.spatial.distance import jensenshannon
from sklearn.neighbors import kneighbors_graph

if TYPE_CHECKING:
    from revise.preprocess.sim2real_pseudospot.workflow import ProposalConfig


Bounds = tuple[float, float, float, float]


@dataclass(frozen=True)
class Candidate:
    candidate_id: str
    role: str
    bounds: Bounds
    total_score: float
    cell_count: int = 0
    composition_score: float = 0.0
    simplicity_score: float = 0.0
    continuity_score: float = 0.0
    coverage_score: float = 0.0
    density: float = 0.0
    composition: tuple[tuple[str, float], ...] = ()


def filter_region_cells(
    cells: AnnData,
    *,
    bounds: Bounds,
    label_key: str,
    unknown_label: str,
    min_cells_per_type: int,
) -> AnnData:
    """Select one rectangle, then apply the final label-retention rule."""
    if label_key not in cells.obs:
        raise KeyError(f"Input cells are missing obs[{label_key!r}].")
    if "spatial" not in cells.obsm:
        raise KeyError("Input cells are missing obsm['spatial'].")
    xmin, xmax, ymin, ymax = bounds
    coordinates = np.asarray(cells.obsm["spatial"])
    inside = (
        (coordinates[:, 0] >= xmin)
        & (coordinates[:, 0] <= xmax)
        & (coordinates[:, 1] >= ymin)
        & (coordinates[:, 1] <= ymax)
    )
    selected = cells[inside].copy()
    labels = selected.obs[label_key].astype(str)
    selected = selected[labels != unknown_label].copy()
    counts = selected.obs[label_key].astype(str).value_counts()
    retained_labels = counts[counts > min_cells_per_type].index
    return selected[selected.obs[label_key].astype(str).isin(retained_labels)].copy()


def _iou(left: Bounds, right: Bounds) -> float:
    left_xmin, left_xmax, left_ymin, left_ymax = left
    right_xmin, right_xmax, right_ymin, right_ymax = right
    width = max(0.0, min(left_xmax, right_xmax) - max(left_xmin, right_xmin))
    height = max(0.0, min(left_ymax, right_ymax) - max(left_ymin, right_ymin))
    intersection = width * height
    left_area = (left_xmax - left_xmin) * (left_ymax - left_ymin)
    right_area = (right_xmax - right_xmin) * (right_ymax - right_ymin)
    union = left_area + right_area - intersection
    return 0.0 if union <= 0 else intersection / union


def suppress_overlaps(
    candidates: list[Candidate], *, max_iou: float, limit: int
) -> list[Candidate]:
    """Keep highest-scoring candidates whose rectangles do not overlap too much."""
    selected: list[Candidate] = []
    ordered = sorted(
        candidates,
        key=lambda candidate: (-candidate.total_score, candidate.candidate_id),
    )
    for candidate in ordered:
        if all(_iou(candidate.bounds, chosen.bounds) <= max_iou for chosen in selected):
            selected.append(candidate)
        if len(selected) == limit:
            break
    return selected


def _candidate_windows(coordinates: np.ndarray, config: "ProposalConfig") -> list[Bounds]:
    xmin, ymin = coordinates.min(axis=0)
    xmax, ymax = coordinates.max(axis=0)
    windows: list[Bounds] = []
    seen: set[Bounds] = set()
    for scale in config.scales:
        for width, height in (
            (config.base_width * scale, config.base_height * scale),
            (config.base_height * scale, config.base_width * scale),
        ):
            x_start = np.floor(xmin / config.step) * config.step
            y_start = np.floor(ymin / config.step) * config.step
            x_stop = np.ceil(xmax / config.step) * config.step
            y_stop = np.ceil(ymax / config.step) * config.step
            for left in np.arange(x_start, x_stop + config.step, config.step):
                for top in np.arange(y_start, y_stop + config.step, config.step):
                    bounds = (float(left), float(left + width), float(top), float(top + height))
                    if bounds not in seen:
                        windows.append(bounds)
                        seen.add(bounds)
    return windows


def _label_counts(
    labels: np.ndarray,
    mask: np.ndarray,
    *,
    unknown_label: str,
    min_cells_per_type: int,
) -> pd.Series:
    selected = labels[mask]
    counts = pd.Series(selected[selected != unknown_label]).value_counts()
    return counts[counts > min_cells_per_type]


def _core_continuity(labels: np.ndarray, coordinates: np.ndarray) -> np.ndarray:
    if labels.size < 2:
        return np.zeros(labels.size, dtype=float)
    graph = kneighbors_graph(
        coordinates,
        n_neighbors=min(8, labels.size - 1),
        mode="connectivity",
        include_self=False,
    ).tocsr()
    scores = np.zeros(labels.size, dtype=float)
    for index in range(labels.size):
        neighbors = graph.indices[graph.indptr[index] : graph.indptr[index + 1]]
        scores[index] = np.mean(labels[neighbors] == labels[index])
    return scores


def _leading_edge_contacts(labels: np.ndarray, coordinates: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    tumor_scores = np.zeros(labels.size, dtype=float)
    fibroblast_scores = np.zeros(labels.size, dtype=float)
    if labels.size < 2:
        return tumor_scores, fibroblast_scores
    graph = kneighbors_graph(
        coordinates,
        n_neighbors=min(8, labels.size - 1),
        mode="connectivity",
        include_self=False,
    ).tocsr()
    for index in range(labels.size):
        neighbors = graph.indices[graph.indptr[index] : graph.indptr[index + 1]]
        if labels[index] == "Tumor":
            tumor_scores[index] = float(np.any(labels[neighbors] == "Fibroblast"))
        elif labels[index] == "Fibroblast":
            fibroblast_scores[index] = float(np.any(labels[neighbors] == "Tumor"))
    return tumor_scores, fibroblast_scores


def _composition(
    counts: pd.Series, categories: list[str]
) -> tuple[dict[str, float], np.ndarray]:
    total = float(counts.sum())
    values = np.array([float(counts.get(category, 0.0)) / total for category in categories])
    return dict(zip(categories, values.tolist())), values


def _continuity_for_window(
    role: str,
    retained: np.ndarray,
    labels: np.ndarray,
    core_scores: np.ndarray,
    tumor_scores: np.ndarray,
    fibroblast_scores: np.ndarray,
) -> float:
    if role != "leading_edge":
        return float(core_scores[retained].mean())
    tumor = retained & (labels == "Tumor")
    fibroblast = retained & (labels == "Fibroblast")
    if not tumor.any() or not fibroblast.any():
        return 0.0
    return float(np.sqrt(tumor_scores[tumor].mean() * fibroblast_scores[fibroblast].mean()))


def propose_candidates(
    cells: AnnData,
    *,
    role: str,
    template_composition: Mapping[str, float],
    config: "ProposalConfig",
    label_key: str,
    unknown_label: str,
    min_cells_per_type: int,
) -> list[Candidate]:
    """Score and return up to three non-overlapping candidate rectangles."""
    if role not in {"leading_edge", "normal_core", "tumor_core"}:
        raise ValueError(f"Unsupported region role: {role}")
    coordinates = np.asarray(cells.obsm["spatial"], dtype=float)
    labels = cells.obs[label_key].astype(str).to_numpy()
    categories = sorted(set(template_composition) | set(labels[labels != unknown_label]))
    if not categories:
        raise ValueError("No retained Level 1 categories are available for proposals.")
    template = np.array([float(template_composition.get(category, 0.0)) for category in categories])
    template /= template.sum()
    core_scores = _core_continuity(labels, coordinates)
    tumor_scores, fibroblast_scores = _leading_edge_contacts(labels, coordinates)

    candidates: list[Candidate] = []
    for index, bounds in enumerate(_candidate_windows(coordinates, config), start=1):
        xmin, xmax, ymin, ymax = bounds
        inside = (
            (coordinates[:, 0] >= xmin)
            & (coordinates[:, 0] <= xmax)
            & (coordinates[:, 1] >= ymin)
            & (coordinates[:, 1] <= ymax)
        )
        counts = _label_counts(
            labels,
            inside,
            unknown_label=unknown_label,
            min_cells_per_type=min_cells_per_type,
        )
        cell_count = int(counts.sum())
        if cell_count < config.min_cells:
            continue
        retained = inside & np.isin(labels, counts.index.to_numpy())
        composition, values = _composition(counts, categories)
        composition_score = float(1.0 - jensenshannon(values, template, base=2.0) ** 2)
        active_categories = int(np.count_nonzero(values))
        simplicity_score = (
            1.0
            if active_categories <= 1
            else float(1.0 + np.sum(values[values > 0] * np.log(values[values > 0])) / np.log(active_categories))
        )
        area = (xmax - xmin) * (ymax - ymin)
        candidates.append(
            Candidate(
                candidate_id=f"{role}-raw-{index}",
                role=role,
                bounds=bounds,
                total_score=0.0,
                cell_count=cell_count,
                composition_score=composition_score,
                simplicity_score=simplicity_score,
                continuity_score=_continuity_for_window(
                    role,
                    retained,
                    labels,
                    core_scores,
                    tumor_scores,
                    fibroblast_scores,
                ),
                density=cell_count / area,
                composition=tuple(composition.items()),
            )
        )
    if not candidates:
        return []

    density_reference = float(np.percentile([candidate.density for candidate in candidates], 95))
    scored = []
    for candidate in candidates:
        coverage_score = min(candidate.density / density_reference, 1.0) if density_reference else 0.0
        total_score = (
            0.5 * candidate.composition_score
            + 0.2 * candidate.simplicity_score
            + 0.2 * candidate.continuity_score
            + 0.1 * coverage_score
        )
        scored.append(
            replace(
                candidate,
                coverage_score=coverage_score,
                total_score=total_score,
            )
        )
    selected = suppress_overlaps(scored, max_iou=config.max_iou, limit=3)
    return [
        replace(candidate, candidate_id=f"{role}-{index}")
        for index, candidate in enumerate(selected, start=1)
    ]


def candidate_table(candidates_by_role: Mapping[str, list[Candidate]]) -> pd.DataFrame:
    """Return the human-review table with geometry, scores, and proportions."""
    rows: list[dict[str, float | int | str]] = []
    for role, candidates in candidates_by_role.items():
        for candidate in candidates:
            xmin, xmax, ymin, ymax = candidate.bounds
            row: dict[str, float | int | str] = {
                "role": role,
                "candidate_id": candidate.candidate_id,
                "xmin": xmin,
                "xmax": xmax,
                "ymin": ymin,
                "ymax": ymax,
                "width": xmax - xmin,
                "height": ymax - ymin,
                "cell_count": candidate.cell_count,
                "composition_score": candidate.composition_score,
                "simplicity_score": candidate.simplicity_score,
                "continuity_score": candidate.continuity_score,
                "coverage_score": candidate.coverage_score,
                "total_score": candidate.total_score,
            }
            row.update({f"proportion_{label}": value for label, value in candidate.composition})
            rows.append(row)
    return pd.DataFrame(rows)


def plot_candidates(
    cells: AnnData,
    *,
    label_key: str,
    candidates_by_role: Mapping[str, list[Candidate]],
    output_path: str | Path,
) -> None:
    """Plot Level 1 spatial coordinates and proposal rectangles for review."""
    coordinates = np.asarray(cells.obsm["spatial"], dtype=float)
    labels = cells.obs[label_key].astype(str)
    roles = ("tumor_core", "normal_core", "leading_edge")
    palette = {
        label: color
        for label, color in zip(sorted(labels.unique()), plt.get_cmap("tab20").colors)
    }
    figure, axes = plt.subplots(1, 3, figsize=(24, 8), sharex=True, sharey=True)
    for axis, role in zip(axes, roles):
        for label, color in palette.items():
            mask = labels == label
            axis.scatter(
                coordinates[mask, 0],
                coordinates[mask, 1],
                s=1,
                color=color,
                label=label,
                rasterized=True,
            )
        for candidate in candidates_by_role[role]:
            xmin, xmax, ymin, ymax = candidate.bounds
            axis.add_patch(
                Rectangle(
                    (xmin, ymin),
                    xmax - xmin,
                    ymax - ymin,
                    fill=False,
                    edgecolor="red",
                    linewidth=2,
                )
            )
            axis.text(xmin, ymin, candidate.candidate_id, color="red", fontsize=9)
        axis.set_title(role.replace("_", " "))
        axis.set_aspect("equal")
        axis.invert_yaxis()
    axes[0].legend(markerscale=5, bbox_to_anchor=(1.02, 1), loc="upper left")
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def plot_selected_region(
    cells: AnnData,
    *,
    label_key: str,
    output_path: str | Path,
) -> None:
    """Save the final confirmed cut in the legacy-compatible cut.png form."""
    coordinates = np.asarray(cells.obsm["spatial"], dtype=float)
    labels = cells.obs[label_key].astype(str)
    figure, axis = plt.subplots(figsize=(10, 10))
    for label, color in zip(sorted(labels.unique()), plt.get_cmap("tab20").colors):
        mask = labels == label
        axis.scatter(coordinates[mask, 0], coordinates[mask, 1], s=2, color=color, label=label)
    axis.set_aspect("equal")
    axis.invert_yaxis()
    axis.legend(markerscale=4)
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)
