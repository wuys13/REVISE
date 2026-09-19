"""Summarize transient Global Anchoring distributions for reference ranking.

Adapted from ReviseSTAtlas ``reference_selection/scoring.py`` at source commit
``fe4ac4a769a0518854e14be6426f1ad0bf327e95``.  The numerical implementation
is preserved from the portable source bundle.
"""

from __future__ import annotations

import numpy as np

from revise.utils.confidence import entropy_certainty
from revise.utils.confidence import max_confidence
from revise.utils.confidence import top2_margin


SCORE_METHODS = {
    "max_median": "max_confidence",
    "certainty_median": "normalized_certainty",
    "margin_median": "top2_margin",
}


def row_metrics(matrix: np.ndarray) -> dict[str, np.ndarray]:
    return {
        "max_confidence": max_confidence(matrix),
        "normalized_certainty": entropy_certainty(matrix),
        "top2_margin": top2_margin(matrix),
    }


def summarize_scores(metrics: dict[str, np.ndarray]) -> dict[str, float]:
    return {
        method: float(np.median(metrics[column]))
        for method, column in SCORE_METHODS.items()
    }
