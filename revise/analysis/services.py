from __future__ import annotations

from typing import Dict

from revise.analysis.metrics import compute_clustering_metrics
from revise.svc import SVC


class SpSVCAnalysisService:
    """Basic downstream metrics for sp-SVC outputs."""

    def __init__(self, svc: SVC) -> None:
        if svc.spatial is None:
            raise ValueError("SpSVCAnalysisService requires svc.spatial")
        self.svc = svc

    def clustering_metrics(self, pred_col: str, ref_col: str) -> Dict[str, float]:
        ari, nmi = compute_clustering_metrics(self.svc.spatial, pred_col, ref_col)
        return {"ari": float(ari), "nmi": float(nmi)}
