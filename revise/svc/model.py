from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Optional

import pandas as pd
from anndata import AnnData


@dataclass
class SVC:
    """Canonical result carrier for all unified REVISE runs."""

    expr: Optional[AnnData]
    spatial: Optional[AnnData]
    svc_kind: str
    cell_type_probs: Optional[pd.DataFrame] = None
    cell_type_label: Optional[pd.Series] = None
    confidence: Optional[pd.Series] = None
    provenance: Dict[str, Any] = field(default_factory=dict)
    quality_metrics: Dict[str, Any] = field(default_factory=dict)
    artifacts: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_obs(self) -> int:
        if self.expr is not None:
            return int(self.expr.n_obs)
        if self.spatial is not None:
            return int(self.spatial.n_obs)
        return 0

    @property
    def n_vars(self) -> int:
        if self.expr is not None:
            return int(self.expr.n_vars)
        if self.spatial is not None:
            return int(self.spatial.n_vars)
        return 0

    def summary(self) -> Dict[str, Any]:
        return {
            "svc_kind": self.svc_kind,
            "n_obs": self.n_obs,
            "n_vars": self.n_vars,
            "quality_metric_keys": sorted(self.quality_metrics.keys()),
            "artifact_keys": sorted(self.artifacts.keys()),
        }
