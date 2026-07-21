from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Optional

from anndata import AnnData


@dataclass
class REVISEDataBundle:
    """Algorithm-facing input bundle.

    REVISE algorithms operate on AnnData matrices. Optional SpatialData objects
    stay attached only as provenance/service-layer context so downstream code
    can write results back to the spatial ecosystem without making SpatialData
    a mandatory runtime dependency.
    """

    st_adata: AnnData
    sc_ref_adata: AnnData
    real_st_adata: Optional[AnnData] = None
    sdata: Any = None
    spatial_unit: Optional[str] = None
    coordinate_system: Optional[str] = None
    source_report: Dict[str, Any] = field(default_factory=dict)
