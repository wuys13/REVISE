from revise.analysis.biological_metrics import compute_asw
from revise.analysis.biological_metrics import compute_conditional_moran_i
from revise.analysis.biological_metrics import compute_identity_metrics
from revise.analysis.biological_metrics import compute_local_label_entropy
from revise.analysis.biological_metrics import compute_tmp_mer
from revise.analysis.biological_metrics import make_cell_type_mean_baseline
from revise.analysis.biological_metrics import shannon_entropy_from_labels
from revise.analysis.biological_metrics import summarize_conditional_moran_i
from revise.analysis.biological_metrics import summarize_tmp_mer
from revise.analysis.metrics import compute_clustering_metrics
from revise.analysis.metrics import compute_metric

__all__ = [
    "compute_metric",
    "compute_clustering_metrics",
    "compute_conditional_moran_i",
    "summarize_conditional_moran_i",
    "compute_tmp_mer",
    "summarize_tmp_mer",
    "compute_local_label_entropy",
    "compute_identity_metrics",
    "compute_asw",
    "make_cell_type_mean_baseline",
    "shannon_entropy_from_labels",
    "ScSVCAnalysisService",
    "SpSVCAnalysisService",
]


def __getattr__(name):
    if name in {"ScSVCAnalysisService", "SpSVCAnalysisService"}:
        from revise.analysis.services import ScSVCAnalysisService, SpSVCAnalysisService

        return {"ScSVCAnalysisService": ScSVCAnalysisService, "SpSVCAnalysisService": SpSVCAnalysisService}[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
