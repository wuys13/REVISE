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
    "SpSVCAnalysisService",
]


_BIOLOGICAL_EXPORTS = {
    "compute_asw",
    "compute_conditional_moran_i",
    "compute_identity_metrics",
    "compute_local_label_entropy",
    "compute_tmp_mer",
    "make_cell_type_mean_baseline",
    "shannon_entropy_from_labels",
    "summarize_conditional_moran_i",
    "summarize_tmp_mer",
}


def __getattr__(name):
    if name in {"compute_metric", "compute_clustering_metrics"}:
        from revise.analysis import metrics

        return getattr(metrics, name)
    if name in _BIOLOGICAL_EXPORTS:
        from revise.analysis import biological_metrics

        return getattr(biological_metrics, name)
    if name == "SpSVCAnalysisService":
        from revise.analysis.services import SpSVCAnalysisService

        return SpSVCAnalysisService
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
