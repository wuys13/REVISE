import os
from importlib import import_module
from typing import TYPE_CHECKING

# Lazily import heavy modules to avoid ImportError during package import when
# optional dependencies are missing in some environments.
if TYPE_CHECKING:  # pragma: no cover - only for type checkers
    from revise.benchmark.sc_svc_impute import ScSVCImpute
    from revise.benchmark.sc_svc_sr import ScSVCSr
    from revise.benchmark.sp_svc import SpSVC


def main(svc):
    from revise.tools.metric import compute_metric

    os.makedirs(svc.config.result_dir, exist_ok=True)
    if svc.st_adata.shape[0] == svc.real_st_adata.shape[0]:
        metric_df = compute_metric(
            svc.st_adata, svc.real_st_adata, svc.logger, adata_process=False, gene_list=None, normalize=True)
        metric_df.to_csv(os.path.join(svc.config.result_dir, "metrics_raw.csv"))
    svc.global_anchoring()
    svc.local_refinement()
    for key, svc_adata in svc.svc.items():
        common_index = svc_adata.obs.index.intersection(svc.real_st_adata.obs.index)
        svc_adata = svc_adata[common_index, :].copy()
        real_st_adata = svc.real_st_adata[common_index, :].copy()
        svc.logger.info(f"compute metric for {key}, imputed data shape {svc_adata.shape}, ground truth data shape {real_st_adata.shape}")
        metrics_df = compute_metric(
            svc_adata, real_st_adata, svc.logger, adata_process=False, gene_list=None, normalize=True)
        metrics_df.to_csv(os.path.join(svc.config.result_dir, f"metrics_normalized.csv"))


def __getattr__(name):
    """Lazily load benchmark classes so imports keep working after packaging."""
    mapping = {
        "SpSVC": "revise.benchmark.sp_svc",
        "ScSVCSr": "revise.benchmark.sc_svc_sr",
        "ScSVCImpute": "revise.benchmark.sc_svc_impute",
    }
    if name in mapping:
        module = import_module(mapping[name])
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = ['SpSVC', 'ScSVCSr', 'ScSVCImpute']
