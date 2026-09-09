from revise.backend.kernels.factory import KERNEL_REGISTRY, build_kernel

_KERNEL_EXPORTS = {
    "GeneImputeKernel": ("revise.backend.kernels.gene_impute", "GeneImputeKernel"),
    "GeneUncertaintyKernel": ("revise.backend.kernels.gene_uncertainty", "GeneUncertaintyKernel"),
    "GlobalAnchoringKernel": ("revise.backend.kernels.global_anchoring", "GlobalAnchoringKernel"),
    "LocalAnchoringKernel": ("revise.backend.kernels.local_anchoring", "LocalAnchoringKernel"),
    "OTKernel": ("revise.backend.kernels.ot", "OTKernel"),
    "GraphAggregateKernel": ("revise.backend.kernels.graph_aggregate", "GraphAggregateKernel"),
    "GraphClusterKernel": ("revise.backend.kernels.graph_cluster", "GraphClusterKernel"),
    "SegEvaluateKernel": ("revise.backend.kernels.seg_evaluate", "SegEvaluateKernel"),
    "SpotSrKernel": ("revise.backend.kernels.spot_sr", "SpotSrKernel"),
}


def __getattr__(name):
    if name not in _KERNEL_EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    from importlib import import_module

    module_name, class_name = _KERNEL_EXPORTS[name]
    obj = getattr(import_module(module_name), class_name)
    globals()[name] = obj
    return obj

__all__ = [
    "KERNEL_REGISTRY",
    "build_kernel",
    "GeneImputeKernel",
    "GeneUncertaintyKernel",
    "GlobalAnchoringKernel",
    "LocalAnchoringKernel",
    "OTKernel",
    "GraphAggregateKernel",
    "GraphClusterKernel",
    "SegEvaluateKernel",
    "SpotSrKernel",
]
