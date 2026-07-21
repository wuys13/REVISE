from __future__ import annotations

from importlib import import_module
from typing import Dict, Tuple, Type

from revise.backend.kernels.base import BaseKernel


KERNEL_REGISTRY: Dict[str, Tuple[str, str]] = {
    "global_anchoring": ("revise.backend.kernels.global_anchoring", "GlobalAnchoringKernel"),
    "graph_aggregate": ("revise.backend.kernels.graph_aggregate", "GraphAggregateKernel"),
    "graph_cluster": ("revise.backend.kernels.graph_cluster", "GraphClusterKernel"),
    "spot_sr": ("revise.backend.kernels.spot_sr", "SpotSrKernel"),
    "gene_impute": ("revise.backend.kernels.gene_impute", "GeneImputeKernel"),
    "gene_uncertainty": ("revise.backend.kernels.gene_uncertainty", "GeneUncertaintyKernel"),
    "seg_evaluate": ("revise.backend.kernels.seg_evaluate", "SegEvaluateKernel"),
}


def _load_kernel_class(kernel_name: str) -> Type[BaseKernel]:
    module_name, class_name = KERNEL_REGISTRY[kernel_name]
    module = import_module(module_name)
    return getattr(module, class_name)


def build_kernel(kernel_name: str, config, logger=None) -> BaseKernel:
    if kernel_name not in KERNEL_REGISTRY:
        raise KeyError(f"Unknown backend kernel: {kernel_name}")
    return _load_kernel_class(kernel_name)(config=config, logger=logger)
