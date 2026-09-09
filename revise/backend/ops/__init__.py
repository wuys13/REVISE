_OPS_EXPORTS = {
    "GlobalAssignment": ("revise.backend.ops.assignment", "GlobalAssignment"),
    "GlobalAssignmentContractError": (
        "revise.backend.ops.assignment",
        "GlobalAssignmentContractError",
    ),
    "validate_global_assignment": (
        "revise.backend.ops.assignment",
        "validate_global_assignment",
    ),
    "bhattacharyya_distance": ("revise.backend.ops.distance", "bhattacharyya_distance"),
    "compute_autocorr_metrics": ("revise.backend.ops.coefficients", "compute_autocorr_metrics"),
    "construct_sc_ref": ("revise.backend.ops.meta", "construct_sc_ref"),
    "get_adjacency_graph": ("revise.backend.ops.topology", "get_adjacency_graph"),
    "get_sc_obs": ("revise.backend.ops.meta", "get_sc_obs"),
    "get_spatial_score": ("revise.backend.ops.coefficients", "get_spatial_score"),
    "get_subcluster": ("revise.backend.ops.meta", "get_subcluster"),
    "get_true_cell_type": ("revise.backend.ops.meta", "get_true_cell_type"),
    "get_weighted_align_score": ("revise.backend.ops.coefficients", "get_weighted_align_score"),
    "merge_subcluster": ("revise.backend.ops.meta", "merge_subcluster"),
    "align_posterior_categories": (
        "revise.backend.ops.posterior_conditioning",
        "align_posterior_categories",
    ),
    "condition_local_ot_cost": (
        "revise.backend.ops.posterior_conditioning",
        "condition_local_ot_cost",
    ),
    "posterior_reference_allocation": (
        "revise.backend.ops.posterior_conditioning",
        "posterior_reference_allocation",
    ),
    "meta_get_prune_adata": ("revise.backend.ops.meta", "get_prune_adata"),
    "shaver_get_prune_adata": ("revise.backend.ops.shaver", "get_prune_adata"),
    "trim_sp_adata": ("revise.backend.ops.shaver", "trim_sp_adata"),
    "trim_spatial_genes": ("revise.backend.ops.shaver", "trim_spatial_genes"),
}

__all__ = sorted(_OPS_EXPORTS)


def __getattr__(name):
    if name not in _OPS_EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    from importlib import import_module

    module_name, attr_name = _OPS_EXPORTS[name]
    obj = getattr(import_module(module_name), attr_name)
    globals()[name] = obj
    return obj
