_OPS_EXPORTS = {
    "AssignmentState": ("revise.backend.ops.assignment", "AssignmentState"),
    "AssignmentStateError": ("revise.backend.ops.assignment", "AssignmentStateError"),
    "aggregate_assignment": ("revise.backend.ops.assignment", "aggregate_assignment"),
    "align_assignment_categories": ("revise.backend.ops.assignment", "align_assignment_categories"),
    "align_assignment_observations": ("revise.backend.ops.assignment", "align_assignment_observations"),
    "argmax_assignment": ("revise.backend.ops.assignment", "argmax_assignment"),
    "assignment_state_evidence": ("revise.backend.ops.assignment", "assignment_state_evidence"),
    "one_hot_assignment": ("revise.backend.ops.assignment", "one_hot_assignment"),
    "project_assignment": ("revise.backend.ops.assignment", "project_assignment"),
    "validate_assignment": ("revise.backend.ops.assignment", "validate_assignment"),
    "AssignmentGuidanceCollector": ("revise.backend.ops.assignment_guidance", "AssignmentGuidanceCollector"),
    "assignment_compatibility": ("revise.backend.ops.assignment_guidance", "assignment_compatibility"),
    "graph_guidance": ("revise.backend.ops.assignment_guidance", "graph_guidance"),
    "ot_cost_guidance": ("revise.backend.ops.assignment_guidance", "ot_cost_guidance"),
    "resolve_assignment_guidance": ("revise.backend.ops.assignment_guidance", "resolve_assignment_guidance"),
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
    "solve_local_ot": ("revise.backend.ops.local_ot", "solve_local_ot"),
    "condition_cost_matrix": ("revise.backend.ops.posterior_conditioning", "condition_cost_matrix"),
    "condition_sparse_graph": ("revise.backend.ops.posterior_conditioning", "condition_sparse_graph"),
    "get_posterior_matrix": ("revise.backend.ops.posterior_conditioning", "get_posterior_matrix"),
    "neighbor_posterior_affinity": ("revise.backend.ops.posterior_conditioning", "neighbor_posterior_affinity"),
    "posterior_affinity": ("revise.backend.ops.posterior_conditioning", "posterior_affinity"),
    "posterior_conditioning_enabled": ("revise.backend.ops.posterior_conditioning", "posterior_conditioning_enabled"),
    "posterior_conditioning_mode": ("revise.backend.ops.posterior_conditioning", "posterior_conditioning_mode"),
    "posterior_conditioning_strict": ("revise.backend.ops.posterior_conditioning", "posterior_conditioning_strict"),
    "posterior_reference_allocation": (
        "revise.backend.ops.posterior_conditioning",
        "posterior_reference_allocation",
    ),
    "reference_measure_from_marginals": (
        "revise.backend.ops.posterior_conditioning",
        "reference_measure_from_marginals",
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
