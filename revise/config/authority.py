"""Typed, package-owned authority for engine defaults and semantic routes."""

from __future__ import annotations

import copy
from dataclasses import dataclass
import hashlib
import json
from typing import Any


@dataclass(frozen=True)
class RouteSpec:
    profile: str
    task: str
    svc_kind: str
    strategy: str
    overrides: dict[str, Any]


ENGINE_DEFAULTS: dict[str, Any] = {
    "runtime": {"seed": 42, "deterministic": True, "compatibility_mode": False},
    "io": {
        "data_root": "raw_data",
        "output_root": "results_unified",
        "sample_name": "P2CRC",
        "st_path": None,
        "sc_ref_path": None,
        "pm_on_cell_path": None,
        "st_file": "Xenium.h5ad",
        "sc_ref_file": "adata_sc_all_reanno.h5ad",
        "gt_svc_file": "selected_xenium.h5ad",
        "seg_method": "seg_1",
        "spot_size": 50,
        "patient_key": "Patient",
        "sample_size": None,
        "save_outputs": True,
        "input_format": "h5ad",
        "spatialdata_path": None,
        "spatialdata_reader": "zarr",
        "spatialdata_table": None,
        "spatialdata_spatial_element": None,
        "spatialdata_coordinate_system": "global",
    },
    "columns": {
        "cell_type_col": "Level1",
        "sub_cell_type_col": "Level2",
        "confidence_col": "Confidence",
        "unknown_key": "Unknown",
    },
    "preprocess": {
        "st_min_counts": 20,
        "st_min_cells": 30,
        "sc_min_counts": 20,
        "sc_min_cells": 50,
        "st_min_transcripts": 60,
    },
    "graph": {
        "method": "joint",
        "alpha": 0.5,
        "n_neighbors": 10,
        "exp_neighbors": 10,
        "spatial_neighbors": 10,
        "random_state": 0,
    },
    "ot": {
        "ga": {"solver": "pot", "pot": {"reg": 0.1, "reg_m": 0.0, "reg_type": "entropy"}},
        "lr": {"solver": "pot", "pot": {"reg": 0.05, "reg_m": 1.0, "reg_type": "kl"}},
        "impute": {"reg": 5.0, "reg_m": 0.0, "reg_type": "kl"},
    },
    "plot": {
        "enabled": False,
        "cluster_resolutions": [0.1, 0.3, 0.5],
        "min_genes": 20,
        "min_cells": 3,
        "sample_size": 10000,
    },
    "reconstruct": {"alpha": 1.0},
    "benchmark": {
        "evaluate": False,
        "dropout_total_counts": 60,
        "swapping_total_counts": 300,
        "lower_ts": 0.2,
        "upper_ts": 0.8,
    },
    "sc": {
        "select_ct": None,
        "resolutions": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        "select_resolution": None,
        "match_spot_sum": False,
        "svc_completeness": True,
        "sr_graph_agg_enabled": False,
        "sr_graph_agg_low_conf_only": False,
        "sr_graph_agg_low_conf_quantile": 0.2,
        "sr_graph_agg_anchor_only": False,
        "sr_graph_agg_anchor_high_conf_quantile": 0.8,
        "sr_graph_agg_confidence_mode": "auto",
        "sr_graph_agg_conf_weighted_alpha": False,
        "sr_graph_agg_conf_alpha_min": 0.0,
        "sr_graph_agg_conf_alpha_max": -1.0,
        "sr_graph_agg_conf_alpha_power": 1.0,
    },
    "impute": {
        "merge_subcluster_method": "mean",
        "subcluster_resolution": 3,
        "in_panel_subcluster_resolution": None,
        "prune": True,
        "n_neighbors": 1,
        "method": "mean",
        "graph_preprocess": True,
        "graph_n_pcs": 50,
    },
}


_SR_BENCHMARK_OVERRIDES = {
    "runtime": {"compatibility_mode": True},
    "benchmark": {"evaluate": True},
    "graph": {"alpha": 0.2},
    "local_refinement": {"strength": 0.0},
    "sc": {
        "sr_graph_agg_enabled": True,
        "sr_graph_agg_low_conf_only": True,
        "sr_graph_agg_low_conf_quantile": 0.1,
        "sr_graph_agg_anchor_only": True,
        "sr_graph_agg_anchor_high_conf_quantile": 0.9,
        "sr_graph_agg_confidence_mode": "auto",
        "sr_graph_agg_conf_weighted_alpha": True,
        "sr_graph_agg_conf_alpha_min": 0.0,
        "sr_graph_agg_conf_alpha_max": 0.25,
        "sr_graph_agg_conf_alpha_power": 1.0,
    },
    "ot": {"ga": {"pot": {"reg": 0.01, "reg_m": 0.0001, "reg_type": "kl"}}},
    "io": {
        "st_file": "xenium_spot.h5ad",
        "gt_svc_file": "selected_xenium.h5ad",
        "sc_ref_file": "real_sc_ref_all.h5ad",
        "spot_size": 50,
    },
}


_IMPUTE_BENCHMARK_OVERRIDES = {
    "runtime": {"compatibility_mode": True},
    "benchmark": {"evaluate": True},
    "graph": {"alpha": 0.2, "n_neighbors": 15},
    "ot": {"ga": {"pot": {"reg": 0.01, "reg_m": 0.0001, "reg_type": "kl"}}},
    "io": {
        "st_file": "selected_xenium.h5ad",
        "gt_svc_file": "selected_xenium.h5ad",
        "sc_ref_file": "real_sc_ref.h5ad",
    },
}


ROUTES: dict[str, dict[str, RouteSpec]] = {
    "application": {
        "sp-SVC": RouteSpec(
            "application_sp", "sp_svc", "sp", "SpSvcApplicationStrategy",
            {
                "reconstruct": {"alpha": 0.5},
                "local_refinement": {"strength": 0.2},
            },
        ),
        "sc-SVC:cluster": RouteSpec(
            "application_sc", "sc_svc", "sc", "ScSvcApplicationStrategy",
            {
                "preprocess": {"st_min_transcripts": 60, "st_min_cells": 100, "sc_min_cells": 100},
                "graph": {"alpha": 0.2, "n_neighbors": 10, "exp_neighbors": 15, "spatial_neighbors": 6},
                "reconstruct": {"alpha": 0.5},
                "ot": {
                    "ga": {"solver": "tacco", "pot": {"reg": 0.1, "reg_m": 0.0, "reg_type": "entropy"}},
                    "lr": {"solver": "tacco", "pot": {"reg": 0.1, "reg_m": 0.0, "reg_type": "entropy"}},
                },
                "sc": {"resolutions": [0.6, 0.7, 0.8], "tacco_annotate": {"multi_center": 1, "lamb": 0.001}},
            },
        ),
        "sc-SVC:sr": RouteSpec(
            "application_sc_super_resolution", "sc_svc_super_resolution", "sc", "ScSvcSuperResolutionApplicationStrategy",
            {
                "graph": {"alpha": 0.2},
                "local_refinement": {"strength": 0.0},
            },
        ),
    },
    "benchmark": {
        "segmentation": RouteSpec(
            "benchmark_seg", "sp_svc", "sp", "SpSvcBenchmarkSegStrategy",
            {
                "runtime": {"compatibility_mode": True},
                "benchmark": {"evaluate": True},
                "graph": {"n_neighbors": 50, "exp_neighbors": 30, "spatial_neighbors": 30, "alpha": 0.8},
                "ot": {"lr": {"pot": {"reg": 1.0, "reg_m": 0.0, "reg_type": "kl"}}},
                "io": {"seg_method": "seg_1"},
                "local_refinement": {"strength": 0.2},
            },
        ),
        "bin2cell": RouteSpec(
            "benchmark_bin2cell", "sp_svc", "sp", "SpSvcBenchmarkSegStrategy",
            {
                "runtime": {"compatibility_mode": True},
                "benchmark": {"evaluate": True},
                "graph": {"n_neighbors": 50, "exp_neighbors": 30, "spatial_neighbors": 30, "alpha": 0.8},
                "ot": {"lr": {"pot": {"reg": 1.0, "reg_m": 0.0, "reg_type": "kl"}}},
                "io": {"seg_method": "bin2cell"},
                "local_refinement": {"strength": 0.2},
            },
        ),
        "batch_effect": RouteSpec(
            "benchmark_sr_batch", "sc_svc_sr", "sc", "ScSvcSrBenchmarkStrategy",
            copy.deepcopy(_SR_BENCHMARK_OVERRIDES),
        ),
        "spot_size": RouteSpec(
            "benchmark_sr_spot_size", "sc_svc_sr", "sc", "ScSvcSrBenchmarkStrategy",
            copy.deepcopy(_SR_BENCHMARK_OVERRIDES),
        ),
        "gene_panel": RouteSpec(
            "benchmark_impute_panel", "sc_svc_impute", "sc", "ScSvcImputeBenchmarkStrategy",
            copy.deepcopy(_IMPUTE_BENCHMARK_OVERRIDES),
        ),
        "gene_dropout": RouteSpec(
            "benchmark_impute_dropout", "sc_svc_impute", "sc", "ScSvcImputeBenchmarkStrategy",
            copy.deepcopy(_IMPUTE_BENCHMARK_OVERRIDES),
        ),
    },
}


LOCKED_KEYS = frozenset(
    {
        "ot.ga.pot.reg",
        "ot.ga.pot.reg_m",
        "ot.ga.pot.reg_type",
        "ot.lr.pot.reg",
        "ot.lr.pot.reg_m",
        "ot.lr.pot.reg_type",
        "sc.tacco_annotate.multi_center",
        "sc.tacco_annotate.lamb",
    }
)


def _hash(payload: Any) -> str:
    serialized = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()


def _authority_document() -> dict[str, Any]:
    profiles = {}
    router = {}
    for namespace, routes in ROUTES.items():
        router[namespace] = {}
        for selector, spec in routes.items():
            profiles[spec.profile] = copy.deepcopy(spec.overrides)
            router[namespace][selector] = {
                "profile": spec.profile,
                "task": spec.task,
                "svc_kind": spec.svc_kind,
                "strategy": spec.strategy,
            }
    return {
        "version": 2,
        "defaults": copy.deepcopy(ENGINE_DEFAULTS),
        "router": router,
        "profiles": profiles,
        "locked_params": {"keys": sorted(LOCKED_KEYS)},
    }


ENGINE_DEFAULTS_HASH = _hash(ENGINE_DEFAULTS)
AUTHORITY_HASH = _hash(_authority_document())


__all__ = [
    "AUTHORITY_HASH",
    "ENGINE_DEFAULTS",
    "ENGINE_DEFAULTS_HASH",
    "LOCKED_KEYS",
    "ROUTES",
    "RouteSpec",
]
