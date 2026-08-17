from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Optional

"""Runner-side configuration contracts.

This module is intentionally separate from `revise.config.authority`:
- `revise.config.authority` owns engine defaults and semantic routes.
- Dataclasses in this module are internal runner contracts consumed by compatibility
  kernels/runners, created in `revise.backend.adapters` from validated merged config.

Keeping this layer explicit makes route resolution and runner compatibility concerns
independent and easier to evolve.
"""


REQUIRED_IO_BY_MODE_TASK = {
    ("application", "sp_svc"): {"st_file", "sc_ref_file"},
    ("application", "sc_svc"): {"st_file", "sc_ref_file"},
    ("application", "sc_svc_super_resolution"): {"st_file", "sc_ref_file"},
    ("benchmark", "sp_svc"): {"st_file", "sc_ref_file", "gt_svc_file"},
    ("benchmark", "sc_svc_sr"): {
        "st_file",
        "sc_ref_file",
        "gt_svc_file",
        "spot_size",
    },
    ("benchmark", "sc_svc_impute"): {
        "st_file",
        "sc_ref_file",
        "gt_svc_file",
    },
}


@dataclass(frozen=True)
class InputSpec:
    role: str
    path: str


def application_st_path(
    data_root: str,
    sample_name: str,
    filename: str,
) -> str:
    return os.path.join(data_root, f"{sample_name}_{filename}")


def application_sc_ref_path(data_root: str, filename: str) -> str:
    return os.path.join(data_root, filename)


def benchmark_role_path(
    data_root: str,
    sample_name: str,
    filename: str,
) -> str:
    return os.path.join(data_root, sample_name, filename)


def benchmark_st_path(
    data_root: str,
    sample_name: str,
    filename: str,
    *,
    task: str,
    seg_method: Optional[str] = None,
    spot_size: Optional[int] = None,
) -> str:
    if task == "sp_svc":
        if seg_method is None:
            raise ValueError("seg_method is required for benchmark sp_svc")
        return os.path.join(
            data_root,
            sample_name,
            str(seg_method),
            filename,
        )
    if task == "sc_svc_sr":
        if spot_size is None:
            raise ValueError("spot_size is required for benchmark sc_svc_sr")
        return os.path.join(
            data_root,
            sample_name,
            f"spot_{int(spot_size)}",
            filename,
        )
    if task == "sc_svc_impute":
        return benchmark_role_path(data_root, sample_name, filename)
    raise ValueError(f"Unsupported benchmark task for ST path: {task}")


def pm_on_cell_path_from_data_root(data_root: str) -> str:
    return os.path.join(data_root, "PM_on_cell.csv")


def configured_st_source_path(io_config, h5ad_path: str) -> str:
    input_format = str(io_config["input_format"]).lower()
    spatialdata_path = io_config.get("spatialdata_path")
    if input_format in {"spatialdata", "auto"} and spatialdata_path:
        return str(spatialdata_path)
    return h5ad_path


def resolve_input_specs(runtime, io_config) -> tuple[InputSpec, ...]:
    """Resolve the one role/path contract shared by preflight and loading."""
    mode = str(runtime.get("mode"))
    task = str(runtime.get("task"))
    key = (mode, task)
    if key not in REQUIRED_IO_BY_MODE_TASK:
        raise ValueError(f"Unsupported mode/task input layout: {mode}:{task}")

    if mode == "application" and io_config.get("st_path") and io_config.get("sc_ref_path"):
        return (
            InputSpec(
                "st",
                configured_st_source_path(io_config, str(io_config["st_path"])),
            ),
            InputSpec("sc_ref", str(io_config["sc_ref_path"])),
        )

    required = REQUIRED_IO_BY_MODE_TASK[key]
    missing = sorted(
        name for name in required if io_config.get(name) in (None, "")
    )
    if missing:
        raise ValueError(f"Missing input path keys for {mode}:{task}: {missing}")

    data_root = str(io_config["data_root"])
    sample_name = str(io_config["sample_name"])
    st_file = str(io_config["st_file"])
    sc_ref_file = str(io_config["sc_ref_file"])
    if mode == "application":
        return (
            InputSpec(
                "st",
                configured_st_source_path(
                    io_config,
                    application_st_path(data_root, sample_name, st_file),
                ),
            ),
            InputSpec("sc_ref", application_sc_ref_path(data_root, sc_ref_file)),
        )

    gt_file = str(io_config["gt_svc_file"])
    return (
        InputSpec(
            "st",
            configured_st_source_path(
                io_config,
                benchmark_st_path(
                    data_root,
                    sample_name,
                    st_file,
                    task=task,
                    seg_method=io_config["seg_method"],
                    spot_size=io_config["spot_size"],
                ),
            ),
        ),
        InputSpec(
            "sc_ref",
            benchmark_role_path(data_root, sample_name, sc_ref_file),
        ),
        InputSpec(
            "gt",
            benchmark_role_path(data_root, sample_name, gt_file),
        ),
    )


def resolved_input_path(specs, role: str) -> str:
    """Return one preflight-resolved role path."""
    for spec in specs:
        if spec.role == role:
            return spec.path
    raise KeyError(f"Missing resolved input role: {role}")


@dataclass
class BaseConf:
    # runtime parameters
    sample_name: str
    raw_data_path: str
    result_root_path: str

    # annotate column keys
    cell_type_col: str
    confidence_col: str
    unknown_key: str


@dataclass(kw_only=True)
class ApplicationSpConf(BaseConf):
    st_file: str
    sc_ref_file: str
    annotate_mode: str

    # annotate parameters
    annotate_pot_reg: float
    annotate_pot_reg_m: float
    annotate_pot_reg_type: str

    # plot parameters
    plot_flag: bool
    plot_cluster_resolution: list
    plot_min_genes: int
    plot_min_cells: int
    plot_sample_size: int

    # reconstruct parameters
    rec_graph_n_neighbors: int
    rec_graph_exp_neighbor_num: int
    rec_graph_spatial_neighbor_num: int
    rec_graph_method: str
    rec_graph_alpha: float

    # reconstruct ot
    rec_pot_reg: float
    rec_pot_reg_m: float
    rec_pot_reg_type: str
    rec_ot_method: str
    rec_alpha: float
    local_refinement_strength: float

    @property
    def result_dir(self):
        return os.path.join(self.result_root_path, self.sample_name)

    @property
    def st_file_path(self):
        return application_st_path(
            self.raw_data_path,
            self.sample_name,
            self.st_file,
        )

    @property
    def sc_ref_file_path(self):
        return application_sc_ref_path(self.raw_data_path, self.sc_ref_file)


@dataclass(kw_only=True)
class ApplicationScConf(BaseConf):
    st_file: str
    sc_ref_file: str
    annotate_mode: Optional[str]

    # annotate parameters
    annotate_pot_reg: float
    annotate_pot_reg_m: float
    annotate_pot_reg_type: str
    tacco_annotate_multi_center: Optional[int]
    tacco_annotate_lamb: Optional[float]

    # reconstruct parameters
    rec_graph_n_neighbors: int
    rec_graph_exp_neighbor_num: int
    rec_graph_spatial_neighbor_num: int
    rec_graph_method: str
    rec_graph_alpha: float
    rec_random_state: int

    # reconstruct ot
    rec_pot_reg: float
    rec_pot_reg_m: float
    rec_pot_reg_type: str
    rec_ot_method: str
    rec_alpha: float
    rec_match_spot_sum: bool

    @property
    def result_dir(self):
        return os.path.join(self.result_root_path, self.sample_name)

    @property
    def st_file_path(self):
        return application_st_path(
            self.raw_data_path,
            self.sample_name,
            self.st_file,
        )

    @property
    def sc_ref_file_path(self):
        return application_sc_ref_path(self.raw_data_path, self.sc_ref_file)


@dataclass(kw_only=True)
class ApplicationScSuperResolutionConf(BaseConf):
    st_file: str
    sc_ref_file: str
    annotate_mode: Optional[str]

    # annotate parameters
    annotate_pot_reg: float
    annotate_pot_reg_m: float
    annotate_pot_reg_type: str

    # graph parameters
    rec_graph_n_neighbors: int
    rec_graph_method: str
    rec_graph_alpha: float
    rec_graph_exp_neighbor_num: int
    rec_graph_spatial_neighbor_num: int

    # pot parameters
    rec_pot_reg: float
    rec_pot_reg_m: float
    rec_pot_reg_type: str
    rec_ot_method: str
    rec_alpha: float
    rec_match_spot_sum: bool
    rec_graph_agg_enabled: bool

    # svc parameters
    svc_completeness: bool
    sr_assignment_seed: int
    local_refinement_strength: float

    @property
    def result_dir(self):
        return os.path.join(self.result_root_path, self.sample_name)

    @property
    def st_file_path(self):
        return application_st_path(
            self.raw_data_path,
            self.sample_name,
            self.st_file,
        )

    @property
    def sc_ref_file_path(self):
        return application_sc_ref_path(self.raw_data_path, self.sc_ref_file)

@dataclass(kw_only=True)
class BenchmarkSegConf(BaseConf):
    st_file: str
    gt_svc_file: str
    sc_ref_file: str
    seg_method: str
    annotate_mode: str
    case_subdir: Optional[str] = None

    # annotate parameters
    annotate_pot_reg: float
    annotate_pot_reg_m: float
    annotate_pot_reg_type: str

    # segmentation evaluation policy
    dropout_total_counts: int
    swapping_total_counts: int
    lower_ts: float
    upper_ts: float

    # reconstruct graph
    rec_graph_n_neighbors: int
    rec_graph_exp_neighbor_num: int
    rec_graph_spatial_neighbor_num: int
    rec_graph_method: str
    rec_graph_alpha: float

    # reconstruct ot
    rec_pot_reg: float
    rec_pot_reg_m: float
    rec_pot_reg_type: str
    rec_ot_method: str
    rec_alpha: float
    local_refinement_strength: float

    @property
    def result_dir(self):
        leaf = self.case_subdir or self.seg_method
        return os.path.join(self.result_root_path, self.sample_name, leaf)

    @property
    def st_file_path(self):
        return benchmark_st_path(
            self.raw_data_path,
            self.sample_name,
            self.st_file,
            task="sp_svc",
            seg_method=self.seg_method,
        )

    @property
    def gt_svc_file_path(self):
        return benchmark_role_path(
            self.raw_data_path,
            self.sample_name,
            self.gt_svc_file,
        )

    @property
    def sc_ref_file_path(self):
        return benchmark_role_path(
            self.raw_data_path,
            self.sample_name,
            self.sc_ref_file,
        )


@dataclass(kw_only=True)
class BenchmarkSrConf(BaseConf):
    st_file: str
    gt_svc_file: str
    sc_ref_file: str
    spot_size: int
    annotate_mode: str
    case_subdir: Optional[str] = None

    # annotate parameters
    annotate_pot_reg: float
    annotate_pot_reg_m: float
    annotate_pot_reg_type: str

    # svc parameters
    svc_completeness: bool
    sr_assignment_seed: int

    # optional graph aggregation (SR robustness benchmark)
    rec_graph_n_neighbors: int
    rec_graph_exp_neighbor_num: int
    rec_graph_spatial_neighbor_num: int
    rec_graph_method: str
    rec_graph_alpha: float
    rec_pot_reg: float
    rec_pot_reg_m: float
    rec_pot_reg_type: str
    rec_ot_method: str
    rec_alpha: float
    rec_graph_agg_enabled: bool
    rec_graph_agg_low_conf_only: bool
    rec_graph_agg_low_conf_quantile: float
    rec_graph_agg_anchor_only: bool
    rec_graph_agg_anchor_high_conf_quantile: float
    rec_graph_agg_confidence_mode: str
    rec_graph_agg_conf_weighted_alpha: bool
    rec_graph_agg_conf_alpha_min: float
    rec_graph_agg_conf_alpha_max: float
    rec_graph_agg_conf_alpha_power: float

    local_refinement_strength: float

    @property
    def result_dir(self):
        leaf = self.case_subdir or f"spot_{self.spot_size}"
        return os.path.join(self.result_root_path, self.sample_name, leaf)

    @property
    def st_file_path(self):
        return benchmark_st_path(
            self.raw_data_path,
            self.sample_name,
            self.st_file,
            task="sc_svc_sr",
            spot_size=self.spot_size,
        )

    @property
    def gt_svc_file_path(self):
        return benchmark_role_path(
            self.raw_data_path,
            self.sample_name,
            self.gt_svc_file,
        )

    @property
    def sc_ref_file_path(self):
        return benchmark_role_path(
            self.raw_data_path,
            self.sample_name,
            self.sc_ref_file,
        )

@dataclass(kw_only=True)
class BenchmarkImputeConf(BaseConf):
    st_file: str
    gt_svc_file: str
    sc_ref_file: str
    annotate_mode: str
    case_subdir: Optional[str] = None

    # preprocess parameters
    prep_min_cells: int
    prep_min_counts: int

    # annotate parameters
    annotate_pot_reg: float
    annotate_pot_reg_m: float
    annotate_pot_reg_type: str

    # reconstruct graph
    rec_graph_preprocess: bool
    rec_graph_n_pcs: int
    rec_graph_n_neighbors: int

    # reconstruct ot
    rec_impute_pot_reg: float
    rec_impute_pot_reg_m: float
    rec_impute_pot_reg_type: str
    rec_ot_method: str

    # reconstruct impute
    rec_merge_subcluster_method: str
    rec_subcluster_resolution: int
    rec_in_panel_subcluster_resolution: Optional[int]
    rec_impute_prune_flag: bool
    rec_impute_n_neighbors: int
    rec_impute_method: str

    @property
    def result_dir(self):
        if self.case_subdir:
            return os.path.join(self.result_root_path, self.sample_name, self.case_subdir)
        return os.path.join(self.result_root_path, self.sample_name)

    @property
    def st_file_path(self):
        return benchmark_st_path(
            self.raw_data_path,
            self.sample_name,
            self.st_file,
            task="sc_svc_impute",
        )

    @property
    def sc_ref_file_path(self):
        return benchmark_role_path(
            self.raw_data_path,
            self.sample_name,
            self.sc_ref_file,
        )

    @property
    def gt_svc_file_path(self):
        return benchmark_role_path(
            self.raw_data_path,
            self.sample_name,
            self.gt_svc_file,
        )
