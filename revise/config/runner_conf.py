from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Optional

"""Runner-side configuration contracts.

This module is intentionally separate from `revise/revise.yaml`:
- `revise/revise.yaml` is the external source-of-truth config and routing surface.
- Dataclasses in this module are internal runner contracts consumed by compatibility
  kernels/runners, created in `revise.backend.adapters` from merged YAML config.

Keeping this layer explicit makes route resolution and runner compatibility concerns
independent and easier to evolve.
"""


REQUIRED_IO_BY_MODE_TASK = {
    ("application", "sp_svc"): {"st_file", "sc_ref_file"},
    ("application", "sc_svc"): {"st_file", "sc_ref_file"},
    ("application", "sc_svc_sr"): {"st_file", "sc_ref_file"},
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
        return os.path.join(
            data_root,
            sample_name,
            str(seg_method or "seg_1"),
            filename,
        )
    if task == "sc_svc_sr":
        return os.path.join(
            data_root,
            sample_name,
            f"spot_{int(spot_size if spot_size is not None else 50)}",
            filename,
        )
    if task == "sc_svc_impute":
        return benchmark_role_path(data_root, sample_name, filename)
    raise ValueError(f"Unsupported benchmark task for ST path: {task}")


def pm_on_cell_path_from_data_root(data_root: str) -> str:
    return os.path.join(data_root, "PM_on_cell.csv")


def configured_st_source_path(io_config, h5ad_path: str) -> str:
    input_format = str(io_config.get("input_format", "h5ad")).lower()
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
                    seg_method=io_config.get("seg_method"),
                    spot_size=io_config.get("spot_size"),
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


def resolved_input_path(specs, role: str, fallback: str) -> str:
    """Return one preflight-resolved role path, with a test-helper fallback."""
    for spec in specs or ():
        if spec.role == role:
            return spec.path
    return fallback


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


@dataclass
class ApplicationSpConf(BaseConf):
    st_file: str
    sc_ref_file: str
    annotate_mode: str

    # annotate parameters
    annotate_pot_reg: float = 0.1
    annotate_pot_reg_m: float = 0.0
    annotate_pot_reg_type: str = "entropy"

    # preprocess parameters
    prep_st_min_counts: int = 20
    prep_st_min_cells: int = 30
    prep_sc_min_counts: int = 20
    prep_sc_min_cells: int = 50

    # plot parameters
    plot_flag: bool = True
    plot_cluster_resolution: list = field(default_factory=lambda: [0.1, 0.3, 0.5])
    plot_min_genes: int = 20
    plot_min_cells: int = 3
    plot_sample_size: int = 10000

    # reconstruct parameters
    rec_graph_n_neighbors: int = 10
    rec_graph_exp_neighbor_num: int = 10
    rec_graph_spatial_neighbor_num: int = 10
    rec_graph_method: str = "joint"
    rec_graph_alpha: float = 0.4

    # reconstruct ot
    rec_pot_reg: float = 0.05
    rec_pot_reg_m: float = 1.0
    rec_pot_reg_type: str = "kl"
    rec_ot_method: str = "pot"
    rec_alpha = 0.5

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


@dataclass
class ApplicationScConf(BaseConf):
    st_file: str
    sc_ref_file: str
    annotate_mode: Optional[str] = None

    # annotate parameters
    annotate_pot_reg: float = 0.06
    annotate_pot_reg_m: float = 0.015
    annotate_pot_reg_type: str = "entropy"
    tacco_annotate_multi_center: Optional[int] = None
    tacco_annotate_lamb: Optional[float] = None

    # preprocess parameters
    prep_st_min_counts: int = 60
    prep_st_min_cells: int = 100
    prep_sc_min_counts: int = 0
    prep_sc_min_cells: int = 100

    # reconstruct parameters
    rec_graph_n_neighbors: int = 10
    rec_graph_exp_neighbor_num: int = 15
    rec_graph_spatial_neighbor_num: int = 6
    rec_graph_method: str = "joint"
    rec_graph_alpha: float = 0.2

    # reconstruct ot
    rec_pot_reg: float = 0.06
    rec_pot_reg_m: float = 0.015
    rec_pot_reg_type: str = "entropy"
    rec_ot_method: str = "pot"
    rec_alpha = 0.5
    rec_match_spot_sum: bool = False

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


@dataclass
class ApplicationScSrConf(BaseConf):
    st_file: str
    sc_ref_file: str
    annotate_mode: Optional[str] = None

    # annotate parameters
    annotate_pot_reg: float = 0.01
    annotate_pot_reg_m: float = 0.0001
    annotate_pot_reg_type: str = "kl"

    # preprocess parameters
    prep_st_min_counts: int = 60
    prep_st_min_cells: int = 100
    prep_sc_min_counts: int = 0
    prep_sc_min_cells: int = 100

    # graph parameters
    rec_graph_n_neighbors: int = 20
    rec_graph_method: str = "joint"
    rec_graph_alpha: float = 0.2
    rec_graph_exp_neighbor_num: int = 10
    rec_graph_spatial_neighbor_num: int = 20

    # pot parameters
    rec_pot_reg: float = 0.05
    rec_pot_reg_m: float = 1.0
    rec_pot_reg_type: str = "kl"
    rec_ot_method: str = "pot"
    rec_alpha: float = 1.0
    rec_match_spot_sum: bool = False

    # svc parameters
    svc_completeness: bool = True
    sr_assignment_seed: int = 42

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

@dataclass
class BenchmarkSegConf(BaseConf):
    st_file: str
    gt_svc_file: str
    sc_ref_file: str
    seg_method: str
    annotate_mode: str
    case_subdir: Optional[str] = None

    # annotate parameters
    annotate_pot_reg: float = 0.1
    annotate_pot_reg_m: float = 0.0
    annotate_pot_reg_type: str = "entropy"

    # segmentation effect parameters
    dropout_total_counts: int = 60
    swapping_total_counts: int = 300
    lower_ts: float = 0.2
    upper_ts: float = 0.8

    # reconstruct graph
    rec_graph_n_neighbors: int = 50
    rec_graph_exp_neighbor_num: int = 30
    rec_graph_spatial_neighbor_num: int = 30
    rec_graph_method: str = "joint"
    rec_graph_alpha: float = 0.8

    # reconstruct ot
    rec_pot_reg: float = 1.0
    rec_pot_reg_m: float = 0.0
    rec_pot_reg_type: str = "kl"
    rec_ot_method: str = "pot"
    rec_alpha: float = 1.0

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


@dataclass
class BenchmarkSrConf(BaseConf):
    st_file: str
    gt_svc_file: str
    sc_ref_file: str
    spot_size: int
    annotate_mode: str
    case_subdir: Optional[str] = None

    # annotate parameters
    annotate_pot_reg: float = 0.01
    annotate_pot_reg_m: float = 0.0001
    annotate_pot_reg_type: str = "kl"

    # svc parameters
    svc_completeness: bool = True
    sr_assignment_seed: int = 42

    # optional graph aggregation (SR robustness benchmark)
    rec_graph_n_neighbors: int = 20
    rec_graph_exp_neighbor_num: int = 10
    rec_graph_spatial_neighbor_num: int = 20
    rec_graph_method: str = "joint"
    rec_graph_alpha: float = 0.2
    rec_pot_reg: float = 0.05
    rec_pot_reg_m: float = 1.0
    rec_pot_reg_type: str = "kl"
    rec_ot_method: str = "pot"
    rec_alpha: float = 1.0
    rec_graph_agg_enabled: bool = False
    rec_graph_agg_low_conf_only: bool = False
    rec_graph_agg_low_conf_quantile: float = 0.2
    rec_graph_agg_anchor_only: bool = False
    rec_graph_agg_anchor_high_conf_quantile: float = 0.8
    rec_graph_agg_confidence_mode: str = "auto"
    rec_graph_agg_conf_weighted_alpha: bool = False
    rec_graph_agg_conf_alpha_min: float = 0.0
    rec_graph_agg_conf_alpha_max: float = -1.0
    rec_graph_agg_conf_alpha_power: float = 1.0

    # spot-level spatial leakage noise (benchmark stress test)
    sr_noise_enabled: bool = False
    sr_noise_lambda: float = 0.0
    sr_noise_k: int = 4
    sr_noise_weight: str = "distance"
    sr_noise_preserve_total_counts: bool = True
    sr_noise_seed: int = 42

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

@dataclass
class BenchmarkImputeConf(BaseConf):
    st_file: str
    gt_svc_file: str
    sc_ref_file: str
    annotate_mode: str
    case_subdir: Optional[str] = None

    # preprocess parameters
    prep_min_cells: int = 30
    prep_min_counts: int = 60

    # annotate parameters
    annotate_pot_reg: float = 0.01
    annotate_pot_reg_m: float = 0.0001
    annotate_pot_reg_type: str = "kl"

    # reconstruct graph
    rec_graph_preprocess: bool = True
    rec_graph_n_pcs: int = 50
    rec_graph_n_neighbors: int = 15

    # reconstruct ot
    rec_impute_pot_reg: float = 5.0
    rec_impute_pot_reg_m: float = 0.0
    rec_impute_pot_reg_type: str = "kl"
    rec_ot_method: str = "pot"

    # reconstruct impute
    rec_merge_subcluster_method: str = "mean"
    rec_subcluster_resolution: int = 3
    rec_in_panel_subcluster_resolution: Optional[int] = None
    rec_impute_prune_flag: bool = True
    rec_impute_n_neighbors: int = 1
    rec_impute_method: str = "mean"

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
