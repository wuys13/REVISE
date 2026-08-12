import os
from dataclasses import dataclass
from typing import Optional

from revise.conf.base_conf import BaseConf


@dataclass
class ApplicationScSrConf(BaseConf):
    """
    sc-SVC super-resolution configuration class for application usage.

    """
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
    rec_graph_n_neighbors = 20
    rec_graph_method = "joint"
    rec_graph_alpha = 0.2
    rec_graph_exp_neighbor_num = 10
    rec_graph_spatial_neighbor_num = 20

    # pot parameters
    rec_pot_reg = 0.05
    rec_pot_reg_m = 1.0
    rec_pot_reg_type = "kl"
    rec_alpha = 1.0
    rec_match_spot_sum = False

    # svc parameters
    svc_completeness: bool = True

    @property
    def result_dir(self):
        return os.path.join(self.result_root_path, self.sample_name)

    @property
    def st_file_path(self):
        return os.path.join(self.raw_data_path, f"{self.sample_name}_{self.st_file}")

    @property
    def sc_ref_file_path(self):
        return os.path.join(self.raw_data_path, self.sc_ref_file)

    @property
    def pm_on_cell_file(self):
        return os.path.join(self.raw_data_path, "PM_on_cell.csv")
