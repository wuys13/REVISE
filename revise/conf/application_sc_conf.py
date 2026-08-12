import os
from dataclasses import dataclass
from typing import Optional

from revise.conf.base_conf import BaseConf


@dataclass
class ApplicationScConf(BaseConf):
    """
    sc-SVC configuration class of  for application usage.

    """
    st_file: str
    sc_ref_file: str
    annotate_mode: Optional[str] = None

    # annotate parameters
    annotate_pot_reg: float = 0.06
    annotate_pot_reg_m: float = 0.015
    annotate_pot_reg_type: str = "entropy"

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
    rec_alpha = 0.5
    rec_match_spot_sum: bool = False

    @property
    def result_dir(self):
        return os.path.join(self.result_root_path, self.sample_name)

    @property
    def st_file_path(self):
        return os.path.join(self.raw_data_path, f"{self.sample_name}_{self.st_file}")

    @property
    def sc_ref_file_path(self):
        return os.path.join(self.raw_data_path, self.sc_ref_file)
