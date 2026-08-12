import os
from dataclasses import dataclass
from dataclasses import field

from revise.conf.base_conf import BaseConf


@dataclass
class ApplicationSpConf(BaseConf):
    """
    sp-SVC configuration class of  for application usage.

    """
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
    rec_alpha = 0.5

    @property
    def result_dir(self):
        return os.path.join(self.result_root_path, self.sample_name)

    @property
    def st_file_path(self):
        return os.path.join(self.raw_data_path, f"{self.sample_name}_{self.st_file}")

    @property
    def sc_ref_file_path(self):
        return os.path.join(self.raw_data_path, self.sc_ref_file)
