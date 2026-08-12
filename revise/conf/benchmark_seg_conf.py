import os
from dataclasses import dataclass

from revise.conf.base_conf import BaseConf


@dataclass
class BenchmarkSegConf(BaseConf):
    """
    sp-SVC configuration class of  for segmentation/bin2cell benchmark.

    """
    st_file: str
    gt_svc_file: str
    sc_ref_file: str
    seg_method: str
    annotate_mode: str

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
    rec_alpha: float = 1.0

    @property
    def result_dir(self):
        return os.path.join(self.result_root_path, self.sample_name, self.seg_method)

    @property
    def st_file_path(self):
        return os.path.join(self.raw_data_path, self.sample_name, self.seg_method, self.st_file)

    @property
    def gt_svc_file_path(self):
        return os.path.join(self.raw_data_path, self.sample_name, self.gt_svc_file)

    @property
    def sc_ref_file_path(self):
        return os.path.join(self.raw_data_path, self.sample_name, self.sc_ref_file)
