import os
from dataclasses import dataclass

from revise.conf.base_conf import BaseConf


@dataclass
class BenchmarkImputeConf(BaseConf):
    """
    sc-SVC configuration class of  for gene_dropout/gene_panel benchmark.

    """
    st_file: str
    gt_svc_file: str
    sc_ref_file: str
    annotate_mode: str

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

    # reconstruct impute
    rec_merge_subcluster_method: str = "mean"
    rec_subcluster_resolution: int = 3
    rec_impute_prune_flag: bool = True
    rec_impute_n_neighbors: int = 1
    rec_impute_method: str = "mean"

    @property
    def result_dir(self):
        return os.path.join(self.result_root_path, self.sample_name)

    @property
    def st_file_path(self):
        return os.path.join(self.raw_data_path, self.sample_name, self.st_file)

    @property
    def sc_ref_file_path(self):
        return os.path.join(self.raw_data_path, self.sample_name, self.sc_ref_file)

    @property
    def gt_svc_file_path(self):
        return os.path.join(self.raw_data_path, self.sample_name, self.gt_svc_file)
