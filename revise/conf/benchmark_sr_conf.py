import os
from dataclasses import dataclass

from revise.conf.base_conf import BaseConf


@dataclass
class BenchmarkSrConf(BaseConf):
    """
    sc-SVC configuration class of  for spot_size/batch_effect benchmark.

    """
    st_file: str
    gt_svc_file: str
    sc_ref_file: str
    spot_size: int
    annotate_mode: str

    # annotate parameters
    annotate_pot_reg: float = 0.01
    annotate_pot_reg_m: float = 0.0001
    annotate_pot_reg_type: str = "kl"

    # svc parameters
    svc_completeness: bool = True

    @property
    def result_dir(self):
        return os.path.join(self.result_root_path, self.sample_name, f"spot_{self.spot_size}")

    @property
    def st_file_path(self):
        return os.path.join(self.raw_data_path, self.sample_name, f"spot_{self.spot_size}", self.st_file)

    @property
    def gt_svc_file_path(self):
        return os.path.join(self.raw_data_path, self.sample_name, self.gt_svc_file)

    @property
    def sc_ref_file_path(self):
        return os.path.join(self.raw_data_path, self.sample_name, self.sc_ref_file)

    @property
    def pm_on_cell_file(self):
        return os.path.join(os.path.dirname(os.path.join(self.raw_data_path, self.sample_name)), "PM_on_cell.csv")
