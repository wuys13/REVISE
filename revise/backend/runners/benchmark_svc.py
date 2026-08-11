import numpy as np
from scipy.sparse import issparse

from revise.backend.runners.base_svc import BaseSVC
from revise.utils.format import warn_if_processed


class BenchmarkSVC(BaseSVC):
    """
    Base class for SVC methods in benchmark scenarios.
    
    This class provides data processing methods suitable for benchmark
    evaluation, including handling of ground truth data (real_st_adata).
    """

    def _adata_processing(self):
        """
        Process data for benchmark scenarios.
        
        Makes variable names unique and computes total_counts for
        real_st_adata if not already present.
        """
        self.st_adata.var_names_make_unique()
        if 'total_counts' not in self.st_adata.obs.columns:
            s = self.real_st_adata.X.sum(axis=1)
            # If scipy sparse matrix row sum returns sparse/matrix object, try using .A1; otherwise convert to ndarray and flatten
            if issparse(self.real_st_adata.X) and hasattr(s, "A1"):
                total = s.A1
            else:
                total = np.asarray(s).ravel()
            self.real_st_adata.obs['total_counts'] = total

    def _adata_validate(self):
        if self.real_st_adata is None:
            raise ValueError(
                "Invalid reconstruction input: "
                f"gt_svc_file_path={getattr(self.config, 'gt_svc_file_path', '<unavailable>')}; "
                "field=real_st_adata; expected=loaded; actual=None"
            )
        overlap_count = len(
            self.st_adata.var_names.intersection(self.sc_ref_adata.var_names)
        )
        if overlap_count == 0:
            raise ValueError(
                "Invalid reconstruction input: "
                f"st_file_path={getattr(self.config, 'st_file_path', '<unavailable>')}; "
                f"sc_ref_file_path={getattr(self.config, 'sc_ref_file_path', '<unavailable>')}; "
                "field=var_names_overlap; expected=>=1; "
                f"actual={overlap_count}"
            )
        warn_if_processed(self.st_adata, self.logger, "st_adata")
        warn_if_processed(self.sc_ref_adata, self.logger, "sc_ref_adata")
        if self.real_st_adata is not None:
            warn_if_processed(self.real_st_adata, self.logger, "real_st_adata")
