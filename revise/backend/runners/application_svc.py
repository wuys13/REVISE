from revise.backend.runners.base_svc import BaseSVC
from revise.utils.format import warn_if_processed


class ApplicationSVC(BaseSVC):
    """SVC class for application scenarios (real data analysis).
    
    This class provides data validation methods
    suitable for real-world spatial transcriptomics data analysis.
    """

    def _adata_validate(self):
        """Validate that spatial and reference data have overlapping genes.
        
        Raises:
            ValueError: If no common genes are found between
                st_adata and sc_ref_adata.
        """
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
