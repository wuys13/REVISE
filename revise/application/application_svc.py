import scanpy as sc

from revise.base_svc_anchor import BaseSVCAnchor
from revise.tools.format import warn_if_processed


class ApplicationSVC(BaseSVCAnchor):
    """SVC class for application scenarios (real data analysis).
    
    This class provides data preprocessing and validation methods
    suitable for real-world spatial transcriptomics data analysis.
    """

    def _adata_processing(self):
        """Preprocess spatial and single-cell reference data.
        
        This method filters cells and genes, and normalizes cell type
        column names by replacing '/' with '_' to avoid issues in
        downstream processing.
        """
        sc.pp.filter_cells(self.st_adata, min_counts=self.config.prep_st_min_counts)
        sc.pp.filter_genes(self.st_adata, min_counts=self.config.prep_st_min_cells)

        sc.pp.filter_cells(self.sc_ref_adata, min_counts=self.config.prep_sc_min_counts)
        sc.pp.filter_genes(self.sc_ref_adata, min_counts=self.config.prep_sc_min_cells)
        replace_columns = {k: k.replace("/", "_") for k in self.sc_ref_adata.obs['Level1'].unique().tolist() if '/' in k}
        self.sc_ref_adata.obs['Level1'].replace(replace_columns, inplace=True)
        replace_columns = {k: k.replace("/", "_") for k in self.sc_ref_adata.obs['Level2'].unique().tolist() if '/' in k}
        self.sc_ref_adata.obs['Level2'].replace(replace_columns, inplace=True)

    def _adata_validate(self):
        """Validate that spatial and reference data have overlapping genes.
        
        Raises:
            AssertionError: If no common genes are found between
                st_adata and sc_ref_adata.
        """
        assert len(self.st_adata.var_names.intersection(self.sc_ref_adata.var_names)) > 0, "st_adata and sc_ref_adata have no common genes"
        warn_if_processed(self.st_adata, self.logger, "st_adata")
        warn_if_processed(self.sc_ref_adata, self.logger, "sc_ref_adata")
