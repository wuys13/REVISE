from revise.base_svc import BaseSVC
from revise.methods.global_anchoring import GlobalAnchoring


class BaseSVCAnchor(BaseSVC):
    """Base class for SVC methods that include annotation functionality.
    
    This class provides a default implementation of the annotate method
    using the Annotate method class. Subclasses should implement the
    reconstruct method.
    """
    def __init__(self, st_adata, sc_ref_adata, config, real_st_adata, logger):
        """Initialize BaseSVCAnchor.
        
        Args:
            st_adata: Spatial transcriptomics AnnData object
            sc_ref_adata: Single-cell reference AnnData object
            config: Configuration object containing method parameters
            real_st_adata: Ground truth spatial data (for benchmarking, can be None)
            logger: Logger instance for logging
        """
        super().__init__(st_adata, sc_ref_adata, config, real_st_adata, logger)
        self.annotate_method = GlobalAnchoring(config=config, logger=logger)

    def global_anchoring(self, *args, **kwargs):
        """Annotate spatial spots using the configured annotation method.
        
        This method uses the Annotate class to assign cell type labels
        to spatial spots based on the single-cell reference.
        """
        self.st_adata = self.annotate_method.run(
            self.st_adata,
            self.sc_ref_adata,
            **self.config.__dict__
        )

    def local_refinement(self, *args, **kwargs):
        """Reconstruct single-cell resolution expression profiles.
        
        This method must be implemented by subclasses.
        """
        raise NotImplementedError("Reconstruct method not implemented.")
