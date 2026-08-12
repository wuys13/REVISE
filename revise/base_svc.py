from abc import ABC
from abc import abstractmethod

from anndata import AnnData

from revise.tools.log import ensure_logger


class BaseSVC(ABC):
    """Base class for Spatial Virtual Cell (SVC) reconstruction methods.
    
    This abstract base class defines the interface for SVC methods that
    annotate spatial transcriptomics data and reconstruct single-cell
    resolution expression profiles.
    """
    def __init__(self, st_adata: AnnData, sc_ref_adata: AnnData, config, real_st_adata, logger=None):
        """Initialize BaseSVC.
        
        Args:
            st_adata: Spatial transcriptomics AnnData object
            sc_ref_adata: Single-cell reference AnnData object
            config: Configuration object containing method parameters
            real_st_adata: Ground truth spatial data (for benchmarking, can be None)
            logger: Logger instance for logging; falls back to the global default when None
        """
        self.st_adata = st_adata
        self.sc_ref_adata = sc_ref_adata
        self.config = config
        self.real_st_adata = real_st_adata
        self.logger = ensure_logger(logger)

    def _validate_inputs(self):
        """
        Validate the format of input AnnData objects.
        
        This method should be overridden by subclasses to implement
        specific validation logic for their use cases.
        """
        pass

    @abstractmethod
    def global_anchoring(self, *args, **kwargs):
        """
        Annotate spatial spots with cell type labels via global anchoring.
        
        This method assigns cell type probabilities or labels to each spot
        in the spatial transcriptomics data based on the single-cell reference.
        
        Args:
            *args: Variable positional arguments (implementation-specific)
            **kwargs: Variable keyword arguments (implementation-specific)
            
        Returns:
            Implementation-specific return value (typically None, modifies self.st_adata)
        """
        raise NotImplementedError("Annotate method not implemented.")

    @abstractmethod
    def local_refinement(self, *args, **kwargs):
        """
        Reconstruct single-cell resolution expression profiles.
        
        This method generates virtual single-cell expression profiles from
        spatial transcriptomics data, typically stored in self.svc dictionary.
        
        Args:
            *args: Variable positional arguments (implementation-specific)
            **kwargs: Variable keyword arguments (implementation-specific)
            
        Returns:
            Implementation-specific return value (typically None, stores results in self.svc)
        """
        raise NotImplementedError("Reconstruct method not implemented.")
