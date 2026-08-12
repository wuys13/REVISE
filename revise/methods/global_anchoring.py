import logging

import numpy as np
import ot
import pandas as pd
from anndata import AnnData

from revise.conf.base_conf import BaseConf
from revise.methods.base_method import BaseMethod
from revise.tools.distance import bhattacharyya_distance


class GlobalAnchoring(BaseMethod):
    """
    Cell type annotation method for spatial transcriptomics data.
    
    This class provides annotation functionality using either tacco
    or optimal transport (POT) based methods.
    """
    def __init__(self, config: BaseConf, logger: logging.Logger):
        """Initialize Annotate method.
        
        Args:
            config: Configuration object containing annotation parameters
            logger: Logger instance for logging
        """
        super().__init__(config, logger)
        self.mode = self.config.annotate_mode
        self.cell_type_col = self.config.cell_type_col
        self.confidence_col = self.config.confidence_col
        self.unknown_key = self.config.unknown_key

    def run(self, st_adata: AnnData, sc_ref_adata: AnnData, **kwargs):
        """
        Annotate spatial transcriptomics data using single-cell reference.
        
        This method assigns cell type probabilities to each spot in the
        spatial data based on the single-cell reference. Supports two modes:
        - 'tacco': Uses tacco library for annotation
        - 'pot': Uses optimal transport (POT) for annotation
        
        Args:
            st_adata: Spatial transcriptomics AnnData object to annotate
            sc_ref_adata: Single-cell reference AnnData object
            **kwargs: Additional parameters including:
                - annotate_pot_reg: Regularization parameter for POT (if mode='pot')
                - annotate_pot_reg_m: Marginal regularization for POT (if mode='pot')
                - annotate_pot_reg_type: Regularization type for POT (if mode='pot')
                - cell_type_col: Column name for cell type (optional override)
        
        Returns:
            AnnData: Annotated spatial data with cell type probabilities
                stored in obsm[self.cell_type_col] and labels in obs[self.cell_type_col]
        
        Raises:
            ValueError: If required parameters are missing for POT mode
            NotImplementedError: If annotation mode is not 'tacco' or 'pot'
        """
        if "cell_type_col" in kwargs:
            self.cell_type_col = kwargs["cell_type_col"]
        st_adata = st_adata.copy()
        if self.mode == "pot":
            if "annotate_pot_reg" not in kwargs or "annotate_pot_reg_m" not in kwargs or "annotate_pot_reg_type" not in kwargs:
                raise ValueError(
                    f"mode {self.mode} requires 'annotate_pot_reg', 'annotate_pot_reg_m' and 'annotate_pot_reg_type'")

            # calculate cost matrix
            overlap_genes = list(st_adata.var_names.intersection(sc_ref_adata.var_names))
            sc_ref_adata_overlap = sc_ref_adata[:, overlap_genes]
            st_adata_overlap = st_adata[:, overlap_genes]

            dums = pd.get_dummies(sc_ref_adata_overlap.obs[self.cell_type_col], dtype=sc_ref_adata_overlap.X.dtype)
            ncats = dums.sum(axis=0)
            dums /= ncats.to_numpy()
            profiles = sc_ref_adata_overlap.X.T @ dums.to_numpy()
            profiles = pd.DataFrame(profiles, index=sc_ref_adata_overlap.var.index, columns=dums.columns)
            sc_ref_adata_overlap.varm[self.cell_type_col] = profiles
            dist = bhattacharyya_distance(profiles.values.T, st_adata_overlap.X.toarray())

            # calculate margin distribution
            cell_profile_mapping = pd.get_dummies(sc_ref_adata_overlap.obs[self.cell_type_col])
            cell_profile_mapping /= cell_profile_mapping.sum(axis=1).to_numpy()[:, None]
            type_prior = np.array(sc_ref_adata_overlap.X.sum(axis=1)).flatten() @ cell_profile_mapping
            spot_prior = pd.Series(np.array(st_adata_overlap.X.sum(axis=1)).flatten(), index=st_adata_overlap.obs.index)

            # normalize
            spot_prior_sum = spot_prior.sum()
            type_prior_sum = type_prior.sum()
            if spot_prior_sum > 0:
                spot_prior /= spot_prior_sum
            else:
                spot_prior = spot_prior / len(spot_prior)  # uniform distribution if sum is 0
            if type_prior_sum > 0:
                type_prior /= type_prior_sum
            else:
                type_prior = type_prior / len(type_prior)  # uniform distribution if sum is 0

            # Handle NaN values in distance matrix (same as original implementation)
            dist = np.nan_to_num(dist)
            dist_max = dist.max()
            if dist_max <= 0:
                dist_max = 1.0
            dist = dist.T / dist_max

            self.logger.info(f'pot params: reg: {kwargs["annotate_pot_reg"]}, reg_m: {kwargs["annotate_pot_reg_m"]}, reg_type: {kwargs["annotate_pot_reg_type"]}')

            T_transform = ot.unbalanced.sinkhorn_unbalanced(
                spot_prior.values,
                type_prior.values,
                dist,
                reg=kwargs["annotate_pot_reg"],
                reg_m=kwargs["annotate_pot_reg_m"],
                reg_type=kwargs["annotate_pot_reg_type"],
                verbose=True,
                numItermax=5000
            )
            cell_type_ot = pd.DataFrame(T_transform, index=spot_prior.index, columns=type_prior.index)

            cell_type_ot *= (cell_type_ot > 0)
            cell_type_ot /= cell_type_ot.sum(axis=1).to_numpy()[:, None]
            st_adata.obsm[self.cell_type_col] = cell_type_ot
            self.logger.info(st_adata.obsm[self.cell_type_col].value_counts())
            # Assign result to adata.obs
            st_adata = self._assign_result(cell_type_ot, st_adata)
        else:
            import tacco as tc
            st_adata_raw = st_adata.copy()
            tc.tl.annotate(st_adata_raw, sc_ref_adata,
                                   self.cell_type_col, result_key=self.cell_type_col,
                                   multi_center=1, lamb=1e-3)

            # st_adata = tacco_anno(st_adata, sc_ref_adata, self.cell_type_col)
            st_adata.obsm = st_adata_raw.obsm.copy()
            cell_type_ot = st_adata_raw.obsm[self.cell_type_col]
            st_adata = self._assign_result(cell_type_ot, st_adata)

        return st_adata

    def _assign_result(self, cell_type_probability: pd.DataFrame, st_adata: AnnData):
        """
        Assign cell type labels and confidence scores to AnnData.
        
        Args:
            cell_type_probability: DataFrame with cell type probabilities
                (rows: spots, columns: cell types)
            st_adata: Spatial AnnData object to update
        
        Returns:
            AnnData: Updated AnnData with cell type labels and confidence scores
        """
        max_columns = cell_type_probability.idxmax(axis=1)
        max_values = cell_type_probability.max(axis=1)
        st_adata.obs[self.cell_type_col] = max_columns.values.astype(object).copy()
        st_adata.obs[self.cell_type_col].replace({np.nan: self.unknown_key}, inplace=True)
        st_adata.obs[self.confidence_col] = max_values.values.copy()
        return st_adata
