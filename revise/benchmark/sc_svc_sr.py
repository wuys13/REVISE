import numpy as np
import scanpy as sc

from revise.benchmark.benchmark_svc import BenchmarkSVC
from revise.methods.spot_sr import SpotSr
from revise.tools.meta import construct_sc_ref
from revise.tools.meta import get_sc_obs
from revise.tools.meta import get_true_cell_type


class ScSVCSr(BenchmarkSVC):
    """
    sc-SVC super-resolution for benchmark CFs: spot size/ batch effect.
    
    This class reconstructs single-cell resolution expression profiles
    from spatial transcriptomics data by redistributing spot-level
    expressions to virtual cells using cell type contributions.
    """
    def __init__(self, st_adata, sc_ref_adata, config, real_st_adata, logger):
        super().__init__(st_adata, sc_ref_adata, config, real_st_adata, logger)
        self._adata_validate()
        self._adata_validate_dec()
        self._adata_processing()
        self.svc_obs = self._get_svc_obs()
        self.spot_sr = SpotSr(self.config, self.logger)
        self.svc = {}

    def _adata_validate_dec(self):

        assert "all_cells_in_spot" in self.st_adata.uns, "spot-sc mapping is not in st_adata.uns"

    def _get_svc_obs(self):
        svc_obs = get_sc_obs(self.st_adata.obs.index, self.st_adata.uns['all_cells_in_spot'], self.st_adata.obsm["spatial"])
        svc_obs = get_true_cell_type(svc_obs, self.real_st_adata)
        return svc_obs

    def local_refinement(self, *args, **kwargs):
        """Reconstruct single-cell expression profiles from spot-level data.

        1. Assigns cell types to each virtual cell using SpotSr
        2. Constructs cell type reference profiles
        3. Calculates gene expression for each cell based on spot contributions
        4. Normalizes expressions to 10,000 counts per cell
        
        The reconstructed data is stored in self.svc["sc_svc_dec"].
        """
        overlap_genes = list(self.st_adata.var_names.intersection(self.sc_ref_adata.var_names))
        st_adata_common = self.st_adata[:, overlap_genes]
        sc.pp.normalize_total(st_adata_common, target_sum=1e4)
        cell_contributions = st_adata_common.obsm["Level1"].values if hasattr(st_adata_common.obsm["Level1"], 'values') else st_adata_common.obsm["Level1"]

        sc.pp.normalize_total(self.st_adata, target_sum=1e4)
        sc.pp.normalize_total(self.sc_ref_adata, target_sum=1e4)
        self.spot_sr.run(self)
        key_type = "clusters"
        if key_type not in self.sc_ref_adata.obs.columns:
            self.sc_ref_adata.obs[key_type] = self.sc_ref_adata.obs["Level1"].astype(str)
        type_list = sorted(list(self.sc_ref_adata.obs[key_type].unique().astype(str)))
        self.logger.info(f'There are {len(type_list)} cell types: {type_list}')
        sc_ref_all = construct_sc_ref(self.sc_ref_adata, key_type=key_type, type_list=type_list)
        sc_ref_all = sc_ref_all.loc[:, overlap_genes]
        type_list = sorted(list(sc_ref_all.index))
        spots = self.svc_obs['spot_name'].unique()

        spot_to_idx = {spot: idx for idx, spot in enumerate(spots)}
        self.logger.info("Using simple allocation method...")

        adata_spot = st_adata_common.copy()
        X = adata_spot.X if type(adata_spot.X) is np.ndarray else adata_spot.X.toarray()

        Y = cell_contributions[:, np.newaxis, :] * sc_ref_all.values.T
        Y = Y / (np.sum(Y, axis=2, keepdims=True) + 1e-10)
        Y = Y * X[:, :, np.newaxis]

        spot_indices = np.array([spot_to_idx[spot] for spot in self.svc_obs['spot_name']])
        type_indices = np.array([type_list.index(t) for t in self.svc_obs['cell_type']])

        SVC_X = Y[spot_indices, :, type_indices]
        self.logger.info(f"Extracted SVC expressions using simple allocation method")

        SVC_X = SVC_X / (np.sum(SVC_X, axis=1, keepdims=True) + 1e-10) * 1e4

        self.logger.info(f"Number of cells processed: {len(self.svc_obs)}")
        self.logger.info(f"Number of unique spots: {len(spots)}")
        self.logger.info(f"Shape of SVC_X: {SVC_X.shape}")

        svc_adata = sc.AnnData(SVC_X)
        svc_adata.var_names = st_adata_common.var_names
        svc_adata.obs = self.svc_obs.copy()
        svc_adata.obs.set_index("cell_id", inplace=True)
        self.svc["sc_svc_dec"] = svc_adata
