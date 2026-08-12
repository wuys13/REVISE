import os.path

import numpy as np
import ot
import pandas as pd
import scanpy as sc
from tqdm import tqdm

from revise.benchmark.benchmark_svc import BenchmarkSVC
from revise.methods.gene_impute import GeneImpute
from revise.methods.gene_uncertainty import GeneUncertainty
from revise.tools.distance import bhattacharyya_distance
from revise.tools.meta import get_subcluster
from revise.tools.meta import merge_subcluster
from revise.tools.shaver import get_prune_adata


class ScSVCImpute(BenchmarkSVC):
    """
    Single-cell SVC imputation for benchmark CFs: gene panel/gene dropout.
    
    This class performs gene imputation by comparing in-panel vs all-panel
    HVG selection strategies and using optimal transport for imputation.
    """
    def __init__(self, st_adata, sc_ref_adata, config, real_st_adata, logger):
        super().__init__(st_adata, sc_ref_adata, config, real_st_adata, logger)
        self._adata_validate()
        self._adata_processing()
        self._adata_processing_impute()
        self.gene_uncertainty = GeneUncertainty(self.config, self.logger)
        self.gene_impute = GeneImpute(self.config, self.logger)
        self.svc = {}

    def _adata_processing_impute(self):
        """
        Process data with transcript count filtering.
        """
        if "cell_id" in self.st_adata.obs.columns:
            self.st_adata.obs_names = self.st_adata.obs["cell_id"]
        self.st_adata = self.st_adata[self.st_adata.obs['transcript_counts'] >= self.config.prep_min_counts, :]
        sc.pp.filter_genes(self.st_adata, min_cells=self.config.prep_min_cells)

        self.sc_ref_adata.obs = self.sc_ref_adata.obs[[self.config.cell_type_col]]
        sc.pp.filter_genes(self.sc_ref_adata, min_cells=self.config.prep_min_cells)
        self.sc_ref_adata.obs[self.config.cell_type_col].replace({"Mono/Macro": "Mono_Macro"}, inplace=True)

    def local_refinement(self, *args, **kwargs):
        """
        Reconstruct expression profiles using gene imputation.

        1. Evaluates gene uncertainty comparing in-panel vs all-panel strategies
        2. Generates subclustered single-cell data for both strategies
        3. Performs local imputation for each cell type using optimal transport
        4. Optionally prunes imputed data
        
        Results are stored in:
        - self.svc["sc_svc_impute_all_panel"]: Imputation using all-panel strategy
        - self.svc["sc_svc_impute_in_panel"]: Imputation using in-panel strategy
        """
        overlap_genes = list(self.st_adata.var_names.intersection(self.sc_ref_adata.var_names))
        assert len(overlap_genes) > 0, "overlap genes not found"
        gene_compare_file = os.path.join(self.config.result_dir, "compare_in_vs_all_panel_moranI.csv")
        if not os.path.exists(gene_compare_file):
            compare_df = self.gene_uncertainty.run(self.sc_ref_adata, overlap_genes)
            compare_df.to_csv(gene_compare_file)
        else:
            compare_df = pd.read_csv(gene_compare_file, index_col=0)
        compare_df = compare_df[~compare_df['test']]

        in_panel_file = os.path.join(self.config.result_dir, "adata_sc_in_panel.h5ad")
        all_panel_file = os.path.join(self.config.result_dir, "adata_sc_all_panel.h5ad")
        if os.path.exists(in_panel_file) and os.path.exists(all_panel_file):
            self.logger.info(f"Load {in_panel_file} and {all_panel_file}")
            adata_sc_in_panel = sc.read(in_panel_file)
            adata_sc_all_panel = sc.read(all_panel_file)
        else:
            self.logger.info(f"Build {in_panel_file} and {all_panel_file}")
            adata_sc_all_panel, adata_sc_in_panel = get_subcluster(
                self.sc_ref_adata,
                compare_df,
                celltype_col=self.config.cell_type_col)
            adata_sc_in_panel.write(in_panel_file)
            adata_sc_all_panel.write(all_panel_file)

        self.svc["sc_svc_impute_all_panel"] = self.local_impute(
            adata_sc_all_panel, f"leiden_{self.config.rec_subcluster_resolution}"
        )

        # impute in panel
        self.svc["sc_svc_impute_in_panel"] = self.local_impute(
            adata_sc_in_panel, f"leiden_{self.config.rec_subcluster_resolution}"
        )
        # self.st_adata = self.st_adata[self.svc["sc_svc_impute_in_panel"].obs_name, :]

        # metrics_in_panel = compute_metric(
        #     adata_to_metric, adata_sp_impute_in_panel, self.logger,
        #     adata_process=False,
        #     gene_list=gene_list,
        #     normalize=True
        # )
        # metrics_in_panel.to_csv(os.path.join(self.config.metric_dir, f"metrics_in_panel.csv"))

    def local_impute(
            self,
            adata_sc,
            sc_subcluster
    ):
        """
        Perform local imputation for each cell type using subclustered reference.
        
        Args:
            adata_sc: Subclustered single-cell reference AnnData
            sc_subcluster: Column name in adata_sc.obs containing subcluster labels
            
        Returns:
            AnnData: Imputed spatial data with reconstructed expressions

        1. Processes each cell type separately
        2. Computes subcluster profiles and distances
        3. Uses optimal transport to find spot-subcluster mappings
        4. Imputes gene expressions using OT coupling weights
        """
        adata_sp = self.st_adata.copy()
        adata_sc = adata_sc.copy()
        cts = list(adata_sc.obs[self.config.cell_type_col].unique())
        adata_sp_cts = []
        for select_ct in tqdm(cts, "Imputation by cell type"):
            self.logger.info(f"Conducting cell type: {select_ct} ........")
            ct_adata_sc = adata_sc[adata_sc.obs[self.config.cell_type_col] == select_ct].copy()
            ct_adata_sp = adata_sp[adata_sp.obs[self.config.cell_type_col] == select_ct].copy()

            overlap_genes = ct_adata_sc.var_names.intersection(ct_adata_sp.var_names)
            ct_adata_sc_overlap = ct_adata_sc[:, overlap_genes].copy()
            ct_adata_sp_overlap = ct_adata_sp[:, overlap_genes].copy()
            dums = pd.get_dummies(ct_adata_sc_overlap.obs[sc_subcluster],
                                  dtype=ct_adata_sc_overlap.X.dtype)
            ncats = dums.sum(axis=0)
            dums /= ncats.to_numpy()
            profiles = ct_adata_sc_overlap.X.T @ dums.to_numpy()
            profiles = pd.DataFrame(profiles, index=ct_adata_sc_overlap.var.index, columns=dums.columns)
            ct_adata_sc_overlap.varm[sc_subcluster] = profiles

            dist = bhattacharyya_distance(profiles.values.T, ct_adata_sp_overlap.X.toarray())

            cell_profile_mapping = pd.get_dummies(ct_adata_sc_overlap.obs[sc_subcluster])
            cell_profile_mapping /= cell_profile_mapping.sum(axis=1).to_numpy()[:, None]
            type_prior = np.array(ct_adata_sc_overlap.X.sum(axis=1)).flatten() @ cell_profile_mapping
            spot_prior = pd.Series(np.array(ct_adata_sp_overlap.X.sum(axis=1)).flatten(),
                                   index=ct_adata_sp_overlap.obs.index)

            spot_prior /= spot_prior.sum()
            type_prior /= type_prior.sum()

            T_matrix = ot.unbalanced.sinkhorn_unbalanced(
                spot_prior.values,
                type_prior.values,
                dist.T / dist.max(),
                reg=self.config.rec_impute_pot_reg,
                reg_m=self.config.rec_impute_pot_reg_m,
                reg_type=self.config.rec_impute_pot_reg_type,
                verbose=True,
                numItermax=5000
            )
            T_matrix = pd.DataFrame(T_matrix, index=spot_prior.index, columns=type_prior.index)

            ct_adata_sc = merge_subcluster(
                ct_adata_sc,
                subcluster=sc_subcluster,
                mode=self.config.rec_merge_subcluster_method
            )

            overlap_genes = ct_adata_sp.var_names.intersection(ct_adata_sc.var_names)
            adata_sp_impute = self.gene_impute.run(
                ct_adata_sp, ct_adata_sc,
                genes_to_predict=overlap_genes,
                neighbor_weights=T_matrix,
            )
            adata_sp_cts.append(adata_sp_impute)

        adata_sp_impute = sc.concat(adata_sp_cts)
        if self.config.rec_impute_prune_flag:
            adata_sp_impute = get_prune_adata(adata_sp_impute)

        return adata_sp_impute
