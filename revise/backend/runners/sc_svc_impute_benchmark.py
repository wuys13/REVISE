import os
import os.path
import shutil

import numpy as np
import pandas as pd
import scanpy as sc
from tqdm import tqdm

from revise.backend.runners.benchmark_svc import BenchmarkSVC
from revise.backend.kernels import GeneImputeKernel as GeneImpute
from revise.backend.kernels import GeneUncertaintyKernel as GeneUncertainty
from revise.backend.ops.distance import bhattacharyya_distance
from revise.backend.ops.local_ot import solve_local_ot
from revise.backend.ops.meta import get_subcluster
from revise.backend.ops.meta import merge_subcluster
from revise.backend.ops.shaver import get_prune_adata


def _align_coupling_to_merged_subclusters(
    coupling,
    merged_ref,
    subcluster_key,
):
    """Align coupling columns to merged-reference rows by subcluster label."""
    if not isinstance(coupling, pd.DataFrame):
        raise ValueError("coupling must be a pandas DataFrame with labeled columns")
    if subcluster_key not in merged_ref.obs.columns:
        raise ValueError(
            f"merged reference is missing obs[{subcluster_key!r}]"
        )

    coupling_labels = pd.Index(coupling.columns)
    merged_labels = pd.Index(merged_ref.obs[subcluster_key].to_numpy())
    for name, labels in (
        ("coupling columns", coupling_labels),
        (f"merged obs[{subcluster_key!r}]", merged_labels),
    ):
        if bool(pd.isna(labels).any()):
            raise ValueError(f"{name} contain null subcluster labels")
        duplicated = labels[labels.duplicated()].tolist()
        if duplicated:
            raise ValueError(f"{name} contain duplicate subcluster labels: {duplicated}")

    coupling_set = set(coupling_labels.tolist())
    merged_set = set(merged_labels.tolist())
    missing = [label for label in merged_labels if label not in coupling_set]
    extra = [label for label in coupling_labels if label not in merged_set]
    if missing or extra:
        raise ValueError(
            "coupling and merged subcluster labels differ: "
            f"missing={missing}, extra={extra}"
        )

    aligned = coupling.reindex(columns=merged_labels.tolist()).copy()
    aligned.columns = pd.Index(merged_ref.obs_names)
    return aligned


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
        self._materialize_cached_gene_compare(gene_compare_file)
        if not os.path.exists(gene_compare_file):
            compare_df = self.gene_uncertainty.run(self.sc_ref_adata, overlap_genes)
            compare_df.to_csv(gene_compare_file)
        else:
            compare_df = pd.read_csv(gene_compare_file, index_col=0)
        compare_df = compare_df[~compare_df['test']]

        in_panel_file = os.path.join(self.config.result_dir, "adata_sc_in_panel.h5ad")
        all_panel_file = os.path.join(self.config.result_dir, "adata_sc_all_panel.h5ad")
        self._materialize_cached_subcluster(in_panel_file, all_panel_file)
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
            adata_sc_all_panel,
            f"leiden_{self.config.rec_subcluster_resolution}",
        )

        # impute in panel
        in_panel_resolution = getattr(self.config, "rec_in_panel_subcluster_resolution", None)
        if in_panel_resolution is None:
            in_panel_resolution = self.config.rec_subcluster_resolution
        self.svc["sc_svc_impute_in_panel"] = self.local_impute(
            adata_sc_in_panel,
            f"leiden_{in_panel_resolution}",
        )
        # self.st_adata = self.st_adata[self.svc["sc_svc_impute_in_panel"].obs_name, :]

        # metrics_in_panel = compute_metric(
        #     adata_to_metric, adata_sp_impute_in_panel, self.logger,
        #     adata_process=False,
        #     gene_list=gene_list,
        #     normalize=True
        # )
        # metrics_in_panel.to_csv(os.path.join(self.config.metric_dir, f"metrics_in_panel.csv"))

    def _materialize_cached_gene_compare(self, target_file: str) -> None:
        """
        Optionally hydrate compare CSV from cache to avoid re-running uncertainty.

        This path is only enabled when REVISE_GENE_COMPARE_CACHE is set and
        target_file does not already exist.
        """
        if os.path.exists(target_file):
            return
        cache_file = os.environ.get("REVISE_GENE_COMPARE_CACHE")
        if not cache_file:
            return
        if not os.path.exists(cache_file):
            self.logger.warning("REVISE_GENE_COMPARE_CACHE does not exist: %s", cache_file)
            return
        shutil.copyfile(cache_file, target_file)
        self.logger.info("Loaded cached compare file for impute benchmark: %s -> %s", cache_file, target_file)

    def _materialize_cached_subcluster(self, in_panel_file: str, all_panel_file: str) -> None:
        """
        Optionally hydrate subcluster AnnData files from the compare cache directory.
        """
        if os.path.exists(in_panel_file) and os.path.exists(all_panel_file):
            return
        cache_file = os.environ.get("REVISE_GENE_COMPARE_CACHE")
        if not cache_file:
            return
        cache_dir = os.path.dirname(cache_file)
        src_in_panel = os.path.join(cache_dir, "adata_sc_in_panel.h5ad")
        src_all_panel = os.path.join(cache_dir, "adata_sc_all_panel.h5ad")
        if not (os.path.exists(src_in_panel) and os.path.exists(src_all_panel)):
            return
        if not os.path.exists(in_panel_file):
            shutil.copyfile(src_in_panel, in_panel_file)
        if not os.path.exists(all_panel_file):
            shutil.copyfile(src_all_panel, all_panel_file)
        self.logger.info(
            "Loaded cached subcluster files for impute benchmark: %s , %s -> %s",
            src_in_panel,
            src_all_panel,
            self.config.result_dir,
        )

    def local_impute(
            self,
            adata_sc,
            sc_subcluster,
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
            if isinstance(
                ct_adata_sc.obs[sc_subcluster].dtype,
                pd.CategoricalDtype,
            ):
                ct_adata_sc.obs[sc_subcluster] = (
                    ct_adata_sc.obs[sc_subcluster]
                    .cat.remove_unused_categories()
                )
            if ct_adata_sc.n_obs == 0:
                continue
            if ct_adata_sp.n_obs == 0:
                continue
            overlap_genes = ct_adata_sc.var_names.intersection(ct_adata_sp.var_names)
            if len(overlap_genes) == 0:
                continue
            ct_adata_sc_overlap = ct_adata_sc[:, overlap_genes].copy()
            ct_adata_sp_overlap = ct_adata_sp[:, overlap_genes].copy()
            dums = pd.get_dummies(ct_adata_sc_overlap.obs[sc_subcluster],
                                  dtype=ct_adata_sc_overlap.X.dtype)
            ncats = dums.sum(axis=0)
            if dums.shape[1] == 0 or np.any(
                np.asarray(ncats, dtype=np.float64) <= 0
            ):
                continue
            dums /= ncats.to_numpy()
            profiles = ct_adata_sc_overlap.X.T @ dums.to_numpy()
            profiles = pd.DataFrame(profiles, index=ct_adata_sc_overlap.var.index, columns=dums.columns)
            ct_adata_sc_overlap.varm[sc_subcluster] = profiles

            cell_profile_mapping = pd.get_dummies(ct_adata_sc_overlap.obs[sc_subcluster])
            cell_profile_mapping /= cell_profile_mapping.sum(axis=1).to_numpy()[:, None]
            target_mass = np.asarray(
                np.array(ct_adata_sc_overlap.X.sum(axis=1)).flatten()
                @ cell_profile_mapping,
                dtype=np.float64,
            )
            source_mass = np.asarray(
                np.array(ct_adata_sp_overlap.X.sum(axis=1)).flatten(),
                dtype=np.float64,
            )
            if not np.all(np.isfinite(source_mass)):
                continue
            if not np.all(np.isfinite(target_mass)):
                continue
            if np.any(source_mass <= 0):
                continue
            if np.any(target_mass <= 0):
                continue
            type_prior = pd.Series(
                target_mass,
                index=cell_profile_mapping.columns,
            )
            spot_prior = pd.Series(
                source_mass,
                index=ct_adata_sp_overlap.obs.index,
            )

            spot_prior /= spot_prior.sum()
            type_prior /= type_prior.sum()

            dist = bhattacharyya_distance(
                profiles.values.T,
                ct_adata_sp_overlap.X.toarray(),
            )
            ot_cost = np.asarray(dist.T, dtype=np.float64)
            ot_cost = np.nan_to_num(ot_cost, nan=0.0, posinf=0.0, neginf=0.0)

            dist_max = float(np.nanmax(ot_cost)) if ot_cost.size else 1.0
            if not np.isfinite(dist_max) or dist_max <= 0:
                dist_max = 1.0
            T_values = solve_local_ot(
                spot_prior.values,
                type_prior.values,
                ot_cost / dist_max,
                method=self.config.rec_ot_method,
                pot_reg=self.config.rec_impute_pot_reg,
                pot_reg_m=self.config.rec_impute_pot_reg_m,
                pot_reg_type=self.config.rec_impute_pot_reg_type,
                pot_verbose=True,
                pot_num_iter_max=5000,
                reference_measure=None,
                event_callback=getattr(
                    self.config,
                    "ot_event_callback",
                    None,
                ),
            )

            T_matrix = pd.DataFrame(
                T_values,
                index=spot_prior.index,
                columns=type_prior.index,
            )

            ct_adata_sc = merge_subcluster(
                ct_adata_sc,
                subcluster=sc_subcluster,
                mode=self.config.rec_merge_subcluster_method
            )
            T_matrix = _align_coupling_to_merged_subclusters(
                T_matrix,
                ct_adata_sc,
                sc_subcluster,
            )

            overlap_genes = ct_adata_sp.var_names.intersection(
                ct_adata_sc.var_names
            )
            adata_sp_impute = self.gene_impute.run(
                ct_adata_sp,
                ct_adata_sc,
                genes_to_predict=overlap_genes,
                neighbor_weights=T_matrix,
            )
            adata_sp_cts.append(adata_sp_impute)

        if not adata_sp_cts:
            return adata_sp[:0].copy()
        adata_sp_impute = sc.concat(adata_sp_cts)
        if self.config.rec_impute_prune_flag:
            adata_sp_impute = get_prune_adata(adata_sp_impute)

        return adata_sp_impute
