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
from revise.backend.ops.posterior_conditioning import (
    align_posterior_categories,
    condition_cost_matrix,
    get_posterior_matrix,
    posterior_affinity,
    posterior_conditioning_enabled,
    posterior_conditioning_mode,
    posterior_conditioning_strict,
    reference_measure_from_marginals,
)
from revise.backend.ops.shaver import get_prune_adata


def _posterior_columns(adata, key):
    if adata is None or key not in adata.obsm:
        return None
    values = adata.obsm[key]
    if hasattr(values, "columns"):
        return [str(col) for col in values.columns]
    return None


def _label_column_index(columns, label):
    if columns is None:
        return None
    label = str(label)
    candidates = {label, label.replace("/", "_"), label.replace("_", "/")}
    normalized = {col: idx for idx, col in enumerate(columns)}
    normalized.update({col.replace("/", "_"): idx for idx, col in enumerate(columns)})
    normalized.update({col.replace("_", "/"): idx for idx, col in enumerate(columns)})
    for candidate in candidates:
        if candidate in normalized:
            return normalized[candidate]
    return None


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
            adata_sc_all_panel, f"leiden_{self.config.rec_subcluster_resolution}"
        )

        # impute in panel
        in_panel_resolution = getattr(self.config, "rec_in_panel_subcluster_resolution", None)
        if in_panel_resolution is None:
            in_panel_resolution = self.config.rec_subcluster_resolution
        self.svc["sc_svc_impute_in_panel"] = self.local_impute(
            adata_sc_in_panel, f"leiden_{in_panel_resolution}"
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

    def _build_impute_posterior_affinity(
            self,
            ct_adata_sp,
            ct_adata_sc,
            sc_subcluster,
            subcluster_order,
            select_ct,
    ):
        posterior_key = getattr(
            self.config,
            "posterior_conditioning_key",
            self.config.cell_type_col,
        )
        q_spot = get_posterior_matrix(ct_adata_sp, posterior_key)
        if q_spot is None and posterior_key != self.config.cell_type_col:
            posterior_key = self.config.cell_type_col
            q_spot = get_posterior_matrix(ct_adata_sp, posterior_key)
        if q_spot is None:
            msg = (
                "Posterior conditioning requested for impute cell type "
                f"{select_ct} but obsm[{posterior_key!r}] is unavailable"
            )
            if posterior_conditioning_strict(self.config):
                raise ValueError(msg)
            self.logger.warning("%s; falling back to the unconditioned OT objective.", msg)
            return None

        q_type = None
        q_sc = get_posterior_matrix(ct_adata_sc, posterior_key)
        spot_posterior = ct_adata_sp.obsm.get(posterior_key)
        sc_posterior = ct_adata_sc.obsm.get(posterior_key)
        spot_has_labels = isinstance(spot_posterior, pd.DataFrame)
        sc_has_labels = isinstance(sc_posterior, pd.DataFrame)
        if q_sc is not None:
            if spot_has_labels and sc_has_labels:
                q_sc = align_posterior_categories(
                    pd.DataFrame(
                        q_sc,
                        index=sc_posterior.index,
                        columns=sc_posterior.columns,
                    ),
                    spot_posterior.columns,
                    posterior_name="reference posterior columns",
                    reference_name="spot posterior columns",
                ).to_numpy(dtype=np.float64, copy=False)
            elif spot_has_labels:
                # A labeled spot axis plus select_ct can derive the existing
                # one-hot target fallback without trusting reference positions.
                q_sc = None
            elif sc_has_labels:
                msg = (
                    "Posterior conditioning requested for impute cell type "
                    f"{select_ct} but only the reference posterior has labeled "
                    "category axes; spot category identity cannot be aligned"
                )
                if posterior_conditioning_strict(self.config):
                    raise ValueError(msg)
                self.logger.warning(
                    "%s; falling back to the unconditioned OT objective.", msg
                )
                return None
        if q_sc is not None and q_sc.shape[1] == q_spot.shape[1]:
            groups = ct_adata_sc.obs[sc_subcluster].astype(str).to_numpy()
            q_type = np.zeros((len(subcluster_order), q_spot.shape[1]), dtype=np.float64)
            for i, group in enumerate(subcluster_order):
                mask = groups == str(group)
                if np.any(mask):
                    q_type[i] = q_sc[mask].mean(axis=0)
            zero_rows = q_type.sum(axis=1) <= 1e-12
            if np.any(zero_rows):
                q_type[zero_rows] = 1.0 / q_type.shape[1]
        else:
            if spot_has_labels:
                spot_posterior = align_posterior_categories(
                    spot_posterior,
                    spot_posterior.columns,
                    posterior_name="spot posterior columns",
                    reference_name="spot posterior columns",
                )
                columns = list(spot_posterior.columns)
            else:
                columns = _posterior_columns(ct_adata_sp, posterior_key)
            label_idx = _label_column_index(columns, select_ct)
            if label_idx is not None:
                q_type = np.zeros((len(subcluster_order), q_spot.shape[1]), dtype=np.float64)
                q_type[:, label_idx] = 1.0

        if q_type is None:
            msg = (
                "Posterior conditioning requested for impute cell type "
                f"{select_ct} but no compatible target-side posterior can be derived"
            )
            if posterior_conditioning_strict(self.config):
                raise ValueError(msg)
            self.logger.warning("%s; falling back to the unconditioned OT objective.", msg)
            return None

        return posterior_affinity(
            q_spot,
            q_type,
            beta=getattr(self.config, "posterior_conditioning_beta", 1.0),
            min_affinity=getattr(self.config, "posterior_conditioning_min_affinity", 0.05),
        )

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

            ot_cost = np.asarray(dist.T, dtype=np.float64)
            ot_cost = np.nan_to_num(ot_cost, nan=0.0, posinf=0.0, neginf=0.0)
            posterior_affinity_matrix = None
            pc_mode = posterior_conditioning_mode(self.config)
            if pc_mode != "off":
                posterior_affinity_matrix = self._build_impute_posterior_affinity(
                    ct_adata_sp_overlap,
                    ct_adata_sc_overlap,
                    sc_subcluster,
                    type_prior.index,
                    select_ct,
                )
                if (
                    posterior_affinity_matrix is not None
                    and posterior_conditioning_enabled(self.config, "cost")
                ):
                    ot_cost = condition_cost_matrix(
                        ot_cost,
                        posterior_affinity_matrix,
                        getattr(self.config, "posterior_conditioning_cost_strength", 0.2),
                    )

            dist_max = float(np.nanmax(ot_cost)) if ot_cost.size else 1.0
            if not np.isfinite(dist_max) or dist_max <= 0:
                dist_max = 1.0
            sinkhorn_kwargs = {}
            if (
                posterior_affinity_matrix is not None
                and posterior_conditioning_enabled(self.config, "reference")
            ):
                sinkhorn_kwargs["c"] = reference_measure_from_marginals(
                    spot_prior.values,
                    type_prior.values,
                    posterior_affinity_matrix,
                )

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
                reference_measure=sinkhorn_kwargs.get("c"),
                event_callback=getattr(self.config, "ot_event_callback", None),
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
