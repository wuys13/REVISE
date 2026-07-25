import json
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
from revise.backend.ops.assignment import (
    AssignmentState,
    AssignmentStateError,
    aggregate_assignment,
    align_assignment_categories,
    align_assignment_observations,
    one_hot_assignment,
    validate_assignment,
)
from revise.backend.ops.assignment_guidance import (
    assignment_guidance_mode,
    assignment_compatibility,
    ot_cost_guidance,
    resolve_assignment_guidance,
)
from revise.backend.ops.distance import bhattacharyya_distance
from revise.backend.ops.local_ot import solve_local_ot
from revise.backend.ops.meta import get_subcluster
from revise.backend.ops.meta import merge_subcluster
from revise.backend.ops.posterior_conditioning import (
    posterior_conditioning_mode,
    reference_measure_from_marginals,
)
from revise.backend.ops.shaver import get_prune_adata


def _guidance_mode(config):
    return assignment_guidance_mode(config)


def _assignment_categories(*adatas, key):
    categories = []
    for adata in adatas:
        if adata is None:
            continue
        posterior = adata.obsm.get(key)
        if isinstance(posterior, pd.DataFrame):
            candidates = list(posterior.columns)
        elif key in adata.obs:
            candidates = adata.obs[key].tolist()
        else:
            candidates = []
        for value in candidates:
            if value not in categories:
                categories.append(value)
    return categories


def _assignment_state_from_adata(adata, *, key, category_labels):
    if key in adata.obsm:
        raw = adata.obsm[key]
        if not isinstance(raw, pd.DataFrame):
            raise AssignmentStateError("category_labels_missing")
        state = validate_assignment(
            AssignmentState(
                values=raw.to_numpy(dtype=np.float64),
                observation_labels=raw.index,
                category_labels=raw.columns,
                source=f"obsm[{key}]",
                level=str(key),
                value_semantics="soft",
                lineage=[
                    {
                        "operation": "load",
                        "container": "obsm",
                        "key": str(key),
                    }
                ],
            )
        )
        return align_assignment_observations(state, adata.obs_names)
    if key not in adata.obs:
        return None
    if not category_labels:
        raise AssignmentStateError("category_labels_missing")
    return one_hot_assignment(
        adata.obs[key].tolist(),
        observation_labels=adata.obs_names,
        category_labels=category_labels,
        source=f"obs[{key}]",
        level=str(key),
        lineage=[
            {
                "operation": "fallback",
                "reason": "broad_label_only",
                "container": "obs",
                "key": str(key),
            }
        ],
    )


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
            guidance_scope="all_panel",
        )

        # impute in panel
        in_panel_resolution = getattr(self.config, "rec_in_panel_subcluster_resolution", None)
        if in_panel_resolution is None:
            in_panel_resolution = self.config.rec_subcluster_resolution
        self.svc["sc_svc_impute_in_panel"] = self.local_impute(
            adata_sc_in_panel,
            f"leiden_{in_panel_resolution}",
            guidance_scope="in_panel",
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

    def _start_impute_guidance_event(
            self,
            *,
            problem_key,
            applicability="applicable",
    ):
        callback = getattr(
            self.config,
            "assignment_guidance_callback",
            None,
        )
        if callback is not None:
            callback(
                "start",
                problem_key=problem_key,
                route=str(
                    getattr(
                        self.config,
                        "assignment_guidance_route",
                        "sim2real:gene_panel",
                    )
                ),
                operator="imputation_ot",
                phase="lr",
                mode=_guidance_mode(self.config),
                applicability=applicability,
                numerics={
                    "beta": float(
                        getattr(
                            self.config,
                            "posterior_conditioning_beta",
                            1.0,
                        )
                    ),
                    "min_affinity": float(
                        getattr(
                            self.config,
                            "posterior_conditioning_min_affinity",
                            0.05,
                        )
                    ),
                    "operator_strength": float(
                        getattr(
                            self.config,
                            "posterior_conditioning_cost_strength",
                            0.2,
                        )
                    ),
                },
                solver=str(self.config.rec_ot_method),
            )
        return callback

    def _record_impute_not_applicable(
            self,
            *,
            problem_key,
            reason,
    ):
        callback = self._start_impute_guidance_event(
            problem_key=problem_key,
            applicability="not_applicable",
        )
        if callback is not None:
            callback(
                "terminal",
                problem_key=problem_key,
                outcome="not_applicable",
                reason=reason,
            )

    def _prepare_impute_assignment_guidance(
            self,
            *,
            ct_adata_sp,
            ct_adata_sc,
            sc_subcluster,
            subcluster_order,
            ot_cost,
            spot_prior,
            type_prior,
            problem_key,
    ):
        callback = self._start_impute_guidance_event(
            problem_key=problem_key,
        )
        mode = _guidance_mode(self.config)
        if mode == "off":
            if callback is not None:
                callback(
                    "terminal",
                    problem_key=problem_key,
                    outcome="off",
                    reason="guidance_off",
                )
            return ot_cost, None, False

        loaded = {}

        def load_state():
            key = self.config.cell_type_col
            categories = _assignment_categories(
                self.st_adata,
                self.sc_ref_adata,
                key=key,
            )
            left = _assignment_state_from_adata(
                ct_adata_sp,
                key=key,
                category_labels=categories,
            )
            right_cells = _assignment_state_from_adata(
                ct_adata_sc,
                key=key,
                category_labels=categories,
            )
            if left is None or right_cells is None:
                raise KeyError(key)
            right = aggregate_assignment(
                right_cells,
                ct_adata_sc.obs[sc_subcluster].tolist(),
                source=f"aggregate({right_cells.source})",
                level=str(sc_subcluster),
            )
            right = align_assignment_observations(
                right,
                subcluster_order,
            )
            loaded["left"] = left
            loaded["right"] = align_assignment_categories(
                right,
                left.category_labels,
            )
            return left

        resolution = resolve_assignment_guidance(mode, load_state)
        if resolution.availability != "available":
            if callback is not None:
                callback(
                    "terminal",
                    problem_key=problem_key,
                    outcome=resolution.outcome,
                    availability=resolution.availability,
                    reason=resolution.reason,
                )
            if resolution.outcome == "failed":
                raise ValueError(
                    f"assignment guidance unavailable: {resolution.reason}"
                )
            return ot_cost, None, False

        affinity = assignment_compatibility(
            loaded["left"],
            loaded["right"],
            beta=getattr(self.config, "posterior_conditioning_beta", 1.0),
            min_affinity=getattr(
                self.config,
                "posterior_conditioning_min_affinity",
                0.05,
            ),
        )
        compatibility_mode = posterior_conditioning_mode(self.config)
        reference_measure = None
        if compatibility_mode == "cost":
            ot_cost = ot_cost_guidance(
                ot_cost,
                affinity,
                getattr(
                    self.config,
                    "posterior_conditioning_cost_strength",
                    0.2,
                ),
            )
        elif compatibility_mode == "reference":
            reference_measure = reference_measure_from_marginals(
                spot_prior,
                type_prior,
                affinity,
            )
        if callback is not None:
            callback(
                "attempt",
                problem_key=problem_key,
                availability="available",
                left_assignment=loaded["left"],
                right_assignment=loaded["right"],
            )
        return ot_cost, reference_measure, True

    def _record_impute_guidance_terminal(
            self,
            *,
            problem_key,
            attempted,
            outcome,
            reason=None,
    ):
        callback = getattr(
            self.config,
            "assignment_guidance_callback",
            None,
        )
        if not attempted or callback is None:
            return
        fields = {"outcome": outcome}
        if reason is not None:
            fields["reason"] = reason
        callback("terminal", problem_key=problem_key, **fields)

    def local_impute(
            self,
            adata_sc,
            sc_subcluster,
            *,
            guidance_scope="impute",
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
        for candidate_ordinal, select_ct in enumerate(
            tqdm(cts, "Imputation by cell type"),
            start=1,
        ):
            self.logger.info(f"Conducting cell type: {select_ct} ........")
            label_text = str(select_ct)
            encoded_label = json.dumps(
                label_text,
                ensure_ascii=False,
                separators=(",", ":"),
            )
            problem_key = (
                f"sc_svc_impute:{guidance_scope}:"
                f"{candidate_ordinal}:"
                f"{len(label_text.encode('utf-8'))}:{encoded_label}"
            )
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
                self._record_impute_not_applicable(
                    problem_key=problem_key,
                    reason="empty_reference_support",
                )
                continue
            if ct_adata_sp.n_obs == 0:
                self._record_impute_not_applicable(
                    problem_key=problem_key,
                    reason="empty_spatial_support",
                )
                continue
            overlap_genes = ct_adata_sc.var_names.intersection(ct_adata_sp.var_names)
            if len(overlap_genes) == 0:
                self._record_impute_not_applicable(
                    problem_key=problem_key,
                    reason="overlap_genes_empty",
                )
                continue
            ct_adata_sc_overlap = ct_adata_sc[:, overlap_genes].copy()
            ct_adata_sp_overlap = ct_adata_sp[:, overlap_genes].copy()
            dums = pd.get_dummies(ct_adata_sc_overlap.obs[sc_subcluster],
                                  dtype=ct_adata_sc_overlap.X.dtype)
            ncats = dums.sum(axis=0)
            if dums.shape[1] == 0 or np.any(
                np.asarray(ncats, dtype=np.float64) <= 0
            ):
                self._record_impute_not_applicable(
                    problem_key=problem_key,
                    reason="zero_target_mass",
                )
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
                self._record_impute_not_applicable(
                    problem_key=problem_key,
                    reason="invalid_source_mass",
                )
                continue
            if not np.all(np.isfinite(target_mass)):
                self._record_impute_not_applicable(
                    problem_key=problem_key,
                    reason="invalid_target_mass",
                )
                continue
            if np.any(source_mass <= 0):
                self._record_impute_not_applicable(
                    problem_key=problem_key,
                    reason="zero_source_mass",
                )
                continue
            if np.any(target_mass <= 0):
                self._record_impute_not_applicable(
                    problem_key=problem_key,
                    reason="zero_target_mass",
                )
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
            (
                ot_cost,
                reference_measure,
                guidance_attempted,
            ) = self._prepare_impute_assignment_guidance(
                ct_adata_sp=ct_adata_sp_overlap,
                ct_adata_sc=ct_adata_sc_overlap,
                sc_subcluster=sc_subcluster,
                subcluster_order=type_prior.index,
                ot_cost=ot_cost,
                spot_prior=spot_prior.values,
                type_prior=type_prior.values,
                problem_key=problem_key,
            )

            dist_max = float(np.nanmax(ot_cost)) if ot_cost.size else 1.0
            if not np.isfinite(dist_max) or dist_max <= 0:
                dist_max = 1.0
            try:
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
                    reference_measure=reference_measure,
                    event_callback=getattr(
                        self.config,
                        "ot_event_callback",
                        None,
                    ),
                )
            except Exception:
                self._record_impute_guidance_terminal(
                    problem_key=problem_key,
                    attempted=guidance_attempted,
                    outcome="failed",
                    reason="solver_failed",
                )
                raise

            try:
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
            except KeyboardInterrupt:
                self._record_impute_guidance_terminal(
                    problem_key=problem_key,
                    attempted=guidance_attempted,
                    outcome="interrupted",
                    reason="result_assembly_interrupted",
                )
                raise
            except Exception:
                self._record_impute_guidance_terminal(
                    problem_key=problem_key,
                    attempted=guidance_attempted,
                    outcome="failed",
                    reason="result_assembly_failed",
                )
                raise
            adata_sp_cts.append(adata_sp_impute)
            self._record_impute_guidance_terminal(
                problem_key=problem_key,
                attempted=guidance_attempted,
                outcome="applied",
            )

        if not adata_sp_cts:
            return adata_sp[:0].copy()
        adata_sp_impute = sc.concat(adata_sp_cts)
        if self.config.rec_impute_prune_flag:
            adata_sp_impute = get_prune_adata(adata_sp_impute)

        return adata_sp_impute
