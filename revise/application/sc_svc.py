import scanpy as sc

from revise.application.application_svc import ApplicationSVC
from revise.methods.graph_cluster import GraphCluster
from revise.tools.bio import get_degs
from revise.tools.bio import conclusions_write
from revise.tools.bio import plot_volcano


"""
The sc-SVC application is organized into two classes:
1. ScSVC: reconstructs sc-SVC. `local_refinement` returns a cell-type-specific sc-SVC,
   including `sc_SVC_adata_spatial` and `sc_SVC_adata_expr` for that cell type (h5ad format).
2. ScSVCAnalysis: performs downstream bioinformatics analyses based on
   `sc_SVC_adata_spatial` and `sc_SVC_adata_expr`.
"""


class ScSVC(ApplicationSVC):
    """
    sc-SVC class for application usage.
    
    This class handles single-cell resolution spatial transcriptomics data,
    filtering cells and genes based on transcript counts and preparing
    data for downstream annotation and reconstruction.
    """
    def __init__(self, st_adata, sc_ref_adata, config, logger):
        super().__init__(st_adata, sc_ref_adata, config, None, logger)
        self._adata_validate()
        self.sc_ref_adata_raw = self.sc_ref_adata.copy()
        self.graph_cluster = GraphCluster(self.config, self.logger)
        self.cluster_col = "SVC_cluster"

    def local_refinement(self, select_ct, sub_cell_type_col, resolutions, select_res=None):

        ct_adata_sp = self.st_adata[self.st_adata.obs['Level1'] == select_ct]
        ct_adata_sc = self.sc_ref_adata[self.sc_ref_adata.obs['Level1'] == select_ct]
        annotate_kwargs = dict(self.config.__dict__)
        annotate_kwargs["cell_type_col"] = sub_cell_type_col
        ct_adata_sp = self.annotate_method.run(ct_adata_sp, ct_adata_sc, **annotate_kwargs)
        sc_SVC_adata, merge_df, best_res = self.graph_cluster.run(ct_adata_sp, resolutions, sub_cell_type_col)
        if select_res is None:
            self.logger.info(f"User does not input select_res, use best_res {best_res} based on spatial alignment score")
            select_res = best_res
        else:
            self.logger.info(f"Use resolution {select_res} from user input")

        sc_SVC_adata.obs[self.cluster_col] = sc_SVC_adata.obs[f'leiden_{select_res}'].astype('category')
        sp_cluster_num = merge_df.loc[merge_df['resolution'] == best_res, 'cluster_num'].values[0]
        self.logger.info(f"resolution {select_res} got cluster number {sp_cluster_num}")

        annotate_kwargs = dict(self.config.__dict__)
        annotate_kwargs["cell_type_col"] = self.cluster_col
        ct_adata_sc = self.annotate_method.run(ct_adata_sc, sc_SVC_adata, **annotate_kwargs)
        return sc_SVC_adata, ct_adata_sc


class ScSVCAnalysis:
    def __init__(self, sc_SVC_adata_spatial, sc_SVC_adata_expr, cluster_col):
        self.sc_SVC_adata_spatial = sc_SVC_adata_spatial
        self.sc_SVC_adata_expr = sc_SVC_adata_expr
        self.cluster_col = cluster_col
        self.sc_SVC_degs = get_degs(self.sc_SVC_adata_expr, groupby=self.cluster_col, method='t-test', fc_threshold=None)
        # self.sc_SVC_degs.to_csv(f"{self.config.result_root_path}/degs_all.csv")

    def get_cm_df(self, sub_cell_type_col):
        grouped = self.sc_SVC_adata_expr.obs.groupby([self.cluster_col, sub_cell_type_col]).size()
        cm_df = grouped.unstack(fill_value=0)
        return cm_df

    def get_svc_degs(self, target_cluster=None, fc_threshold=None):
        if not target_cluster:
            return get_degs(self.sc_SVC_adata_expr, groupby=self.cluster_col, method='t-test', fc_threshold=fc_threshold)
        else:
            assert isinstance(target_cluster, (list, tuple, set)), f"target_cluster {target_cluster} should be a list, tuple, or set"
            assert len(target_cluster) > 0, "target_cluster should not be empty"
            select_adata = self.sc_SVC_adata_expr[self.sc_SVC_adata_expr.obs[self.cluster_col].isin(target_cluster)]
            return get_degs(select_adata, groupby=self.cluster_col, method='t-test', fc_threshold=fc_threshold)

    def get_dot_plot(self, cluster_nums, marker_dict, normalize=False):
        if cluster_nums is None:
            print("Using all clusters")
            cluster_nums = self.sc_SVC_adata_expr.obs['SVC_cluster'].cat.categories.tolist()
        assert isinstance(cluster_nums, (list, tuple, set)), f"target_cluster {cluster_nums} should be a list, tuple, or set"
        assert len(cluster_nums) > 0, "target_cluster should not be empty"
        select_adata = self.sc_SVC_adata_expr[self.sc_SVC_adata_expr.obs[self.cluster_col].isin(cluster_nums)]
        select_adata = select_adata.copy()
        if normalize:
            sc.pp.normalize_total(select_adata, target_sum=1e4)
            sc.pp.log1p(select_adata)
        sc.pl.dotplot(select_adata, marker_dict, groupby=self.cluster_col, dendrogram=False)

    def get_violin_plot(self, cluster_nums, normalize=False):
        if cluster_nums is None:
            print("Using all clusters")
            cluster_nums = self.sc_SVC_adata_expr.obs['SVC_cluster'].cat.categories.tolist()
        assert isinstance(cluster_nums, (list, tuple, set)), f"target_cluster {cluster_nums} should be a list, tuple, or set"
        assert len(cluster_nums) > 0, "target_cluster should not be empty"
        select_adata = self.sc_SVC_adata_expr[self.sc_SVC_adata_expr.obs[self.cluster_col].isin(cluster_nums)]
        if normalize:
            sc.pp.normalize_total(select_adata, target_sum=1e4)
            sc.pp.log1p(select_adata)
        sc.pl.violin(select_adata, "LPL", groupby=self.cluster_col)
    
    def get_pathway_conclusion(
        self,
        cluster_nums,
        fc_threshold,
        pathway_num,
        gene_num,
        geneset_file=["MSigDB_Hallmark_2020", "KEGG_2021_Human", "GO_Biological_Process_2025"],
        normalize=False
    ):
        if cluster_nums is None:
            print("Using all clusters")
            cluster_nums = self.sc_SVC_adata_expr.obs['SVC_cluster'].cat.categories.tolist()
        if len(geneset_file) == 0:
            raise NotImplementedError(f"Got empty geneset_file!")
        select_adata = self.sc_SVC_adata_expr[self.sc_SVC_adata_expr.obs['SVC_cluster'].isin(cluster_nums)]
        if normalize:
            sc.pp.normalize_total(select_adata, target_sum=1e4)
            sc.pp.log1p(select_adata)
        deg_df = self.get_svc_degs(cluster_nums, fc_threshold)
        all_pathway = conclusions_write(
            deg_df, geneset_file, gene_num=gene_num, pathway_num=pathway_num,
            print_flag=True, conclusion_file_name=None
        )
        return all_pathway

    def get_volcano_plot(self, cluster_nums, target_group, replace_cols=None, fc_threshold=None, log_fold_changes=10, logfc_threshold=1, padj_threshold=1e-6, top_k=10):
        indices = self.sc_SVC_adata_expr.obs['SVC_cluster'].isin(cluster_nums)
        select_adata = self.sc_SVC_adata_expr[indices, :]
        if replace_cols is not None:
            select_adata.obs['SVC_cluster'].replace(replace_cols, inplace=True)

        select_deg_df = get_degs(select_adata, groupby='SVC_cluster', method='t-test', fc_threshold=fc_threshold)
        select_deg_df = select_deg_df[select_deg_df['group'] == target_group]
        select_deg_df.reset_index(drop = True, inplace = True)
        select_deg_df = select_deg_df[select_deg_df['logfoldchanges'].abs() <= log_fold_changes]

        plot_volcano(select_deg_df, logfc_threshold=logfc_threshold, padj_threshold=padj_threshold, top_k=top_k, save_file_name=None)
        return select_deg_df
