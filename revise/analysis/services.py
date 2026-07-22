from __future__ import annotations

from typing import Dict

import scanpy as sc

from revise.analysis.bio import conclusions_write
from revise.analysis.bio import get_degs
from revise.analysis.bio import plot_volcano
from revise.analysis.metrics import compute_clustering_metrics
from revise.svc import SVC


class ScSVCAnalysisService:
    """Downstream iST-SVC analysis that consumes unified SVC only."""

    def __init__(self, svc: SVC, cluster_col: str = "SVC_cluster") -> None:
        if svc.expr is None or svc.spatial is None:
            raise ValueError("ScSVCAnalysisService requires both expr and spatial AnnData")
        if cluster_col not in svc.expr.obs:
            raise ValueError(f"cluster_col '{cluster_col}' not found in svc.expr.obs")
        self.svc = svc
        self.cluster_col = cluster_col
        self.sc_SVC_adata_spatial = svc.spatial
        self.sc_SVC_adata_expr = svc.expr
        # Keep notebook contract: precomputed DEG table on init.
        self.sc_SVC_degs = get_degs(self.sc_SVC_adata_expr, groupby=self.cluster_col, method="t-test", fc_threshold=None)

    def get_cm_df(self, sub_cell_type_col: str):
        grouped = self.sc_SVC_adata_expr.obs.groupby([self.cluster_col, sub_cell_type_col]).size()
        return grouped.unstack(fill_value=0)

    def get_svc_degs(self, target_cluster=None, fc_threshold=None):
        if target_cluster is None:
            select_adata = self.sc_SVC_adata_expr
        else:
            if not isinstance(target_cluster, (list, tuple, set)):
                raise ValueError("target_cluster must be list/tuple/set when provided")
            if len(target_cluster) == 0:
                raise ValueError("target_cluster should not be empty")
            select_adata = self.sc_SVC_adata_expr[self.sc_SVC_adata_expr.obs[self.cluster_col].isin(target_cluster)]
        return get_degs(select_adata, groupby=self.cluster_col, method="t-test", fc_threshold=fc_threshold)

    def get_dot_plot(self, cluster_nums, marker_dict, normalize: bool = False):
        if cluster_nums is None:
            series = self.sc_SVC_adata_expr.obs[self.cluster_col]
            if hasattr(series, "cat"):
                cluster_nums = series.cat.categories.tolist()
            else:
                cluster_nums = series.astype(str).unique().tolist()
        if not isinstance(cluster_nums, (list, tuple, set)):
            raise ValueError("cluster_nums should be list/tuple/set")
        if len(cluster_nums) == 0:
            raise ValueError("cluster_nums should not be empty")
        select_adata = self.sc_SVC_adata_expr[self.sc_SVC_adata_expr.obs[self.cluster_col].isin(cluster_nums)].copy()
        if normalize:
            sc.pp.normalize_total(select_adata, target_sum=1e4)
            sc.pp.log1p(select_adata)
        sc.pl.dotplot(select_adata, marker_dict, groupby=self.cluster_col, dendrogram=False)

    def get_pathway_conclusion(
        self,
        cluster_nums,
        fc_threshold,
        pathway_num,
        gene_num,
        geneset_file,
        normalize: bool = False,
    ):
        if cluster_nums is None:
            series = self.sc_SVC_adata_expr.obs[self.cluster_col]
            if hasattr(series, "cat"):
                cluster_nums = series.cat.categories.tolist()
            else:
                cluster_nums = series.astype(str).unique().tolist()
        if not geneset_file:
            raise ValueError("geneset_file cannot be empty")
        select_adata = self.sc_SVC_adata_expr[self.sc_SVC_adata_expr.obs[self.cluster_col].isin(cluster_nums)]
        if normalize:
            select_adata = select_adata.copy()
            sc.pp.normalize_total(select_adata, target_sum=1e4)
            sc.pp.log1p(select_adata)
        deg_df = get_degs(select_adata, groupby=self.cluster_col, method="t-test", fc_threshold=fc_threshold)
        return conclusions_write(
            deg_df,
            geneset_file,
            gene_num=gene_num,
            pathway_num=pathway_num,
            print_flag=True,
            conclusion_file_name=None,
        )

    def get_volcano_plot(
        self,
        cluster_nums,
        target_group,
        replace_cols=None,
        fc_threshold=None,
        log_fold_changes=10,
        logfc_threshold=1,
        padj_threshold=1e-6,
        top_k=10,
    ):
        indices = self.sc_SVC_adata_expr.obs[self.cluster_col].isin(cluster_nums)
        select_adata = self.sc_SVC_adata_expr[indices, :].copy()
        if replace_cols is not None:
            select_adata.obs[self.cluster_col].replace(replace_cols, inplace=True)

        select_deg_df = get_degs(select_adata, groupby=self.cluster_col, method="t-test", fc_threshold=fc_threshold)
        select_deg_df = select_deg_df[select_deg_df["group"] == target_group].copy()
        select_deg_df = select_deg_df[select_deg_df["logfoldchanges"].abs() <= log_fold_changes].copy()
        plot_volcano(
            select_deg_df,
            logfc_threshold=logfc_threshold,
            padj_threshold=padj_threshold,
            top_k=top_k,
            save_file_name=None,
        )
        return select_deg_df


class SpSVCAnalysisService:
    """Basic downstream metrics for hST-SVC outputs."""

    def __init__(self, svc: SVC) -> None:
        if svc.spatial is None:
            raise ValueError("SpSVCAnalysisService requires svc.spatial")
        self.svc = svc

    def clustering_metrics(self, pred_col: str, ref_col: str) -> Dict[str, float]:
        ari, nmi = compute_clustering_metrics(self.svc.spatial, pred_col, ref_col)
        return {"ari": float(ari), "nmi": float(nmi)}
