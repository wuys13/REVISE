import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

from revise.analysis.advanced import enrichment
from revise.analysis.basic import differential_expression

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.useafm"] = False


def _require_gseapy():
    try:
        import gseapy
    except ImportError as exc:
        raise ImportError(
            "gseapy is required for pathway enrichment features; "
            "install it with `python -m pip install \"revise-svc[pathway]\"`."
        ) from exc
    return gseapy


def _require_networkx():
    try:
        import networkx
    except ImportError as exc:
        raise ImportError(
            "networkx is required for pathway network plots; install it with "
            "`python -m pip install \"revise-svc[pathway]\"`."
        ) from exc
    return networkx


def get_degs(adata, groupby, method="t-test", fc_threshold=None, reference="rest"):
    """Run differential expression and return a tidy DataFrame of DEG statistics."""
    return differential_expression.get_degs(
        adata,
        groupby=groupby,
        method=method,
        fc_threshold=fc_threshold,
        reference=reference,
    )


# get enrichment pathway
def get_enrichment(deg_genes, geneset_file, cutoff=0.05):
    """Run Enrichr on a gene list and return the results DataFrame."""
    try:
        return enrichment.get_enrichment(
            deg_genes,
            geneset_file,
            cutoff=cutoff,
        )
    except enrichment.EnrichmentExecutionError as exc:
        cause = exc.__cause__ or exc
        print(
            "Skipping enrichment analysis: "
            f"{type(cause).__name__}: {cause}"
        )
        return pd.DataFrame(columns=enrichment.ENRICHMENT_RESULT_COLUMNS)


# get conclusions used for input to LLM
def conclusions_write(
    deg_df,
    geneset_file,
    gene_num=60,
    pathway_num=20,
    print_flag=False,
    conclusion_file_name=None,
):
    """Summarize top genes/pathways per group and optionally write to a file."""
    groups = deg_df["group"].unique()
    all_pathway = pd.DataFrame()
    for g in tqdm(groups, desc="Enrichment analysis"):
        deg_df_group = deg_df.loc[deg_df["group"] == g]
        deg_df_group.reset_index(drop=True, inplace=True)
        deg_genes = deg_df_group["gene"][:gene_num].tolist()
        pathway = get_enrichment(deg_genes, geneset_file)
        pathway = pathway.copy()
        pathway["group"] = g
        if "Adjusted P-value" in pathway.columns:
            pathway.sort_values(by="Adjusted P-value", ascending=True, inplace=True)
        all_pathway = pd.concat([all_pathway, pathway], axis=0)

        if print_flag:
            print(f"Cluster {g}: ")
            print(f"Up-regulated genes: {deg_genes}")
            print(f"Enriched pathways: {pathway['Term'].values[:pathway_num]}\n")

        if conclusion_file_name:
            with open(conclusion_file_name, "a") as f:
                f.write(f"Cluster {g}:\n")
                f.write(f"Up-regulated genes: {deg_genes}\n")
                f.write(
                    f"Enriched pathways: {pathway['Term'].values[:pathway_num]}\n\n"
                )
    return all_pathway


def plot_volcano(
    select_deg_df,
    logfc_threshold=0.5,
    padj_threshold=0.05,
    top_k=10,
    figsize=(10, 8),
    title=None,
    save_file_name=None,
):
    """Plot a volcano plot for a DEG DataFrame."""
    df = select_deg_df.copy()

    df["-log10_padj"] = -np.log10(df["pvals_adj"] + 1e-300)

    colors = []
    sig_up = []
    sig_down = []

    for i, row in df.iterrows():
        if row["pvals_adj"] < padj_threshold and row["logfoldchanges"] > logfc_threshold:
            colors.append("#E64B35")
            sig_up.append(i)
        elif (
            row["pvals_adj"] < padj_threshold
            and row["logfoldchanges"] < -logfc_threshold
        ):
            colors.append("#4DBBD5")
            sig_down.append(i)
        else:
            colors.append("#BDBDBD")

    df["color"] = colors

    plt.figure(figsize=figsize)

    plt.scatter(
        df["logfoldchanges"],
        df["-log10_padj"],
        c=df["color"],
        alpha=1,
        s=40,
        edgecolors="none",
    )

    plt.axhline(y=-np.log10(padj_threshold), color="grey", linestyle="--", alpha=0.8)
    plt.axvline(x=logfc_threshold, color="grey", linestyle="--", alpha=0.8)
    plt.axvline(x=-logfc_threshold, color="grey", linestyle="--", alpha=0.8)

    if len(sig_up) > 0:
        top_up = df.loc[sig_up].nlargest(top_k, "logfoldchanges")
        for i, row in top_up.iterrows():
            plt.annotate(
                row["gene"],
                (row["logfoldchanges"], row["-log10_padj"]),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=8,
                alpha=0.8,
            )

    if len(sig_down) > 0:
        top_down = df.loc[sig_down].nsmallest(top_k, "logfoldchanges")
        for i, row in top_down.iterrows():
            plt.annotate(
                row["gene"],
                (row["logfoldchanges"], row["-log10_padj"]),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=8,
                alpha=0.8,
            )

    plt.xlabel("log2(Fold Change)")
    plt.ylabel("-log10(Adjusted p-value)")
    if title:
        plt.title(title)
    else:
        plt.title("Volcano Plot")

    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    if save_file_name is not None:
        plt.savefig(save_file_name, dpi=300, bbox_inches="tight")
    else:
        plt.show()


def pathway_barplot(pathway, pathway_num=5, save_file_name=None):
    """Plot a horizontal bar chart of top enriched pathways."""
    if pathway is None or pathway.empty or "Adjusted P-value" not in pathway.columns:
        print("Skipping pathway barplot: no enriched terms.")
        return
    pathway = pathway.head(pathway_num).copy()

    pathway["-log10(Adjusted P-value)"] = -np.log10(pathway["Adjusted P-value"])

    pathway = pathway.iloc[::-1].reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.barh(
        range(len(pathway)),
        pathway["-log10(Adjusted P-value)"],
        color="skyblue",
        edgecolor="none",
    )

    ax.set_yticks(range(len(pathway)))
    ax.set_yticklabels(pathway["Term"])

    ax.set_xlabel("-log10(Adjusted P-value)")
    ax.set_title("Pathway Enrichment Analysis")

    ax.grid(axis="x", linestyle="--", alpha=0.7)

    plt.tight_layout()

    if save_file_name is not None:
        plt.savefig(save_file_name, dpi=300, bbox_inches="tight")
    else:
        plt.show()


def pathway_network_plot(pathway, top_term=6, save_file_name=None):
    """Plot a pathway similarity network using gseapy enrichment_map."""
    if pathway is None or pathway.empty or "Term" not in pathway.columns:
        print("Skipping pathway network plot: no enriched terms.")
        return
    gp = _require_gseapy()
    nx = _require_networkx()

    nodes, edges = gp.enrichment_map(pathway, top_term=top_term)
    nodes["gene_num"] = (
        nodes["Genes"]
        .fillna("")
        .str.strip()
        .pipe(lambda s: s.where(s != "", np.nan))
        .str.split(r"\s*,\s*")
        .str.len()
        .fillna(0)
        .astype(int)
    )

    G = nx.from_pandas_edgelist(
        edges,
        source="src_idx",
        target="targ_idx",
        edge_attr=["jaccard_coef", "overlap_coef", "overlap_genes"],
    )

    fig, ax = plt.subplots(figsize=(8, 8))

    # init node cooridnates
    pos = nx.layout.spiral_layout(G)

    node_list = list(G.nodes())
    node_sizes = []
    node_colors = []
    node_labels = {}

    for node in node_list:
        if node < len(nodes):
            node_sizes.append(nodes.iloc[node]["gene_num"] * 1000)
            node_colors.append(nodes.iloc[node]["p_inv"])
            node_labels[node] = nodes.iloc[node]["Term"]
        else:
            node_sizes.append(1000)
            node_colors.append(0.0)
            node_labels[node] = f"Term_{node}"

    # draw node
    nx.draw_networkx_nodes(
        G,
        pos=pos,
        cmap=plt.cm.Oranges,
        node_color=node_colors,
        node_size=node_sizes,
    )

    nx.draw_networkx_labels(G, pos=pos, labels=node_labels, font_size=8)
    # draw edge
    edge_weights = list(nx.get_edge_attributes(G, "jaccard_coef").values())
    if edge_weights:
        edge_widths = [x * 30 for x in edge_weights]
        nx.draw_networkx_edges(
            G, pos=pos, width=edge_widths, edge_color="#CDDBD4"
        )
    if save_file_name is not None:
        plt.savefig(save_file_name, dpi=300, bbox_inches="tight")
    else:
        plt.show()


__all__ = [
    "conclusions_write",
    "get_degs",
    "get_enrichment",
    "pathway_barplot",
    "pathway_network_plot",
    "plot_volcano",
]
