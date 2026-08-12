import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from tqdm import tqdm


def get_prune_adata(adata):
    """
    Prune zero-variance genes by converting them to integer type.
    
    Args:
        adata: AnnData object to prune
        
    Returns:
        AnnData: Copy of adata with zero-variance genes converted to integers
    """
    adata = adata.copy()
    adata_int = adata.copy()
    adata_int.X = adata_int.X.astype(int)

    zero_var_mask = np.asarray(adata_int.X.var(axis=0) == 0).ravel()
    adata.X[:, zero_var_mask] = adata_int.X[:, zero_var_mask]

    return adata


def trim_sp_adata(adata_sp, adata_sc, celltype_col, logfc_threshold=-3):
    """
    Trim low-expression genes from spatial data based on cell type-specific analysis.
    
    This function identifies and removes genes that are lowly expressed in spatial
    data compared to their expected expression based on single-cell reference.
    The process involves:
    1. Finding cell type with highest expression for each gene
    2. Performing differential expression analysis
    3. Identifying low-expression genes
    4. Setting their expression to zero in spatial data
    
    Args:
        adata_sp: Spatial transcriptomics AnnData object (will be modified)
        adata_sc: Single-cell reference AnnData object
        celltype_col: Column name in adata_sc.obs containing cell type labels
        logfc_threshold: Log fold change threshold for identifying low-expression genes
            
    Returns:
        tuple: (adata_sp_filtered, low_expr_genes_dict)
            - adata_sp_filtered: Filtered spatial data with low-expression genes set to zero
            - low_expr_genes_dict: Dictionary mapping cell types to lists of trimmed genes
    """
    # Step 1: Find cell type with highest expression for each gene
    celltype_genes = get_max_expr_celltypes(adata_sc, celltype_col)

    # Step 2: Perform differential expression analysis using each cell type as reference
    de_results = calculate_celltype_differential_expression(adata_sc, celltype_col, celltype_genes)

    # Step 3: Identify low-expression genes
    low_expr_genes_dict = identify_low_expr_genes(
        de_results, logfc_threshold = logfc_threshold,
        pval_threshold = 0.001, percentile_threshold = None
    )

    adata_sp_filtered = trim_spatial_genes(adata_sp, low_expr_genes_dict, celltype_col)
    return adata_sp_filtered, low_expr_genes_dict


def get_max_expr_celltypes(adata, groupby):
    """
    Identify cell type with highest average expression for each gene.
    
    This function groups genes by the cell type that expresses them most highly,
    which is useful for cell type-specific gene analysis.
    
    Args:
        adata: AnnData object with normalized and log-transformed expression
        groupby: Column name in adata.obs containing cell type labels
        
    Returns:
        dict: Dictionary mapping cell type names to lists of gene names
            that have highest expression in that cell type
    """
    from revise.tools.log import Logger
    logger = Logger().get_logger()
    logger.info("Calculating maximum expression cell types for each gene...")

    adata.obs[groupby] = adata.obs[groupby].astype('category')
    cell_types = adata.obs[groupby].cat.categories

    # Preprocessing
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)

    # Calculate average expression matrix for each cell type
    mean_expr_matrix = np.zeros((len(cell_types), adata_norm.n_vars))
    for i, cell_type in enumerate(cell_types):
        cell_type_mask = adata_norm.obs[groupby] == cell_type
        if sparse.issparse(adata_norm.X):
            mean_expr_matrix[i] = np.array(adata_norm[cell_type_mask, :].X.mean(axis=0)).flatten()
        else:
            mean_expr_matrix[i] = adata_norm[cell_type_mask, :].X.mean(axis=0)

    # Find cell type with highest expression for each gene
    max_expr_indices = np.argmax(mean_expr_matrix, axis=0)

    # Group genes by cell type with highest expression
    celltype_genes = {cell_type: [] for cell_type in cell_types}
    for gene_idx, celltype_idx in enumerate(max_expr_indices):
        gene = adata_norm.var_names[gene_idx]
        max_celltype = cell_types[celltype_idx]
        celltype_genes[max_celltype].append(gene)

    # Print number of specific genes for each cell type
    for cell_type, genes in celltype_genes.items():
        logger.info(f"{cell_type}: {len(genes)} genes with maximum expression")

    return celltype_genes


def calculate_celltype_differential_expression(adata, groupby, celltype_genes):
    """
    Perform differential expression analysis for cell type-specific genes.
    
    For each cell type, this function uses its specific genes (genes with
    highest expression in that type) to perform differential expression
    analysis against all other cell types.
    
    Args:
        adata: Single-cell AnnData object with normalized expression
        groupby: Column name in adata.obs containing cell type labels
        celltype_genes: Dictionary mapping cell types to lists of specific genes
        
    Returns:
        dict: Nested dictionary structure:
            {ref_celltype: {test_celltype: DataFrame with DE results}}
            Each DataFrame contains columns: 'gene', 'logfoldchanges', 'pvals_adj', 'ref_celltype'
    """
    from revise.tools.log import Logger
    logger = Logger().get_logger()
    logger.info("Performing differential expression analysis for each cell type's specific genes...")

    adata.obs[groupby] = adata.obs[groupby].astype('category')
    cell_types = adata.obs[groupby].cat.categories

    # Preprocessing
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)

    de_results = {}

    for ref_celltype, specific_genes in tqdm(celltype_genes.items(), desc="Cell types"):
        if not specific_genes:
            logger.info(f"Skipping {ref_celltype}: no specific genes")
            continue

        logger.info(f"\nAnalyzing {ref_celltype} with {len(specific_genes)} specific genes")

        # Create adata subset containing only this cell type's specific genes
        specific_gene_mask = adata_norm.var_names.isin(specific_genes)
        adata_subset = adata_norm[:, specific_gene_mask].copy()

        # Ensure reference cell type has enough cells in the data
        if (adata_subset.obs[groupby] == ref_celltype).sum() < 3:
            logger.info(f"Skipping {ref_celltype}: not enough cells")
            continue

        # Use current cell type as reference for differential expression analysis
        try:
            sc.tl.rank_genes_groups(
                adata_subset,
                groupby=groupby,
                reference=ref_celltype,  # Should be a string, not a dict
                use_raw=False
            )
        except Exception as e:
            logger.warn(f"Error analyzing {ref_celltype}: {e}")
            continue

        # Extract results
        celltype_results = {}
        for test_celltype in cell_types:
            if test_celltype == ref_celltype:
                continue
            try:
                gene_names = adata_subset.uns['rank_genes_groups']['names'][test_celltype]
                logfoldchanges = adata_subset.uns['rank_genes_groups']['logfoldchanges'][test_celltype]
                pvals_adj = adata_subset.uns['rank_genes_groups']['pvals_adj'][test_celltype]

                results_df = pd.DataFrame({
                    'gene': gene_names,
                    'logfoldchanges': logfoldchanges,
                    'pvals_adj': pvals_adj,
                    'ref_celltype': ref_celltype
                })

                celltype_results[test_celltype] = results_df
            except (KeyError, ValueError, IndexError) as e:
                logger.warn(f"Error extracting results for {test_celltype} vs {ref_celltype}: {e}")
                continue

        de_results[ref_celltype] = celltype_results

    return de_results


def identify_low_expr_genes(de_results, logfc_threshold=0, pval_threshold=0.05, percentile_threshold=None):
    """
    Identify low-expression genes based on differential expression analysis results
    """
    from revise.tools.log import Logger
    logger = Logger().get_logger()

    low_expr_genes_dict = {}

    # Initialize low-expression gene list for each cell type
    all_celltypes = set()
    for ref_results in de_results.values():
        all_celltypes.update(ref_results.keys())

    for cell_type in all_celltypes:
        low_expr_genes_dict[cell_type] = []

    # Iterate through all reference cell type results
    for ref_celltype, celltype_results in de_results.items():
        for test_celltype, de_df in celltype_results.items():
            # Filter low-expression genes
            low_expr_mask = (de_df['logfoldchanges'] < logfc_threshold) & (de_df['pvals_adj'] < pval_threshold)
            low_expr_genes = de_df[low_expr_mask]['gene'].tolist()

            low_expr_genes_dict[test_celltype].extend(low_expr_genes)

    # Remove duplicates
    for cell_type in low_expr_genes_dict:
        low_expr_genes_dict[cell_type] = list(set(low_expr_genes_dict[cell_type]))

    for cell_type, genes in low_expr_genes_dict.items():
        logger.info(f"{cell_type}: {len(genes)} low-expression genes.")

    if percentile_threshold is not None:
        low_expr_genes_dict = apply_percentile_threshold(
            low_expr_genes_dict, de_results, percentile_threshold
        )

    return low_expr_genes_dict


def apply_percentile_threshold(low_expr_genes_dict, de_results, percentile_threshold=10):
    """
    Apply percentile threshold to further filter low-expression genes.
    
    This function refines the low-expression gene list by applying a
    percentile threshold on log fold changes, keeping only genes below
    the threshold.
    
    Args:
        low_expr_genes_dict: Dictionary mapping cell types to lists of low-expression genes
        de_results: Differential expression results from calculate_celltype_differential_expression
        percentile_threshold: Percentile threshold (0-100) for log fold change filtering
            
    Returns:
        dict: Refined dictionary with filtered gene lists per cell type
    """
    from revise.tools.log import Logger
    logger = Logger().get_logger()
    logger.info(f"Applying percentile threshold: {percentile_threshold}%")

    refined_dict = {}

    for cell_type, genes in low_expr_genes_dict.items():
        # Find reference corresponding to this cell type (need to infer from de_results)
        ref_celltype = None
        for ref, results in de_results.items():
            if cell_type in results:
                ref_celltype = ref
                break

        if ref_celltype is None:
            refined_dict[cell_type] = genes
            continue

        # Get logfc for all genes of this cell type
        de_df = de_results[ref_celltype][cell_type]

        # Calculate percentile threshold for logfc
        logfc_threshold = np.percentile(de_df['logfoldchanges'], percentile_threshold)

        # Filter genes: logfc below percentile threshold and pval significant
        filtered_genes = []
        for gene in genes:
            gene_data = de_df[de_df['gene'] == gene]
            if len(gene_data) > 0 and gene_data['logfoldchanges'].values[0] < logfc_threshold:
                filtered_genes.append(gene)

        refined_dict[cell_type] = filtered_genes
        logger.info(f"{cell_type}: {len(filtered_genes)} genes after percentile filtering")

    return refined_dict


def trim_spatial_genes(adata_sp, low_expr_genes_dict, celltype_col, inplace=True, drop_stored_zeros=True):
    """
    Set gene expression to 0 for genes specified in low_expr_genes_dict for each cell type, and update adata_sp.X

    Parameters
    ----------
    adata_sp : AnnData
    low_expr_genes_dict : dict[str, list[str]]
        {cell_type: [gene1, gene2, ...]}
    celltype_col : str
        Cell type column name in adata_sp.obs
    inplace : bool
        True to modify in place; False to return a copy
    drop_stored_zeros : bool
        Whether to clean stored zeros (eliminate_zeros) after setting sparse matrix to zero

    Returns
    -------
    AnnData
    """
    from revise.tools.log import Logger
    logger = Logger().get_logger()
    if not inplace:
        adata_sp = adata_sp.copy()

    if celltype_col not in adata_sp.obs:
        raise KeyError(f"Column not found in obs: {celltype_col}")

    # Warn about non-unique var_names (affects gene name -> index mapping)
    if hasattr(adata_sp.var_names, "is_unique") and not adata_sp.var_names.is_unique:
        logger.warn("[WARN] var_names not unique, mapping will use the index of the last occurrence. Recommend calling adata_sp.var_names_make_unique() first.")

    # Gene name to index mapping
    gene_to_idx = {g: i for i, g in enumerate(adata_sp.var_names)}
    is_sparse = sparse.issparse(adata_sp.X)

    if is_sparse:
        orig_fmt = adata_sp.X.getformat()
        X = adata_sp.X.tolil(copy=True)  # Row assignment friendly
    else:
        # Ensure it's a writable ndarray
        X = np.asarray(adata_sp.X)
        orig_fmt = None
        if not X.flags.writeable:
            X = X.copy()

    total_cells_touched = 0
    total_gene_hits = 0

    # Process each cell type
    ct_values = adata_sp.obs[celltype_col].values
    for ct, gene_list in low_expr_genes_dict.items():
        # Row indices
        rows = np.where(ct_values == ct)[0]
        if rows.size == 0:
            # No cells of this cell type
            continue

        # Column indices (filter out non-existent gene names)
        cols = np.array([gene_to_idx[g] for g in gene_list if g in gene_to_idx], dtype=int)
        if cols.size == 0:
            # None of this cell type's genes found in var_names
            continue

        if is_sparse:
            # For LIL, set zero row by row for safety
            for r in rows:
                X[r, cols] = 0
        else:
            # Dense: block assignment
            X[np.ix_(rows, cols)] = 0

        total_cells_touched += rows.size
        total_gene_hits += cols.size

    # Write back & cleanup
    if is_sparse:
        # Convert back to original format
        if orig_fmt == "csr":
            X = X.tocsr()
        elif orig_fmt == "csc":
            X = X.tocsc()
        elif orig_fmt == "coo":
            X = X.tocoo().tocsr()  # Back to CSR is more general
        else:
            X = X.tocsr()

        if drop_stored_zeros:
            X.eliminate_zeros()

    adata_sp.X = X

    logger.info(f"[trim_spatial_genes] touched cells: {total_cells_touched}, genes per-ct (avg not shown). is_sparse={is_sparse}, cleaned_zeros={is_sparse and drop_stored_zeros}")
    return adata_sp
