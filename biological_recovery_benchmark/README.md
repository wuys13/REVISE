# Biological recovery benchmark

This directory contains a standalone, single-`AnnData` evaluation layer for
the Figure 3 real-data biological-recovery analysis. It does not import or
modify `revise`.

## Metrics

- Cell-state agreement: ARI, AMI, NMI and ASW.
- Marker specificity: cell-type-level TMP, MER and log2MER.
- Conditional spatial coherence: per-gene MISC and MIDC.
- Ordinary spatial coherence: per-gene global Moran's I and independently
  graphed per-cell-type Moran's I.

TMP/MER uses `adata.X` exactly as supplied. Spatial metrics operate on a copy,
apply total-count normalization to 10,000 counts followed by `log1p`, and do
not mutate the input. The caller supplies clustering labels; this module does
not run Leiden or select a resolution.

## Example

```python
from biological_recovery_benchmark import evaluate_adata, save_evaluation_results

results = evaluate_adata(
    adata,
    cell_type_col="Level1",
    pred_label_col="leiden",
    embedding_key="X_pca",
    marker_map=marker_map,
    min_cell_type_size=10,
)
save_evaluation_results(
    results,
    "output/P1CRC/raw",
    metadata={
        "dataset": "P1CRC",
        "platform": "Visium HD",
        "method": "Raw",
        "input_kind": "raw",
        "input_source": "P1CRC_HD.h5ad",
    },
)
```

Global and cell-type-specific Moran results are intentionally saved separately.
Cell types below `min_cell_type_size` are omitted only from the latter and are
listed in `run_metadata.json`. No hashes or Git identifiers are recorded.
