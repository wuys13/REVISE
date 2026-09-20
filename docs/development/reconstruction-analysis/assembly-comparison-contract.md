# Four-assembly comparison contract

**Current status (2026-09-20): real mini execution completed.** [Acceptance](acceptance.md) records results and limitations; full-scale and scientific acceptance remain open.

## Fixed roles and scope

The methods are `mean`, `random`, `within_cluster`, `outside_cluster`. The explicit types are `T`, `Mono_Macro`, `Fibroblast`; broad column is explicitly `Level1`. Only cell-type `/` spelling becomes `_`. There is no type alias dictionary or silent column fallback. Original IDs, missing values and historical `SVC_cluster` identity remain intact.

The zero-gene baseline contains the three historical carriers' own Level1, SVC_cluster and coordinates. Its expression and historical Leiden columns are never used. Historical SVC_cluster is the designated comparison baseline, not independent biological truth; cluster identity is local to each broad type.

For each type, take exact observation IDs shared by all methods and baseline, plus exact gene names shared by all methods. Do not pair by position or fill missing genes. On these IDs, compare each method's own spatial coordinates against baseline by ID, requiring matching shape, finite values and exact equality. An invalid common coordinate blocks that type's clustering/metrics/plot; excluded rows remain outside the comparison. No rescaling or automatic tolerance relaxation.

Each method gets its own working copy, normalize_total(1e4), log1p, PCA, neighbours and Leiden. Retain seed42 and resolutions0.6/0.7/0.8. Baseline labels do not enter expression clustering. Preserve coverage/exclusions, all metrics and contingencies; do not select a winner automatically.

## Executable inputs and outputs

- [Baseline builder](../../../scripts/prepare_assembly_baseline.py): backed obs/spatial reads, duplicate/blank ID protection, zero gene axis, source SHA256 and coordinate provenance.
- [Real mini configuration](../../../examples/assembly-comparison-real.json): formal SVC.h5ad paths, baseline, explicit types/columns, pixel scale0.2125 and output location.
- [Notebook](../../../reproduce/case/assembly_comparison.ipynb) and [comparison module](../../../revise/analysis/assembly_comparison.py): input hashes, common scopes, independent clustering and saving.

The current baseline is49279×0. Four mini SVCs are each501×13088. Common IDs are T160, Mono_Macro157, Fibroblast68; all share13088 genes. All coordinate gates passed. [Executed Notebook](../../../output/mini-acceptance/20260920/assembly/comparison/assembly_comparison.executed.ipynb) saves coverage, scope summary, effective parameters, input hashes, baseline provenance,36 metric rows,36 contingency tables and12 spatial figures. It does not overwrite inputs.

```bash
python scripts/prepare_assembly_baseline.py --output output/assembly-new/original_spatial.h5ad
python scripts/execute_acceptance_notebook.py --repo . --source reproduce/case/assembly_comparison.ipynb --comparison-config examples/assembly-comparison-real.json --output output/assembly-new/notebook
```

For a new run, copy the JSON and explicitly set its baseline_path, four method_paths and output_dir to the new artifacts; the checked-in real JSON describes the completed mini run. Existing baseline/Notebook outputs are not silently overwritten by the helper scripts.

## Historical Notebook source and path evidence

The four preserved case Notebooks establish where the real reconstruction
inputs and earlier spatial carriers came from. They are path/provenance
evidence, not evidence that the four new comparison methods have been
generated.

| Historical Notebook | Raw/reference inputs named by the Notebook | Historical spatial carrier | Handoff meaning |
|---|---|---|---|
| [`Xenium_sc_SVC_T.ipynb`](../../../reproduce/case/Xenium_sc_SVC_T.ipynb) | `../../raw_data/Real_application/P2CRC_Xenium.h5ad` and `../../raw_data/Real_application/adata_sc_all_reanno.h5ad` | `../../results/sc_SVC_case/P2CRC_Xenium/T/spatial.h5ad` | T-cell carrier with `SVC_cluster`; expression-side `expr.h5ad` is not the spatial baseline. |
| [`Xenium_sc_SVC_Monocyte.ipynb`](../../../reproduce/case/Xenium_sc_SVC_Monocyte.ipynb) | The same P2CRC Xenium and reference H5ADs | `../../results/sc_SVC_case/P2CRC_Xenium/Mono_Macro/spatial.h5ad` | Macro/TAM carrier with `SVC_cluster`; expression-side `expr.h5ad` is not the spatial baseline. |
| [`Xenium_sc_SVC_Fibroblast.ipynb`](../../../reproduce/case/Xenium_sc_SVC_Fibroblast.ipynb) | The same P2CRC Xenium and reference H5ADs | The Notebook guide names `sc_SVC_spatial.h5ad`, while the local file is `../../results/sc_SVC_case/P2CRC_Xenium/Fibroblast/spatial.h5ad` | CAF carrier with `SVC_cluster`; the filename drift must be resolved explicitly. |
| [`VisiumHD_sp_SVC.ipynb`](../../../reproduce/case/VisiumHD_sp_SVC.ipynb) | `../../raw_data/Real_application/P1CRC_HD.h5ad` and the matched reference | `../../results/sp_SVC_case/P1CRC/sp_SVC.h5ad` | hST/VisiumHD historical path; it is not the P2CRC Xenium T/Macro/CAF baseline. |

## Metadata-only audit of the local snapshot

The audit below used `h5py` to read HDF5 metadata, `/obs/_index`, selected
categorical codes, and `/obsm/spatial`. The `/X` matrix was not loaded or
examined. Counts are metadata evidence only; they do not establish biological
truth or method quality.

| H5AD | Native shape | IDs / coordinates | Broad and subtype metadata |
|---|---:|---|---|
| `raw_data/Real_application/P2CRC_Xenium.h5ad` | `340837 × 422` | Unique string IDs (`0` … `340836`); `obsm["spatial"]` is `(340837, 2)`, finite | `obs["Level1"]`: T `23392`, Mono/Macro `31675`, Fibroblast `26152`; `SVC_cluster` is absent |
| `raw_data/Real_application/adata_sc_all_reanno.h5ad` | `222815 × 18071` | Unique string IDs; no spatial coordinate key | Reference-only `Level1`/`Level2` labels include T, Mono/Macro, and Fibroblast; `SVC_cluster` is absent |
| `results/sc_SVC_case/P2CRC_Xenium/T/spatial.h5ad` | `14219 × 343` | Unique IDs, all present in raw P2CRC IDs; `spatial` is `(14219, 2)` and exactly matches raw coordinates by ID | `Level1=T`; `SVC_cluster` has 7 classes: `0:3480, 1:2463, 2:2303, 3:2003, 4:1550, 5:1403, 6:1017` |
| `results/sc_SVC_case/P2CRC_Xenium/Mono_Macro/spatial.h5ad` | `17605 × 343` | Unique IDs, all present in raw P2CRC IDs; `spatial` is `(17605, 2)` and exactly matches raw coordinates by ID | `Level1=Mono_Macro`; `SVC_cluster` has 8 classes: `0:3981, 1:3795, 2:2635, 3:2504, 4:1954, 5:1608, 6:715, 7:413` |
| `results/sc_SVC_case/P2CRC_Xenium/Fibroblast/spatial.h5ad` | `17455 × 343` | Unique IDs, all present in raw P2CRC IDs; `spatial` is `(17455, 2)` and exactly matches raw coordinates by ID | `Level1=Fibroblast`; `SVC_cluster` has 10 classes: `0:3806, 1:2754, 2:2529, 3:2453, 4:2142, 5:1404, 6:800, 7:776, 8:621, 9:170` |

The three historical spatial carriers have disjoint IDs and contain `49279`
combined observations. Their IDs are all raw P2CRC IDs and their coordinate
rows match the raw object exactly, but they do not cover every raw broad-label
cell: the carrier-to-raw-label matches are T `13987/14219`, Macro
`17357/17605`, and CAF `17236/17455`. The metadata audit cannot establish why
the remaining IDs have a different raw `Level1` label, so it does not infer a
selection rule from that discrepancy. The carrier's own `Level1` field remains
the historical routing field.

The existing spatial files also contain historical `leiden_0.6`, `leiden_0.7`,
and `leiden_0.8` columns. The comparison contract deliberately ignores those
columns and reruns expression Leiden for each of the four generated method
matrices.

## 后续执行任务与验收

H2–H4 are implemented and exercised on the mini inputs. Full-scale execution must retain the same upstream configuration across methods and use new output directories. Coordinate perturbation regressions must prevent scores for the affected common scope. Review coverage before interpreting ARI/NMI; smaller mini intersections do not establish full-sample biological validity or a final assembly default. The consumer links these results and does not duplicate this experiment.
