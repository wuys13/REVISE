# Real-data handoff contract for assembly comparison

**Status (2026-09-19): metadata handoff only.** This document records the
inputs and roles required by the existing comparison Notebook. No real
reconstruction, four-method assembly run, or biological interpretation was
performed for this handoff.

The executable entrypoints are
[`reproduce/case/assembly_comparison.ipynb`](../../../reproduce/case/assembly_comparison.ipynb)
and
[`revise/analysis/assembly_comparison.py`](../../../revise/analysis/assembly_comparison.py).
The Notebook reads already generated H5AD files; it does not reconstruct or
publish them.

## Fixed roles

The comparison has exactly four generated expression methods:

```text
mean, random, within_cluster, outside_cluster
```

Each method is copied, normalized to a per-observation total of `1e4`,
transformed with `log1p`, and clustered independently with the same PCA,
neighbor, and Leiden settings. The generated method matrices are declared to
be finite, nonnegative, unlogged linear expression before that working-copy
step.

The existing spatial `SVC_cluster` field is the designated **Ground Truth
spatial subtype baseline** for this comparison's label role. It is used for
the baseline label panel, contingency tables, ARI/NMI, and no expression
clustering. The spatial
coordinates are used only for plotting after expression clustering. The old
per-type sc-SVC outputs are therefore baseline carriers and historical input
evidence; **old sc-SVC is not a fifth comparison method**.

The implementation fixes the other comparison rules:

- correspondence uses exact observation IDs and exact gene names; row
  position is never used and missing values are never filled with zero;
- each broad type is filtered independently, then all methods and the baseline
  are intersected on real IDs and generated methods on common genes;
- missing types, IDs, genes, and baseline labels remain visible in the
  coverage/metrics tables;
- resolutions `0.6`, `0.7`, and `0.8` are retained with seed `42`; no best
  resolution or winning method is selected automatically;
- the native H5AD files are hashed before and after the Notebook run and must
  remain unchanged.

The implementation source for these rules is the module's
[`DEFAULT_METHODS` and `INPUT_ASSUMPTIONS`](../../../revise/analysis/assembly_comparison.py#L25-L33),
the four-method path check and read-only loader
[`load_assembly_inputs`](../../../revise/analysis/assembly_comparison.py#L82-L115),
and the independent Leiden/evaluation path
[`compare_assembly_methods`](../../../revise/analysis/assembly_comparison.py#L472-L533).

## Required H5AD contract

The four method paths must resolve to four generated H5AD files. Each method
object must provide unique, non-empty `obs_names` and unique `var_names`; its
`X` must be numeric, finite, nonnegative, and declared unlogged linear
expression. Its broad-label column must be shared with the baseline, or the
configuration must be changed before running.

The one baseline path must resolve to a single whole-sample spatial H5AD that
contains:

```text
obs[SVC_cluster]       designated baseline subtype labels
obs[<broad column>]    T/Macro/CAF routing labels
obsm["spatial"]        finite coordinates with at least two columns
```

The baseline does not need to be used as an expression method, but its
observation IDs must overlap the generated objects for a type to be scored.
The real-data configuration in
[`examples/assembly-comparison-real.json`](../../../examples/assembly-comparison-real.json)
uses `Level1` because that is the broad-label field present in the local
historical spatial H5ADs, and includes both observed Macro spellings
(`Mono/Macro` and `Mono_Macro`). If a future whole-sample handoff exposes only
`revise_Level1`, update both the column setting and the metadata check before
running; do not rely on the module's default silently.

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

The current comparison Notebook itself uses the four method placeholders
`output/assembly/{mean,random,within_cluster,outside_cluster}.h5ad` and
`output/assembly/original_spatial.h5ad`; it does not point at the historical
per-type files. Its current default broad column is `revise_Level1` and its
current Macro aliases omit both `Mono/Macro` and `Mono_Macro`. The example
configuration supplies the observed local fields explicitly.

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

## Confirmed baseline sources and later run preparation

The user confirmed that Ground Truth refers to the already reconstructed data, not the raw input. The three historical spatial carriers above contain those labels and are the available baseline sources. There is no single preassembled whole-sample baseline H5AD in the local snapshot:
the raw P2CRC input lacks `SVC_cluster`, while the existing `SVC_cluster`
carriers are three separate per-type spatial files. There are also no four
generated method H5ADs under the Notebook's `output/assembly/` placeholders.
The metadata-only scan of historical local `*.h5ad` inputs (before this round’s synthetic evidence generation) found only five objects with
both `obs/SVC_cluster` and `obsm/spatial`: the three per-type carriers above and
two parity copies of the T carrier. It found no whole-sample baseline.
Consequently, the JSON example is a runnable **path/configuration template
once** those five files exist, but it is intentionally not runnable against
the current checkout. It does not silently substitute the raw H5AD, one
per-type carrier, or old sc-SVC output for a missing method.

Before a real comparison run, the executor can assemble a labels-and-coordinates baseline from the three confirmed carriers, preserving IDs, their own Level1, SVC_cluster and spatial coordinates. Duplicate IDs must be rejected; cluster labels are interpreted within each broad type. This is input preparation from existing data, not missing labels or a request for the user to supply new annotation. Generate the four method files separately before comparison. Then run the existing Notebook and review
the coverage table before interpreting any ARI/NMI or spatial figure. This
handoff does not perform that generation, real reconstruction, or scientific
acceptance.


## 后续执行任务与验收

本节为已确定的执行约定，尚未实现汇集工具、生成真实四方法产物或新增坐标校验。双库审阅统一从 [cross-repo-review.md](cross-repo-review.md) 进入。

| 步骤 | 执行规则 | 验收证据 |
|---|---|---|
| H2 基线汇集 | 从上述 T、Mono_Macro、Fibroblast 三份 carrier 提取原 ID、Level1、SVC_cluster、spatial；拒绝空或重复 ID，包括跨文件重复；保留各文件来源与 SHA256，不改写源文件 | 汇集前后行数、每类型每 cluster 计数、ID/坐标逐行一致；旧快照合计 49,279，仅作核对参照 |
| 基线载体 | 单独标签/坐标 AnnData，零基因轴；不读取或拼接历史 X 作为比较表达。保留原 cluster 值，cluster 身份按 (broad type, SVC_cluster) 解释 | 新分群流程完全不使用基线表达和历史 Leiden；不同大类的数字 cluster 不被混作同一类 |
| H3 四方法生成 | 使用现有 `reconstruct.py --config` 或 `run_application`；四份配置固定原始输入、同一 reference（使用外部配置时固定同一 `--reference-config`）、QC/GA/LR 与 seed，仅修改 output.ist_mapping 和独立输出位置；不承诺已有 GA 缓存 | 保存四份有效配置、输入摘要、版本和完成状态；各目录正式 SVC.h5ad 传入 Notebook，不用旧 expr.h5ad 代替 |
| H4 共同范围 | 沿用每类型真实 ID 交集与四方法基因名交集；保留原 ID，不拼造映射、不补零；展示类型路由不一致和各侧排除数量 | coverage 含原始数量、共同数量和排除原因；无共同范围的类型明确不可比较，不产生误导分数 |
| H4 坐标门槛 | 当前绘图只读取基线坐标，因此需在绘图前补校验：按真实 ID 对齐后逐方法比较 spatial 前两维，形状/有限性及坐标来源单位一致；原样复制坐标默认要求相等，发现差异先调查，不自动平移或放宽容差 | 不一致案例应阻止该类型比较；通过后才使用基线坐标作共同绘图框架 |
| Notebook | 使用同一预处理与 Leiden seed/参数，保留所有预定分辨率；按类型解释 ARI/NMI、列联和空间图 | 输入哈希不变；输出参数、共同范围、图表和全部分辨率结果；不自动选赢家 |

这些任务由 REVISE 侧负责准备；分析库可复核或消费比较结果，不要求复制一套实现。现有 JSON 为路径模板，尚不是一键生成四方法的执行器。H3 允许分别重跑现有主流程，不新增缓存复用协议；若相同配置仍产生上游标签或范围差异，先记录并解释，不能把其影响全部归因于 assembly。
