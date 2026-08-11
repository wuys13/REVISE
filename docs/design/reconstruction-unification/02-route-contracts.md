# Expected Runtime Route Contracts

Parent index: [Reconstruction Unification Design Package](README.md)

本文定义实施和测试共同使用的预期线路。Current fact 描述 baseline 已有算法步骤；Decided target 描述统一后必须观察到的调用边界。任何 route trace 测试只能观察这些真实调用，不能另造一套 production state machine。

## Shared lifecycle

Decided target：所有 runtime 路线共享：

```text
validated effective config
→ mode-specific preparation
→ GA facade
→ route-owned LR
→ finalize SVC
→ mode-specific publication/evaluation
```

GA facade 只调用 `OTKernel.annotate`。LR 可以调用 `LocalAnchoringKernel → OTKernel.annotate`，或自行构造 masses/cost/support 后调用 `OTKernel.couple`。LR 中不得出现 `GlobalAnchoringKernel`。

## Application routes

| Route | Expected stages | Conditions and outputs |
| --- | --- | --- |
| `sp-SVC` | 三个 Application preprocessing → slash normalization → spatial/reference overlap contract → GA Level1 annotation → per-cell-type graph/top-k → masses/cost/support → matrix coupling → graph aggregation → finalize | LR 可以按既有小样本/空 support 条件保留原表达；输出一个 primary AnnData |
| `sc-SVC` | reference filter → spatial preprocess → reference preprocess → reference obs 投影到 Level1/Level2 → slash normalization → spatial-only overlap subset → GA Level1 annotation → select broad cohort → LR Level2 annotation → gene/spatial/joint graph → Leiden/align score/best resolution → `SVC_cluster` → reverse LR annotation to reference → finalize pair | `alpha` 与 `resolutions` 来自 Application YAML；输出顺序固定为 `(spatial_adata, expression_adata)` |
| `sc-SVC-sr` | `ensure_all_cells_in_spot` → 三个 Application preprocessing → slash normalization → GA Level1 annotation → spot posterior → mandatory virtual-cell allocation → optional per-cell-type graph/matrix coupling → optional graph aggregation → finalize | mandatory allocation 不受 optional refinement 开关控制；graph aggregation 启用时 primary 是 graph-aggregated AnnData |

### Standard sc-SVC axis contract

Decided target：

- GA 和第一轮 LR annotation 的 target observation 轴保持完整输入/选中 cohort，不预先删掉零质量行来隐藏问题。
- spatial 输出使用 spatial/reference overlap gene axis。
- expression 输出保留过滤后 reference 的完整 gene axis。
- expression obs 最终公开列为 `Level1`、`Level2`、`SVC_cluster`、`Confidence`；`Patient` 仅参与过滤，随后移除。
- best resolution 在最大 alignment score 并列时取列表中最后一个，与 baseline `values[-1]` 一致。

## Benchmark routes

Benchmark preprocessing 保持 route-owned；不改造成 Application 三函数。六类 YAML 是唯一 route selector。

| Route YAML | Case expansion | Expected stages | Primary outputs/evaluation |
| --- | ---: | --- | --- |
| `segmentation` | 4 (`seg_1..seg_4`) | Benchmark prep → GA annotation → segmentation evaluation when metadata exists → replace/candidate partition → graph/matrix coupling → merge | reconstructed sp-SVC + existing segmentation metrics |
| `bin2cell` | 1 | Benchmark prep → GA annotation → replace/candidate partition → graph/matrix coupling → merge | reconstructed sp-SVC + existing evaluation |
| `batch_effect` | 16 (4 discovered spot sizes × 4 batches in the frozen benchmark fixture) | input/GT mapping → Benchmark prep → GA annotation → mandatory virtual-cell allocation → optional graph aggregation/matrix coupling | raw allocation and graph-aggregated output when enabled + existing evaluation |
| `spot_size` | 4 (20/50/100/200) | input/GT mapping → Benchmark prep → GA annotation → mandatory virtual-cell allocation → optional graph aggregation/matrix coupling | raw allocation and graph-aggregated output when enabled + existing evaluation |
| `gene_panel` | 1 | benchmark/imputation prep → GA annotation → gene uncertainty → in-panel and all-panel subclustering → per-cell-type matrix coupling → gene imputation/prune | both in-panel and all-panel imputed outputs + existing evaluation |
| `gene_dropout` | 1 | benchmark/imputation prep → GA annotation → gene uncertainty → in-panel and all-panel subclustering → per-cell-type matrix coupling → gene imputation/prune | both in-panel and all-panel imputed outputs + existing evaluation |

Case cardinality contract is therefore `4 / 1 / 16 / 4 / 1 / 1`. A dataset that cannot supply the frozen batch-effect spot-size fixture must fail or report the actual preflight discrepancy; it must not silently change the documented benchmark definition.

## Benchmark noise deletion boundary

Current fact：baseline SR Benchmark strategy contains `_inject_sr_spatial_leakage_noise`, `sr_noise_*`, `sr_spatial_noise`, `st_input_noisy` and a noise-specific GA override.

Decided target：上述旧污染全部删除，SR Benchmark 直接使用共享 `RunnerBackedStrategy.global_anchoring`。以下内容不是 noise，必须保留：

- `sr_graph_agg_*` 与 confidence/anchor graph aggregation；
- mandatory virtual-cell allocation；
- gene panel/dropout imputation；
- in-panel/all-panel outputs；
- route evaluation 和普通 GA。

## Trace assertions

每条路线至少断言：

1. strategy resolution 与 YAML route 一致；
2. preprocessing/GA/LR/finalize 的顶层顺序；
3. annotation 或 coupling 的内层调用次数与顺序；
4. 条件分支未执行时有明确条件，而不是另一路径 fallback；
5. output role、primary object 和 evaluation 数量；
6. Benchmark case expansion 数量与每个 case 的 effective seed/provenance。

继续阅读：[OT boundaries and contracts](03-ot-boundaries-and-contracts.md)。
