# Moran 与 EMT 功能方法

框架节点：[2.4 Moran](analysis-framework.md#node-2-4)、[2.5 EMT](analysis-framework.md#node-2-5)、[3.3 EMT 空间场](analysis-framework.md#node-3-3)。

下文保留既有方法与约束，供 Moran/EMT 两个 overview 行及按需审计引用；不因本轮文档
更新改变计算定义，也不把已有执行记录写成本轮验收。

本页负责 overall 中的全基因 Moran 与 EMT coverage/score，以及 localization 使用的
EMT spatial field。它不新增 marker、DEG、trajectory、interaction 或应用基因地图。
输入和 projection 见 [input-views.md](input-views.md)，partition/window/Region 见
[impact-analysis.md](impact-analysis.md)，记录和图的保存见
[outputs-and-test-plan.md](outputs-and-test-plan.md)。

Moran 的 full-side 与 shared-valid、EMT coverage 与 EMT score 是四个独立小问题；
它们在 Notebook、record 和 HTML 中分别保存、分别判读。

正式Moran/EMT观察由同一结论记录提供；本页维护完整方法与解释边界，短目的与节点对应由本地静态内容映射供两种展示复用。

## 1. Full-gene space and analysis states

**Purpose.** 保留 Raw 与 Recon 各自完整的可用 gene space，区分测量能力与统计量
是否可计算，使 coverage/score 或 Moran 的 NA 不被误读为零。

**Computation.** 以两侧 gene ID 并集登记 `analysis/inputs/gene_availability.csv`，
再为每个 aspect 写自己的 availability 表。每个 gene/侧至少记录 provided or
unmeasured、表达状态、实际使用/排除原因、expression view、unit/edge support 和
状态。Raw-derived QC/HVG 只属于 partition，不改变 Moran/AUCell 的 gene universe。

状态至少区分 `computed`、`unmeasured`、`low_coverage`、`insufficient_support`、
`not_computable` 和执行 `blocked`。Zero expression、constant expression、无边、
未测量和缺资源不能互换。

**Meaning.** 完整 gene space 描述“每侧输入能测到什么”；shared-valid 只描述两侧
共同有效 gene 的对应变化。它不增加 Raw 的未测量能力，也不把 projection 变成原生
unit-level measurement。

**Direction and conditions.** Recon 可以有 Raw 没有的 gene；Raw 未测量时不补零、
不裁剪 Recon。单侧有效结果可以发布，delta 只有两侧都有效时生成；缺少软件、损坏
资源、轴违反契约或代码异常是执行失败，不能伪装成 scientific unavailability。

<a id="moran"></a>

## 2. Moran：full-side distributions

**Purpose.** 在固定 cell type spatial units 内，描述每侧完整有效基因的空间
autocorrelation 分布，以及重建前后各自的有效 gene 数和摘要变化。

**Computation.** 每个 `gene × scope × comparison_id` 保存 Raw/Recon Moran、两侧
状态、`n_units`、有效 `n_edges`、graph ID、expression view 和 reason。按 Raw exact
Level1 建立 group；VisiumHD 保留 `All` 及三个 parent，Xenium 只保留三个 parent。
同一 group 的两侧使用同一 shared observation/coordinate graph；最低支持
`n_units >= 51`。当前 HD graph 可沿用历史 6-neighbor 语义，但实际参数、坐标来源、
权重变换、components、isolates、版本和 observation 顺序必须写入 graph audit。

两侧各自执行 `normalize_total(target_sum=1e4) → log1p`，随后在各自完整 gene axis
上计算描述性 Moran I；不做 permutation inference、p 值或多重检验。图只构建一次，
`compare_moran` 只消费同一 weights 和两侧表达，不隐式构图或归一化。

**Meaning.** `moran_all_valid` 的 median、Q1、Q75、有效数和分布回答“该侧完整
输入的空间结构如何”；它不回答哪些 unit 发生变化，也不自动表示 biological
correctness。Raw 与 Recon 的有效 gene 数可以不同，不能强行合并。

**Direction and conditions.** Moran summary 可以是 higher、lower、mixed 或 unknown。
只有图、unit/edge support、表达变换和 full-side availability 有效时才解释方向；
常量、无边、未测量或低支持写状态，不能用 NaN→0。结果表和展示表同时保留：

- `analysis/moran/moran_summary.csv`：Raw/Recon full-side `n_valid`、median、Q1、
  Q75 与状态；
- `analysis/moran/moran_distribution_summary.csv`：`gene_set=all_valid` 的逐侧
  compact distribution，以及 `shared_valid` 的对应摘要；
- `analysis/moran/gene_availability.csv`、graph weights 和有序 observations：
  重新加载与图审计所需证据。

## 3. Moran：shared-valid paired delta

**Purpose.** 在两侧均有效的同名 gene 上，保存逐 gene `Recon − Raw`，作为与
full-side 分布并行的配对视角。

**Computation.** 先按每侧 availability、同一 group graph 和同一 expression rule
筛选 shared-valid gene，再逐 gene 计算 Moran delta；保存 shared count、两侧值、
delta、状态、`n_units/n_edges`、graph ID 和 gene IDs。不能因为 shared gene 较少就
替代 full-side 表，也不能把两侧的 full gene union 改成 overlap panel。

**Meaning.** delta 描述当前两侧在共同有效 gene 上的自相关变化；与 full-side
分布一致或不一致都应如实报告。它不是 paired spatial observation 的 biological
effect，也不支持把 Moran 变化解释为 pathway activation。

**Direction and conditions.** 只有 shared-valid gene 数、图和两侧值均有效时才写
paired direction；shared 数不足或一侧缺失时保留 `insufficient_support`/
`unavailable`。Moran overview 行内部必须把 `moran_all_valid` 与
`moran_shared_valid` 作为两个子问题保留，但不增加 overview 行；正文是否展开由本批
记录和审计需要决定。

<a id="emt-coverage"></a>

## 4. EMT resource coverage

**Purpose.** 回答固定 EMT resource 在两侧表达载体中有多少基因可供评分，避免把
coverage 扩大写成新增实测基因或生物学改善。

**Computation.** 首轮固定本地 Hallmark 2025.1 GMT 中的
`HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`，记录资源路径、版本、digest、物种、
gene ID 类型和实际 `resource_total`。当前批准目标为 200 genes；如果 digest 的
实际总数不同，保存真实总数和原因，不静默改写为 200。按 scope/side 保存资源内
available/detected 数、coverage、provider 参数和状态，另写 resource gene rows。

**Meaning.** coverage 是表达载体对该资源的支持，不能解释成“EMT 活性提高”，也
不能证明某个 gene 在 Raw 中未测量却被恢复。Xenium cluster-mean projection 的
coverage 仍带有 reference-derived 条件。

**Direction and conditions.** coverage 可上升、下降或保持不变；无命中、资源缺失、
gene ID 不匹配或单侧不可用写状态与 reason。coverage 与 score 不共用一条结论或
一个总分；HTML/Notebook 单独链接 `availability.csv` 和 `resource.json`。

<a id="emt-score"></a>

## 5. EMT score / AUCell

**Purpose.** 在完整 Raw/Recon expression space 中描述固定 EMT resource 的 unit-level
score，并为 localization 提供可审计的 score/delta 空间字段。

**Computation.** 首轮 scorer 为 AUCell，复用现有 [aucell.py](../../../revise/analysis/advanced/aucell.py)
wrapper 的能力；参数固定记录 `AUC_threshold=0.01`、`seed=42`、provider 版本和
实际 per-side rank cutoff。`0.01` 按 provider 的检出基因 quantile 规则记录，不写
成固定前 1% gene。Raw/Recon 在各自完整 gene space 评分，不裁剪共同 panel。

iST `paired` 必须先按 [输入视图](input-views.md) 求 cluster mean，再在 projected
expression 上评分，即 `AUCell(mean expression)`；不能先给 donor cells 分别评分后
再平均却沿用同一 expression view。score observation 是 spatial unit；同一 cluster
共享 projected expression，不增加独立 reference-cell observation。

保存 `analysis/pathway_activity/scores.csv.gz` 的每个 `unit × EMT pathway` 两侧
score、paired delta 和状态；`summary.csv` 按冻结 Level1/指定 anatomy 汇总有效 n、
mean、median、Q1、Q3、coverage 与 paired delta 的 mean/median；`resource.json`、
`cutoff_audit.csv`、`resource_genes.csv` 保存 resource、rank length、ties/zero
ranks、seed、provider、顺序和逐侧动态 cutoff。

**Meaning.** score delta 是当前完整 gene space、resource coverage、rank cutoff 和
expression carrier 条件下的描述性变化。`median(Recon) − median(Raw)` 不能替代
逐 unit paired delta 的 median。它不单独证明 biological EMT activity 增幅；尤其
Xenium 结果要与 cluster-mean、coverage 和 cutoff audit 一起阅读。

**Direction and conditions.** score 可升高、降低、混合或无法判定；只有两侧 score
和 paired IDs 有效时才写 delta。coverage 受限时可以保留有效 Recon score，但结论
必须写 Raw limitation。没有新增人为 coverage threshold；缺少 AUCell 依赖时不自动
用 `score_genes` 替代。

<a id="emt-spatial-fields"></a>

## 6. EMT spatial field

**Purpose.** 把已保存的 EMT score/delta 定位到同一 spatial unit 和坐标，供局部双 baseline
空间阅读，而不是增加新的功能分析。

**Computation.** 每行至少保存 `scope`、`unit_id`、`feature`、坐标、Raw/Recon/delta、
两侧状态、expression view 和坐标单位；使用 input audit 与 graph/window 的同一
observation basis。局部 spatial 入口只消费 `scores.csv.gz`/空间字段，不从 Moran 表或 H5AD
重新计算 score。

**Meaning.** spatial field 显示 score/delta 在当前 carrier 下的位置。iST 的字段是
cluster-mean expression projection；hST 使用其已登记的 reconstructed carrier；
二者都不是追加的 spatial validation。

**Direction and conditions.** 只有 spatial IDs、坐标、score condition 和表达载体
有效时才绘图。parent-specific image 归属于对应 scope row；common score distribution
或资源图只嵌入一次，其他 row 链接同一原图。缺失一侧保留单侧 field 和状态，delta
留空。

## 7. Gene/functional Notebook 块

Notebook 只展示核心参数、短调用、检查、保存、图和一段 brief batch observation；
长实现由本地 helper 承载。现有分析块可以引用 [notebook_analysis.py](../../../reproduce/case/reconstruction_impact/notebook_analysis.py)
和 [notebook_helpers.py](../../../reproduce/case/reconstruction_impact/notebook_helpers.py)。

| 块 | staged short calls | audit/table | plot / brief observation |
| --- | --- | --- | --- |
| Moran full/shared | register gene union → build group graph → normalize/log → compute both views → save states | gene availability、moran summary/distribution、graph weights/observations | distribution/Q75 与 shared scatter；分别观察 full-side 与 shared-valid |
| EMT coverage | load resource → resolve IDs → count available/detected genes → save availability | resource、resource genes、coverage/reason | coverage table/summary；说明资源覆盖和 expression carrier |
| EMT score | score Raw/Recon with fixed params → pair IDs → summarize → save fields | scores、summary、cutoff audit、status | score distribution/delta；说明 paired n、cutoff 与是否有 Raw limitation |
| EMT spatial field | read saved scores → join coordinates → save feature field | field CSV、coordinate/expression provenance | parent spatial fields；一句定位观察，不升级为机制 |

HTML 的事实入口是 [report_records.py](../../../reproduce/case/reconstruction_impact/report_records.py)
和 [report_html.py](../../../reproduce/case/reconstruction_impact/report_html.py)；它们
消费保存表，不重算 Moran/AUCell。若后续需要重新验证，沿用保存的 source/executed
identity、错误/警告流和 CSV/JSON 可重载证据；仅刷新阅读产物时保留原 execution identity。
