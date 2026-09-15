# Partition、空间多样性与 Region 方法

框架节点：[2.1–2.3 分群](analysis-framework.md#node-2-1)、[3.1–3.9 空间与局部状态](analysis-framework.md#node-3-1)、[4.1–4.5 Region](analysis-framework.md#node-4-1)。

下文保留既有方法与约束，供五行报告和按需审计引用；不因本轮文档更新改变计算定义，
也不把已有执行记录写成本轮验收。

本页负责五行中的分群、局部双 baseline 和 Region，以及 localization 需要的空间条件中由 reconstruction-impact 直接计算的
方法：partition complexity、matched-K、Raw Level2、anatomy、window、局部多样性、
State/Gain Region 和 scale sensitivity。它不定义输入载体或 Moran/AUCell；对应
边界见 [input-views.md](input-views.md) 和 [gene-and-function.md](gene-and-function.md)。
保存格式与报告消费见 [outputs-and-test-plan.md](outputs-and-test-plan.md)。

所有方法段用 Purpose → Computation → Meaning → Direction and conditions。`K`、
`Neff`、evenness、区域面积和 changed fraction 都是描述性事实；它们不会被合成
一个 improvement score。

## 1. Shared cohort and partition prerequisites (on-demand audit)

**Purpose.** 固定一个 scope 的 Raw-defined spatial units、feature selection 和
seed，使 partition、changed units、window 和 Region 共享可审计的 comparison basis。

**Computation.** 对 global 或 parent scope 先验证 reconstruction IDs 是 Raw 子集，
再选择 scope，最后执行 Raw-defined QC。VisiumHD 当前 partition QC 为
`min_genes=50`、`min_cells=3`，排除线粒体基因；Raw 的 canonical normalize/log 和
Seurat-v3 HVG 使用 `n_top_genes=2000`，feature set 共享给重建表达。两侧仍在各自
graph preparation 中完成其实际 normalize/log，不能写成只有 Raw 被预处理。partition
backend 为 igraph，iterations 为 2，seed 由配置传入（当前 route 为 42）。

Raw-derived HVG 只约束 partition graph，不是 Moran/AUCell 的共同 gene panel。所有
scope 的抽样、实际 ID、QC 保留与排除原因写入 observations 和 partition audit；不能
从历史表或文件顺序恢复 cohort。

**Meaning.** 后续的 assignment、Moran spatial graph 和 window denominator 只在
同一保存的 observation basis 上解释；`n_units`、`n_changed` 和 valid windows
不能混用 full tissue、paired parent 和 donor-cell axis。

**Direction and conditions.** QC/HVG 或轴检查失败时状态是 `blocked`；支持不足可
发布 `insufficient_support`。上游 partition 限制只传递给依赖 assignment 的局部
结论；Moran/EMT 根据各自输入条件独立判读。

<a id="partition-complexity"></a>

## 2. Partition complexity diagnostic

**Purpose.** 在保留各侧自身 cluster 语义的前提下，描述 K、cluster size 和 split
pattern 的变化；它不是 changed-unit headline。

**Computation.** VisiumHD 用 Raw 侧选择的 resolution（global 在候选
`0.3/0.5/0.8` 中按 Level1 ARI 选择，parent 使用 `0.5`），在同一 resolution
运行另一侧，并保存 candidate resolution、Raw/Recon K、cluster sizes、ARI、
contingency（absolute 与 normalized 及分母）和选择理由。Xenium 不做重建侧
resolution sweep；它比较 reference-resolution Raw Leiden 与 fixed final clusters，
`comparison_kind=raw_reference_vs_fixed_final_clusters`，表头写 `Fixed final-cluster K`。

**Meaning.** K 增加只表示在当前 graph/feature 条件下观察到更多 partition，K 减少
只表示更少；ARI 和 size distribution 是结构诊断。它们不能单独证明 resolution
更真实，也不等于 unit identity change。

**Direction and conditions.** 复杂度方向可以是 higher、lower 或 unchanged；不
预设“更多就是更好”。Xenium fixed-final 与 HD same-resolution 是不同比较，不能
混在同一 headline。该方法不需要 matched-K 成立，但只支持复杂度诊断。

<a id="matched-k"></a>

## 3. Matched-K assignment

**Purpose.** 在 cluster 数量尽量可比时，描述哪些 spatial units 的 assignment 改变，
并明确复杂度不匹配何时限制解释。

**Computation.** 独立扫描另一侧（sc-SVC 时扫描 Raw）以接近目标 K，选择一次最大
重叠 Hungarian mapping，保存 Raw/Recon K、target/actual K、absolute/normalized
contingency、mapping、unit assignments、ARI、changed fraction、balanced change
和 Wilson interval。若实际 K 与目标差距超过 1，标记
`unmatched_cluster_complexity`；保留诊断表，但不进入 headline change 估计。

Global Level1 与 parent-internal assignment 使用不同 scope。parent 内不重新匹配；
global mapping 只用于它的 global comparison。task 名仍按 Raw exact `Level1` 映射。

**Meaning.** `changed_units / paired_units` 是 assignment difference，不是细胞身份
变化；Wilson interval 只描述该 scope 分数的不确定性，不是 biological confidence。

**Direction and conditions.** change 可以上升、下降或零；只有 matched-K 状态满足、
ID/坐标 basis 有效且 denominator 明确时才可作 headline。unmatched 时记录数值供
审计，但下游依赖它的局部方向标为受限。`matched_k_resolution_sweep.csv` 与
`complexity_raw_resolution_sweep.csv` 的语义不能因字段名相似而交换。

## 4. Raw Level2 baseline

**Purpose.** 给每个 parent 一个与 Raw Leiden 平行的 reference-derived baseline，
用于检验局部状态方向是否依赖 Raw Level1 的粗粒度 partition。

**Computation.** 对原始 Raw expression 做 parent selection，再读取配置的
`raw_level2_mapping.reference_h5ad`：VisiumHD 使用 route-native POT，Xenium 使用
TACCO 和 P2CRC reference subset。记录 labels、posterior、mapping method、reference
filter/count、path/digest 和 seed。Raw Level1 在 parent 内是 uniform label，描述性
baseline 为 `Kobs=1`、`Neff=1`、`evenness=1`。

**Meaning.** Raw Level2 是 reference mapping 的比较基线，不是 ground truth，也不
是 Recon-derived Level2。它回答“Raw 的参考映射层级怎样呈现局部 diversity”，不能
被写成重建真值或与 iST `expr.h5ad` 混为一个表达载体。

**Direction and conditions.** `local_vs_raw_level2` 只在 Raw Level2 labels/posterior
和同一 window/draws 有效时解释；reference 缺失、映射失败或支持不足保留状态。
不能因为重建 expression output 存在就跳过配置指定 reference。

<a id="window-support"></a>

## 5. Physical windows and common rarefaction

**Purpose.** 在同一 physical window、同一 paired units 和同一 deterministic draws
下比较 Raw Leiden、Recon 和 Raw Level2 的局部状态，尽量把 local density 的影响
与 label diversity 分开。

**Computation.** 原始坐标先按 route 配置转换为微米。当前候选 square side 为
`16, 24, 32, 40, 56, 80 um`；full anatomy context 和 parent grid 共用 full
Raw Level1 坐标最小值作为 origin，但各自根据 occupancy 选择 scale。保留至少
4 个 parent units 的窗口，按 retained-unit fraction 对 `log(window_side)` 的
chord-distance knee 选择主尺度；保存完整 support curve、candidate scales、origin、
转换系数、有效窗口和 selected scale。

<a id="local-diversity"></a>

### Local diversity

每个有效窗口以 4 个 units 为 rarefaction size、200 个 draws；每个 draw 计算
`Kobs`、entropy、`Neff=exp(entropy)` 和 `evenness=Neff/Kobs`，再按相同 draws 汇总
Raw Leiden、Recon、Raw Level2。结果按两条 baseline 写出：

```text
Raw Leiden  → Recon 及 delta
Raw Level2  → Recon 及 delta
```

seed 由配置传入（当前为 42），并绑定 cohort subsampling、rarefaction、State/Gain
bootstrap 和 scale sensitivity。每个 scope 使用自己的 deterministic RNG；实际
IDs 和 draws 写入 parameters/audit。

**Meaning.** `Kobs` 是观察到的 label 数；`Neff` 是 abundance-aware effective
cluster count；evenness 是条件于 richness 的 balance。三者不是同一个指标。两个
baseline 的 Recon 行可以在阅读上重复，但统计不能重复累计。
摘要表的 Raw/Recon 列分别是有效窗口侧别值的中位数，delta 列是逐窗口
`Recon − baseline` 的中位数；后者通常不等于前两个中位数相减。保留这一区别，
不能为了让摘要表看起来可直接相减而改写已保存的 paired delta。

**Direction and conditions.** 在 `baseline_status` 有效、window support 满足、并且
使用已选主尺度时，Raw Leiden 的 Neff 降低才支持“更集中”的方向；Raw Level2 的
Neff 升高才支持“更细内部状态”的方向。Kobs/evenness 和 scale sensitivity 只作
辅助。条件不满足时方向写 `direction unknown`/受限，不能从另一条 baseline 推断。

## 6. Anatomy context and spatial localization (detail referenced by the five rows)

<a id="anatomy"></a>

### Anatomy context

**Purpose.** 提供与 partition、window scale 和 threshold 调参独立的组织背景，并把
changed units 与 EMT score/delta 定位到可重载的空间单位。

**Computation.** anatomy 在 full Raw Level1 context 上独立计算，当前 source labels
为 `Tumor` 与 `Intestinal Epithelial`：前者存在而后者不存在为 `Tumor`，反之为
`Normal`，二者都存在为 `Interface`，其余为 `Other`。窗口可包含其他 cell types；
分类按两种 source label 是否出现。parent windows 使用实际 observed anatomy units
映射到 context，不能用空 anatomy cell 的几何中心替代。changed-unit 图使用 matched-K
已保存 assignments；EMT spatial field 使用已保存的 unit-level score/delta。

**Meaning.** anatomy 是 context 和分层分母，不是结果调参变量。changed-unit fraction
描述 assignment difference；EMT spatial field 描述当前表达载体和评分条件下的
localization，不能升级为独立 biological validation。

**Direction and conditions.** anatomy/coordinate audit 有效时可报告区域分层；缺少
context 或空间 field 时保留 `unavailable`。parent-specific image 只归属于该
parent row；full/shared common image 只嵌入一次并从其他 row 链接。

<a id="changed-units"></a>

### Changed units and EMT localization

changed-unit 图使用 matched-K 已保存的 assignment；它的 `changed_units/paired_units`
只描述该 assignment edge 下的差异，不能改写为 biological identity change。EMT
spatial field 只消费已保存的 score/delta；parent-specific 图留在对应 scope row。

<a id="region-reliability"></a>

## 7. State Region 与 Gain Region

**Purpose.** 区分“重建后局部状态达到高多样性”与“相对 Raw Leiden 的局部增益”，
避免把 Region mask 写成变化最大的区域或统一 improvement score。

**Computation.** State threshold 作用于 rarefied `Neff_recon`；Gain threshold 作用
于 `delta_neff > 0` 的分布。两者分别使用 survival-curve segmented breakpoint 和
500 次 window bootstrap，保存 threshold、valid/total bootstrap、CI width、有效
Recon Neff range、valid windows/units、region windows/units、area、area fraction
和 `region_windows/valid_windows`。结果表中 `region_type` 区分 State/Gain；State 是本层的主要状态展示，
不是重建收益的总headline；Gain 只写 audit/extent。

### Threshold reliability

阈值可靠性在 mask 之前评估：有效 bootstrap 至少 80%，且 CI width 不超过有效
Recon Neff 观测范围的 25%。这些是既有阈值可靠性规则，不是重建效果好坏的门槛。阈值状态和 `threshold_reason` 进入结论记录，同时链接
原始 threshold/extent 表；科学 CSV 的既有字段不能被报告重绘静默改名。

<a id="region-extent"></a>

### Region extent

**Meaning.** State area/fraction 表示当前 route、parent、主尺度和有效窗口下的
可识别高多样性覆盖；它不是跨平台 effect size。Gain 说明局部增益分布的 audit，
不是另一个 headline。

**Direction and conditions.** 阈值可靠时才可展示 mask 与 area。无稳定 breakpoint、
有效窗口不足、bootstrap 有效率不足或 CI 过宽时状态为 `no_stable_threshold`，
Region extent 的 area、area fraction、unit fraction 和 `region_windows` 写 `NA`/`null`，
而不是 0；`valid_windows` 和 `valid_units` 等支持计数仍保留。overview 只显示有效窗口空间中的 `area_fraction`；面积、有效窗口和 unit 分母在正文小表和具名数据链接中展开，并链接完整
bootstrap/extent 表。State 与 Gain 的 threshold 不能互换。

<a id="scale-sensitivity"></a>

## 8. Scale sensitivity

**Purpose.** 检查局部 diversity 方向和摘要是否依赖 window scale，为主尺度结论
提供条件性证据。

**Computation.** 对每个候选 scale 重复同一 paired three-assignment rarefaction、
seed、draws、baseline edges 和有效窗口条件，保存每个 scale 的 support、Neff/delta
摘要和主尺度标记。

**Meaning.** 曲线可以说明某个方向在若干尺度上保持或改变；它只支持主 Neff/delta
的 scale sensitivity，不验证 State/Gain mask 跨尺度不变。

**Direction and conditions.** 只有各尺度 support 和 baseline 语义一致时才可比较；
低支持尺度保留状态。主尺度以预先定义的 support rule 选择，不读取 cluster label、
diversity 或 Region 结果反向调参。

## 9. Notebook 块、产物和 brief batch observation

Notebook 保留目的与方法链接、短调用、必要检查及既有科学图；正式批次观察读取本次records，
不在绘图cell中重复生成判读。长实现使用 [notebook_helpers.py](../../../reproduce/case/reconstruction_impact/notebook_helpers.py)
和 [notebook_analysis.py](../../../reproduce/case/reconstruction_impact/notebook_analysis.py)；
不要把长循环和绘图实现重新塞入 cell。推荐块如下：

| 块 | staged short calls | 必须保存 | 图与观察 |
| --- | --- | --- | --- |
| Partition complexity | prepare Raw features → run resolution diagnostic → save summary/contingency | K、size、ARI、resolution、absolute/normalized denominator、status | complexity/cluster/contingency 图；一句说明复杂度方向及限制 |
| Matched-K | scan target K → choose once → Hungarian once → save assignments | K、mapping、unit change、Wilson denominator、unmatched state | mapping/change 图；一句说明 assignment difference |
| Raw Level2 | load configured reference → map labels/posterior → register baseline | method/filter/count、labels/posterior、reference digest | baseline audit；说明 reference-derived boundary |
| Anatomy/windows | build full context → choose scale → register windows/support | coordinates、anatomy、scale curve、window denominator | anatomy/support 图；说明 cohort/scale conditions |
| Local diversity | run three assignments on same draws → summarize metrics | Kobs/Neff/evenness、draw count/seed、both baselines、anatomy | 2×3 matrix and spatial fields；一句 baseline-specific observation |
| Region | bootstrap threshold → plot reliability → emit extent/mask | threshold state, CI width, valid bootstrap, valid/region windows and units, area/fraction | reliability before State mask；一句 stability/NA observation |

输出文件由 [outputs-and-test-plan.md](outputs-and-test-plan.md) 定义。当前科学
函数参考为 [reconstruction_impact.py](../../../revise/analysis/reconstruction_impact.py)、
[partition_change.py](../../../revise/analysis/basic/partition_change.py) 和
[spatial_region.py](../../../revise/analysis/basic/spatial_region.py)；这些链接是实现
入口，不改变本页的方法定义。

## 10. 解释边界和验收

必须验证 Raw QC/HVG 由 Raw 定义、paired ID/坐标严格一致、same-resolution 与
matched-K 不混淆、Hungarian 只执行一次、Wilson denominator 正确、三个 assignment
共享 draws、anatomy 不参与调参、State/Gain threshold 可重放。结构上输入审计和
overall 结果必须在 localization/Region 之前，threshold reliability 必须在 mask
之前；HTML/Notebook 不通过重排标题来改变计算顺序。

负 delta、低支持、unmatched cluster、`no_stable_threshold` 和 NA 都是可复核结果。
报告只能消费保存表和记录；renderer 不重新读取 H5AD 或重算指标。
