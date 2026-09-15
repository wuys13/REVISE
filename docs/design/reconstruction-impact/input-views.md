# 输入视图与表达载体

框架节点：[1.1 载体](analysis-framework.md#node-1-1)、[1.2 范围](analysis-framework.md#node-1-2)、[1.3 配对](analysis-framework.md#node-1-3)、[1.4 表达](analysis-framework.md#node-1-4)、[1.5 预处理](analysis-framework.md#node-1-5)、[1.6 条件总结](analysis-framework.md#node-1-6)。

下文保留既有方法与约束，作为报告五行入口的按需审计依据。这里的定义、计算和状态
仍然有效；本页不把它们强制放进主文，也不把已有执行记录写成本轮验收。

本页负责 foundation layer 的输入方法。它定义 Raw、spatial、reference/donor
载体、cohort、坐标、route 和全基因可用性；不定义 partition、Moran 或 Region
的科学指标。总览和作者归属见 [README](README.md)，保存字段和状态见
[outputs-and-test-plan.md](outputs-and-test-plan.md)。

每个方法段都用 Purpose → Computation → Meaning → Direction and conditions。这里
的“视图”是带有来源、轴和映射说明的只读分析载体，不是要新增的 H5AD；cluster
mean 可以在内存构造，但不得修改 Raw、reference 或已发布 reconstruction。

<a id="comparison-foundation"></a>

## 1. Carrier identity and provenance

<a id="carrier-presentation"></a>
### 1.1 内容组织（按需审计）

载体和路线语义从主文移到按需 audit。读者需要核对来源时，使用
**目的说明 → 表达来源流程 → 所需字段与具名证据入口**；不规定主文必须出现的表格
数量、列数或固定图。空间范围和纳入/排除数量由 1.2 audit 负责，基因范围和可用性
由 1.4 audit 负责。这里的角色字段只是审计字段定义，不能被解释为主文的基础表要求。

**本节逻辑与预期作用。** 重建后的统计差异只有在表达来源和观察单位明确时才能解释。
因此先展示表达形成过程，再列载体角色：HD直接消费空间重建表达；Xenium先由表达侧
reference cells求cluster均值，再按空间单位的cluster映射回空间轴。这样的组织让读者
在看到Moran或EMT变化之前，知道比较针对什么表达对象，避免把投影字段误认为实测表达，
或把reference cell数量当作空间样本分母。本节提供后续评估的语义基础，本身不证明改善。

```text
HD
Raw 空间表达 ──────────────────┐
                              ├─ 后续 Raw / Recon 比较
Recon 空间重建表达 ────────────┘

Xenium
Raw 空间表达 ──────────────────────────────────┐
                                              ├─ 后续 Raw / Recon 比较
表达侧 reference cells → 按 cluster 求均值 ─┐   │
空间 units 的 cluster ──────────────────────┴─→ 投影表达
```

流程说明表达来源，不替代 1.2 的 cohort 选择或 1.3 的 ID/坐标配对检查。
Raw Level2 reference 在后续角色表中单独列出，说明其对第二 baseline 的标签映射作用；
不与 Xenium 表达侧输出混为一个角色。

若审计需要结构化查看，可按载体、观察单位、提供内容、后续用途、解释边界记录下表
字段。下表是字段示例，不是主文固定表，也不表示本批结论或支持条件已经验收。

**Purpose.** 说明每个文件承担的角色，并把 spatial observation 轴与表达侧轴分开，
使后续的 paired、partition、空间图和功能评分不会误用 carrier。

**Computation.** foundation 首先登记 `sample_id`、`task_name`、`route`、
`source_role`、source path/digest、observation ID 轴、gene ID 轴、Level1 列名、
坐标来源/单位、表达层、映射方式和未覆盖数量。至少区分：

| 载体 | 观察单位 | 提供内容 | 后续用途 | 解释边界 |
| --- | --- | --- | --- | --- |
| Raw（两路线） | 原始空间 unit | 测量表达、空间ID、坐标和Raw标签 | 定义比较对象、Raw baseline与完整组织背景 | 基因未测量不等于零；具体范围与配对另行审计 |
| HD Recon | 重建空间 unit | 空间重建表达及其ID/坐标 | 提供Recon表达与空间比较输入 | 重建估计不等于真实测量；配对资格由1.3检查 |
| Xenium spatial carrier | 空间 unit | 坐标、ID、SVC_cluster | 空间配对、分群标签和局部窗口；接收投影表达 | 不能把表达侧reference cells当成这些空间units |
| Xenium expression carrier | 表达侧reference cell | 按SVC_cluster分组的表达 | 求cluster mean后映射回空间unit，供Moran/EMT比较 | 未投影的表达侧轴不能进入空间分母；不是逐空间unit实测值 |
| Xenium projected expression | 映射后的空间 unit | 对应cluster的均值表达 | Recon侧Moran与EMT输入 | 同cluster空间units共享表达；该视图由内存构造，不要求新增H5AD |
| Raw Level2 reference（两路线） | 参考单细胞 | 参考表达及Level2标签 | 对Raw进行Level2映射，提供第二baseline | 与Xenium expression carrier是不同角色，不混为同一参考用途 |

**Meaning.** `spatial carrier` 与 `reference/donor carrier` 即使来自同一次
reconstruction，也不因此拥有相同 observation 轴。Xenium `expr.h5ad` 只能进入
表达侧分析；它不进入 spatial pairing、Leiden、坐标、window 或 anatomy 分母。
Raw Level2 reference 与 iST 的 `expr.h5ad` 是两个独立角色。

**Direction and conditions.** 视图只有在轴唯一、非空、来源可追溯时才可标记
`available`。缺少 digest、重复轴、未知 role 或错误表达层是 `blocked`/执行错误；
不能把这些问题写成数据低 coverage。任何重新排序都写入 audit，不能依赖文件行号。

## 2. Cohort、ID 和坐标配对

<a id="cohort-presentation"></a>
### 1.2 内容组织（按需审计）

范围和 cohort 从主文移到按需 audit。需要追溯时，按
**目的说明 → 选择流程 → 需要的计数/分母 → 具名证据入口**读取；不规定主文必须
出现来源图、选择范围表或分母索引，也不改变实际 cohort 或 QC。对象选择和各分析
分母仍需分开记录，防止把并行分析误读为累积筛选流程。

**本节逻辑与预期作用。** 重建前后的差异也可能受到比较对象范围不同的影响。因此先说明
Raw标签与重建覆盖如何决定候选对象，再展示各parent及分析的实际范围。这样既能追溯
纳入和排除过程，又能横向识别不同结果使用的分母，避免将Partition专用QC范围套到
Moran/EMT，或将HD global与parent结果累计。本节建立效果评估的范围条件，不以数量多寡判改善。

```text
Raw 标签＋重建覆盖
├─ 各 parent 独立选取 → 分析输入
│                     ├─ Partition 专用 QC
│                     └─ Moran/EMT 各自使用范围
└─ HD global 独立选取 → 独立记录和分母（Xenium无此分支）
```

完整 Raw 的 anatomy 背景另行登记，不从 Partition QC 后队列推出。审计记录按
Fibroblast、Mono/Macro、T 并列，HD global 单独记录；按路线实际步骤填写，不为 Xenium
复制不存在的 HD 筛选。数量相同不证明 ID 相同，配对证据由 1.3 负责。

审计所需的范围记录应能区分来源范围、重建覆盖、候选对象、实际纳入对象、排除原因和
各分析分母；每个计数注明对象、统计范围与分母。Partition QC 后数量不作为 Moran 或
EMT 的共同前置条件，输入数量与实际有效/配对数量按分析语义分别记录。基因数和窗口
数不混入 spatial unit 计数。
基因数与窗口数不混入空间unit计数；其有效性在对应节点展开。

审计填写遵循“实际范围 → 对本项判断的限定 → 后续条件入口”。例如 Partition 写明
QC 后对象范围；Moran 另外指向有效基因条件；EMT 另外指向资源与评分条件。该审计
只说明范围条件，不提前判定后续结果有效或重建获益。

HD Fibroblast 的 30,000 个纳入 unit、23,151 个 Partition QC 保留 unit 以及 Moran/EMT
各自使用的范围，可作为记录示例；具体批次仍以保存的 observations 和对应分析分母为准。


<a id="pairing-presentation"></a>
### 1.3 内容组织（按需审计）

配对和坐标从主文移到按需 audit。需要追溯时，按
**目的说明 → 必要核查字段/状态 → 本批证据入口**读取；不要求主文渲染全量核查表，
也不增加 ID 教学示意。

**本节逻辑与预期作用。** ID与坐标是确认同对象、同位置比较成立的技术条件。
用简洁核查表说明实际检查范围、结果与证据，使读者知道后续配对差异是否有可靠对应关系；
不将技术检查扩展成独立科学展示，也不以配对通过证明重建改善。

审计字段覆盖 ID 唯一及对应关系、按 ID 对齐、坐标来源与单位、有限性及配对坐标一致性，
并注明实际范围与证据入口。本批解释说明支持哪些比较、还有哪些条件未获证据支持；
缺少审计证据时保留 pending/blocked，不从数量相等推断配对成立。

**Purpose.** 规定哪些 Raw spatial units 可以进入一个 scope，避免用 reconstruction
label 反向定义 cohort，也避免把 donor cells 当成 spatial units。

**Computation.** 对声明为 paired 的任务按以下顺序：

1. 验证 reconstruction observation IDs 是 Raw IDs 的子集；
2. 用 ID 显式对齐，不依赖行号；
3. 验证保留 ID 都有有限二维坐标，并在同一坐标系/单位下相等；
4. 在覆盖范围内选择 global 或 parent cohort，再应用 Raw-defined QC；
5. 记录纳入、排除和未覆盖 ID 的 reason、numerator/denominator、坐标转换和
   `observation_basis`。

Global 与 parent cohort 独立选择；parent 选择使用 full Raw 的 exact `Level1`，
不受 global 抽样子集限制。VisiumHD 的当前上限是每个 scope 30,000 个覆盖范围
内的 ID；`USE_FULL_VISIUMHD_COHORT=True` 只取消这个上限，仍执行 Raw QC。配置和
实际保留 ID 由 [VisiumHD config](../../../configs/analysis/reconstruction_impact_visiumhd_p1crc.yaml)
和结果根中的 observations/assignment 表共同证明。

**Meaning.** `n_units`、`n_changed`、valid windows 和 anatomy 分母只在声明的
scope 内解释。它们不能混用 full tissue、paired parent 和 expression-side donor
cell count。

**Direction and conditions.** ID 缺失、坐标不一致、空标签、重复轴或重建输出包含
Raw 不认识的 unit 时，视图不能标记 `paired`；保留 `unavailable`/`blocked` 与
原因，不能静默取交集后继续。parent-internal change 只能描述该 parent 内的
assignment difference。

### Exact task label mapping

安全 task 名、文件夹名和 biological label 不互相猜测。当前显式映射为：

| task key | Raw exact `Level1` |
| --- | --- |
| `Fibroblast` | `Fibroblast` |
| `Mono_Macro` | `Mono/Macro` |
| `T` | `T` |

`Mono_Macro` 只是安全 ID；内部对 `/` 的 canonical normalization 不能变成 task
discovery 规则。未知 label、近似字符串或自动新增 `Mono` 别名必须失败。实现参考
是 [PARENT_SOURCE](../../../reproduce/case/reconstruction_impact/build_notebooks.py)。

## 3. Route 和表达 projection

**Purpose.** 区分 hST、iST 和 sST 的 observation 与 expression 语义，让 Moran、
AUCell 和局部空间方法只使用它们真正支持的载体。

**Computation.** foundation 保存 route、输入文件和 projection 方式：

| route | spatial 轴 | expression 轴 | 当前允许的比较 |
| --- | --- | --- | --- |
| hST | Raw spatial 与 native reconstructed spatial units | 通常与 spatial units 同轴的 native expression | 只有 handoff 证明 ID/坐标关系时才称 paired |
| iST `paired` | Raw spatial 与 reconstructed `spatial.h5ad` ST units | `expr.h5ad` reference/donor cells，经 cluster mean 投影 | partition/window/Moran 使用 spatial 轴；表达方法使用 projected expression |
| iST `mean` | native `SVC.h5ad` output units | 同一 native carrier | 保留 native mean 语义，不自动声称逐 observation paired |
| iST `random` | native `SVC.h5ad` output units | native donor assignment | 保存 donor ID 和 seed；只有显式 handoff 才进入 paired 比较 |
| sST | generated spatial units | parent spot 关系生成的 baseline/output expression | 首轮只验证 parent-spot 映射与守恒，不宣称 native pairing |

iST `paired` 的 cluster mean 按固定顺序构造：

1. 从 `expr.h5ad` 读取 reference/donor cells 的 `SVC_cluster`；
2. 在实际 `X`/layer 上按 cluster 求 sparse mean，保持 counts、log 或 normalized
   语义并记录；
3. 用完全相同的 key 映射到 `spatial.h5ad` 的 ST units；
4. 标记 `expression_view=cluster_mean_projection`，并保存 donor cluster 数、
   每 cluster donor cell 数、空间覆盖数和实际表达语义。

当前 Notebook 检查每个 retained spatial cluster 都能找到 expression mean，
即 retained spatial keys 是 expression keys 的子集。它没有证明两个完整输入的
cluster 集合完全相同；更强的集合一致要求仍是待讨论事项，本轮不增加该检查，
也不把现有运行标成已满足全等条件。一个 cluster 内的 spatial units 共享
projected expression；这是一种 reference-derived field，
不是 donor cell 的逐细胞真值。`expr.h5ad` observation ID 不能进入 spatial 坐标
或 window 分母。sST 的 `parent_spot_equal_split` 必须验证每个 generated unit 的
Raw counts 和守恒，并保存 `spot_name` 映射。

**Meaning.** projection 说明“该表达字段如何得到”，不增加 measured-gene 或
ground-truth 证据。iST donor 是 reference cell ID，不是患者、实验供体或 sample
ID；记录中同时写 `donor_id_source=reference_cell` 和来源文件。

**Direction and conditions.** 只有 cluster key、轴、来源和表达层均通过 audit 才
可评分或进入 Moran。不能因为 `expr.h5ad` 存在就把它当作 Raw Level2 reference，
也不能给 donor cells 分别评分后再冒充 `cluster_mean_projection`。

## 4. Full-gene availability

<a id="gene-space-presentation"></a>
### 1.4 内容组织（按需审计）

完整基因空间从主文移到按需 audit。需要追溯时，按**目的说明 → 两侧输入摘要 → 状态/限制
与具名数据入口**读取；不要求主文显示基因总数表、集合构成图或正文集合拆分表。逐基因
明细继续保留在已有 availability 数据中。

本节说明后续分析使用什么表达输入，避免将输入范围与分析有效性混为一谈。
审计摘要保留 Raw/Recon、输入提供基因数及表达来源；不同 scope 确有差异时注明 scope。
基因更多不直接判作改善，输入提供也不等于统计量可计算。有效基因数与不可计算原因
在Moran中展开，EMT资源覆盖与评分条件在EMT中展开，不在基础层重复展示。
以下完整方法与状态要求继续由本节维护，简洁展示不等于删除底层审计。


**Purpose.** 让每侧保留自己的完整 gene space，并区分未测量、零表达、低覆盖和
不可计算，为 Moran full/shared 与 EMT coverage/score 提供可重载的能力状态。

**Computation.** 以两侧 gene ID 的并集登记 foundation 能力；对每个分析另外记录：

- Raw 是否提供/测量该 gene；
- 测量存在但全零、低覆盖或不满足计算条件；
- Recon 是否提供该 gene；
- 实际计算 gene 数、排除数、资源版本和 expression view；
- `unmeasured`、`low_coverage`、`insufficient_support`、`not_computable` 或
  `computed` 状态与 reason。

Raw-derived QC/HVG 只属于 partition 的 feature selection，不是 Moran/AUCell 的
共同 panel。Raw 缺少某个 gene 时不补零、不删掉 Recon 的 gene，也不以 overlap
裁剪最终结果。大型矩阵引用原 H5AD 及 digest，不另写 unit × full-gene 巨型 CSV。

**Meaning.** availability 回答“当前输入允许开展什么分析”，不预先把 Raw 判成
失败。full-side Moran 需要各自完整有效基因，shared-valid Moran 才需要共同有效
基因；EMT coverage 与 score 也各自登记两侧条件。

**Direction and conditions.** `NA` 只表示该项在当前条件下没有可用数值；它不是
zero。缺失软件/损坏资源/invalid axes 属于执行失败，不能伪装成 `unmeasured`。
任何一侧不可用仍保留另一侧有效结果和状态；delta 只有两侧有效时生成。

<a id="preprocessing-presentation"></a>
### 1.5–1.6 预处理与比较条件的按需审计

本节把不同分析的输入要求分清，使后续变化可在各自条件下解释；不把全部分析塞进
初始化步骤。主文不强制展示预处理表；需要复核时从对应记录和方法入口查看，不增加教学图
或重复的基础总结矩阵。

| 分析 | 基础层说明什么 | 在对应分析块完成什么 | 限制如何传递 |
| --- | --- | --- | --- |
| Partition | 对应路线的输入、Raw专用QC/HVG与scope | 聚类、复杂度诊断和matched-K检查 | 仅传给依赖assignment的比较；不限制Moran/EMT的完整输入 |
| Moran | 两侧完整基因、同组空间对象/图、normalize/log原则 | 构图、统计量、各侧与共同有效数及原因 | 按本项图、表达及有效基因条件解释 |
| EMT | 固定资源、两侧表达来源、Xenium先projection后评分 | 资源coverage、实际cutoff、评分及配对delta | coverage与score分开解释，投影条件保留 |
| 局部状态与Region | 引用后续窗口、baseline及阈值方法入口 | localization 完成窗口/同draws比较，Region 完成阈值与覆盖 | 不在基础 audit 预先判定窗口或Region可用 |

具体参数分别由[Partition方法](impact-analysis.md)、[Moran/EMT方法](gene-and-function.md)
和[窗口方法](impact-analysis.md#window-support)维护，不在这里另建一套参数规范。

1.6 作为按需 audit 的收尾短段，概括已知输入支持条件及未解决事项，引用 1.2 的分母审计、
1.3 的配对审计与本表。尚未执行的分析写为待检验，不预写可计算或有效；不新增整体
“通过/失败”总评，也不重新抄写前述数值。这里的组织规则不代表正式报告、Notebook 或
审计实现已经按新规范对齐。

## 5. 坐标、单位与 provenance

**Purpose.** 固定空间图、window 和 anatomy 使用的坐标来源，避免从图像尺寸或
数据范围反推尺度。

**Computation.** 优先读取 `adata.obsm['spatial']`；generated sST 只有在 audit 明确
时才读取 `obs[['x','y']]`。坐标必须为有限二维数值，并登记原坐标、转换系数、
转换后单位、origin、source digest 和 provenance：

| route | 原坐标语义 | 当前转换 |
| --- | --- | --- |
| P1CRC VisiumHD | spatial output coordinate | `0.27380817798463214 um / coordinate`，来自分析配置 |
| P2CRC Xenium | morphology pixel | `0.2125 um / pixel`，来自分析配置与 Xenium 空间语义 |

full anatomy context 与 parent window 使用 full Raw Level1 坐标最小值作为 origin，
但各自依据 occupancy 选择 scale；两者不是同一个 square grid。

**Meaning.** 同一坐标系和单位支持空间图、shared Moran graph 和 window 分母；它
不意味着不同 route 的 observation 是 biological matched。

**Direction and conditions.** 坐标转换缺失、非有限、维度错误或来源不明时状态为
`blocked`。任何未来 route 必须提供显式变换和单位，不能仅因两个数组都叫
`spatial` 就直接建图。

## 6. Foundation audit 的 Notebook 入口

基础 Notebook 块只在需要复核时提供短调用、关键参数、中间检查、audit/table 输出和一段
brief batch observation；它们不要求出现在主文。长实现由 [notebook_helpers.py](../../../reproduce/case/reconstruction_impact/notebook_helpers.py)
与 [notebook_analysis.py](../../../reproduce/case/reconstruction_impact/notebook_analysis.py)
承载。推荐的可审阅块如下：

| 块 | staged short calls | audit/table | plot / brief observation |
| --- | --- | --- | --- |
| Route and carrier | load → validate axes → register source/digest | `analysis/inputs/carriers.csv`、`source_files.csv` | route/axis summary；说明哪些载体进入哪些任务 |
| Cohort and coordinates | align IDs → select scope → apply Raw QC → record exclusions | `observations.csv.gz`、partition audit、coordinate/provenance fields | cohort/coordinate coverage；说明 paired basis 与排除分母 |
| Expression views | register gene union → build cluster mean when required → validate cluster map | `gene_availability.csv`、mapping/reference clusters | expression-view summary；说明 projection/原生表达条件 |
| Analysis prerequisites | prepare partition view → register Moran graph/AUCell resources | parameters/audit、graph/resource metadata | 必要的输入或坐标图；说明本批次可计算/受限项 |

上述调用名表示按需审计的可见阶段，不要求创建同名公共 API。实际 source
notebook、执行 notebook 和其 digest 由 [build_notebooks.py](../../../reproduce/case/reconstruction_impact/build_notebooks.py)
及 runner 负责。分析只有在实际依赖某项 foundation 条件时才读取对应审计；不为填充
主文而提前执行所有基础检查。

## 7. 失败与边界

输入视图通过检查后才进入科学计算：轴唯一性、标签列、cluster 集合、坐标形状、
paired ID/坐标关系、gene axis、donor/parent provenance 和 source digest。支持不足
返回 `insufficient_support`；单侧未测量保留 `unmeasured`/`unavailable`；依赖缺失、
资源损坏、轴违反契约和代码错误是 `blocked`/失败。状态、reason 和原始证据必须
写进 [输出契约](outputs-and-test-plan.md) 规定的 manifest 或长表。
