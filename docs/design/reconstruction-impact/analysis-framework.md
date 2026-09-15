# Reconstruction-impact 分析框架与逐项确认

本页负责分析目的、完整内容树、节点职责、确认进度与实现差距；是这些内容的唯一归属。
总入口仍为 [README](README.md)。方法细则留在对应方法文档，存储和呈现契约留在
[输出契约](outputs-and-test-plan.md)，不在这里复制公式、CSV schema 或另一份报告。

## 1. 目的与证据逻辑

最终评估：**重建相较 Raw 在哪些方面获得更好效果，证据是否充分，结论在什么条件下成立。**
描述变化是证据步骤，不是分析终点；评估改善不等于预设每项结果都更好。

```text
预期收益 → 比较问题 → 科学证据 → 条件下的判断 → 空间解释
```

- 区分改善证据、变化诊断、解释背景、可靠性条件；不把所有指标规定为越高或越低越好。
- 保留支持改善、未支持改善、相反证据和无法判断的情况。这是讨论语义，尚不是新增 record enum。
- range、窗口及 high-diversity Region 都须说明如何直接或间接服务于 impact 判断。
  单独识别出 Region 不构成改善证明；需结合对应 Raw baseline、局部差异与可靠性解释。
- 不新增未经确认的分析、门槛或综合分数。现有数值和科学定义不能因改写规范而被静默替换。
- 主树共用，HD/Xenium 差异在节点内说明；固定 parent 顺序为 Fibroblast、Mono/Macro、T。
  HD global 使用独立 cohort 与分母，Xenium 不增设 global。
- 本轮组织方案的主文 overview 固定五行：分群、基因信息 + Moran、EMT、双 baseline、State。
  载体、范围/配对、完整基因和预处理属于按需审计；它们仍有节点编号和方法链接，但不要求在主文显示。
- 27 个节点是科学追溯注册表，不是 HTML 左侧目录或 27 个必须渲染的模块。一个 overview
  行可以引用多个节点，一个节点也可以只在审计或证据索引中出现；HTML 的私有阅读顺序另行固定。
- 下列顺序是阅读顺序，不是全串行依赖：matched-K 受限不自动使 Moran、EMT 或 Region 无效。
  上游限制只传递到实际依赖该条件的判断。

<a id="content-tree"></a>
## 2. 共用内容树（27 节点追溯注册表）

下表和后面的节点表保留完整 27 节点、方法锚点和旧 question 映射。它们用于把一条
科学结论追溯到输入条件、计算和记录，不规定 HTML 必须逐项展开。基础节点 1.1–1.6
继续作为按需 audit；HTML 左侧目录及长页顺序见 [README 的 HTML 私有阅读树](README.md#html-私有阅读树与长页)，
不由本批发现、记录数量或图是否存在决定，也不为了凑齐节点而编造文字、表或图。

```text
Reconstruction-impact scientific traceability registry (27 nodes; not HTML navigation)
├─ Foundation audit (on demand; 1.1–1.6)
│  ├─ 1.1 载体与路线语义
│  ├─ 1.2 比较范围与细胞类型
│  ├─ 1.3 配对与空间坐标
│  ├─ 1.4 表达视图与完整基因空间
│  ├─ 1.5 分析专用预处理
│  └─ 1.6 基础条件收尾
├─ Overall evidence nodes (2.1–2.6)
│  ├─ 2.1–2.3 分群结构、K 与 assignment 诊断
│  ├─ 2.4 基因信息与 Moran
│  └─ 2.5–2.6 EMT 与总体条件收尾
├─ Local-neighborhood evidence nodes (3.4–3.8)
│  ├─ 3.4 窗口支持、共同 draws 与尺度
│  └─ 3.5–3.8 双 baseline 的 Kobs、Neff、evenness 与 anatomy 分层
├─ Spatial evidence nodes (3.1–3.3, 3.9, 4.1–4.5)
│  ├─ 3.1–3.3 anatomy、分群/assignment、EMT 空间场
│  ├─ 3.9 空间定位与局部状态收尾
│  └─ 4.1–4.5 State 可靠性、连续场/mask、覆盖、尺度与 Gain audit
└─ summary 跨 parent 综合解释（不合成 improvement score）
```

**明确边界：** 2.3 复用已有结果，不重新计算一套 cell-type change。4.4 当前检查局部指标
对尺度的敏感性，不能称为 Region 边界跨尺度稳定；它与 Region 的可识别性和覆盖分开记录。
科学注册表的节点分组不改变旧 question、方法定义或记录粒度；HTML 只按固定阅读树把这些
节点归入总体、局部邻域和空间的页面位置。

## 3. 节点职责与方法归属

下表是组织骨架中的问题边界，不是对现有结果已经支持改善的判定。各节点确认及自主对齐状态见第5节；既有方法继续有效，未验收的实现另行登记。

| 节点 | 接收什么 → 回答什么 → 向后提供什么 | 详细规范归属 |
| --- | --- | --- |
| 开篇 | 已审阅记录 → 本批整体结果与条件 → 正文导航；不提前推断 | [输出契约](outputs-and-test-plan.md) |
| <a id="node-1-1"></a>1.1 | 输入与重建来源 → 比较载体的真实语义 → 各方法可解释的对象 | [输入：载体](input-views.md#comparison-foundation) |
| <a id="node-1-2"></a>1.2 | Raw 标签与覆盖 → 哪些对象纳入各 scope → 明确 cohort 与抽样边界 | [输入](input-views.md) |
| <a id="node-1-3"></a>1.3 | ID/坐标与 cohort → 是否支持配对和空间比较 → 配对轴、单位及分母 | [输入](input-views.md) |
| <a id="node-1-4"></a>1.4 | 表达来源与基因轴 → 测量/重建与可计算范围 → 完整基因及缺失语义 | [输入](input-views.md)、[基因方法](gene-and-function.md) |
| <a id="node-1-5"></a>1.5 | 各分析视图 → 专用预处理是否合适 → Partition/Moran/EMT 各自条件 | [输入](input-views.md)、[impact](impact-analysis.md)、[基因方法](gene-and-function.md) |
| <a id="node-1-6"></a>1.6 | 上述审计 → 哪些比较成立/受限 → 逐分析传递条件 | [输出契约](outputs-and-test-plan.md) |
| <a id="node-2-1"></a>2.1 | Partition 结果 → 第一组 K/ARI 如何变 → 控制 K 比较的诊断背景 | [复杂度](impact-analysis.md#partition-complexity) |
| <a id="node-2-2"></a>2.2 | 第二组 matched-K 状态及 assignment → K、ARI、unit_change 与 balanced 如何变 → 一次 mapping 与限制 | [Matched-K](impact-analysis.md#matched-k) |
| <a id="node-2-3"></a>2.3 | 各 scope 结果 → 变化涉及哪些 cell type → 不混分母的总体解释 | [Matched-K](impact-analysis.md#matched-k) |
| <a id="node-2-4"></a>2.4 | Raw/Recon 的 provided_n、完整表达与同组同图 → full-side 与 shared-valid 两种基因范围的空间自相关变化 → 基因维度证据及限制 | [Moran](gene-and-function.md#moran) |
| <a id="node-2-5"></a>2.5 | 同资源及两侧表达 → 覆盖与活动分别如何变 → EMT 总体证据及空间 score | [Coverage](gene-and-function.md#emt-coverage)、[score](gene-and-function.md#emt-score) |
| <a id="node-2-6"></a>2.6 | 结构/Moran/EMT 结果 → 各维度支持什么 → 待定位现象；不合总分 | [输出契约](outputs-and-test-plan.md) |
| <a id="node-3-1"></a>3.1 | 完整 Raw context → 位置属于什么背景 → 固定 anatomy 解释上下文 | [Anatomy](impact-analysis.md#anatomy) |
| <a id="node-3-2"></a>3.2 | 已保存 mapping/change → 归属变化在哪里 → 带控制条件的定位证据 | [Changed units](impact-analysis.md#changed-units) |
| <a id="node-3-3"></a>3.3 | 已保存 EMT score/delta → EMT 变化在哪里 → 功能空间证据 | [EMT fields](gene-and-function.md#emt-spatial-fields) |
| <a id="node-3-4"></a>3.4 | 三 assignment/坐标 → baseline 和窗口是否支持比较 → 同窗口同 draws 的局部基础 | [窗口](impact-analysis.md#window-support) |
| <a id="node-3-5"></a>3.5 | 共同窗口/draws → 状态数量如何变 → 双 baseline richness 证据 | [局部多样性](impact-analysis.md#local-diversity) |
| <a id="node-3-6"></a>3.6 | 共同窗口/draws → 有效状态多样性如何变 → 双 baseline Neff 证据 | [局部多样性](impact-analysis.md#local-diversity) |
| <a id="node-3-7"></a>3.7 | 同窗口状态组成 → 均衡程度如何变 → 辅助解释 Kobs/Neff | [局部多样性](impact-analysis.md#local-diversity) |
| <a id="node-3-8"></a>3.8 | 局部指标与 anatomy → 哪些背景出现哪些差异 → 分层证据与分母 | [局部多样性](impact-analysis.md#local-diversity) |
| <a id="node-3-9"></a>3.9 | 定位和局部证据 → 效果或差异发生在哪里 → Region 解释基础 | [输出契约](outputs-and-test-plan.md) |
| <a id="node-4-1"></a>4.1 | Neff/窗口/阈值审计 → 能否稳定识别 → 阈值状态与适用条件 | [可靠性](impact-analysis.md#region-reliability) |
| <a id="node-4-2"></a>4.2 | 连续场及阈值状态 → 状态与可识别区域在哪 → 连续场、mask 或不可识别 | [Region](impact-analysis.md#region-extent) |
| <a id="node-4-3"></a>4.3 | Region/anatomy/有效范围 → 覆盖多少 → 不混淆窗口与 unit 的覆盖证据 | [覆盖](impact-analysis.md#region-extent) |
| <a id="node-4-4"></a>4.4 | 已有候选尺度结果 → 局部结论是否依赖尺度 → 敏感性限制 | [尺度](impact-analysis.md#scale-sensitivity) |
| <a id="node-4-5"></a>4.5 | 可靠性、状态、覆盖 → Region 如何服务 impact → 不夸大的区域解释 | [输出契约](outputs-and-test-plan.md) |
| <a id="node-summary"></a>结尾 | 五行证据及限制 → 哪里支持改善/不支持/无法判断 → 可追溯综合解释 | [输出契约](outputs-and-test-plan.md) |

<a id="node-writing"></a>
## 4. 每个节点的填写与确认结构

**统一说明段（规范要求）。** 每个固定阅读位置在图表或方法细节之前，用一段连贯文字说明：
本节在重建前后评估中的逻辑位置、为什么采用这些比较/组织方式、希望帮助读者作出什么判断。
预期效果指分析希望解决的问题，不预写本批结果为正向。该段连接前后节点，不能只是重复标题。
完整说明在方法文档维护；网页在固定问题位置显示简短结果、表/图和解读，Notebook 保留
简洁英文目的及方法链接。文字模板与填写职责见[输出契约](outputs-and-test-plan.md#narrative-templates)。
科学节点编号继续服务追溯；HTML 的固定左侧目录和单页长页顺序由 README/输出契约维护，
按“总体 → 局部邻域 → 空间”阅读。固定问题标题是“分群”“基因信息 + Moran”“EMT”“双
baseline”“State”，三列 parent 始终为 `Fibroblast`、`Mono/Macro`、`T`。这棵展示树不
随本批发现、record 数量或图可用性改变，缺失事实只显示 pending/NA。

**选择方式（已确认）。** 存在实质取舍时，先用流程图、树、并列表或具名图表示例展示差异，
明确真实数据与模拟内容，再给选项卡并说明优劣。无需选择的事实和常规衔接直接推进，
不询问是否继续；每次确认后先落盘再展开下一事项。

**展示取舍（规范要求）。** 先判断是否服务于实际科学对象及重建收益判断，再决定是否作图。
细胞类型、基因、EMT、局部状态及 Region 的比较可以使用有信息价值的图示；ID 等技术
配对检查默认放在按需 audit 中，不为每个节点机械添加教学图，也不把常规排版细节
反复拆成选择题。对比示例用于帮助决定有实质影响的内容或判读方案。
默认采用能支持判断的最简形式；仅仅能画出图、能细分集合，不构成增加展示的理由。
网页正文按“固定问题标题 → 简短结果 → 表/图 → 1–3 条解读”组织，基础 audit 不因节点注册
而强制显示。保存的 `node_narratives.next` 继续保留在结果包并供 Notebook 使用；网页不沿用
旧的 `next` 顺序组织 module 或导航，网页导航服从固定左侧目录。

“全方位”要求已约定的结构、基因、EMT、空间局部与 Region 问题及支持条件形成完整证据链，
不自动授权增加分析类别。

```text
节点（需要生成 narrative 时）
├─ 结论标题：节点双语标题，指明 scope 或 baseline
├─ 本批 lead：只引用本批已保存事实，说明该节点在五行中的位置
├─ 小表：数字、单位、分母、条件、状态和具名 evidence link
├─ 图：只有能改变判断的已登记图，写清 scope/单位/方向
├─ 1–3 条 interpretations：事实 → 有条件的判断 → 限制（双语同源）
├─ next：保存在 narrative/Notebook，网页不以它排序或导航
└─ record_ids/review：绑定来源记录与审阅身份；缺证据就保持 pending
```

先读取真实证据，按已确认原则自主确定常规组织；仅有实质影响且无法依原则确定的选择，
才提供具体对比示例并通过选项卡讨论。
方法细则写回方法文档；字段/结论填写规则写回输出契约；本页仅保存确认摘要与链接。
确认后立即落盘再进入下一节点。数据不足时先记缺口，是否新增分析另行确认。

## 5. 规范决策与状态边界

**协作与决策授权。** 用户授权主Agent先依既有原则判断其余节点，不逐项询问常规组织。
子代理仅使用内部collaboration通道或子任务最终返回；禁止使用跨任务对话消息工具
（包括`codex_app.send_message_to_thread`）向主任务传递结果。主Agent整合后统一汇报。
“依原则确定”指本轮文档组织决定，不冒充用户逐项确认或实现已验收。


已批准的阅读安排与实际科学结果分别记录；本轮实现及检查状态见
[规范入口](README.md#report-display-status)。

| 对象 | 本轮组织方案 | 完成后应（检查条件） |
| --- | --- | --- |
| 目的原则、路线差异和五行 overview | 固定标题为分群、基因信息 + Moran、EMT、双 baseline、State；parent 顺序固定 | 固定标题和 parent 顺序与正式页面一致 |
| 1.1–1.5 基础定义、范围/配对、基因、预处理 | 收在固定目录的按需 audit | 可从固定入口展开，方法锚点、分母和缺失语义可追溯；不占五行 |
| 1.6 基础条件收尾 | 收在按需 audit，只引用已有记录 | 没有记录时保持 pending；不得为填导航编造结果 |
| 2.1–4.5 27 节点追溯 | 保留科学注册表，按 HTML 私有树归入总体、局部邻域或空间 | 27 节点、方法链接和记录关联可追溯；不要求逐节点渲染 |
| 文字模板与 `node_narratives` | 作为同源内容契约继续保留 | 页面与 Notebook 可引用同源事实，但是否已填充、审阅需按实际结果核对 |
| 旧 16 类 question | 保持兼容 | 不改 question ID、记录粒度或科学定义 |

**本轮组织方案的中心报告规则。** 1.1–1.5 的基础定义、范围/配对、完整基因和预处理
仍是方法与审计依据，但从主文移到固定目录的按需入口；overview 只保留分群、基因信息 +
Moran、EMT、双 baseline、State 五行。27 节点和旧 16 类 question 都继续用于追溯，不能
为了五行而删除 record、method anchor 或计算定义。真实页面仍需核对固定左侧目录、长页顺序和
每行内容；当前浏览器检查尚未确认用户阅读效果，不能写成完整验收通过。

可见节点使用输出契约中的文字模板；只有存在本批事实时才填入 lead、表、图和
interpretations。缺失事实写明 pending/状态，不用空 shell 凑齐节点。

**保留的证据差距。** 当前foundation `gene_availability.csv`按gene并集登记
`provided`与表达来源，代码依据`var_names`决定是否提供。该检查不能单独证明
缺失原因是平台未测量，也不能证明被提供的基因在当前scope非零或可计算Moran。
这项审计语义差距留待对齐，不因选择简洁展示而删除，也不自动更改现有结果状态。
1.5/1.6采用[预处理短表及条件收尾](input-views.md#preprocessing-presentation)。
其他节点按以下表自主对齐，不为已批准的方法和常规展示重复发起选择。

<a id="remaining-node-alignment"></a>
### 固定阅读树中的节点职责

以下表索引每组科学节点在固定 HTML 阅读树中的逻辑、填写职责与写法；科学计算和详细判读
仍以对应方法章节为唯一归属。页面位置不由发现、record 数量或图是否存在决定，不要求每个
节点在主文中可见，也不要求每个节点强行生成一张图、一张新表或一种新结论记录。

| 节点 | 为什么看、希望判断什么 | 保留的核心展示 | 本批结论与衔接 |
| --- | --- | --- | --- |
| 1.5–1.6 | 分清各分析条件，避免误用共同预处理或预先宣称支持 | 按需 audit 的参数、状态和数据入口 | 有记录才写基础条件；进入五行结果，不在主文预填 |
| 2.1–2.3 | 分清两组 K、成员改变与 ARI/balanced 的分群诊断，定位涉及的 cell type | 网页 overview 只预览两组 K 与成员改变；正文保留完整 K、ARI、unit_change、balanced 表；HD global 独立 | 两组按已有 record 分开填；`balanced = 1 - macroF1`，缺字段保留 NA/状态；变化量不等于改善量 |
| 2.4 | 分别检查已有 provided_n、完整有效基因的空间结构与共同有效基因的对应变化 | 原有 provided_n、full/shared 分布、Q75、有效数、逐基因对应与 delta | provided_n 来自 gene availability/foundation；两种 Moran 视角分别给出方向和幅度；不把自相关升高直接判为生物学恢复 |
| 2.5 | 区分EMT可供评分的资源范围与实际评分变化 | coverage、cutoff、score 分布和配对摘要 | coverage与score各有解释，保留表达来源条件；空间场按需引用 |
| 2.6 | 把已观察到的结构、基因、EMT现象带到空间层 | 短段引用前述结论，不复制完整表或合总分 | 明确接下来要在空间检查的现象 |
| 3.1–3.3 | 在空间顺序中判断已有变化出现在何处及何种组织背景 | anatomy、分群图/assignment、EMT spatial field 的已登记证据 | 复用同一 mapping/score；按固定顺序展开，不另算一种空间总效果 |
| 3.4 | 确定局部邻域比较的窗口支持和尺度 | 既有 parent 及 anatomy 窗口支持图、所选尺度与 baseline 摘要 | 说明同窗口同 draws、两条 Raw baseline 的各自角色 |
| 3.5–3.8 | 并列检验相对 Raw Leiden 与 Raw Level2 的三项局部指标 | Kobs/Neff/evenness 的完整共享 2×3 图组；anatomy 分层摘要 | 两 baseline 并列，保留窗口支持/尺度；Neff 方向依既有条件解释，Kobs/evenness 辅助 |
| 3.9 | 归纳变化位置及局部状态证据 | 短段引用既有空间证据 | 为Region解释建立上下文，不重复整套矩阵 |
| 4.1–4.3 | 在空间末段先确定 State 是否可识别，再说明连续场、mask 与覆盖 | 三 parent 可靠性比较、连续 Neff 场、State mask/extent 表 | 保留已有阈值规则；不稳定时覆盖 NA；识别到 State 不等于证明改善 |
| 4.4–4.5 | 说明尺度条件和 State 能支持的解释，给 Gain 独立 audit 入口 | 已有同定义 rarefaction 敏感性图/表、Gain audit | 不冒称 State 边界稳定；Gain 只保留 audit，无新主图 |
| 五行 overview 与结尾 | 按固定标题预览并归纳哪些方面支持收益、仅见变化或证据不足 | 五行摘要；结尾简短综合解释及证据入口 | 基于同一记录逐项对应结果及条件；无新增门槛或综合分 |

报告沿用科学问题与 cell type 的证据；Notebook 保留短调用、核心检查与既有主图。
共用图只完整展示一次；专属空间图按 parent 归属。基础 audit 没有科学图时保留
状态和数据入口，不补图。记录、五行总览、正文顺序和文字模板见[输出文档](outputs-and-test-plan.md)。

## 6. 实现对照与当前检查

科学节点顺序和目的以 `content_contract.py` 为程序侧唯一映射；正文锚点使用
`node-<编号>`，Notebook 使用同名 `impact_node`。HTML 左侧目录和长页顺序是独立的私有
展示契约，按 README 的总体 → 局部邻域 → 空间组织，不从科学树、文件目录或本批发现推导。
旧 `module-<question>` 锚点保留为兼容入口。
`content_contract.py` 保留 27 个科学节点、16 类 question 和方法追溯；网页阅读顺序
由 renderer 的轻量映射维护。已保存的 `node_narratives.next` 留在结果包与 Notebook，
网页不沿用旧的排序指引。两份正式报告及检查状态见[规范入口](README.md#report-display-status)。

HD All 的成员改变按 Raw Level1 类型分层，保留独立分母与 Wilson 区间；其侧别 `provided_n`
未保存时明确显示缺失，完整载体的基因范围单独命名。既有 Gain audit 只读展示。

## 7. 推进与验收

- 每轮检查节点如何服务重建收益、证据是否够、复杂度是否值得、是否误把诊断当改善、是否重复规范。
- 常规排版、节点挂载和页面细节按既定阅读树自主推进，不逐模块发起审批。
- 仅在排除或新增分析，或改变比较含义/比较条件时单独讨论；其余节点沿用既有定义和记录契约。
- 本轮不修改算法、cohort、baseline、阈值或科学结果定义；实现与展示检查已完成，用户阅读效果尚未确认。
- 展示修改按实际受影响的数值、分母、图和导航检查；科学执行身份变化才重新检查其依赖。
