# Reconstruction-impact 规范入口

本目录是 reconstruction-impact 的唯一规范目录。README 提供总入口，各文档按职责维护分析逻辑、输入与
比较条件、保存的事实和阅读报告；不创建另一份 plan，也不把静态旧 CSV 当作新
执行证据。科学定义在方法文档中，Notebook 负责执行调用，HTML 只读取已
保存的记录、表和图。

现有正式科学结果根只有：

```text
output/reconstruction_impact/P1CRC_VisiumHD/
output/reconstruction_impact/P2CRC_Xenium_ThreeParents/
```

每个结果根由 schema 1 配置的 `output.dir`（或已登记的环境覆盖）确定，并包含
`analysis/`、`figures/`、`notebook/` 和同一结果集的 `report.html`。这些结果和已审阅的
中英解读作为只读科学依据；本轮只更新两份 HTML 的阅读组织。
实现与检查状态见[下方说明](#report-display-status)。

结果入口：

| 样本 | 阅读报告 | 执行/阅读证据 |
| --- | --- | --- |
| VisiumHD | [报告](../../../output/reconstruction_impact/P1CRC_VisiumHD/report.html) | [validation.json](../../../output/reconstruction_impact/P1CRC_VisiumHD/notebook/validation.json) |
| Xenium | [报告](../../../output/reconstruction_impact/P2CRC_Xenium_ThreeParents/report.html) | [validation.json](../../../output/reconstruction_impact/P2CRC_Xenium_ThreeParents/notebook/validation.json) |

## 1. 先回答什么

Impact 的最终目的是评估：**重建相较 Raw 在哪些方面获得更好效果，证据是否充分，结论在什么条件下成立。**
描述变化是证据步骤，不预设改善；保留未支持改善、相反证据和无法判断的情况。
range、窗口和 high-diversity Region 必须服务这一证据链，不能单独作为改善证明。

本轮组织方案把报告总览固定为五行：**分群、基因信息 + Moran、EMT、双 baseline、State**。
基础定义、范围/配对、完整基因空间和分析专用预处理属于按需审计入口，不占用总览行，
也不要求在主文正文重复展开。完整的 27 个科学节点仍作为追溯注册表，节点职责和旧
question 映射见 [分析框架](analysis-framework.md#content-tree)；它们不是 HTML 左侧目录，
也不是必须逐项渲染的报告清单。HTML 私有阅读树另行固定为“总体 → 局部邻域 → 空间”，
阅读顺序不等于全串行计算依赖，限制只传递到实际依赖该条件的结论。现有方法、正式结果
和已审阅解读是对齐依据。

### Scope 规则

报告的 parent 列固定为 `Fibroblast`、`Mono/Macro`、`T`。任务安全名
`Mono_Macro` 只用于路径和记录键，必须显式映射到 Raw exact label `Mono/Macro`。
所有 overview 行和正文 scope 都按这一顺序排列，不按结果大小、状态或发现动态排序。

VisiumHD 的 global `All` 是单独的 HD scope，有自己的 cohort、分母、图和结论；它
不进入三个 parent 的分母。Xenium 不创建 global 比较，只报告三个 parent。parent
内部的 assignment change 不能写成 Level1 identity change。`All` 只使用已有的六类
global record：`foundation`、`complexity`、`matched_k`、`moran_all_valid`、
`moran_shared_valid`、`changed_units`；不为 EMT、局部邻域、State 或 Gain 补造 All record。
All 开头显示独立抽样、分群和 Moran 的各自分母。`changed_units` 专表已按
Raw Level1 类型口径保留 Wilson 区间；All 的 `provided_n` 归档缺失时保持 `NA` 并显式说明，另列
已保存完整载体计数 `18085 → 12926`，不把它填入 `provided_n`。

Moran 的 full-side 与 shared-valid 是两个独立小问题；EMT coverage 与 EMT score
也是两个独立小问题。Raw Leiden 与 Raw Level2 是两条独立 baseline，分别保存
`local_vs_raw_leiden` 和 `local_vs_raw_level2`；同一个 Recon 数值在阅读上可以
出现两行，但不能重复累计。本轮组织方案以 State 行的 headline、Gain 只作为
audit；这规定呈现分工，不表示 State 本身证明改善，解释边界见框架节点 4.1–4.5。
Moran 的 full-side/shared-valid 和 EMT 的 coverage/score 在对应五行中保持子问题分开，
但不因此增加总览行。

State 的 overview 只放一个 headline metric：`area_fraction`。area、
`region_windows/valid_windows` 和有效窗口/单位分母留在正文条件列和 named data link。
阈值不稳定时，headline `area_fraction` 为 `NA`/`null`，并写
`no_stable_threshold` 与原因；这表示当前条件下不可识别，不表示零面积。

## 2. 谁写什么

作者先写方法文档，再由 Notebook 和记录引用同一语义。下面的归属是规范的一部分：

| 内容 | 唯一归属 | 允许出现的位置 |
| --- | --- | --- |
| 完整内容树、节点职责、逐项确认与实现差距 | [分析框架](analysis-framework.md) | 各详细规范链接对应节点，不复制整树 |
| Purpose、Computation、Meaning、可能方向与成立条件 | `input-views.md`、`impact-analysis.md`、`gene-and-function.md` 的对应方法段 | HTML module 的简短链接和 Notebook 的简短提示；不在多处重新定义 |
| 节点顺序、双语短目的、方法链接及既有question映射 | [content_contract.py](../../../reproduce/case/reconstruction_impact/content_contract.py) 静态元数据 | 生成器、records与renderer共同引用；不含科学算法或判读规则 |
| 核心参数、分阶段短调用、audit 输出、图调用和一段批次观察 | source Notebook；由 [build_notebooks.py](../../../reproduce/case/reconstruction_impact/build_notebooks.py) 生成 | Notebook 英文 Markdown/代码 cell；长实现留在同目录 helper |
| 长实现和绘图 helper | `reproduce/case/reconstruction_impact/notebook_helpers.py`、`notebook_analysis.py` 以及现有分析模块 | Notebook 只保留可审阅的调用和检查，不复制长函数 |
| 数值、分母、状态、结论、named data links | `analysis/conclusions/`、科学 CSV/JSON 和 artifact index | Notebook/HTML 只消费这些已保存事实 |
| HTML 分组、矩阵、图和证据链接 | [report_records.py](../../../reproduce/case/reconstruction_impact/report_records.py) 与 [report_html.py](../../../reproduce/case/reconstruction_impact/report_html.py) | `report.html`；不得在 renderer 里重新计算科学指标 |
| 稳定的 HTML shell 和 Notebook 英文块格式 | [HTML 模板](templates/report.html)、[Notebook 块模板](templates/notebook-block.md) | 只能各有这一份；不创建副本或四套模板 |

方法段统一采用以下写作入口；它是文档格式，不是新的分析类别：

```markdown
### <method>
**Purpose.** 要回答的科学问题和比较单位。
**Computation.** 输入、公式/算法、参数、保存字段。
**Meaning.** 结果能支持什么，不能支持什么。
**Direction and conditions.** 可能方向、headline 条件、低支持/不可计算状态。
```

逐部分文字模板、字段职责和真实 HD 填写例子统一放在
[输出契约](outputs-and-test-plan.md#narrative-templates)；不另建平行的文字规范。
HTML/Notebook 的既有 shell 仍由各自路径提供，不能反向改变本页的科学定义。Notebook
的 Markdown 和注释使用英文；本目录的设计说明可使用中文。

## 3. HTML 私有阅读树与长页

### 本轮组织方案

科学 27 节点树和 HTML 私有阅读树各自承担不同职责：前者保存科学问题、方法锚点和
旧 question 的追溯关系，后者只规定读者在页面中的移动顺序。HTML 左侧层级目录固定，
正文使用一张可连续滚动的长页；目录不由当前发现、记录数量、图是否存在或结果状态动态
重排。缺少事实时在固定位置显示 `pending`/`NA` 和原因，按需审计可以折叠，但不改变主层级。
Notebook 保留自己的执行/authoring 顺序，不要求与这棵 HTML 阅读树同排版。

```text
HTML private reading tree (fixed left navigation + one long page)
├─ 00 总览（五行固定问题标题）
├─ 01 总体
│  ├─ 分群
│  ├─ 基因信息 + Moran
│  └─ EMT
├─ 02 局部邻域
│  └─ 双 baseline
├─ 03 空间
│  ├─ anatomy
│  ├─ 分群图 / assignment
│  ├─ EMT spatial
│  ├─ local 分层
│  ├─ State 可靠性 → 连续场 → mask / coverage
│  └─ Gain audit（独立入口）
├─ 04 按需审计（foundation、范围/配对、完整基因、预处理）
└─ 05 HD All（仅已有六类 global record；仅 VisiumHD）
```

### 总览五行

| 行 | 固定问题标题 | 本轮组织方案中的可见事实 | 细节入口 |
| --- | --- | --- | --- |
| 1 | 分群 | 两组 K（`complexity`、`matched_k`）与成员改变预览（`unit_change`）；ARI 和 `balanced`（`1 - macroF1`）在正文分群表完整呈现 | 2.1–2.3、按需 partition audit |
| 2 | 基因信息 + Moran | Raw/Recon 的已有 `provided_n`，并列 Moran full-side 与 shared-valid；只有 shared-valid 使用配对 delta | 1.4、2.4、gene/graph audit |
| 3 | EMT | coverage 与 score 分开 | 2.5、3.3、resource/cutoff audit |
| 4 | 双 baseline | Raw Leiden 与 Raw Level2 并列；Kobs、Neff、evenness 共用完整 2×3 图组，并显示窗口支持与尺度 | 3.4–3.8、window/anatomy audit |
| 5 | State | reliability → 连续场 → mask/coverage；Gain 仅从独立 audit 入口进入 | 4.1–4.5、threshold/extent audit |

两组分群记录缺少某字段时保留 `NA`/状态，不从另一组推算；`provided_n` 取已有 gene
availability/foundation 事实，不用 Moran valid n 代替。每个 cell 的 parent 顺序固定为
`Fibroblast`、`Mono/Macro`、`T`，HD `All` 仍按 scope 规则独立处理。

<a id="report-display-status"></a>
### Renderer 实现与检查状态

两份正式报告已按“总体 → 局部邻域 → 空间”重排，具有固定左侧层级目录、五行总览、
双 baseline 对照、State/Gain 独立入口和 HD All 独立分母。图旁补充简短结果文字；
已有分析、记录身份与已审阅中英解读继续保留，未运行科学重计算。

2026-09-15 展示检查已完成：两份报告各 38 张登记图完整展示一次，内部锚点无重复或缺失，
所有本地证据链接存在；1440px 与 390px 下无整页横向溢出，图片加载、目录点击、折叠展开、
浏览历史和旧行锚点跳转正常。展示层测试文件 26 项通过。**用户对阅读效果的最终评价尚未确认。**
历史科学执行与 R0–R4 记录保持原有身份，不改标为本轮重新通过。

正文 module 的 renderer anchor 是 `module-<question>`，例如
`#module-moran_all_valid`；HD All 在末尾加 `-all`。旧 scope-row 锚点保留或映射到对应新问题。overview 和正文通过 renderer-owned IDs 互链，
不依赖中文标题自动 slug。record 的 `evidence.tables` 保存 named CSV paths，
renderer 将它们作为具名数据链接消费；不要另造 `evidence_tables` 字段或在 HTML 中猜路径。

网页正文 module 采用“**固定问题标题 → 简短结果 → 表/图 → 1–3 条同源解读**”，并按
“总体 → 局部邻域 → 空间”的长页顺序排列。保存的 `node_narratives.next` 继续保留在
结果包并供 Notebook 使用；网页不沿用旧的 `next` 顺序组织 module 或导航，网页导航服从固定
左侧目录。
小表只放该问题需要的数字、分母、条件和具名数据链接；有科学作用的图紧随小表，
没有图时保留文字和链接，不为填充版式补图。基础节点不作为主文 module，只有读者请求
或需要审计时才展开。`next` 只作为保存的叙述字段，不新增网页科学结论。

同一张 common figure 只在拥有它的 module/row 完整嵌入一次，其他 scope 只给链接。
parent-specific spatial images 随对应 scope row 展示；不能把 Fibroblast 的图
当作三个 parent 的共同证据。原图路径只能取自 record 的 evidence；固定的展示映射规定已登记图的所属模块和图注，不扫描或猜测未登记文件。

`report.html` 不读 H5AD、不调用科学计算、不补齐缺失值。完成后重绘只改变展示；如果
记录、参数、实现或资源改变，必须产生新的执行身份或待审阅状态。

正式示例见 [输出契约的示例索引](outputs-and-test-plan.md#report-examples)；链接指向两个唯一正式结果集，不复制报告或数据。

## 4. Notebook 入口（不要求同 HTML 排版）

两个 canonical source notebook 是：

- [VisiumHD_sp_SVC_Reconstruction_Impact.ipynb](../../../reproduce/case/reconstruction_impact/VisiumHD_sp_SVC_Reconstruction_Impact.ipynb)
- [Xenium_sc_SVC_Fibroblast_Reconstruction_Impact.ipynb](../../../reproduce/case/reconstruction_impact/Xenium_sc_SVC_Fibroblast_Reconstruction_Impact.ipynb)

五行与科学节点职责见 [共用内容树](analysis-framework.md#content-tree)，HTML 私有阅读树见
本页第 3 节；具体 authoring 契约见 [输出文档](outputs-and-test-plan.md)。Notebook 只需
保留可审阅的执行顺序、核心参数、短调用、检查、图和同源观察，不要求复制左侧目录或
长页排版。

Notebook 按“目的与方法链接 → 短调用 → 必要检查与保存 → 有科学作用的图 → 同源批次观察”组织。
技术审计不强制配图，也不强制在主文显示；正式观察由同一批 records 在相应节点填入，
不在 cell 中另写一套判断。实现细节由 [notebook helpers](../../../reproduce/case/reconstruction_impact/notebook_helpers.py)
和 [notebook analysis](../../../reproduce/case/reconstruction_impact/notebook_analysis.py)
承载；helper 绿色测试或旧执行结果不能替代 source/executed parity 和真实网页展示检查。

## 5. 保留的科学边界

- Raw 和 Recon 各自保留完整可用基因空间。共同有效基因只用于 shared-valid 配对；
  不把 missing、zero、low coverage 和 not computable 混成零。
- Raw-derived QC/HVG 只约束 partition。Moran 使用同组同图，最低支持
  `n_units >= 51`；EMT 使用固定资源、实际命中数和 provider cutoff；二者按各自
  availability 记录。
- State 是重建后的高多样性状态，不是变化最大区域。Gain 保留为
  audit-only。State 阈值至少需要 80% 有效 bootstrap 且 CI width 不超过有效 Recon
  Neff 范围的 25%；否则不输出可用 mask。
- 负 delta、低支持、unmatched cluster、资源覆盖不足和 no-stable-threshold 都是
  可以发布的科学状态。缺依赖、损坏资源、invalid axes 和代码异常是执行失败。
- 本轮不新增 marker/DEG/trajectory/interaction 或应用基因空间图，也不修改结果
  目录外的数据和方法代码。报告呈现检查与 route 科学结果验收分别记录。

## 6. 文档索引

| 文档 | 只负责什么 | 何时读取 |
| --- | --- | --- |
| [README](README.md) | 总入口、scope、内容归属与导航 | 开始任务、写记录或审阅报告 |
| [analysis-framework.md](analysis-framework.md) | 目的、完整内容树、节点职责、确认状态和实现差距 | 逐层讨论与内容审查 |
| [input-views.md](input-views.md) | Raw/spatial/reference 载体、cohort、坐标、projection、完整 gene availability | foundation 或路线改变 |
| [impact-analysis.md](impact-analysis.md) | partition、matched-K、Raw Level2、window、anatomy、local diversity、State/Gain Region | overall/localization/region 方法实现 |
| [gene-and-function.md](gene-and-function.md) | Moran full/shared、EMT coverage/score 和空间字段 | Moran、AUCell 或 EMT 空间图 |
| [outputs-and-test-plan.md](outputs-and-test-plan.md) | 记录、表、状态、HTML/Notebook 模板契约和 R0–R4 验收 | 保存结果、渲染、重载和验收 |
| [batch-integration.md](batch-integration.md) | 2.0/历史 handoff 的未来接入边界 | 独立 batch 设计；不是当前执行入口 |

跨文档链接指向上述规范和实际 renderer/builder；不能再创建平行的 reconstruction-impact
计划。任何新判断先更新相应唯一归属，再由 Notebook/record/HTML 引用。
