# REVISE 重建与分析衔接：事项总览

本目录组织 REVISE 与 `REVISE_Analysis_Agent` 的协议迭代：为什么做、已经讨论到哪里、后续需要什么，以及如何进入执行计划。核对日期：2026-09-20。

**当前已完成 R1–R8 已定范围的工程实现（R8b 限本地来源），以及 D1/C1；本轮完成双仓修改及约1%连续ROI真实联测。** 真实mini四方法比较已执行；全量、生物学解释与最终默认值后置；D2 保持现状，D3 和在线来源不进入近期开发。见[双库审阅入口](cross-repo-review.md)和[当前验收表](acceptance.md)。

## 从问题进入

| 想了解什么 | 入口 |
|---|---|
| GPT 同时审阅双库 | [协议、证据与审阅说明](cross-repo-review.md) |
| 四方法对照原空间亚型 | [Notebook 比较协议](assembly-comparison-contract.md) |
| 全部事项及当前状态 | [完整事项表](#items) |
| 每一项落实到哪个代码文件、测试如何检查 | [落实文件索引](implementation-index.md) |
| 第一批能先做什么、为什么 | [一：范围已定的小改动](categories/01-small-changes.md) |
| 如何交付完整 Raw、统一 SVC 并接入分析 | [二：统一样本交付](categories/02-sample-delivery.md) |
| 第二批如何一次实施、如何验收 | [第二批总计划](plans/whole-sample-delivery.md) · [验收清单](acceptance.md) |
| reference 如何独立准备并接入单样本/batch | [Reference Preparation](plans/reference-preparation.md) |
| 两种 OT assembly 与比较 Notebook | [三：assembly 方法与比较](categories/03-assembly-research.md) |
| 等待用户方案、外部材料和后续研究的事项 | [四：后续扩展](categories/04-future-extensions.md) |
| 已确认什么，历史事项是否遗漏 | [共识与来源索引](decisions-and-sources.md) |

阅读顺序为 **总览 → 类别中的具体事项 → 落实文件索引／对应执行计划**。本页是状态的唯一维护位置；类别页维护原因、共识、依赖和讨论问题；落实索引维护源码与测试证据；计划页维护实施细节。事项编号用于本目录导航，不是已发布的 issue ID。

<a id="items"></a>
## 完整事项表

“讨论已确认”指行为方向；“计划已有草案”仍需逐项审阅，不能视为“执行计划已对齐”。R1 已完成；R2/R3/R4 的生产端实现和工程验收已完成，但mini正式交付与消费已完成，hST unknown表达和sST校正缺口单独保留；全量与科学解释未完成；R7/R8 本地准备、R5/R6/D1/C1 均已完成工程验证；当前真实结果见验收表。

| 类别 | 事项与细节入口 | 讨论 | 执行计划 | 实施 | 下一步 |
|---|---|---|---|---|---|
| 一 | [R1 sST 最终逐细胞归一化删除](categories/01-small-changes.md#r1) | 已确认 | [已对齐](plans/sr-final-scaling.md) | 已有实现；385 项针对测试通过 | mini已核查，parent校正存在有效零项；见验收表 |
| 二 | [R2 全类型重建与 sample 级 SVC](categories/02-sample-delivery.md#r2) | 已确认：批次行为已对齐 | [第二批总计划](plans/whole-sample-delivery.md#r2) · [详细计划](plans/whole-sample-svc.md) | 生产端实现；工程验收完成；mini已联测；全量待验收 | P2 mini已核对；后续服务器全量 |
| 二 | [R3 完整 Raw 与已有推断](categories/02-sample-delivery.md#r3) | 已确认：Raw、推断列和冲突规则已对齐 | [第二批总计划](plans/whole-sample-delivery.md#r3) · [详细计划](plans/raw-publication.md) | 生产端实现；工程验收完成；mini已联测；全量待验收 | 使用 P2CRC_Xenium 完成真实 Raw/provenance 核对 |
| 二 | [R4 分析配置与联合验收](categories/02-sample-delivery.md#r4) | 已确认：交接边界已对齐 | [第二批总计划](plans/whole-sample-delivery.md#r4-c1) · [详细计划](plans/analysis-handoff.md) | fixture 消费已验证；mini已联测；全量待验收；真实 sST 表达验收未完成 | 使用 P2CRC_Xenium 完成后续 loader 与标签/空间消费检查 |
| 二 | [C1 分析端消费协议协同](categories/02-sample-delivery.md#c1) | 已确认线性输入契约与责任边界 | [已对齐](plans/ot-assembly.md) | 生产端及iST/sST mini表达联测通过 | hST unknown、sST校正缺口和全量/科学解释后置 |
| 三 | [R5 OT assembly](categories/03-assembly-research.md#r5) | 已确认两种模式 | [已对齐](plans/ot-assembly.md) | 已有实现；工程验证完成 | 用户逐项审阅；真实验收另定 |
| 三 | [R6 assembly 比较与最终默认值](categories/03-assembly-research.md#r6) | 已确认 Notebook 目标；科学结论后置 | [已对齐](plans/ot-assembly.md#r6比较-notebook) | 已有实现；工程验证完成 | mini Notebook已执行；科学解释后置 |
| 四 | [D1 confidence](categories/04-future-extensions.md#d1) | 已确认与 R7 合并定义 | [已对齐](plans/ot-assembly.md) | 已有实现；工程验证完成 | 已共用函数并记录阶段来源 |
| 四 | [R7 Reference 评估、排序、选择](categories/04-future-extensions.md#r7) | 已确认：本地候选、三种分数、max_median top-1 | [已对齐](plans/reference-preparation.md) | 已实现并通过工程验收 | 真实规模与科学质量后续评估 |
| 四 | [R8a Patient ID 候选来源](categories/04-future-extensions.md#r8a) | 已确认：显式配对列和值，不推断 Patient | [已对齐](plans/reference-preparation.md) | paired 提取及消费已验证 | 按实际 pool 提供明确配对列和值 |
| 四 | [R8b CELLxGENE 候选来源](categories/04-future-extensions.md#r8b) | 已确认本地候选；在线发现仍待讨论 | [本地部分已对齐](plans/reference-preparation.md) | 本地筛选已验证；在线获取未开始 | 公开检索继续后置 |
| 四 | [D2 gene-wise uncertainty](categories/04-future-extensions.md#d2) | 已确认保持现状 | 本批不修改 | 保留已有实验能力 | 不列近期开发 |
| 四 | [D3a out-of-spot](categories/04-future-extensions.md#d3a) | 部分确认：不伪造 parent spot | 尚未撰写 | 未开始 | 讨论生成与来源语义 |
| 四 | [D3b 复杂 correspondence](categories/04-future-extensions.md#d3b) | 待真实需求再讨论 | 尚未撰写 | 不新增通用映射层 | 有关系实例后重新判断 |

以下内容同样保留，避免误列成新开发或以后重复提出：

| 性质 | 事项 | 当前依据与处理 |
|---|---|---|
| 已有能力 | [mean/random assembly](categories/02-sample-delivery.md#existing-assembly) | 实现已静态核对；不等于 whole-sample 已实现，本批 application 针对套件通过，非真实样本验收 |
| 已有能力 | [分析库 backend 解耦](categories/02-sample-delivery.md#backend-independence) | 已核对包依赖与源码；不替代新协议的联合验收 |
| 明确不新增 | [普通 gene-overlap 持久化](categories/02-sample-delivery.md#no-gene-overlap) | 现场可推导，保留原因，不创建删除任务 |

<a id="cadence"></a>
## 先后关系与讨论节奏

1. 第一批 R1 已完成；第二批已完成 R2/R3/R4 的生产端实现和工程验收，内部按“全类型计算 → Raw 与发布 → 分析交接 → 整体验证”执行。P2CRC_Xenium mini与sST真实表达消费已完成；全量、sST校正缺口及科学解释仍待后续，不为每个小事项重新启动流程。其余第三、四类始终保留；R7/R8 已完成本地准备与接入。
2. R2 是 iST 全类型输出主线；R3 的原始对象保护可并行推进，但完整 iST annotation 验收依赖 R2。R4 联合验收依赖 R2/R3 的同次交付。
3. sST 表达分析额外依赖 C1 的真实消费能力。数据产品交付、成功加载、标签/空间分析和具体科学解释分开判断，统一记录在[验收清单](acceptance.md)。
4. R5 两种 OT 模式与 R6 比较 Notebook 已完成工程验证；T/Mono_Macro/Fibroblast mini重分群对照已运行，最终默认值后置。R7/R8 不以真实比较完成为前提。
5. D1 与 reference 评分共用数值定义，保留现有列；不新建 posterior 协议。R5/R6/C1 已按[本批计划](plans/ot-assembly.md)完成工程验证，本轮mini验收见当前验收表，全量随后在服务器执行。

## 后续如何更新

有新讨论时，更新对应类别中的“共识／待讨论”，在[决策索引](decisions-and-sources.md#chronology)记录实质变化，再更新本页三种状态。第二批实施过程中的工程处理写入对应详细计划，最终测试、真实运行和联合验收写入[验收清单](acceptance.md)；“已落实”必须附相应代码或验证依据，静态检查不能代替科学验收。

新增计划只放到已有事项下面，并从类别页索引。未知问题保留具体下一步，不用空的实现计划占位。需要发布 tracker issues 时再沿用本目录的事项与依赖；本轮没有创建外部 issues。
