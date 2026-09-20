# REVISE 事项落实文件索引

[双库审阅入口](cross-repo-review.md) · [Notebook 比较协议](assembly-comparison-contract.md) · [返回总览](README.md) · [完整事项表](README.md#items) · [验收清单](acceptance.md)

本页是每个事项的“落实证据入口”。它把三件容易混淆的事情分开：

1. **落实状态**：当前是否已经有本仓库的生产实现、工程测试，或只是讨论/规划。
2. **核心设计逻辑**：为什么这些文件应当一起出现，以及实现要保持什么边界。
3. **具体检查位置**：源码负责什么，测试能证明什么；计划或现有入口若没有新实现，会明确标为“参考，不是实现证据”。

2026-09-20 更新：双仓mini真实联测已完成，任务编号、产物、局部限制与服务器入口统一见[当前验收](acceptance.md)。下表沿用历史R1–R8编号；历史工程数量不作为本轮通过数。

行号是 2026-09-19 核对时的入口，代码继续变化后应重新核对。链接到源码的行号用于定位，不把静态阅读误报为真实样本或科学验收。

## 状态与证据等级

| 表达 | 含义 |
|---|---|
| **已落实：工程通过；mini 已联测；全量/科学待验收** | 生产代码和针对性测试已存在，但仍需要真实输入、科学解释或用户验收。 |
| **部分落实** | 只有部分协议或 fixture/交接路径已验证，尚不能宣称完整消费协议或真实样本完成。 |
| **已有能力／明确不新增** | 已有代码能力或经过讨论决定不新增；这不是本轮新开发任务。 |
| **未落实** | 没有本事项的新生产实现；下面的文件只是现有入口、实验代码或计划参考。 |

## 全部事项的落实与检查表

| 事项 | 当前是否落实 | 核心设计逻辑 | 具体实现／证据文件 | 大致验收方式 |
|---|---|---|---|---|
| **R1 sST 最终逐生成细胞缩放删除** | **已落实：工程通过** | 保留内部归一化、建图副本 normalization/log1p 和 parent-spot 逐基因校正；删除最后把每个生成细胞单独缩放到 10,000 的分支，不恢复 raw counts。 | 生产校正：[sc_svc_super_resolution_application.py#L295-L299](../../../revise/backend/runners/sc_svc_super_resolution_application.py#L295-L299)；旧键拒绝：[application/config.py#L173-L177](../../../revise/application/config.py#L173-L177)、[config/loader.py#L212-L219](../../../revise/config/loader.py#L212-L219)；测试：[test_sc_sr_guidance.py#L267-L297](../../../tests/backend/test_sc_sr_guidance.py#L267-L297)、[test_request.py#L299-L332](../../../tests/application/test_request.py#L299-L332)。 | 检查每个 parent spot 汇总后逐基因匹配内部 normalized spot；生成细胞自身不再都为 10,000；旧键 `true/false` 都报迁移错误，删除后配置通过。 |
| **R2 iST 全类型重建与 sample 级 SVC** | **已落实：工程通过；mini 已联测；全量/科学待验收** | 一次 GA 结果服务整个样本；遍历实际 broad type；有效 Level2 超过一个才做 LR；不合格类型记录跳过，合格类型失败使样本失败；严格合并基因轴并为 cluster 加类型命名空间。 | GA 复用和全类型编排：[adapters.py#L516-L521](../../../revise/backend/adapters.py#L516-L521)、[adapters.py#L523-L616](../../../revise/backend/adapters.py#L523-L616)；命名与严格合并：[adapters.py#L206-L245](../../../revise/backend/adapters.py#L206-L245)；工程测试：[test_sc_local_ot.py#L626-L708](../../../tests/backend/test_sc_local_ot.py#L626-L708)、[test_sc_local_ot.py#L754-L790](../../../tests/backend/test_sc_local_ot.py#L754-L790)、[test_whole_sample_delivery.py#L10-L93](../../../tests/integration/batch/test_whole_sample_delivery.py#L10-L93)。 | 用多 eligible 类型、单 Level2、部分/空白 Level2、全跳过和 eligible 失败 fixture 检查 GA 次数、跳过原因、失败边界、cluster 不碰撞、direct/batch 一致；再用 P2CRC_Xenium 做真实类型覆盖核对。 |
| **R3 完整 Raw 与已有推断** | **已落实：工程通过；mini 已联测；全量/科学待验收** | 从未被预处理污染的原始对象发布；保留 X、obs/var 轴、坐标和人工标签；真实 ID 回填 `revise_Level1/2`，冲突统计独立记录，不覆盖原始列，不为缺失单位补标签。 | Raw 规则：[delivery.py#L14-L82](../../../revise/application/delivery.py#L14-L82)；原始快照和来源校验：[reconstruct.py#L81-L118](../../../reconstruct.py#L81-L118)；发布事务：[publication.py#L186-L256](../../../revise/application/publication.py#L186-L256)、[publication.py#L266-L323](../../../revise/application/publication.py#L266-L323)；测试：[test_sample_delivery.py#L43-L100](../../../tests/application/test_sample_delivery.py#L43-L100)、[test_sample_delivery.py#L128-L169](../../../tests/application/test_sample_delivery.py#L128-L169)。 | 检查 QC 排除单位仍在 Raw，轴/X/坐标不变，推断按真实 ID 回填，冲突与缺失正确；注入计算、写入、安装失败，确认旧完整产物不变。 |
| **R4 sample.yaml 与分析联合验收** | **生产端及 fixture 联测完成；mini 已联测；全量/科学待验收** | Raw、SVC、`sample.yaml` 同次交付；配置只声明真实知道的列、表达尺度和坐标单位；分析端通过独立 loader 消费已发布文件，不依赖 REVISE backend。 | 配置生成：[delivery.py#L85-L118](../../../revise/application/delivery.py#L85-L118)；交付路径：[publication.py#L17-L21](../../../revise/application/publication.py#L17-L21)；真实消费者入口 fixture：[test_analysis_delivery.py#L19-L72](../../../tests/integration/test_analysis_delivery.py#L19-L72)。 | 在独立 cwd 用分析库真实 `load_sample` 和已有标签/空间分析入口；检查结果状态、报告文件和输入 hash；随后用 P2CRC_Xenium 重复，不能以 fixture 代替真实验收。 |
| **C1 分析库消费协议协同** | **生产端对接及表达 fixture 验证完成** | 消费者已定有限非负非 log 线性 X；有明确来源才声明可用，未知不升级。 | [源表达声明](../../../revise/application/expression.py)、[配置生成](../../../revise/application/delivery.py)、[真实消费者联测](../../../tests/integration/test_analysis_delivery.py)。 | 已知线性浮点表达可运行正式表达分析；unknown 保持不可用，输入哈希不变；不修改分析库。 |
| **已有能力：mean/random assembly** | **已有能力；已回归** | 在统一 spatial observation axis 上使用 expression carrier；mean 取 cluster 均值，random 按稳定 seed 选择 donor；它不是 cluster assembly。 | 实现：[ist_assembly.py#L59-L105](../../../revise/application/ist_assembly.py#L59-L105)；测试：[test_ist_publication.py#L29-L67](../../../tests/application/test_ist_publication.py#L29-L67)。 | 检查 mean 数值、random seed 可复现、gene/cluster 轴一致；真实科学优劣留给 R6，不把工程通过写成默认值科学结论。 |
| **已有能力：分析库 backend 解耦** | **已有能力；交接仍需联合验收** | 消费者从 Raw/SVC/sample.yaml 读取，不把 REVISE backend 作为运行前提；独立 loader 的真实消费测试仍是协议证据。 | 本仓库的独立交接测试：[test_analysis_delivery.py#L19-L72](../../../tests/integration/test_analysis_delivery.py#L19-L72)；发布输入读取：[batch/inputs.py#L213-L282](../../../revise/batch/inputs.py#L213-L282)。 | 在没有 REVISE backend 运行环境的独立 cwd 加载并分析，确认文件和输入未被改写；不能仅凭 import 解耦推断方法学支持。 |
| **明确不新增：普通 gene-overlap 持久化** | **明确不新增** | Raw/SVC 的共同基因和各自 var 轴可由实际对象推导；只有 gene-wise uncertainty 另有科学定义时才讨论保存量。 | 当前对象按自身轴读取和校验：[batch/inputs.py#L255-L282](../../../revise/batch/inputs.py#L255-L282)；决策说明：[02-sample-delivery.md#no-gene-overlap](categories/02-sample-delivery.md#no-gene-overlap)。这里没有待实现代码位置。 | 不创建新的 overlap 字段；检查两侧 var_names 与分析需要时现场推导，避免把“未新增”误报成遗漏。 |
| **R5 两种 OT assembly** | **工程实现及验证完成** | 同群 ST 邻居 0.2/自身 0.8；within 匹配群内 donor，outside 匹配同 broad type 全部 cluster profiles；TACCO 权重重建全基因表达。 | [assembly](../../../revise/application/ist_assembly.py)、[OT 实现](../../../revise/application/ist_ot.py)、[测试](../../../tests/application/test_ist_ot.py)、[计划](plans/ot-assembly.md)。 | 验证空间边界、权重轴、数值、资源上限、真实 TACCO fixture 与失败保护；不代表科学效果。 |
| **R6 比较 Notebook 与最终默认值** | **工程与mini真实比较完成；全量/科学解释后置** | T/Mono_Macro/Fibroblast 重建后重新 Leiden，与原空间 subtype 比较；不把旧标签当新结果。 | [Notebook](../../../reproduce/case/assembly_comparison.ipynb)、[计算模块](../../../revise/analysis/assembly_comparison.py)、[测试](../../../tests/analysis/test_assembly_comparison.py)。 | 小 fixture 执行 Notebook、共同轴/缺失、列联表/ARI/NMI/图表；mini真实运行已完成，全量与默认值另定。 |
| **D1 confidence 与 reference 评分** | **工程实现及验证完成** | 复用现有 confidence；与 R7 共用最大值、certainty/margin 计算；不混淆阶段，不新增 posterior 协议。 | [公共函数](../../../revise/utils/confidence.py)、[reference 评分](../../../revise/reference_preparation/scoring.py)、[annotation](../../../revise/backend/kernels/ot.py)、[测试](../../../tests/reference_preparation/test_confidence_shared.py)。 | 旧数值/排名不变；阶段与类别元数据正确；assembly 不覆盖 annotation confidence。 |
| **R7 Reference 评估与选择** | **工程实现及验证完成** | 独立候选 GA、三种评分、max_median top-1、严格概率与共同 ST 轴、完整证据；科学质量另验收。 | [准备核心](../../../revise/reference_preparation/api.py)、[GA 接口](../../../revise/reference_preparation/host.py)、[覆盖入口](../../../reconstruct.py)、[batch](../../../revise/batch/runner.py)；[测试](../../../tests/reference_preparation)、[真实小型 GA](../../../tests/application/test_reference_ga.py)。 | [两步指南及验收](plans/reference-preparation.md)；真实规模与科学评估后置。 |
| **R8a 显式配对来源** | **工程实现及验证完成** | 指定列和值精确提取，不推断患者关系；主流程消费清空旧 filter。 | [paired](../../../revise/reference_preparation/paired.py)、[提取测试](../../../tests/reference_preparation/test_paired.py)、[direct/batch 交付](../../../tests/integration/batch/test_whole_sample_delivery.py)。 | 运行时提供真实 pool 与明确配对信息；[计划](plans/reference-preparation.md)。 |
| **R8b 本地候选 / CELLxGENE 来源** | **本地筛选完成；在线发现后置** | 支持目录直接 H5AD 与清单，不宣称在线获取完成。 | [配置](../../../revise/reference_preparation/config.py)、[筛选](../../../revise/reference_preparation/screening.py)、[命令入口验证](../../../tests/application/test_reference_screen_entry.py)。 | [后续来源问题](categories/04-future-extensions.md#r8b)。 |
| **D2 gene-wise uncertainty** | **保持现状，本批不修改** | 普通 gene overlap 与 gene-wise uncertainty 不同；要先定义算法、轴、数值含义、保存量和消费者。 | 仓库有实验代码：[gene_uncertainty.py#L88-L119](../../../revise/backend/kernels/gene_uncertainty.py#L88-L119)，**仅是实验线索，不是已确认公共协议**。规划入口：[04-future-extensions.md#d2](categories/04-future-extensions.md#d2)。 | 本批不改方法；未来有扩展需求时再核对数值定义、轴与消费方。 |
| **D3a out-of-spot 与 provenance** | **未落实** | 未来没有真实 parent spot 的生成单位不能伪造 `spot_name`；算法、空间范围、表达来源和标识必须一起定义。 | 当前只实现真实 parent spot 的 sST baseline：[batch/inputs.py#L330-L359](../../../revise/batch/inputs.py#L330-L359)；这不是 out-of-spot 实现。规划入口：[04-future-extensions.md#d3a](categories/04-future-extensions.md#d3a)。 | 有具体生成案例后，验收来源可追溯、无伪造 parent、ID/空间/表达语义一致；此前不新增占位接口。 |
| **D3b 复杂 Raw↔SVC correspondence** | **未落实通用映射层** | 先保留真实 native ID、parent spot 或 cell ID；只有真实存在多对多等关系且分析需要时，才设计 sidecar 或独立映射协议。 | 当前是简单 ID 对齐：[batch/inputs.py#L307-L328](../../../revise/batch/inputs.py#L307-L328) 与真实 sST parent 映射：[batch/inputs.py#L330-L359](../../../revise/batch/inputs.py#L330-L359)；没有通用 correspondence 实现。规划入口：[04-future-extensions.md#d3b](categories/04-future-extensions.md#d3b)。 | 提供实际关系例子和消费问题后，先验证现有字段是否足够，再决定是否要 sidecar；不把当前简单对齐扩大成通用协议。 |

历史工程结果：**964 项测试 + 24 项子测试通过**，真实样本与科学验收未执行。[验证记录](verification-ot-assembly-2026-09-19.json) · [配置与使用](plans/ot-assembly.md)。

## 应如何使用这份索引

- 验收一个已经实施的事项时，先看同一行的生产文件，再看测试文件，最后对照[验收清单](acceptance.md)区分工程、真实样本和科学结论。
- 若一行写的是“现有入口”“实验代码”或“规划入口”，它只说明当前边界，不代表事项已经实现。
- 代码位置和状态发生变化时，先更新本页，再更新总览的实施状态；不要只改计划页而留下无法核对的“已落实”。

历史双库审阅前浅验收：**75 项通过，0 失败、0 跳过**；sST 坐标与 broad 别名交接修复已覆盖。执行 Notebook、4 张 PNG 和消费者结果见[审阅入口](cross-repo-review.md#本次实际验收结果)。旧 964 项为上一批快照，不与本轮相加。


2026-09-19 协同说明更新（仅文档）：双库审阅统一见 [交接入口及 H1–H6 任务](cross-repo-review.md#协同任务与责任清单)。hST/sST 本来不涉及 SVC_cluster，不列为待补功能。75 项为此前快照。2026-09-20 已完成基线汇集、四方法mini真实生成及严格共同坐标校验，当前证据见验收表。
