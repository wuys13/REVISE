# 二、统一样本交付与分析衔接

> **协议更新（2026-09-19 OT assembly 批次）：** 当前消费者已固定线性非负、未取 log 的 X；旧 log 声明已拒绝。现行声明和处理以[新计划](../plans/ot-assembly.md)为准。真实验收继续后置。
[返回总览与状态](../README.md) · [落实文件索引](../implementation-index.md) · [共同决策](../decisions-and-sources.md#confirmed-decisions)

本类解决一个完整路径：重建完成后，分析者拿到原始测量、可分析的统一空间表达及真实语义，不必理解内部 cell-type 目录与 paired carriers。R2、R3、R4 已获用户合并授权，生产端实现和工程验收已经完成；真实 P2CRC_Xenium、科学解释和 sST 表达联合验收仍按[第二批总计划](../plans/whole-sample-delivery.md)记录。

目标产品是 `raw.h5ad + SVC.h5ad + sample.yaml`。重建 manifest、日志和内部载体由 REVISE 管理；不会全部变成分析输入要求。CLI、输出名称、返回值与旧配置迁移已按相应计划落实并进行工程回归；真实样本仍单独验收。

<a id="r2"></a>
## R2：全类型重建、合并与 sample 级发布恢复

**为什么做。** 原先 iST/sc-SVC cluster 路径要求一个具体 broad type，batch 按 sample × cell type 组织。单类型结果不能直接支持整张组织的跨类型分析。whole-sample 在本期指所有符合条件类型的合并结果，不代表所有 Raw 单位均被重建。

**现状。** GA 与 LR 已有统一生命周期；mean/random 已实现。生产端现已支持一次 GA 后遍历类型、合并结果、sample 级失败/恢复及公共交付。输入校验、跳过规则、direct/batch 一致性和发布恢复已由工程套件覆盖；真实样本仍需单独运行核对。

**已确认行为。** GA 后遍历实际出现的全部 broad types；对应有效 reference 中有超过一个不同的有效 Level2 标签才重建，否则提示该类型未纳入、需要提前完善 Level2 注释。只合并重建类型；random 默认、mean 保留。符合条件类型的计算错误使样本失败、不发布新结果，不变成静默跳过。

**本批已对齐。** 未指定单类型选择时使用 whole-sample；显式 `select_ct`、`select_cell_type`、batch `cell_type` 和旧 `cell_types` 列表保留单类型兼容语义，但不能把兼容调用标成完整样本交付。whole-sample 中 reference 部分细胞缺失 Level2 时，LR 只使用有效、非空 Level2 行且要求超过一个类别；显式旧单类型允许一个有效 Level2，零个仍失败。缺失不能转成字符串类别。跨类型合并校验和 cluster 标识必须避免碰撞，不增加未讨论的最低细胞数阈值。全跳过不发布空成功。

**依赖与节奏。** R2 与 R3 的 Raw source snapshot 可并行推进，但完整 iST 回填和统一发布在二者一致后验收。资格记录、发布和恢复随本批统一验收，不单独留下只能在内存运行的半条路径。

**需要用户参与。** 本批无需提前补充材料；旧遍历代码/Notebook 仅作可选参考。资格、random 默认、失败规则和旧入口兼容范围已纳入批次计划。执行中若无法取得可信来源或必须改变公共语义，再上提具体问题。

**下一步与证据。** 工程证据、测试计数和剩余真实样本边界见[验收清单](../acceptance.md)及[验证记录](../verification-2026-09-19.json)；后续使用 P2CRC_Xenium 完成真实输入核对。详细算法边界见[全类型 SVC 计划](../plans/whole-sample-svc.md)。当前落点见 [ScSvcApplicationStrategy](../../../../revise/backend/adapters.py)、[单类型 LR](../../../../revise/backend/runners/sc_svc_application.py)、[batch 任务与恢复](../../../../revise/batch/runner.py)。这里针对 iST，不把 sST spot 分配算法改成同一循环。

<a id="r3"></a>
## R3：完整 Raw 与本次已有推断

**为什么做。** 原始测量不能因 QC、参考基因交集或为了 pairwise 比较而缩小。重建已经计算出的重要 annotation 应交付给分析，避免分析端再次运行 OT/TACCO。

**现状与目标。** 生产端现已完成统一 Raw finalization。新 Raw 保留原始 X、完整 obs/var 轴及坐标，按真实 ID 回填已有 Level1/Level2。hST/sST 与 iST 的可用推断不同；缺失不伪造。sST 工作对象有原地 normalization，不能将其当成 Raw 发布。

**共识与边界。** confidence 复用现有列；最新 D1 与 reference 评分共用定义，不新增 GA/LR 拆列或完整 posterior 保存协议。Raw/SVC 默认独立；hST/iST 保留真实 ID，sST 保留真实 spot_name/cell_id。普通 gene-overlap 不持久化。coverage boolean 是历史候选，不作为本期已确认公共字段。

**本批已对齐。** 原始输入已有同名人工标签时保留原始列，新推断写入独立的 `revise_Level1`、`revise_Level2`；两侧冲突记录统计，不覆盖原始值。没有推断的单位保留已有注释并维持本次推断缺失；不得以“保存 inference”为名覆盖原始文件或补造标签。各 route 能提供的标签和 confidence 仍按实际能力声明。

**依赖与节奏。** 先建立不受预处理污染的原始对象边界，再接各 route 的实际推断；同次 Raw/SVC 成功交付一起验收。完整 iST 路径依赖 R2，原始矩阵保护可先做，但不单独发布半成功 Raw。

**需要用户参与。** D1 已与 R7 合并，无需另交新方案。标签冲突处理已按独立列和冲突记录纳入本批；真实验收优先使用现有 P2CRC_Xenium，确实缺标签来源或原始矩阵证据时再请用户补充。

**下一步与证据。** 工程完整性、transaction、命名空间冲突和剩余真实样本边界见[验收清单](../acceptance.md)、[验证记录](../verification-2026-09-19.json)和[Raw 计划](../plans/raw-publication.md)。来源见 [输入与预处理编排](../../../../reconstruct.py)、[GA/LR 对象](../../../../revise/backend/adapters.py)、[当前 Raw 读取视图](../../../../revise/batch/inputs.py)。

<a id="r4"></a>
## R4：样本配置与联合验收

**为什么做。** 统一文件只是交付的一部分；分析还需要知道哪个样本、哪些列、什么坐标单位和表达来源。交接应让独立分析库直接消费，而非建立旧 carrier 兼容层。

**现状。** 分析库支持 `sample.yaml` 和独立 Raw/SVC；REVISE 已生成与该接口配套的完整产品，并通过合成 fixture 的真实 loader/标签空间分析入口。消费者处于并行修改中，不能把本次工作树快照称为稳定发布契约；真实 P2CRC_Xenium 仍未运行。

**当前工程进展。** 合成交付 fixture 已通过消费者正式入口的 loader、标签/空间
diversity 和报告生成检查，输入未被改写；整体状态因 expression identity unknown
为 `partial`；后续已知线性浮点 fixture 的正式表达分析也已通过。这不代表 P2CRC_Xenium 真实样本或用户科学验收已完成。

**已有共识。** 配置只声明真实已知信息；缺失前提如实报告；只有具体比较需要时才配对。SVC_cluster 必须是真实已有标签，不为 hST/sST 造一列。分析端不依赖 reconstruction backend。

**本批已对齐。** batch 的层级 sample ID 使用确定、可逆的单目录段编码并保留原始映射；单样本坐标、各 route 的矩阵语义和部分分析不可用时的边界按实际来源声明。代表验收样本优先使用 P2CRC_Xenium，分析结果状态与工程交付分开记录。

**依赖与节奏。** R2/R3 提供配套文件后完成联合验收；sST 表达消费依赖 C1。不把 loader 成功、空间标签分析完成、表达分析完成混成同一个结论。

**需要用户参与。** 本批先按已对齐的样本身份映射和 P2CRC_Xenium 联合验收范围执行；表达/坐标来源能由现有材料证明的由 Agent 核对，不能证明时才补材料。无需要求用户先安装新平台或重做原始数据。

**下一步与证据。** fixture 工程交接已记录在[验收清单](../acceptance.md)、[验证记录](../verification-2026-09-19.json)和[分析交接计划](../plans/analysis-handoff.md)；后续使用 P2CRC_Xenium 检查真实 loader、标签/空间步骤和实际 `status/error`。REVISE 当前交接实现见 [batch runner](../../../../revise/batch/runner.py)；它的 reconstruction manifest 不是分析库当前入口。

<a id="c1"></a>
## C1：分析库消费协议的独立协同

**为什么做。** sST 保留算法内部 normalization，结果不是 raw counts，也未 log。分析端不能靠错误的 scale 声明自动进入 normalize/log 流程。

**最新协议。** 当前 Analysis Agent 固定 finite、nonnegative、unlogged linear X，允许小数和可信 normalized 非 log 表达。旧 log 声明已拒绝；来源未知仍不可用。生产端只对明确声明来源的数据放行，不由数值猜测历史。

**责任已确认。** 本批更新 source declaration、sample.yaml 和工程联测，不修改消费者。两库版本须记录；真实科学验收后置。

**依赖与用户参与。** 本批暂无必交材料。源处理历史不明的输入保留 unknown，不要求现在整理全部真实样本。

**下一步与证据。** 已按[OT assembly 配套计划](../plans/ot-assembly.md)完成 producer 更新；实现位置和实际测试见[落实索引](../implementation-index.md)。

<a id="existing-assembly"></a>
## 已有能力：mean/random assembly

现有 [assemble_ist](../../../../revise/application/ist_assembly.py) 支持 cluster mean 与按 seed 选择 donor，使用 expression carrier 的 X，空间轴来自 spatial carrier。已有[针对测试](../../../../tests/application/test_ist_publication.py)，已纳入本批回归。R2 复用该能力并完成多类型合并；无需用户补材料，后续真实样本验收见[验收清单](../acceptance.md)。

<a id="backend-independence"></a>
## 已有能力：分析库 backend 解耦

此前核对分析库包依赖、源码与迁移说明，已经不把 REVISE backend 作为运行依赖，因此历史 A22 不重复立为新开发。这样保持分析库从两个独立对象出发。R4 已完成新交付的小型正式消费验证，真实样本后置；暂无新增用户材料需求，见[证据与边界](../decisions-and-sources.md#evidence)。

<a id="no-gene-overlap"></a>
## 明确不新增：普通 gene-overlap 元数据

Raw/SVC 的共同、实测或重建独有基因可从各自 var_names 低成本推导，无须为了协议完整新增 measured_in_raw 等字段。它不依赖额外材料，也不是要删除现有任意字段的授权。后续只在[昂贵 gene-wise uncertainty](04-future-extensions.md#d2)确有定义时另行讨论保存；当前不新增任务或实现计划。
