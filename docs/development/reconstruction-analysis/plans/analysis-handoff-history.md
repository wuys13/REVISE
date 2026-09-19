> 历史计划快照，旧 scale 支持描述已失效。当前说明见 [交接计划](analysis-handoff.md)。

# R4：分析交接与联合验收计划

> **协议更新（2026-09-19 OT assembly 批次）：** 当前消费者已固定线性非负、未取 log 的 X；旧 log 声明已拒绝。下文此前的表达 scale 支持范围属于历史快照，现行声明和本批处理以[新计划](../plans/ot-assembly.md)为准。真实验收继续后置。
**状态：生产端实现及 fixture 工程消费验收完成；真实样本联合验收待完成。** 用户已批准与 R2/R3 合并执行；实施状态统一见[总览](../README.md#items)，批次级边界和剩余真实输入入口见[第二批总计划](whole-sample-delivery.md)。

[事项 R4](../categories/02-sample-delivery.md#r4) · [协同 C1](../categories/02-sample-delivery.md#c1) · [共同决策](../decisions-and-sources.md#confirmed-decisions)

## 目标、责任与依赖

将 R2/R3 的同次 `raw.h5ad + SVC.h5ad` 配上分析库可读取的 `sample.yaml`，并验证真实消费路径。文件存在、能加载、完成具体分析是不同验收层次。

R4 只负责 REVISE 生产交接与联合验收安排，不重复实现 R1 的 sST 算法修改或 R3 的 Raw 发布。分析库消费调整由 C1 独立推进；不修改其并行工作树，不替它决定表达分析方法。

完整文件交付依赖 R2/R3；sST 表达计算的完成额外依赖 C1。缺少某项表达能力不应阻断与表达无关的标签/空间检查，也不应把受阻表达分析报成成功。

当前已完成合成 fixture 的消费者工程交接：正式入口完成 loader、标签/空间 diversity 和报告生成，输入哈希未改变；由于表达 identity unknown，结果整体为 `partial`。生产端全样本交付、metadata/alias 和 batch/direct 路径已纳入最终工程套件。这不等于 P2CRC_Xenium 真实样本验收，也不等于 sST 表达分析已获得支持。

## 当前接口与剩余差异

2026-09-19 静态核对：REVISE batch 仍按 modality/cell type 记录 reconstruction manifest，见 [runner](../../../../revise/batch/runner.py) 和 [inputs](../../../../revise/batch/inputs.py)。独立分析库不用该内部 loader，而使用 `load_sample(sample_yaml)`、`run_analysis(sample_yaml, analysis, output_dir, parameters)`。

分析库当前未提交工作树已支持 Raw/SVC 分侧 `expression.identity/scale/matrix`，撤掉全局未归一化门槛；矩阵缺失也可加载标签/空间对象。表达消费者仍主要支持 untransformed 或 log1p。此处是工作树快照，不是本任务已验证的稳定发布契约；交接前必须重新核对。

生产端三种 route 不能都贴相同表达声明：iST assembly 使用 donor carrier 的 X，不主动归一化；sST 保留内部 normalization 与 spot 校正，输出非 log、非 raw counts。其他 route/输入历史也需据实际来源声明，不能由非负数值推断尺度。

## 配置与接口变化

以下字段来自当前消费者，配置生成是本计划拟新增的 producer 能力。未知保留未知；新 scale 名称或额外元数据不擅自升级为共同协议。

| 内容 | 本批交接行为 | 实施验证限制 |
|---|---|---|
| `schema_version`、`sample_id` | 使用消费者支持的版本和样本身份；层级 batch ID 编码为确定、可逆的单目录段并保留原始映射 | 映射实现需用碰撞 fixture 验证；不简单取末级目录 |
| `files.raw`、`files.svc` | 指向同次成功发布的两个对象，路径相对样本 YAML | 不扫描文件名猜结果，不要求消费者读内部 carriers |
| `expression.raw` / `expression.svc` | 分别声明 `matrix: X`、真实来源 identity、已知 scale | 当前消费者只支持 X；unknown 可加载但相关表达计算不可用 |
| 表达 scale | 已支持的 untransformed/log1p 必须有来源证据 | normalized nonlog 仍是 C1 缺口；`normalized_nonnegative` 只能作候选名称，不能宣称已支持 |
| `columns.broad/subtype/reconstruction` | broad/subtype 指向两侧实际存在的 `revise_Level1`/`revise_Level2`；reconstruction 只指向真实存在的重建标签 | `SVC_cluster` 不能为 hST/sST 伪造；缺失按分析需求报告 |
| `label_aliases` | 必要时使用已确认的标签映射 | 不悄悄改写 H5AD 或合并不同标签 |
| `spatial.key/unit/microns_per_coordinate` | 采用实际输入配置或数据来源 | producer 的 `um` 与消费者模板 `micron` 需明确等价映射；单位字符串转换不等于坐标重投影，比例不能编造 |

Raw 回填重建阶段已经得到的 annotation 由 R3 负责，这是合法的推断交付；R4 不再次运行 OT，也不从不对应的 SVC 行猜 Raw 标签。已有 spot_name/cell_id 留在相应对象内，不额外构造 universal correspondence。

## 实施顺序（生产端与 fixture 已完成，真实样本独立进行）

1. **固定交接依据。** 核对两库版本及相关工作树差异、消费者接口和已知样本，确定本次联合验收针对哪个状态。
2. **检查 R2/R3 产品。** 核对同次运行身份、Raw 原始矩阵/轴、SVC 坐标和实际标签；发现上下游结果不一致交回对应任务，不在配置生成时修补数据。
3. **生成轻量配置。** 采用对齐后的样本 ID 规则，将路径、标签和真实表达/坐标语义写入 sample.yaml。配置作为同次交付的一部分发布，旧配置不能搭配新文件冒充一致。
4. **调用真实 loader。** 验证两个独立轴、相对路径、列名和 unknown/缺失的读取行为；无需消费者导入 revise。
5. **完成实际分析路径。** 在已满足前提的代表样本上调用正式分析入口，检查结果记录及实际产物。sST normalized nonlog 未被 C1 支持时，只验收已具备前提的步骤，明确保留未完成项。
6. **保存联合证据。** 给出调用、样本配置、两库相关版本、结果及缺口。不同 route 分开报告，不能以某一路成功代表全部已完成。

工程证据见[验证记录](../verification-2026-09-19.json)和[验收清单](../acceptance.md)。真实样本命令只验证实际可消费的标签/空间路径；当坐标物理比例或表达 identity 缺少来源时，必须检查消费者返回的实际 `status`/`error`，不预先承诺所有空间或表达步骤只会返回 `partial`。

## 发布失败与结果边界

生产端发布仍遵守 R2/R3 的样本级成功/失败规则，不能因分析端某项不可用就伪造一个新的“重建 partial”契约。

分析端沿用其真实状态：`succeeded`、`partial`、`skipped`、`failed`。某项请求部分未完成时，不能只凭 loader 成功称为联合验收通过。预期缺前提与程序异常分开报告。

| 情况 | 处理与验收含义 |
|---|---|
| Raw 被归一化工作对象污染、文件混用、ID/轴错误、伪造标签或尺度 | 交接失败，修复生产数据或配置，不进入成功结论 |
| 缺少 Level2 或未支持某种表达 scale | 对应分析部分明确 unavailable；有前提的其他部分可继续 |
| 缺物理比例 | 只阻塞需要物理尺度的窗口等分析，不一概阻止加载或所有空间计算 |
| 单位/样本身份存在歧义 | 先补依据或对齐映射；不猜测、静默合并或重命名 |
| 消费者执行异常 | 记录失败证据，不包装成缺前提或成功 |

## 针对验证与验收

- Raw/SVC 行数、顺序、ID 和基因集合不同仍可独立读取；仅具体比较按真实关系局部对齐。
- 配置相对路径在从其他工作目录调用时仍可解析；样本身份映射不碰撞，原始身份可追溯。
- 不存在的 Level2/SVC_cluster 不伪造；与缺失无关的分析不被全局拒绝。
- sST 不被错误声明为 untransformed 或 log1p；C1 未完成时明确对应表达分析缺口。
- 真实调用至少包含 loader 及一个已满足前提的正式分析步骤，结果文件与状态一致。
- 分析前后输入 H5AD 不被改写；合成 fixture 只证明工程组合，真实样本的语义与科学解释另行核对。

以上 fixture、相对路径、双轴加载、输入哈希和正式 consumer 入口已通过工程套件；真实 P2CRC_Xenium 和 sST normalized nonlog 表达仍未完成。

## 用户参与与重新讨论条件

层级 sample ID 的映射、两侧标签列的交接方式、代表样本与联合验收范围已纳入本批执行基线：使用可逆单段映射、原始列保留且推断写入 `revise_Level1`/`revise_Level2`、优先使用现有 P2CRC_Xenium。normalized nonlog 的具体方法学处理仍由用户与 C1 所在分析任务讨论，本计划不预设只 log 或再次 normalize。

必需证据是样本的真实表达和坐标来源；Agent 先查现有配置和资料，只有无法查明时再请用户补充。不要将“提供全部样本元数据”泛化为用户的新任务。

旧遍历代码不是 R4 的前提；confidence 新方案不阻塞结构交接。科学含义、用户接受标准、消费者接口发生变化时，暂停受影响验收并更新本计划；不扩大为未授权的重建算法改造。

## C1 交接内容

生产端提供实际 route、矩阵处理经过、标签/坐标可用性和一组代表性 Raw/SVC；分析端给出各分析的支持条件、变换策略及结果证据。对方修改在本库只登记，不实施。细节和后续讨论入口见[协同事项](../categories/02-sample-delivery.md#c1)，本任务不建立另一套分析计划或通用兼容框架。
