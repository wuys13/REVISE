# Batch reconstruction and analysis: implementation decisions

## Goal and deliverables

在现有标准 ST 目录和分层配置上，让单项与批量调用共用执行逻辑，正确提供 raw 与重建分析输入，并支持按分析方面执行、比较、发布和恢复。

**本文只负责范围、开发顺序和验收入口。** 最终确定的数据约定统一在[协议](protocol.md)，执行约定统一在[框架](framework.md)；两者维护现行实现；本计划保留取舍与实施依据，完成情况集中在验收文档，不再另写第二套协议。[返回总入口](README.md#contract-and-status)。

本轮交付五项能力：

1. 稳定的单项文件调用接口。
2. hST／iST／sST 的最小分析输入辅助。
3. 一个方面一个完整任务的分析接入。
4. 最少但充分的结果与计算条件记录。
5. 分开的重建／分析复用判断。

## Scope and decisions

- 保留标准输入目录、分层 `batch.yaml`、共享 reference 和现有输出布局；不复制共享数据，不新增输入清单。
- iST 本轮仍按单类型执行；支持原生 `paired`，`mean`／`random` 继续正常读取。先运行 paired 检查结果，以后可用独立输出目录保存 random 结果进行比较。
- iST raw 使用原始实测 panel，重建使用自身基因范围；不强制共同 panel，缺测不补零。
- sST raw 按实际输出的父 spot 映射与生成 cell 数守恒均分，构造分析基线；不改变重建算法。
- 每个分析方面内部协调两侧计算与比较；只运行明确配置的方面，不固定五类 API。
- 保存实际结果、必要空间证据及计算条件；不强制保存所有中间矩阵、理论候选组合或大型证据包。
- 本轮不实现 EMT、Moran、CCI 等科学算法，不建设资源注册表、跨方面调度或分布式执行。CCI 前后比较保留为未来模块能力。
- 外部项目适配 REVISE 的文件接口；不反向绑定 Atlas 的 API 命名、数据组织或调度系统。

## Implementation order

| 顺序 | 工作 | 完成标准 | 详细阅读 |
|---|---|---|---|
| 1 | 抽取单任务执行边界 | 单项不误跑其他任务，与批量共用实现 | [调用接口](framework.md#task-api) |
| 2 | 提供分析输入辅助 | paired 对应正确，sST 均分守恒，输入不变 | [分析输入](protocol.md#carriers-and-alignment) |
| 3 | 扩展分析接入与结果交接 | 模块获得正确输入，产物可读取并能解释比较依据 | [模块接入](framework.md#aspect-modules)、[结果交接](protocol.md#analysis-results) |
| 4 | 分离复用判断 | 分析变化不触发重建，无变化结果可复用 | [复用身份](framework.md#reuse) |
| 5 | 完成验收和文档更新 | 回归、调用示例与文档引用一致 | [验收场景](plan-validation.md#acceptance-scenarios) |

每一步都先补充能够区分正确与错误行为的测试，再接入执行链。框架完成意味着调用、数据语义和恢复正确，不等同于具体科学算法已完成或生物学结论有效。

## Reference decisions

Atlas 提供的是消费方需求参考，来源入口为 [development-reference/README.md](../../../../ReviseSTAtlas/preparation/external/revise-api/development-reference/README.md)。该材料记录的 REVISE 来源快照为 `93c20e3bbdb055e9bc296a6db8a6536e7d44fe85`，不是本计划的实现状态。外部链接假定两个仓库并列检出；不复制其契约。

| 参考方向 | 本计划的判断 |
|---|---|
| [批量与单项调用](../../../../ReviseSTAtlas/preparation/external/revise-api/development-reference/batch-interaction.md) | 采纳精确选择任务、失败恢复和可消费记录；接口沿用 REVISE 的配置与文件边界，不采用另一套内存对象入口。 |
| [五类分析与比较](../../../../ReviseSTAtlas/preparation/external/revise-api/development-reference/analysis.md) | 采纳有实际结果的两侧比较及空间定位；方面不固定为五个函数。CCI 比较待算法分支接入，partition／diversity／regions 不机械改名为 niche。 |
| [验证与交付](../../../../ReviseSTAtlas/preparation/external/revise-api/development-reference/validation-and-delivery.md) | 采纳可重载结果、计算条件和小样本对照；暂缓完整候选全集、全部中间矩阵和独立多层报告体系。 |
| [来源与待讨论事项](../../../../ReviseSTAtlas/preparation/external/revise-api/development-reference/sources-and-decisions.md) | 记录影响解释的来源与真实资源摘要；不引入资源注册表或外部冻结 API。 |

只有标签和坐标参与的 Neighbor／Niche，在两者都不变时不应凭重建表达宣称发生变化；表达驱动的亚型重注释、对象匹配与统计前提由具体模块说明。paired 的 donor 均值视图可用于分析，但不证明真实单细胞异质性。重建影响分析分支 `codex/reconstruction-impact-analysis`（参考快照 `3bd0934`）提供背景，算法不在本轮迁入。

## Question index

| 想解决的问题 | 阅读位置 |
|---|---|
| 现在已经支持什么 | [实现状态](README.md#contract-and-status)、[当前执行链](framework.md#execution-and-module-boundaries) |
| 本轮增加什么、不做什么 | [目标](#goal-and-deliverables)、[边界](#scope-and-decisions) |
| paired 怎样参与分析 | [paired 视图](protocol.md#paired-view) |
| sST raw 怎样构造 | [守恒基线](protocol.md#sst-baseline) |
| 实现或外部数据适配参照什么样式 | [独立输入输出示例](input-output-example.md#input-tree) |
| 输入输出目录和配置由谁定义 | [配置继承](protocol.md#configuration-and-inheritance)、[输出角色](protocol.md#output-roles-and-current-result-semantics) |
| 单项、批量分别怎样调用 | [调用接口与执行边界](framework.md#task-api) |
| 新分析模块需要提供什么 | [模块接入](framework.md#aspect-modules) |
| 哪些结果必须保存 | [结果交接](protocol.md#analysis-results) |
| 配置变化后重跑哪些任务 | [复用](framework.md#reuse)、[失败与恢复](framework.md#recovery) |
| 如何证明实现正确 | [验收场景与测试对应](plan-validation.md#acceptance-scenarios)、[交付检查](plan-validation.md#delivery-checks) |
| Atlas 建议的取舍与来源 | [参考取舍](#reference-decisions) |
