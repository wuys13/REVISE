# REVISE × Analysis Agent：GPT 双库审阅入口

[实施索引](implementation-index.md) · [当前验收](acceptance.md) · [Notebook 比较约定](assembly-comparison-contract.md)

## 审阅目标与边界

审阅 REVISE **当前工作树**生成的交付能否被独立 `REVISE_Analysis_Agent` 当前工作树正确消费。两库均可能有未提交修改；不能仅查看 HEAD，也不能把旧批次代码当成本次对象。版本、文件摘要和浅验收记录见 [证据索引](verification-review-handoff.json)。

REVISE 负责 reference 准备、重建、表达来源与标签事实、同次事务发布；分析库负责加载、能力判断、分析和报告。准备与重建分别运行，正式交接入口是 `sample.yaml + raw.h5ad + SVC.h5ad`，不是内部 spatial/expr 双载体。此轮不修改分析库、不运行真实全量重建、不产生科学优劣结论。


## 当前结论与证据时效

核心交付接口已有工程打通证据；尚未证明所有当前路由、真实数据和下游分析均完成端到端验收。文档就绪、产物加载、分析成功、科学结论是四个不同完成层级。

本次只更新文档。当前源码身份见 [本次协同源码快照](coordination-source-snapshot.json)：REVISE HEAD 为 `a26d36fe3b11d652f2cd9eeafb92e990518eb3ad`，有大量既有未提交修改；消费者 HEAD 为 `4118061983f15d0a51b6fdac8c9e4f534cd2c814`，状态为未跟踪 `other_reviews/`。上一轮消费者 HEAD 为 `291cf4f06aad187af0ccbf5bb465bf3b03f4ab88`。本次摘要覆盖直接审视的接口文件；GPT 仍须检查两库全部当前工作树和新增文件，不能只审提交差异。旧证据索引保留，不覆写为新通过结果。

## 协同任务与责任清单

| 任务 | 当前状态／触发条件 | 责任 | 预期行为与验收 |
|---|---|---|---|
| H1 三件套和状态核对 | 有旧版小型联测；消费者已更新 | 双方，REVISE 提供产物，消费者解释状态 | 用正式发布产物调用真实 loader/runner，核对路径、标签、来源、状态和输入哈希；参考 `tests/integration/test_analysis_delivery.py` 的两个测试，不用手写 YAML 替代生产输出 |
| H2 历史基线汇集 | 三份来源已确认；单一基线未生成 | REVISE | 按[比较约定](assembly-comparison-contract.md#后续执行任务与验收)汇集标签与坐标，拒绝重复 ID，记录三份来源及摘要；旧快照预期 49,279 行，变化时先解释来源，不硬截断 |
| H3 四方法真实产物 | 四模式已有实现；真实结果未生成 | REVISE | 固定输入/reference/上游设置/seed，独立输出目录；保存有效配置和文件摘要，确认只改变 assembly 及相应输出位置 |
| H4 比较前校验 | ID/基因交集已有；方法坐标对照尚未实现 | REVISE | 绘图前按 ID 核对四方法与基线坐标；报告类型覆盖、标签路由差异、共同基因和排除；失败不得用基线坐标掩盖差异 |
| H5 无 cluster 路由 | hST/sST 正常不产生 SVC_cluster | 消费者核查，REVISE 提供真实字段 fixture | 独立表达分析不应被缺少 cluster 阻断；依赖它的步骤明确不适用/能力不可用，无异常堆栈；不要求生产补标签 |
| H6 ID/parent 与物理坐标 | iST ID 与 sST parent 语义不同 | 双方 | 仅在指标要求时对齐；禁止用虚拟细胞 ID 猜 spot 一一配对；缺物理尺度时限制相关步骤，不猜单位或坐标变换 |

H2–H4 是已明确待执行任务，不是待用户提供标签或重新确定算法。本轮不生成基线、真实方法产物或增加校验代码。H1/H5/H6 是给双库审阅者的核查任务，不预判消费者存在缺陷。在线发现、D2、D3、posterior 不因此进入本轮。

## 如何解释消费者结果

以当前 `reconstruction_impact.py` 的 `run_stage/result` 为依据：

- 阶段 `completed` 表示该阶段完成；新增缺失能力记录时，有产物为 `partial`，无产物为 `unavailable`。
- 阶段异常记录为 `error` 并进入 `stage_errors`；总体有阶段异常则为 `failed`。
- 总体有缺失能力或仍有 pending 阶段为 `partial`；无这些情况且有输出才为 `succeeded`，否则为 `skipped`。
- 因而总体 `partial` 不等于计算失败，也不能宣称所有要求都通过；需要同时看 `stages`、`unavailable`、`stage_errors` 和具体产物。其他分析由各自结果实现定义，不能机械套用此状态机。

对 hST/sST，缺 cluster 是预期不适用。审阅应检查限制是否局部生效，而非要求总体 Impact 强行 succeeded。对生产端，eligible 计算失败仍失败且事务回滚，不能借消费者 partial 口径发布半个样本。

## 对接时的最小操作顺序

1. 记录双方 HEAD、未提交与未跟踪文件及关键文件摘要；先识别与历史验收快照的差异。
2. 使用本页小型验收命令在新目录运行，明确设置 `REVISE_ANALYSIS_ROOT` 指向被审阅的消费者工作树；缺依赖或跳过不能算当前联合通过。
3. 核查 H1/H5/H6，分别提交已确认缺陷、未验证风险和正常能力限制，注明归属；消费者修改由其仓库独立落实。
4. 后续真实执行先完成 H2/H3，再通过 H4，最后运行 Notebook。真实执行另行安排，不能把本次文档交付当作完成。

## 当前协议对照

| 接口／含义 | REVISE 生产行为 | Analysis Agent 消费行为 | 结论与验证边界 |
|---|---|---|---|
| 文件与相对路径 | schema_version=1，files.raw/svc 与 YAML 同目录 | `io.read_sample_config/load_sample` 按 YAML 所在目录解析 | 小 fixture 正式入口通过；跨 cwd 有测试 |
| 样本 ID | 原始 ID percent 编码为单目录段，provenance 留原值 | `io.safe_segment` 拒绝路径段与不安全名称 | 已有分层 ID fixture；不是把原 ID 末级截断 |
| Broad/subtype | 两侧实际存在才声明 revise_Level1/2；sST 使用已有 cell_type 作 broad 别名来源，原列保留 | 使用 columns 映射，缺失能力由具体步骤处理 | 不要求 Raw 获得伪造 subtype；消费者默认列名不是列存在证据 |
| 重建主标签 | SVC 实际有 SVC_cluster 才声明 columns.reconstruction | Impact 使用 reconstruction_key 构建重建状态；不以 Level2 替代 | iST 对齐；hST/sST 本来不涉及该标签，不列为生产缺陷；仅相关分析能力不适用 |
| 表达 | 显式来源支持有限非负非 log `.X`；允许小数；unknown 不升级 | 当前固定线性 X；identity unknown 或旧 scale unknown 限制表达 | 已知线性 fixture 正式 Moran 分析；unknown 标签空间路径可 partial |
| 坐标 | 校验两侧坐标形状／有限；sST 已补将真实生成 x/y 写入 obsm.spatial；um 转 micron，未知不猜比例 | 各分析按物理单位／比例判断能力 | 有尺度的小 fixture 已验证；未知真实输入不保证物理分析可用 |
| 轴与关系 | 原 Raw 保持；iST 真实 ID，sST 真实 parent 信息；不强造通用配对 | 默认 native objects 独立；成员关系按明确单位交集处理 | 简单 ID 已验证；sST 一对多不等于已支持全部 membership 分析 |
| 失败与输入保护 | 事务发布与回滚，eligible 计算失败不交付半样本 | 输出 status/stage_errors/unavailable；不改输入 H5AD | 失败注入及输入哈希测试；partial 不是全成功 |
| Confidence | 原最大注释权重；reference 共用数值，阶段来源独立 | 不是所有分析必需的统一科学置信度 | 未建立新 posterior 消费协议，不能跨阶段解释为同一概率 |
| 四方法比较 | 独立 Notebook 重跑表达 Leiden，对照既有 SVC_cluster | 交接说明供核对／后续接入，不要求消费者自动运行该 Notebook | fixture 验证；真实方法优劣及默认值未定 |

生产源码：[delivery](../../../revise/application/delivery.py)、[表达声明](../../../revise/application/expression.py)、[发布事务](../../../revise/application/publication.py)、[OT assembly](../../../revise/application/ist_ot.py)。
消费者源码（相邻仓库）：`revise_analysis/io.py`、`revise_analysis/runner.py`、`revise_analysis/analyses/reconstruction_impact.py`、`revise_analysis/analyses/_shared.py`；协议说明 `docs/input-output.md`。具体文件摘要见证据索引。

## 本轮发现与处理

- **已修复：sST 输出坐标未进入标准槽位。** runner 先前只输出 obs.x/y，正式发布要求 obsm.spatial，导致生成后发布失败。现将同一坐标按原行序复制到 obsm.spatial；不插值、不移动细胞。
- **已修复：sST broad 标签未声明。** 虚拟细胞已有 cell_type；现在以它生成 revise_Level1，保留原字段和值，不从 spot 多数标签替代细胞标签。既有 `/`→`_` 命名处理未改变；跨侧别名可按真实数据显式配置，不猜测统一。
- **已确认的产品边界：** hST/sST 本来不涉及 SVC_cluster，不要求补列、新增分群或实现 iST 的 cluster 状态分析。消费者只限制确实依赖该标签的步骤；其他满足表达、坐标等条件的分析应可独立执行。这是正常设计，不是待补功能。
- **尚未证明：** 坐标有限和 shape 正确不证明任意用户输入均处于同一物理坐标框架。真实跨侧分析须核对来源／单位／原点；本轮不新增通用坐标变换。
- **实际基线位置修正：** raw_data/Real_application/P2CRC_Xenium.h5ad 无 SVC_cluster；已找到三个历史类型 spatial.h5ad，见[元数据记录](assembly-comparison-contract.md)。用户已确认基线指这些重建后文件。标签材料已存在；汇集为单一基线与生成四方法产物属于后续真实运行准备，不要求用户补标签。

## 小型验收与复现

从 REVISE 根目录，使用已有科学 Python 环境执行：

```bash
/Users/stephen/miniconda3/envs/python3.10/bin/python scripts/verify_review_handoff.py --output output/review-handoff/new-review-run
```

脚本使用独立新目录，保留测试日志/JUnit、正式消费者结果和执行后的比较 Notebook。完整调用、实际数量和产物路径在证据索引；已有目录不覆盖。依赖独立消费者、真实 TACCO 和可用 Jupyter kernel。缺失依赖或 skip 不能算联合通过。

此前完整回归 964 项与 24 项子测试见 [上一批记录](verification-ot-assembly-2026-09-19.json)，本轮只刷新直接相关检查，不把两次数量累计。

## 给 GPT 的审阅任务

请同时读取两库当前工作树，并先核对证据索引与当前版本是否一致。沿生产配置→实际 H5AD/YAML→消费者 loader→具体分析检查，不仅比较文档或字段名称。

重点：主标签 SVC_cluster 与 broad/subtype 的区别；来源 unknown 与线性声明；原始 Raw 保护；轴、缺失与关系；坐标物理尺度；部分可用状态；事务和复用；Notebook 的表达分群是否独立于 Ground Truth。可复现的小型产物优先于推测。

每个发现请给出严重度、双方文件位置、实际触发条件、影响、复现证据和建议归属（生产端／消费端／共同约定）。分别列出已确认缺陷、未验证风险、后置科学问题。不要要求把明确暂缓的在线来源、D2、D3 或完整 posterior 扩展提前实现。

真实 T/Macro/CAF 对照与 sST 真实消费尚未执行；Notebook 的 Ground Truth 是用户指定的既有重建标签，不额外宣称为独立生物学真值。


## 本次实际验收结果

以下为此前 `20260919-final` 快照的结果，并非本次文档更新后重新执行。分析库 HEAD 已推进，当前匹配性须按下列审阅任务刷新，不能沿用历史通过结论覆盖新版本。

最终 **75 passed、0 failed、0 skipped**，退出码 0；另保存执行后的多分辨率 Notebook、4 张 PNG、指标／覆盖／列联表。正式消费者标签空间 fixture 为 `partial`（表达来源 unknown），已知线性 fixture 的 spatial_autocorrelation 为 `succeeded`。未修改输入 H5AD 或源 Notebook，未修改消费者代码。

- [执行 Notebook 的 HTML](../../../output/review-handoff/20260919-final/comparison.html) · [已执行 ipynb](../../../output/review-handoff/20260919-final/comparison.executed.ipynb) · [示例空间图](../../../output/review-handoff/20260919-final/comparison-1.png)
- [指标](../../../output/review-handoff/20260919-final/metrics.csv) · [覆盖表](../../../output/review-handoff/20260919-final/coverage.csv) · [运行摘要与产物哈希](../../../output/review-handoff/20260919-final/summary.json)
- [消费者标签空间报告](../../../output/review-handoff/20260919-final/fixtures/test_published_sample_is_consu0/analysis/synthetic%2Fwhole-sample/reconstruction_impact/report.html)
- [复现脚本](../../../scripts/verify_review_handoff.py) · [真实比较配置模板](../../../examples/assembly-comparison-real.json)

输出目录被 Git 忽略，避免把 fixture 二进制与结果当成源代码；本目录的证据索引、复现脚本与配置模板可纳入版本管理。迁移机器时需同时携带输出目录，或在新目录重跑脚本。所有示例图均为合成数据，不是 T/Macro/CAF 的真实科学结果。
