# REVISE × Analysis Agent：双库协同入口

[本轮完整验收与服务器命令](acceptance.md) · [实施索引](implementation-index.md) · [比较协议](assembly-comparison-contract.md) · [本轮机器记录](verification-mini-2026-09-20.json)

2026-09-20 已完成双仓 mini 联测及后续收口。REVISE定义并交付SVC；Analysis Agent按字段消费统一对象。正式入口是同次发布的 `sample.yaml + raw.h5ad + SVC.h5ad`。实际完成重建的类型合并；未重建类型不并入。

## 当前结论与证据时效

本轮 Anatomy 改由完整交付 SVC 的 broad 和坐标定义，保留 Raw 共享原点；Notebook 默认正式 P2 project。sST改为稳定的share-then-target校正，571个旧有效零项恢复，剩余250个校正前真零项明确unresolved。三份Notebook与batch的172/14/45份科学表JSON一致，iST的44份State/Gain及窗口表未变，12份报告只读刷新。

P1 HD表达保持unknown；Visium0.73尺度保持暂定；sST真零支持仍partial。iST四assembly、H2/H4和未受影响重建复用已有证据，不宣称全量或科学验收。新结果位于`output/mini-acceptance/20260920-closeout`，旧产物保留。

本轮基线HEAD为REVISE `464ee26ac22978c3994e40af50325413e70cef13`、消费者 `3e70ab364f89f9e24e981e3e9239e7fe563138a6`。当前验收包含未提交修改；用户审阅未改。机器记录原字段为上一阶段快照，当前看`closeout`。详细修改、测试、真实产物及剩余事项见[验收表](acceptance.md)。

## 协同任务与责任清单

| 任务 | 本轮状态 | 责任与边界 |
|---|---|---|
| H1 正式交付/消费 | 三路线真实mini已联通 | 生产发布器直接交付，消费者按真实前提记录partial/skipped，不能伪造全成功 |
| H2 历史baseline | 49,279×0；源ID/cluster/坐标保留 | REVISE；三个carrier各自broad，不加载历史X |
| H3 四assembly | 六次真实重建中的四个iST run已完成 | REVISE；同ROI/reference/QC/GA/LR/seed，仅assembly与输出位置不同 |
| H4 common scope/坐标 | 三type共同ID 160/157/68；全部共同坐标一致 | REVISE；异常只阻止对应类型可比较结论，保留排除范围 |
| H5 无SVC_cluster | hST/sST实际验证，K-control开启时unavailable | Analysis；sST原生分子分析仍运行，hST另受unknown表达身份限制 |
| H6 parent/物理坐标 | iST/hST原ID原坐标、sST真实parent继承坐标 | 双方；sST无同单位membership，Visium尺度0.73为暂定 |

## 如何解释消费者结果

以当前 `reconstruction_impact.py` 的 `run_stage/result` 为依据：

- 阶段 `completed` 表示该阶段完成；新增缺失能力记录时，有产物为 `partial`，无产物为 `unavailable`。
- 阶段异常记录为 `error` 并进入 `stage_errors`；总体有阶段异常则为 `failed`。
- 总体有缺失能力或仍有 pending 阶段为 `partial`；无这些情况且有输出才为 `succeeded`，否则为 `skipped`。
- 因而总体 `partial` 不等于计算失败，也不能宣称所有要求都通过；需要同时看 `stages`、`unavailable`、`stage_errors` 和具体产物。其他分析由各自结果实现定义，不能机械套用此状态机。

对 hST/sST，缺 cluster 是预期不适用。审阅应检查限制是否局部生效，而非要求总体 Impact 强行 succeeded。对生产端，eligible 计算失败仍失败且事务回滚，不能借消费者 partial 口径发布半个样本。

## 当前协议对照

| 接口／含义 | REVISE 生产行为 | Analysis Agent 消费行为 | 结论与验证边界 |
|---|---|---|---|
| 文件与相对路径 | schema_version=1，files.raw/svc 与 YAML 同目录 | `io.read_sample_config/load_sample` 按 YAML 所在目录解析 | 三路线真实 mini 正式入口通过；跨 cwd 有测试 |
| 样本 ID | 原始 ID percent 编码为单目录段，provenance 留原值 | `io.safe_segment` 拒绝路径段与不安全名称 | 已有分层 ID fixture；不是把原 ID 末级截断 |
| Broad/subtype | 两侧实际存在才声明 revise_Level1/2；sST 使用已有 cell_type 作 broad 别名来源，原列保留 | 使用 columns 映射，缺失能力由具体步骤处理 | 不要求 Raw 获得伪造 subtype；消费者默认列名不是列存在证据 |
| 重建主标签 | SVC 实际有 SVC_cluster 才声明 columns.reconstruction | Impact 使用 reconstruction_key 构建重建状态；不以 Level2 替代 | iST 对齐；hST/sST 本来不涉及该标签，不列为生产缺陷；仅相关分析能力不适用 |
| 表达 | 显式来源支持有限非负非 log `.X`；允许小数；unknown 不升级 | 当前固定线性 X；identity unknown 或旧 scale unknown 限制表达 | P2 iST/sST 真实表达消费；P1 HD unknown 明确限制表达 |
| 坐标 | 校验两侧坐标形状／有限；sST 已补将真实生成 x/y 写入 obsm.spatial；um 转 micron，未知不猜比例 | 各分析按物理单位／比例判断能力 | 有尺度的小 fixture 已验证；未知真实输入不保证物理分析可用 |
| 轴与关系 | 原 Raw 保持；iST 真实 ID，sST 真实 parent 信息；不强造通用配对 | 默认 native objects 独立；成员关系按明确单位交集处理 | 真实 mini ID/parent已验证；sST禁用同单位membership |
| Anatomy 来源 | 实际重建完成的类型合并为 SVC；未重建类型不并入 | 完整 SVC broad 与坐标定义 Anatomy；保留共享 Raw 原点，再映射各侧观测点 | Other 仅描述已交付 SVC 窗口未观察到 Tumor/Normal；无 SVC Anatomy 网格为 Unknown，不代表原组织缺少这些类型 |
| 失败与输入保护 | 事务发布与回滚，eligible 计算失败不交付半样本 | 输出 status/stage_errors/unavailable；不改输入 H5AD | 失败注入及输入哈希测试；partial 不是全成功 |
| Confidence | 原最大注释权重；reference 共用数值，阶段来源独立 | 不是所有分析必需的统一科学置信度 | 未建立新 posterior 消费协议，不能跨阶段解释为同一概率 |
| 四方法比较 | 独立 Notebook 重跑表达 Leiden，对照既有 SVC_cluster | 交接说明供核对／后续接入，不要求消费者自动运行该 Notebook | 真实 mini 四方法已比较；科学解释和最终默认值未定 |

生产源码：[delivery](../../../revise/application/delivery.py)、[表达声明](../../../revise/application/expression.py)、[发布事务](../../../revise/application/publication.py)、[OT assembly](../../../revise/application/ist_ot.py)。
消费者源码（相邻仓库）：`revise_analysis/io.py`、`revise_analysis/runner.py`、`revise_analysis/analyses/reconstruction_impact.py`、`revise_analysis/analyses/_shared.py`；协议说明 `docs/input-output.md`。具体文件摘要见证据索引。

## 审阅与复现

从[验收表](acceptance.md)进入真实配置、命令、产物与剩余事项。核查顺序为配置/来源 → 实际H5AD/YAML → loader能力 → Notebook/batch → 保存结果/报告。重点核对缺失分母、unknown表达、主cluster身份、独立轴、物理坐标和失败传播。

工程fixture仍可用 `python scripts/verify_review_handoff.py --output output/review-handoff/<new-run>`，需将 `REVISE_ANALYSIS_ROOT` 指向当前消费者并使用同一科学kernel；旧结果不覆盖。当前真实验证不等同于全量容量、科学改善或默认assembly选定。历史75/964/759测试记录仅作历史。
