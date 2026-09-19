# R4/C1：分析交接与联合验收

[实施索引](../implementation-index.md) · [当前协议对照](../cross-repo-review.md) · [验收入口](../acceptance.md)

**状态：生产端实现与小型正式消费者联测完成；真实全量样本及科学验收未执行。** 已对齐的范围包括可逆样本 ID、原始列保留、独立推断别名、同次事务发布和线性表达来源声明。

## 当前接口

- `sample.yaml` schema_version=1，使用相对 `files.raw/svc`；消费者从 YAML 所在目录解析。
- 只声明实际存在的 broad/subtype/reconstruction。iST 的 `SVC_cluster` 是消费者重建主标签，`revise_Level1/2` 不是它的替代品。
- 当前正式 `.X` 契约为有限非负非 log 线性表达，可含小数。identity 必须有来源依据；unknown 可加载但不放行表达分析。已知线性输出不再携带旧 scale unknown；log 不是当前支持格式。
- sST 保留内部 normalization 与 parent-spot 校正事实；不宣称 raw counts。真实 sST 表达消费仍未验收，但不能沿用“所有 normalized nonlog 都不支持”的旧结论。
- 坐标单位或比例未知时明确保留未知；不同技术的 Raw/SVC 独立读取，不强造配对。具体分析缺失前提时保留真实 unavailable/partial/failed 状态。

## 联合验证与责任

REVISE 负责生产事实和发布完整性；分析库负责方法前提及输出状态。调用真实 load_sample/run_analysis，核对标签空间结果、已知线性表达结果、unknown 限制与输入哈希。测试路径与产物集中在[审阅入口](../cross-repo-review.md)，真实验收清单见[当前验收](../acceptance.md)。

真实样本需单独核对来源、覆盖、表达身份与物理尺度；不以 fixture 代替科学结论，不在本轮修改分析库。Notebook 四方法重分群属于[独立比较协议](../assembly-comparison-contract.md)，不等于 Impact 直接使用主标签的步骤。

原计划与旧协议理由保留在[历史快照](analysis-handoff-history.md)，不得作为当前消费者支持范围。
