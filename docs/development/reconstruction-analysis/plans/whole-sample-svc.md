# R2：iST whole-sample SVC

**状态：生产端实现及工程验收完成；真实样本待验收。** 用户已批准与 R3/R4 合并执行；批次级边界和剩余真实输入入口见[第二批总计划](whole-sample-delivery.md)。

[返回重建分析总览](../README.md) · [R2 样本交付分类](../categories/02-sample-delivery.md#r2) · [共同约束](../decisions-and-sources.md#confirmed-decisions)

依据：[共同决策](../decisions-and-sources.md#confirmed-decisions)及下列当前代码；旧本地草案仅为可选历史线索，本文不依赖其存在。

## 目标与范围

一次 iST 样本运行完成 GA 后，遍历该空间样本实际出现的全部 broad cell types。

whole-sample 中，reference 对应 broad type 若有超过一个不同的有效 Level2，才进入 LR；显式旧单类型调用保留既有兼容语义，允许一个有效 Level2，零个有效 Level2 仍失败。

有效 Level2 排除缺失和空白；不以 reference 细胞数替代类别数。单类型兼容路径不因本批 whole-sample 资格规则而改变。是否存在需要排除的特殊标签须依据实际数据再讨论，不自行定义“无意义标签”。

不满足条件的类型警告并跳过 LR；最终 SVC 只合并真正完成重建的类型。

符合条件但 LR 或 assembly 报错，整个样本失败且不发布新结果。

首版 random 为默认，mean 保留；cluster assembly 不在 R2 实现。

Confidence 暂复用现状，不在本任务中重构。

## 工程验收记录

最终套件 **759 passed, 43 warnings, 0 failed, 0 skipped**，包含 application/config/backend/batch、分析交接和全样本 batch/direct 端到端测试；机器可复核摘要见[验证记录](../verification-2026-09-19.json)，完整日志为 `/tmp/revise-r2-final.log`。全样本端到端测试走真实 application 预处理、pipeline、publication 和 batch 流程，仅用确定性小 stub 隔离 GA/LR kernel；已核对 GA 只调用一次、多个 broad type 的 skip/合并、random/mean、cluster 命名空间、Raw/SVC 产物以及 direct/batch 的 X、obs、var、坐标一致。

这证明生产端工程协议已落地，不替代 P2CRC_Xenium 真实运行或用户科学解释。旧单类型兼容资格、whole-sample 的单 Level2 skip/eligible failure、历史结果保护和 `uns` 来源元数据边界均由套件覆盖；真实样本仍按[验收清单](../acceptance.md)执行。

## 已确认行为与代码证据

- [single application entry](../../../../reconstruct.py) 当前负责编译、加载、预处理和一次 reconstruction。
- [iST assembly](../../../../revise/application/ist_assembly.py) 已支持 mean/random，并保留 spatial obs 轴与坐标。
- [application publication](../../../../revise/application/publication.py) 已校验 spatial/expression carrier 的 required outputs。
- [iST publication tests](../../../../tests/application/test_ist_publication.py) 已覆盖 paired、mean、random 及 transactional failure。
- 当前 [ScSvcApplicationStrategy](../../../../revise/backend/adapters.py) 的 `solve_ot` 只接受一个 `select_ct`，这是 R2 的最小改动落点。
- [sST strategy](../../../../revise/backend/adapters.py) 是另一条 route，不应作为 iST whole-sample 循环的依据。
- [batch parity 测试](../../../../tests/integration/batch/test_real_sample_parity.py)提供 batch 与 direct application 同语义的验证入口，本轮没有运行它。

## 具体模块与流程

1. 保持 `reconstruct.py` 为单样本入口；batch 每个 sample 继续通过已有 application document/runner 路径进入同一 strategy。
2. 沿用 pipeline 在 `solve_ot` 之前调用的 `global_anchoring` 阶段，全样本只做一次 GA；`ScSvcApplicationStrategy.solve_ot` 消费该结果，不再调用 GA。
3. 按空间样本实际出现的 broad type 逐一做 Level2 eligibility。
4. 对 eligible 类型调用现有 runner 的 local refinement；不得每类型重复 GA、输入预处理或样本级 fingerprint。
5. 复用已有 carriers、`assemble_ist` 和 `anndata` 合并边界，把完成重建的类型合并成一个 sample SVC。
6. 用显式 broad/cluster 标识保持跨类型 ID 不碰撞；public 输出仍是一套 sample 结果。

不新造 `SampleReconstructionResult` 横切架构，不把每类型 artifact 契约提前写成已定。

## 接口变化

- 未指定单类型选择时，iST application/batch 走 whole-sample 默认；显式 `select_ct`、`select_cell_type`、batch `cell_type` 保留单类型兼容语义。旧 batch `cell_types` 列表可作为显式兼容模式，但不能把它的成功记录标成完整样本交付。
- whole-sample 只接受 `random`/`mean`；`paired` 给出迁移说明。单类型兼容调用未指定 mapping 时保留旧默认和返回形式。
- batch 的 sample task 是发布与恢复单位；内部仍可保留 type-level carriers、日志和旧目录，不自动删除历史文件。
- `output.ist_mapping` 继续复用现有 `random`/`mean`，whole-sample 默认迁移到 `random`；fingerprint/provenance 必须记录实际选择。

不要为 R2 增加 cluster assignment、Confidence 新列或新的公共 artifact 类型。

## 失败与恢复

无 eligible 类型或全部类型被跳过时，不发布空 SVC 并声称成功。

样本 GA 失败，或单个 eligible 类型的 LR/assembly 失败，均使 sample failed；不能降级为跳过该类型。

发布前先完成全部类型合并、轴校验和 fingerprint；失败时历史完整结果保持不变。

旧 per-cell-type 成功记录不得直接冒充新的 sample-level 成功；不自动删除旧文件。

## 实施步骤（已完成，剩余真实样本验收独立进行）

1. 固定一次 GA 的调用边界和 shared posterior/axis 传递。
2. 把 eligibility、warning 和 eligible type 顺序放进现有 strategy context。
3. 循环现有 local refinement，收集 spatial/expression carriers。
4. 复用现有 assembly/concat 规则生成 sample SVC。
5. 接入现有 publication transaction 和 batch handoff。
6. 最后再迁移 direct/batch API 的 select_ct 兼容层。

## 针对验证与验收

- 两个 eligible broad type 加一个单 Level2 类型：GA 只调用一次；前两者进入 whole-sample LR，后一者有具体 warning 且不进入合并结果；显式选择该单类型时保留旧行为并允许一个有效 Level2。
- 任一 eligible 类型抛错：sample failed，未写入新的完整 SVC。
- random 在固定 seed 下可复现；mean 与 random 均复用已有 assembly 测试。
- 合并结果 gene axis 一致、obs ID 唯一，cluster 标签可在同簇细胞间重复但不同 broad types 的同名 cluster 不混淆，spatial 坐标保留。
- direct application 与 batch sample 输出一致；不把 sST runner 纳入本测试。

上述工程项已经由最终 759 项套件通过；P2CRC_Xenium 的真实类型覆盖、运行时间/资源和科学解释仍未验收。

## 用户参与、材料与必停边界

本批已对齐：保留旧单类型兼容入口；whole-sample 中同一 reference 类型存在部分 Level2 缺失时，LR 只使用有效、非空 Level2 行且要求超过一个类别；显式旧单类型允许一个有效 Level2，零个仍失败。不把全局 Level2 完整性校验作为逐类型跳过之前的门槛；whole-sample 不 eligible 类型只警告并跳过，只合并完成重建的类型。

旧遍历代码/Notebook只是可选核对材料，不是实施阻塞。本批无需用户额外提供材料；常规工程问题在批次内处理并在验收报告汇报。

cluster assembly、Confidence 新协议和 reference discovery 是后续事项。

若 shared GA 与现有 runner 不能在不重复计算的条件下传递，暂停受影响部分并上提接口证据，不新造未经对齐的横切抽象；其他不依赖该决定的工程检查继续进行。
