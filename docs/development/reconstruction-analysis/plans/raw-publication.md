# R3：完整 Raw 与已有 annotation 发布

**状态：生产端实现及工程验收完成；真实样本待验收。** 用户已批准与 R2/R4 合并执行；批次级边界和剩余真实输入入口见[第二批总计划](whole-sample-delivery.md)。

[返回重建分析总览](../README.md) · [R3 样本交付分类](../categories/02-sample-delivery.md#r3) · [共同约束](../decisions-and-sources.md#confirmed-decisions)

依据：[共同决策](../decisions-and-sources.md#confirmed-decisions)及下列当前实现；旧本地草案仅为可选历史线索。

## 目标与范围

同一次成功的 reconstruction 为分析端交付完整原始观测轴 Raw，以及本次真实获得的 annotation。

Raw 必须保留原始 X、obs/var 轴、坐标和可追溯 ID。

Raw 不得来自已经被 preprocess、normalize 或 runner 原地变换的工作对象。

只回填实际获得的 Level1/Level2 推断；没有推断的单位保留缺失。

Confidence 暂复用现有列语义；GA/LR 分层 confidence/posterior 方案后置。

不在 R3 引入 reference discovery、gene-wise uncertainty 或真实 cell pairing。

## 已确认行为与代码证据

- [reconstruct.py](../../../../reconstruct.py) 当前顺序是 load → `preprocess_data` → reconstruct；Raw 快照必须在工作对象改变前取得。
- [application publication](../../../../revise/application/publication.py) 已提供临时写入、commit/rollback 的发布边界。
- [batch AnalysisInputs.raw](../../../../revise/batch/inputs.py) 当前从 manifest 中的原始 spatial 输入路径读取 Raw，尚无新的 Raw 发布角色。
- [batch spatial/reconstructed views](../../../../revise/batch/inputs.py)是既有内部消费路径；分析库的独立入口由 R4 衔接，不要求其使用 AnalysisInputs。
- [batch/direct parity 测试](../../../../tests/integration/batch/test_real_sample_parity.py)提供 obs、var、X 和 spatial 坐标比较入口，本轮没有运行。
- [publication tests](../../../../tests/application/test_publication.py) 已覆盖失败时旧产物保护，可作为 transaction 验收基础。

这些证据说明接口落点；生产端 Raw publication 现已实现并由第二批工程套件验证，真实 P2CRC_Xenium 的内容和来源仍待验收。

## 具体模块与流程

1. `reconstruct.py` 加载后保留未经处理的原始对象或可校验重载来源，记录原始 obs/var 顺序、坐标和输入 identity；不无条件增加一份大型矩阵深拷贝。
2. 现有 preprocess、GA、LR 和 assembly 继续使用工作副本；不得把 Raw 快照传入会原地改写的 runner。
3. 将本次真实推断按稳定 observation ID 对齐回 Raw 的 annotation 列。
4. 缺失推断保留缺失；原始表达和坐标不因 annotation join 被重排或过滤。
5. 将 Raw 作为现有 publication transaction 的一个 sample-level artifact，与同次 SVC 一起安装。
6. 在 batch reconstruction manifest 中记录 Raw 路径、轴摘要、来源 identity 和 annotation provenance。

## 接口变化

扩展 [output_paths](../../../../revise/application/publication.py) 和 `_published_artifacts`，增加 `raw.h5ad` 的 Raw role；旧 direct 调用的返回值保持兼容，Raw 通过 sample-level 发布与运行记录交付。Raw 输入若与 working spatial 使用同一个内存对象，在进入会原地归一化的 runner 前先 copy 冻结；文件输入仍以最终可校验重载为准，不要求长期保留第二份大型矩阵。

在既有 batch fingerprint 与 publication 记录中纳入 Raw 产物身份，使 Raw 与 SVC 必须来自同一输入和同一成功运行；不另建并行的 fingerprint 框架。

扩展 [AnalysisInputs](../../../../revise/batch/inputs.py) 读取 Raw role 的校验，复用现有 `_load`、axis 和 coordinate checks；分析库正式入口仍从 `sample.yaml` 消费交付文件。

不得把 normalized 工作矩阵、重建 SVC 或历史 Raw 文件冒充本次 Raw。

## 标签与冲突：本批行为

原始 annotation 与本次推断使用独立列：新增 `revise_Level1`、`revise_Level2` 保存真实获得的推断，原始列保留，不覆盖同名人工标签，也不拆分 Confidence 或引入 posterior。

通过唯一、非空的真实 observation ID 回填；没有本次推断的单位保持缺失，QC 排除单位仍在 Raw 中。原始标签与新推断冲突时保留两侧并记录冲突统计，不因冲突覆盖原值或使 sample 自动失败。若输入已经占用新增列名，发布前明确报错，避免静默丢失历史信息。

新的 Raw coverage boolean 仍不是本项必需公共字段；参与重建、跳过原因和 annotation provenance 记录在 manifest/运行记录中。

发布元数据也保护原始命名空间：如果输入已经占用 `uns.revise_delivery`，发布明确报错，不覆盖原始交付记录；已有 `uns.revise_reconstruction` 保留，不能被新交付元数据替换。

Level2 缺失不应通过填充 Unknown 来伪造覆盖；分析端按缺失处理。

## 失败与恢复

Raw 快照失败、ID 对齐失败、X/轴/坐标校验失败时，sample failed，不发布 Raw 或 SVC 新版本。

R2 eligible 类型计算失败同样阻止整套 sample publication；不发布只有 Raw 的“半成功”结果。

沿用 publication 临时文件和 commit 后清理；安装失败恢复历史 Raw/SVC 文件。

输入 identity、配置、参与类型或 annotation mapping 改变时，旧 fingerprint 不得复用。

## 实施步骤（已完成，剩余真实样本验收独立进行）

1. 确认 Raw source snapshot 的时点和 copy/backed 行为。
2. 按已对齐的独立推断列和冲突记录规则，定义稳定 ID 对齐与 annotation 写入规则。
3. 在 publication manifest 增加 Raw role 和 fingerprint 输入。
4. 扩展 batch `AnalysisInputs` 的 raw path/metadata 校验。
5. 接入 transaction failure/recovery，保护历史产物。
6. 用小型 fixture 验证原始轴、缺失标签和坐标保存。

## 针对验证与验收

- fixture 含会被 QC 过滤的单位：发布 Raw 仍包含完整原始轴。
- 比较发布前原始 X、obs/var 顺序和 spatial 坐标，证明未被 normalization 污染。
- 只对本次真实获得 ID 回填标签；无推断单位保持缺失。
- R2 中跳过的类型不在 SVC 重建结果中，但不从 Raw 删除。
- 任一对齐或发布错误后，历史 Raw/SVC 字节内容保持不变。
- batch 与 direct 交付的 Raw 都有可校验的产物记录；既有内部 reader 若仍保留则同步适配；外部分析验收由 R4 负责。不运行真实大样本代替针对验证。

上述工程项已由最终套件覆盖，包括内存 Raw alias 防污染回归、QC 排除单位保留、原始轴/X/坐标和发布回滚；真实 P2CRC_Xenium 仍需按[验收清单](../acceptance.md)运行和人工核对。

## 用户参与、材料与必停边界

原始标签保留、新推断独立保存、冲突记录且不覆盖的行为已在第二批计划中对齐；不要求用户先提供新的 confidence 方案。执行时核对一份含缺失 Level2、QC 排除单位和坐标元数据的真实/代表性 fixture；Agent 先查现有 P2CRC_Xenium 输入，只有证据不足才上提材料需求。

Confidence 分层、posterior 持久化、reference 未标注的科学补全均停止在后续讨论。

若只能从 normalized 工作对象恢复 Raw，必须暂停 R3 发布并上提，不能猜测或复制历史文件；R2/R4 中不依赖该来源的工程检查可以继续。
