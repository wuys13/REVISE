# 第二批：统一样本交付与分析衔接

**状态：生产端实现及工程验收完成；真实样本与科学/联合验收待完成。** 用户已批准将 R2、R3、R4 合并为一个实施批次；本页记录批次级边界、工程证据和剩余真实样本入口。工程批次已经完成，不能把未运行的 P2CRC_Xenium 或未支持的 sST 表达分析写成已验收。

[返回总览](../README.md) · [验收清单](../acceptance.md) · [R2 详细计划](whole-sample-svc.md) · [R3 详细计划](raw-publication.md) · [R4 详细计划](analysis-handoff.md) · [统一样本交付分类](../categories/02-sample-delivery.md)

## 交付目标与边界

一次成功的样本运行提供同一运行身份下的：

- 完整 `raw.h5ad`：保留原始 X、obs/var 轴、坐标、原始注释和真实 ID；只回填本次确实获得的推断。
- 统一 `SVC.h5ad`：iST 在一次 GA 后遍历实际 broad type，完成符合条件的 LR 后合并；首版 whole-sample 默认 `random`，保留 `mean`。
- 分析端 `sample.yaml`：使用相对路径，声明两侧各自真实的标签、坐标和表达语义，并保留样本 ID 映射与来源记录。
- 同次运行的 manifest、fingerprint、日志和发布状态，用于恢复、失败保护和审计。

本批包括生产端 R2/R3/R4；C1 只登记与 `REVISE_Analysis_Agent` 的交接和联合验收要求，不修改兄弟仓库的并行工作树。R1 已完成的 sST 最终缩放调整不回退；第三类 assembly 研究和第四类等待方案的事项不进入本批。confidence 继续复用现有列，不能在本批擅自引入新的 posterior 或不确定性协议。

本批不以真实大样本运行代替工程验证，也不因 loader 成功就宣称科学分析或生物学收益已经验收。各层结果分别记录为工程通过、真实输入可消费、分析结果状态和用户科学审阅。

当前工程证据已经完整记录：最终生产端和交接测试套件 **759 passed, 43 warnings, 0 failed, 0 skipped**，日志见 `/tmp/revise-r2-final.log`，机器可复核摘要见[验证记录](../verification-2026-09-19.json)。全样本端到端 fixture 走真实 application 预处理、pipeline、publication 和 batch/direct 对齐流程，仅以确定性小 stub 替代 GA/LR kernel；已验证 GA 单次调用、多类型 A/B/C 处理、Raw 保留 QC 排除单位、SVC 合并和 direct/batch 的 X/obs/var/坐标一致。消费者正式入口完成 loader、标签/空间 diversity 和报告生成，输入哈希保持不变；整体结果因 expression identity unknown 为 `partial`。P2CRC_Xenium 真实样本和用户科学审阅仍待验收，sST normalized nonlog 表达仍受 C1 限制。

执行中处理的工程问题也已记录：CLI 在沙箱内触发 OpenMP 共享内存错误时，验证命令关闭可选 readline 并在沙箱外重跑，CLI 检查通过；lazy kernel 测试缓存改为隔离，避免跨测试污染。43 个 warning 是依赖弃用提示及刻意覆盖 NaN/overflow、AnnData view/index 和重复 ID 的边界 fixture，不构成失败。

## 已对齐的批次行为

| 方面 | 本批固定行为 | 失败或缺失时的处理 |
|---|---|---|
| 重建范围 | iST 全样本 GA 只运行一次，然后按空间样本实际出现的 broad type 固定顺序遍历 | 未满足资格的类型警告并跳过；符合资格但计算失败使整个 sample failed |
| LR 资格 | whole-sample 对应 reference 中有效、非空且互不相同的 Level2 超过一个；显式旧单类型保留既有行为，允许一个有效 Level2，零个仍失败 | whole-sample 中缺失 Level2 的 reference 行不参与资格判断；不得将缺失值转成字符串类别 |
| SVC 组装 | 只合并真正完成重建的类型；默认 `random`，保留 `mean`；cluster assembly 后置 | 全部跳过或无成功类型不发布空成功；合并前检查基因轴和 ID |
| 旧入口 | 保留显式 `select_ct`、`select_cell_type` 和 batch `cell_type` 的单类型语义；旧 `cell_types` 列表可作为显式兼容模式 | 兼容调用不被标作 whole-sample 完整交付；全样本不接受 `paired`，给出迁移提示 |
| Raw | 预处理前建立可信原始来源；保留完整原始矩阵、轴、坐标和原始注释 | 无可信原始对象或可验证重载来源时，不能声称生成完整 Raw |
| 推断标签 | 新推断保存为独立的 `revise_Level1`、`revise_Level2`；原始标签不覆盖；按唯一、非空真实 ID 回填 | 新列名已被输入占用时明确报错；原始标签与推断冲突记录统计，不因冲突覆盖或伪造 |
| 发布与恢复 | Raw、SVC、`sample.yaml`、manifest 和运行记录是同一发布单元；写入/安装失败保护历史完整结果 | 不发布半成功组合；旧单类型成功记录不能冒充新的 sample-level 成功 |
| 分析衔接 | `sample.yaml` 指向同次产物，标签/坐标/表达声明按实际来源填写；消费者从独立输入入口加载 | 不支持的表达尺度标记 unavailable/partial；sST normalized nonlog 不冒充 raw counts 或 log1p |

Raw 发布还保护来源元数据：如果输入已经占用 `uns.revise_delivery` 命名空间，发布明确失败，避免覆盖原始交付记录；既有 `uns.revise_reconstruction` 保留，不被新交付元数据替换。

这些行为的详细接口和测试边界分别维护在 [R2](whole-sample-svc.md)、[R3](raw-publication.md) 和 [R4/C1](analysis-handoff.md)；本页不复制它们的全部实现说明。

## 实施顺序

<a id="r2"></a>
### 1. R2：全类型重建与统一 SVC

在现有 application、strategy 和 batch 路径上扩展：

1. 复核全样本入口、旧单类型入口和 `random`/`mean` 配置迁移。
2. 在共享 context 中完成一次 GA，保留其 axis、posterior 和 provenance，供后续各类型 LR 使用。
3. 计算实际 broad type 的 whole-sample Level2 资格；跳过不合格类型并记录原因，对合格类型运行现有 LR/assembly。显式旧单类型调用走兼容资格规则。
4. 按确定顺序合并 spatial/expression carriers，校验基因轴、ID 和跨类型 cluster 标识。
5. 将 sample 级状态接到现有 fingerprint、publication transaction 和 batch handoff；不删除旧类型目录。

<a id="r3"></a>
### 2. R3：完整 Raw 与真实推断

R2 的 SVC 组装可与 Raw source snapshot 和发布接口工作并行，但最终联合发布要等两者一致：

1. 从未经预处理的输入或可验证重载来源建立 Raw snapshot，不把 normalized 工作对象当作 Raw。
2. 保留完整原始 X、obs/var 顺序、坐标、原始注释和 QC 排除单位。
3. 只将本次真实获得的推断按稳定 ID 回填到独立列；没有推断的位置继续缺失。
4. 将 Raw role、来源 identity、轴摘要和 annotation provenance 纳入同次 manifest/fingerprint。
5. 在计算、对齐或安装失败时回滚本次临时产物，保护历史 Raw/SVC。

<a id="r4-c1"></a>
### 3. R4/C1：样本配置与联合消费

1. 由生产端生成相对路径的 `sample.yaml`，层级 sample ID 使用确定、可逆的单目录段编码，并保留原始 ID 映射。
2. `columns.broad`/`columns.subtype` 指向两侧真实存在的 `revise_Level1`/`revise_Level2`；不存在的 `SVC_cluster` 或物理比例不伪造。
3. 依据真实来源声明每一侧 matrix identity、scale 和坐标单位；不把未知或 sST normalized nonlog 写成已支持的 raw/log1p。
4. 使用分析库真实 loader 加载本批产物，并在前提满足的代表性 fixture 上完成至少一个标签/空间分析；按 `succeeded`、`partial`、`skipped`、`failed` 记录实际状态。
5. 生产端和消费者分别保留版本、命令、输入哈希和结果记录；不修改消费者工作树。

### 4. 统一工程验收与交付

R2/R3/R4 的工程核对已经一次性完成，覆盖多类型、缺失标签、失败恢复、Raw 完整性、SVC 组成、旧入口和分析消费。结果、759 项测试和剩余真实联合验收见[验收清单](../acceptance.md)与[验证记录](../verification-2026-09-19.json)。未完成的真实样本和科学/表达验收不启动第三、四类事项。

## 测试与验收边界

必须覆盖：whole-sample 多类型和单 Level2 类型、显式旧单类型的一/零 Level2 兼容行为、部分/空白 Level2、全跳过、合格类型失败、GA 单次调用、random 可复现、mean 输出、跨类型 cluster 不碰撞、direct/batch 语义一致、Raw QC 排除单位保留、原始轴/坐标/X 不变、标签冲突及缺失保留、计算/写入/安装失败恢复，以及旧结果不被错误复用。

针对分析交接，至少验证相对路径从其他工作目录仍可解析、两侧独立轴可加载、输入 H5AD 不被改写，并使用分析库的正式入口完成一个已有前提的标签/空间步骤。sST 表达分析尚未具备 normalized nonlog 的完整消费方案，必须单独列为未完成项；它不阻塞已有标签或空间检查。

工程测试通过不等于真实 P2CRC_Xenium 科学验收。真实命令和逐项检查见[一页验收清单](../acceptance.md)；本批默认不执行昂贵的全量重建。

## 当前用户参与与上提边界

本批已获得整体实施授权，无需用户提前提供材料。实现过程中遇到常规接口、路径或测试问题，由执行者在本批内处理并在完成报告中说明。

只有以下情况应暂停受影响部分并上提：无法获得可信 Raw 来源、必须覆盖原始标签才能继续、现有消费者接口要求改变已确认的公共语义、或 shared GA 无法在不重复计算的情况下传递。此时继续完成不依赖该决定的工程核对，不把局部失败包装成整批成功。

## 相关计划和记录

- [R2 详细计划：iST whole-sample SVC](whole-sample-svc.md)
- [R3 详细计划：完整 Raw 与已有 annotation](raw-publication.md)
- [R4 详细计划：分析交接与联合验收](analysis-handoff.md)
- [第二类事项背景](../categories/02-sample-delivery.md)
- [决策、来源与执行边界](../decisions-and-sources.md)
- [验收清单与 P2CRC_Xenium 命令](../acceptance.md)
