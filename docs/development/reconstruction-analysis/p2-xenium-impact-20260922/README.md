# P2 Xenium full random → Impact：真实运行与独立审视

Task ID: `p2-xenium-impact-20260922`。状态：`preparing`；科学接受：尚未授予。

## 目标与责任

本任务不以进程成功结束作为完成条件。目标是在真实全量输入上解释重建前后变化，修复可解决的重大问题，并明确尚不能解释的部分。REVISE 负责参考准备、重建与三文件交付；Analysis Agent 负责独立分析、结果和报告。主代理整合证据；执行代理各自负责所属仓库；独立审视代理只读复核实现、底表和解释。

首轮固定 `P2CRC_Xenium / random / seed=42`，参考过滤 `Patient=P2CRC`；Impact 使用 `All/Fibroblast/Mono_Macro/T`、40 μm、EMT 及既有参数。原生 Raw/SVC 独立，仅有明确需要的局部诊断使用共同 ID。首轮不执行其他平台、其他 assembly 或独立 Moran/pathway。

## 执行边界

- 使用启智当前 RUNNING 的 cpu 或 cpu-4；每次任务启动核实资源与占用。
- 两仓分别安装独立 Python 环境，不修改已有共享环境。
- 复用 `raw_data/zenodo/Real_application/P2CRC_Xenium.h5ad` 与 `adata_sc_all_reanno.h5ad`。准备阶段双端 SHA256 分别为 `2e0f31756e6dc2f84d258d48ec02ee38b1c188e815eb1910536f348afabaa840` 和 `c7c277aad7ccc3b0520eeb3c28fb19d4c89314688f35f04dfb3a8e315c01e714`；执行 receipt 再记录实际身份。
- GMT 为 `h.all.v2025.1.Hs.symbols.gmt`，SHA256 `f22066af72e215ccb7b89d88e492c07e1eef17534c2ca7b0f9902cfecbbdd8e9`。
- 正式运行与诊断使用独立轮次目录。科学代码、参数或输入改变后，完整重跑受影响的 Impact；batch 不支持跨进程恢复阶段。报告展示改变可只读刷新。
- 不为得到正面结果而改变阈值、填零、扩大比较口径或选择参数。

## 审视矩阵与完成条件

逐项记录输入/donor 身份、覆盖与观察单位、标签粒度、Anatomy、共同基因与表达尺度、EMT 评分背景、reference/random 影响、稳定性及报告一致性。每项必须有底表、分母、适用结论和限制。见 [问题与判断表](findings.md)。

所有影响核心判断且可解决的重大问题闭环，最终产物经独立审视与主代理复核后，才可结束本次执行。未知若阻断核心判断，则保持 `partial/blocked`；其他科学未知以完整因果和证据留给用户，不自动写为 accepted。

## 证据索引

运行回执须包含 owner/upstream/downstream、两仓 commit/dirty、解释器/依赖、输入/GMT 哈希、命令、执行位置、实例、观测时间、状态、资源峰值和产物入口。各仓自己的 receipt/result 是事实来源，本页只索引。

| 证据 | 当前状态 | 入口 |
|---|---|---|
| 独立启动审视 | GO_WITH_LIMITS | 必须 run-scoped + Impact-only；未发现必须先修改的科学实现阻塞 |
| 重建环境与资源 | preparing / full gate closed | cpu 实测 cgroup 4 核、8 GiB；baseline 约 2.06 GiB。正在核对稀疏 random 路径的真实容量，不能用宿主机内存代替配额 |
| 分析环境与资源 | installing | cpu-4 已部署 clean `sl@da097ab`，GMT SHA已复核；CPU-only环境安装中，尚未启动Impact |
| iter-001 三文件交付 | not_started | `results/p2-xenium-impact-20260922/iter-001/P2CRC_Xenium/random`（预定） |
| iter-001 Impact | not_started | 待交付配置固定 |
| 独立结果审视 | not_started | 待真实结果 |
| 最终判定 | not_started | 尚无全量或科学完成声明 |

## 运行前实现验证

分析端 `da097ab` 已通过118项全库测试与16项聚焦测试（有重叠，不相加），独立代码审视关闭已有问题。重建端74项聚焦测试通过；日志 `output/p2-xenium-impact-20260922/iter-001/logs/producer-focused-tests-final.log`，SHA256 `61bcef6ae5c357019dcc62476165c7757293aac57b7cc9ad0ec87eb4e9144cd1`。以上均为技术验证，不是full运行或科学接受。

真实输入的稀疏存储测量见 [输入存储观察](input-storage-observation.json)，不等于运行峰值。结果前的诊断触发与边界见 [必要诊断协议](diagnostic-protocol.md)。
