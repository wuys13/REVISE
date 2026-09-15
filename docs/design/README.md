# REVISE Design Records

本目录保存尚未完全落实到代码、但已经形成明确语义边界的设计记录。它面向三类读者：

- 后续实现者：不必重新翻查长对话即可开始工作；
- 审阅者：能区分当前实现、目标契约和仍未验证的假设；
- 未来维护者：能追溯一项调整为什么发生、影响哪些路径、如何验证。

## Design packages

| 主题 | 唯一入口 | 主要内容 |
| --- | --- | --- |
| Batch reconstruction 与 analysis | [Batch 协议与调用](batch/README.md) | 分层配置、输入载体、单项/批量调用、发布与恢复；科学比较链接 reconstruction-impact |
| Reconstruction impact 与应用分析（设计整合，待分阶段实施） | [Reconstruction Impact 主计划](reconstruction-impact/README.md) | 全基因分析目标、现有 impact、Moran/AUCell 与应用能力定位；输入视图、CSV/报告、批量接入及验收按需读取 |
| Notebook 与 analysis 适配（已审阅，待实施计划） | [Notebook 与 Analysis 适配：主决策文档](notebook-analysis-adaptation/README.md) | 精简总入口；`basic/advanced` 架构、第一阶段边界、核心修改与执行顺序，详细数据/现状脚本/route 映射由分文档承载 |
| Reconstruction runtime 统一（active） | [Reconstruction Unification Design Package](reconstruction-unification/README.md) | OTKernel、GA/LR 边界、Application preprocessing/return、5+6 YAML、route trace 与 P2CRC parity |
| Assignment guidance 历史方案（已被简化契约取代） | [Assignment Guidance Design Package](assignment-guidance/README.md) | 保留旧 policy/state-machine 设计作为历史决策证据；当前行为见 active source docs |

设计记录不是运行证据。某项行为只有在代码、测试和对应运行证据完成后，才能从“目标设计”升级为“已实现”或“已验证”。
