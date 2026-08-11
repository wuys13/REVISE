# REVISE Design Records

本目录保存尚未完全落实到代码、但已经形成明确语义边界的设计记录。它面向三类读者：

- 后续实现者：不必重新翻查长对话即可开始工作；
- 审阅者：能区分当前实现、目标契约和仍未验证的假设；
- 未来维护者：能追溯一项调整为什么发生、影响哪些路径、如何验证。

## Design packages

| 主题 | 唯一入口 | 主要内容 |
| --- | --- | --- |
| Reconstruction runtime 统一（active） | [Reconstruction Unification Design Package](reconstruction-unification/README.md) | OTKernel、GA/LR 边界、Application preprocessing/return、5+6 YAML、route trace 与 P2CRC parity |
| Assignment guidance 历史方案（已被简化契约取代） | [Assignment Guidance Design Package](assignment-guidance/README.md) | 保留旧 policy/state-machine 设计作为历史决策证据；当前行为见 active source docs |

设计记录不是运行证据。某项行为只有在代码、测试和对应运行证据完成后，才能从“目标设计”升级为“已实现”或“已验证”。
