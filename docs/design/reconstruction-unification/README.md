# Reconstruction Unification Design Package

这是 REVISE reconstruction runtime 统一工作的唯一导航入口。本文档包把“实施前代码事实”“已经定稿的目标”“已经执行过的验证”和“仍需真实数据回答的问题”分开记录，避免把结构重构、测试通过或一次真实运行误写成算法等价或生物学验证。

## Authority and status vocabulary

- **Current fact**：从本轮开始时的代码、配置或测试直接核对到的事实。实施完成后的当前事实以 [09-delivery-status.md](09-delivery-status.md) 为准。
- **Decided target**：用户已经确认、实现不得自行改变的目标契约。
- **Verified evidence**：有具体测试、构建或运行产物支持的结论；必须注明证据范围。
- **Scientific unknown**：需要真实数据、benchmark 或生物学评估才能回答，不能由静态检查或合成测试推断。

当文档和代码冲突时：baseline 的 Current fact 以当时 checkout 为准；Decided target 以本设计包为准；交付后的实现事实只写入 delivery status。正式实施顺序由 [implementation plan](../../plans/2026-08-11-001-refactor-unify-reconstruction-runtime-plan.md) 负责。

## Frozen decisions

1. GA 可以经过 `GlobalAnchoringKernel`；任何 LR Module、runner 或 strategy 都不得出现 `GlobalAnchoringKernel`。
2. 所有 `tacco.tl.annotate`、`tacco.utils.solve_OT` 和 `ot.unbalanced.sinkhorn_unbalanced` 的实际调用只能位于一个底层 `OTKernel` Module。
3. `OTKernel.annotate` 与 `OTKernel.couple` 是两个同级、互不调用的固定契约；它们不读取 runner/config，也不依赖 GA/LR stage kernel。
4. sc-SVC 的两次 LocalAnchoring annotation 属于 LR；其他路线的 LR 负责构造自身的 matrix coupling problem，然后调用 `OTKernel.couple`。
5. TACCO 返回的 category 顺序是结果权威；验证 observation/category 轴、集合、唯一性和数值，但不按 reference 首次出现顺序重排，也不自动归一化 posterior。
6. Application 使用三个独立预处理函数，直接在 `reconstruct.py` 的 lifecycle callback 中按 reference filter → spatial preprocess → reference preprocess 顺序调用；不增加第四个总包装函数。
7. slash label normalization 保留；本轮不扩大 collision、零交集、solver fallback 或 posterior 修复范围。
8. Benchmark YAML 是唯一 route selector；最终只保留 5 个 Application YAML 和 6 个 Benchmark YAML，不保留 `revise.yaml`、`--confounding`、`noise.yaml` 或 `imputation.yaml`。
9. `run_application` 直接返回 AnnData；标准 sc-SVC 返回 `(spatial_adata, expression_adata)`，dry-run 返回 `None`。
10. 测试先固定每条路线的预期 stage/call 顺序，再观察真实调用；不为测试引入第二套 production trace/state system。

## Document responsibilities

| 文件 | 专职功能 | 不负责什么 |
| --- | --- | --- |
| [01-context-and-goals.md](01-context-and-goals.md) | baseline、目标、非目标、旧脚本的基准地位与证据边界 | 不定义逐路线调用序列 |
| [02-route-contracts.md](02-route-contracts.md) | 3 条 Application、6 类 Benchmark 的预期运行线路、条件和输出 | 不展开 OT 数值实现 |
| [03-ot-boundaries-and-contracts.md](03-ot-boundaries-and-contracts.md) | `OTKernel` 两种契约、GA/LR 依赖边界、TACCO/POT 不变量 | 不决定 YAML 形状 |
| [04-application-entry-and-preprocessing.md](04-application-entry-and-preprocessing.md) | 三个预处理函数、`reconstruct.py` 流程、返回/发布/dry-run | 不描述 Benchmark case expansion |
| [05-configuration-and-provenance.md](05-configuration-and-provenance.md) | typed defaults、5+6 YAML、override、hash 与 provenance | 不宣称运行已通过 |
| [06-verification-and-p2crc-parity.md](06-verification-and-p2crc-parity.md) | route-trace 方法、旧 `test.py` 十一点对比、P2CRC 双跑程序 | 不替代 delivery status |
| [07-implementation-and-risk.md](07-implementation-and-risk.md) | 实施依赖、风险、删除白名单、独立审查决定 | 不记录最终测试结果 |
| [08-core-file-index.md](08-core-file-index.md) | 生产文件、目标文件、测试文件与按任务阅读索引 | 不重复设计论证 |
| [09-delivery-status.md](09-delivery-status.md) | 实施后事实、验证证据、失败项和科学未知 | 不修改已冻结目标 |
| [Implementation plan](../../plans/2026-08-11-001-refactor-unify-reconstruction-runtime-plan.md) | 依赖顺序、实施单元和验收门槛 | 不复制全部契约细节 |

## Recommended reading paths

第一次接手：

```text
README → 09 delivery status → 01 → implementation plan → 07 → 相关专项文档
```

修改 OT、GA 或 LR：

```text
README → 03 → 02 中目标路线 → 08
```

修改 Application 或预处理：

```text
README → 04 → 02 中三个 Application 路线 → 08
```

修改配置或 Benchmark：

```text
README → 05 → 02 → 07
```

验证 P2CRC：

```text
README → 06 → 02 的 sc-SVC → 03 → 09
```

## Scope boundary

本设计包不要求：

- 所有路线共享同一上层 GA/LR runner；统一边界是底层 OT 输入输出与禁止依赖。
- 为旧 flag、旧 `revise.yaml` 或旧 runner-owned GA 增加兼容 shim。
- 本轮顺手修复无关 dead code、防御性校验或科学算法。
- 把真实数据、`test_sc/`、`raw_data/` 或运行输出加入 Git。
- 以软件 wiring、合成数组或一次成功运行宣称生物学有效。

返回上级索引：[REVISE Design Records](../README.md)。
