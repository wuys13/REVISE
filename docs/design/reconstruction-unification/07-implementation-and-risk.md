# Implementation Order, Review Decisions, and Risk

Parent index: [Reconstruction Unification Design Package](README.md)

## Dependency order

1. 冻结 route/config/OT/preprocessing baseline 与 characterization tests。
2. 先将 matrix coupling 的 solver call 收进 `OTKernel.couple`，不改变 problem construction。
3. 再将 annotation solver call 收进 `OTKernel.annotate`，让 GA/LR 同时依赖它并删除 LR→GA 边。
4. 搬出三个 Application preprocessing 函数，删除重复 class/adapter 调用点。
5. 简化 `reconstruct.py`、AnnData return 与 publication。
6. 建立 typed config authority，迁移 Application YAML。
7. 迁移六类 Benchmark YAML，删除旧 route selector 与 noise 污染。
8. 运行全路线 trace、distribution、legacy parity；最后更新 delivery status。

每一步先在原调用点做 characterization，再进行机械迁移。若新 wrapper 需要改变算法参数、轴或 fallback 才能通过测试，应停止并报告，而不是把行为变化藏进 refactor。

## Accepted independent-review decisions

先前从 OT architecture、legacy semantics 和 config/adversarial 三个独立角度审查后，已经接受：

- 统一 seam 位于底层 OT 输入输出，而不是强迫所有 LR 使用同一上层实现。
- annotation 与 matrix coupling 必须是同一 Module 的两个同级接口，不能互相实现或用任意 options dict。
- sc-SVC 两次 TACCO annotation 明确属于 LR；`LocalAnchoringKernel` 不得复用 GA kernel。
- TACCO returned category order 是 canonical result order；验证集合但不按 reference 顺序重排。
- Application 三个 preprocessing 函数直接在 `reconstruct.py` 串行调用，不增加 orchestration class/function。
- Benchmark YAML 是唯一 route selector，删除 `--confounding`，不做 old/new auto-detection。
- `output.name` 只影响 filename prefix，不承担内部 identity。
- dry-run 用轻量 preflight 检查 reference filter 匹配，不读取完整矩阵。
- route trace 用 spies/fakes 观察真实调用，不建第二套状态机。

明确拒绝或延期：

- annotation/coupling 两层互相转换的“更统一”抽象；
- GA 与 LR 共用一个 runner-owned kernel instance；
- 自动 posterior reorder/normalization；
- zero-overlap 自动修复和 solver fallback；
- slash collision 新策略；
- Benchmark noise 或额外 imputation YAML；
- 本轮真实数据生物学验收。

## High-risk areas and controls

| Risk | Failure signal | Control |
| --- | --- | --- |
| OT wrapper 改变 dtype/normalization | POT/TACCO coupling 数值漂移 | 迁移前后同输入 characterization；route 保留 row normalization |
| TACCO category order 被旧 strict check 拒绝 | same labels, different order failure | set/uniqueness validation；保留 returned order；回归测试 |
| preprocessing 重复执行 | obs/var 数量二次下降 | route trace 断言三函数各一次；删除 class/adapter duplicates |
| sc-SVC expression var 被 overlap 截断 | expression output gene 数下降 | 单独断言 spatial overlap axis 与 full filtered reference axis |
| SR mandatory allocation 被误删 | raw sc-SVC-sr 不生成 | noise 删除白名单与 mandatory allocation trace |
| graph aggregation primary 选择错误 | 写盘与返回对象不同 | 对 enabled/disabled 两分支做 identity/axis/value 检查 |
| typed defaults 漏迁移 | runner dataclass fallback 悄悄生效 | old→new field inventory；effective config snapshot；禁隐藏 defaults |
| Benchmark selector 双权威 | YAML 与 CLI route 不一致 | 删除 `--confounding`；launcher 只传 YAML；invalid schema fail-fast |
| paired publication 半提交 | spatial/expr 一新一旧 | 保留既有 staging + atomic replace；本轮不新增 backup/rollback 防御逻辑 |
| 误宣称验证完成 | docs 与证据不匹配 | 只有 delivery status 能升级 Verified evidence |

## Surgical deletion lists

只在替代路径已覆盖后删除：

- runner-owned `annotate_method`/`global_anchoring` 与 `BaseSVCAnchor` 层；
- `LocalAnchoringKernel` 对 GA 的 import/delegation；
- adapter `_ensure_transcript_counts` 与 Application filters；
- `ApplicationSVC._adata_processing`、sc-SVC-sr constructor 重复调用；
- `ApplicationExecution` public wrapper/exports/tests/docs；
- `_inject_sr_spatial_leakage_noise`、`sr_noise_*`、`sr_spatial_noise`、`st_input_noisy` 和 noise-only GA override；
- `--confounding`、旧 auto-detection、`revise.yaml` 与 packaging entry。

不能顺手删除既有无关 dead code。受保护列表见 [02-route-contracts.md](02-route-contracts.md)。

## Review gates

- Gate A：OT static boundaries 与 focused equivalence 通过后才能迁移 routes。
- Gate B：三函数调用一次、轴契约通过后才能删除 class/adapter preprocessing。
- Gate C：所有结果参数完成迁移、repo/package templates 一致后才能删除 `revise.yaml`。
- Gate D：9 类 route trace、distribution、diff check 通过后才具备软件交付候选资格。
- Gate E：P2CRC 十一点双跑另行升级 real numerical evidence；不阻塞结构重构的诚实交付，但必须明确未验证。

继续阅读：[Core file index](08-core-file-index.md)。
