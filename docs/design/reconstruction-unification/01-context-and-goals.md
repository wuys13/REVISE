# Context, Goals, and Evidence Boundary

Parent index: [Reconstruction Unification Design Package](README.md)

> 本文记录实施启动时的 baseline。交付后的 Current fact 与 Verified evidence 只从 [09-delivery-status.md](09-delivery-status.md) 读取。

## Why this work exists

Current fact：当前包已经有统一 pipeline 与 route strategy，但 OT 和 preprocessing 的所有权仍然交错：

- `GlobalAnchoringKernel` 自己实现 POT/TACCO annotation；
- `LocalAnchoringKernel` 在 Application sc-SVC + TACCO 时反向委托 `GlobalAnchoringKernel`；
- matrix coupling 主要通过 `solve_local_ot`，但多个 runner 直接持有该函数；
- Application preprocessing 同时存在于 adapters、`ApplicationSVC._adata_processing` 和 sc-SVC-sr 构造流程；
- `reconstruct.py` 同时承载 parser、config mapping、publication 和 `ApplicationExecution` wrapper，实际 AnnData 流不够直接；
- `revise/revise.yaml`、Application YAML 和 Benchmark CLI 同时承担结果参数与路线选择。

这使“形式 unified”与“算法真正落回旧实现”难以独立审查，尤其容易把 sc-SVC 两次 LR TACCO annotation 错算为 GA。

## Legacy scientific baseline

Current fact：`test_sc/test.py` 是 standard sc-SVC Xenium_T 的原始语义基准。其关键顺序是：

```text
read P2CRC Xenium
→ spatial transcript/gene filtering
→ read common sc reference and filter Patient=P2CRC
→ slash label normalization and reference gene filtering
→ Level1 TACCO annotation
→ select T cohort
→ Level2 TACCO annotation
→ gene/spatial/joint graph and resolution selection
→ SVC_cluster
→ reverse TACCO annotation from SVC cluster to sc reference
→ write spatial and expression H5AD
```

Decided target：包可以保留统一接口、typed config、provenance 和 transaction publication 等形式差异，但所有会影响最终结果的数据子集、轴顺序、solver 输入、graph/resolution 选择和返回对象都必须可逐步与旧脚本比较。

## Goals

Decided target：

1. 建立一个底层 `OTKernel`，统一 annotation OT 与 matrix coupling OT 的真实 solver 调用位置。
2. 把 GA facade 和 LR computation 分开，保证 LR 永远不依赖 GA kernel。
3. 让 `reconstruct.py` 可直接读出 YAML → effective config → 三个预处理函数 → pipeline → publication → AnnData return。
4. 让所有 Application/Benchmark 结果参数在 runner 构造前形成完整 effective config，并能重放其来源与 hash。
5. 用 5 个 Application YAML 和 6 个 Benchmark YAML 表达正式路线，移除旧 noise 污染与双路线权威。
6. 以九类预期运行线路和旧脚本十一个检查点作为验证基准。

## Non-goals

Decided target：本轮不做以下工作：

- 不比较或改变 POT 与 TACCO 的科学优劣。
- 不引入自动 solver fallback、posterior normalization/reorder 或 zero-overlap auto filtering。
- 不扩大 slash-normalization collision contract；保留现有失败行为。
- 不删除 Benchmark 正常 preprocessing、mandatory SR allocation、graph aggregation 或 imputation。
- 不将技术一致性写成真实 P2CRC 数值一致，更不写成生物学验证。

## Evidence levels

| Level | 可以证明 | 不能证明 |
| --- | --- | --- |
| Static boundary | import/call graph 满足 OTKernel、GA、LR 依赖约束 | solver 真正返回相同结果 |
| Focused synthetic | 输入输出轴、dtype、support、发布和 route call 顺序符合契约 | 真实 P2CRC 中间量一致 |
| Legacy parity run | 相同环境/seed 下，新旧 11 个检查点达到约定一致性 | 结果具有生物学有效性 |
| Scientific evaluation | 预先定义的真实数据指标和人工审查支持科学结论 | 不自动推广到未测样本或平台 |

## Protected workspace boundary

Current fact：`raw_data/` 与 `test_sc/` 在当前 checkout 中是未跟踪材料。Decided target：任何实施、测试或打包步骤都不得 stage 它们；运行输出同样保持未跟踪或位于明确忽略路径。

继续阅读：[Route contracts](02-route-contracts.md)。
