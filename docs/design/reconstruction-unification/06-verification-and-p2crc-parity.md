# Verification and P2CRC Parity

Parent index: [Reconstruction Unification Design Package](README.md)

## Verification organizing principle

测试以“这条正式路线实际执行了哪些 stage/calls、产生什么对象”为主轴，而不是先按 unit/integration/delivery 分类。先从 [02-route-contracts.md](02-route-contracts.md) 得到期望 trace，再通过 spy/fake 观察 production 调用；不向 production 增加第二套 trace/state system。

## Route-trace matrix

必须覆盖：

- 3 条 Application runtime 路线和 dry-run；
- 6 类 Benchmark 路线及 `4/1/16/4/1/1` case cardinality；
- GA facade → `OTKernel.annotate`；
- standard sc-SVC 两次 `LocalAnchoringKernel → OTKernel.annotate`；
- sp/SR/imputation LR problem construction → `OTKernel.couple`；
- preprocessing、GA、LR、finalize、publication/evaluation 的顺序；
- conditional skip/graph aggregation/imputation pair outputs；
- write object 与 return object 语义一致。

静态 boundary tests 同时检查：

- `OTKernel` 外没有三种直接 solver 调用；
- LR/runners 中没有 `GlobalAnchoringKernel`；
- `OTKernel` 不导入 runner、runner_conf 或 stage kernel；
- 无旧 `ApplicationExecution`、`--confounding`、noise symbols 或 `revise.yaml` packaging entry。

## Legacy `test_sc/test.py` comparison checkpoints

在相同输入、依赖版本和 seed 下，按顺序比较；第一处不一致立即停止并保存左右证据：

1. `Patient == P2CRC` 后的 reference obs index/count；
2. spatial `transcript_counts >= 60` 后的 obs index/count；
3. spatial/reference `filter_genes(min_cells=100)` 后的 var index；
4. slash label normalization 后 Level1/Level2 值；
5. Level1 posterior、category axis、labels，以及旧实现 `obs["max_score"]` 对新包 `obs["Confidence"]`；
6. selected T spatial/reference obs 与各自 var axis；
7. Level2 posterior、category axis、labels，以及该次 annotate 覆盖后的旧 `max_score` 对新 `Confidence`；
8. gene/spatial/joint graph、各 resolution Leiden partition、spatial/alignment score、best resolution；并列取最后一项；
9. spatial `SVC_cluster`；
10. reverse annotation posterior、category axis、labels，以及 expression reference 上旧 `max_score` 对新 `Confidence`；
11. 最终 spatial/expression `X`、obs/var axes、公开 obs columns 和关键 metadata。

## Comparison rules

- obs/var/category axis 与离散 labels：精确比较。
- dense/sparse numeric arrays：转到共同可比表示后使用 dtype-aware `allclose`；容差必须记录，不能为通过测试临时放宽。
- Leiden：比较 partition equivalence；cluster ID 名称可整体置换。
- graph：shape、support/indices 精确比较，weights 使用 dtype-aware `allclose`。
- metadata：只比较契约字段；package provenance 的新增字段属于允许的形式差异。
- confidence 映射：旧脚本每次 `annotate` 写入/覆盖 `obs["max_score"]`，新包在对应 target 和对应时点写入/覆盖 `obs["Confidence"]`。比较前先精确对齐 observation index，再用 dtype-aware `allclose` 比较数值；不能把不同 annotation 时点的 confidence 互相比。新包公开 schema 始终保留名称 `Confidence`，不为 parity 输出 `max_score` alias。
- TACCO target axis 保持完整；不得在比较层自动过滤、重排或重新归一化。
- expression 输出保留完整过滤后 reference var axis；spatial 输出使用 overlap axis。

## Focused synthetic evidence

真实双跑前先证明：

- 三个 preprocessing 函数的默认、override、缺列/零匹配失败和调用顺序；
- TACCO fresh key、reference copy、returned category order 和 posterior no-normalization；
- POT dtype、route-owned row normalization；
- TACCO hard support、feasibility 与 marginal continuation；
- output.name 四种命名和 paired atomic publication；
- effective config/provenance/hash 与 5+6 template distribution。

## P2CRC real run protocol

1. 冻结当前 commit、Python、TACCO/POT/Scanpy/Squidpy/Leiden 版本、线程环境和 random seed。
2. 对 `test_sc/test.py` 增加只写临时 evidence root 的 checkpoint 输出；不修改/跟踪 `test_sc/` 原文件。
3. 旧脚本与 package 使用同一原始 spatial/reference，分别写入全新输出目录。
4. 每完成一个 checkpoint 即比较；第一处 divergence 停止后续科学结论，保留输入轴、shape、dtype、摘要和必要小矩阵。
5. 若十一步通过，再比较最终 H5AD reload 后 contract。
6. 在 [09-delivery-status.md](09-delivery-status.md) 分别记录 structural、synthetic、real numerical 和 scientific 状态。

## Proof boundary

route trace 通过只证明调用线路；synthetic parity 只证明局部契约；P2CRC 十一步通过只证明该环境/样本下的技术与数值对齐。任何生物学有效性、泛化或 solver 优劣仍是 Scientific unknown。

继续阅读：[Implementation and risk](07-implementation-and-risk.md)。
