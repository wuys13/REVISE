# R5/R6/D1/C1：OT assembly 与分析配套

[事项与证据索引](../implementation-index.md) · [总览状态](../README.md#items)

**用户计划已对齐；生产实现与工程验收完成，真实样本及科学比较未运行。** 本批只实现与工程验证，不运行真实全量重建，不宣称科学收益，不修改 Analysis Agent。真实样本由用户逐项检查后另行选择。

## 为什么做与确定边界

在统一 ST 轴上用 OT 混合 SC 表达，比随机 donor 和固定 cluster mean 提供更细的表达分配。SC 的 Level2 监督分群由外部准备；REVISE 遍历实际 broad types，有效不同 Level2 超过一个才重建，否则汇总提示补齐后重跑。合格类型计算失败仍使整个样本失败。

| 阶段 | 本批确定行为 |
|---|---|
| 空间处理 | 匹配用 ST 副本上，同 SVC_cluster 的空间邻居合计 0.2，自身 0.8；无邻居保留自身 |
| within_cluster | ST 单细胞匹配同群 SC donor 单细胞 |
| outside_cluster | ST 单细胞匹配同 broad type 全部 SVC_cluster mean profiles |
| OT 与输出 | TACCO；按行归一化全部有效运输权重，乘 SC 全基因表达，不用 top-1 或反向权重 |
| 旧能力 | 保留 mean/random/paired 和现有 LR；默认仍为 random |
| 科学边界 | 不把空间平滑视为改善证明；不增加两图结构迭代优化 |

“outside”表示候选 cluster 范围，不是 out-of-spot。两条路径都重新计算 ST × donor/profile 的 OT；旧 SC × SVC_cluster 注释矩阵不能代替它。空间混合系数与原建图 alpha 独立。

## 配置与运行入口

在已有 Application YAML 中修改：

```yaml
output:
  dir: results/assembly/within
  ist_mapping: within_cluster  # 或 outside_cluster
  ist_ot:
    spatial_weight: 0.2
    max_cost_entries: 2000000
    gene_block_size: 256
```

同一参数块可用于 batch 的 `output`。`ist_ot` 只接受两个新模式。质量沿用 Benchmark：ST 端取共同基因表达总量，SC donor 取对应总量，outside profile 的质量为该群 donor 共同基因总量之和（包含 reference 群大小影响，不是等 profile 配额，也不是经过校准的组织比例）。

单次完整 OT 超出矩阵元素上限会在分配代价矩阵前失败，提示改用 outside_cluster 或显式调大上限，不自动替换方法，也不拆行独立求解改变质量约束。表达乘法直接按基因块写入最终矩阵，避免完整群输出副本。代价上限不限制最终 dense SVC（float64）的大小；运行记录保存 `ot_output_bytes`，真实规模运行前仍需按细胞×基因数评估内存。新 assembly 的 OT 固定 TACCO，上游 GA/LR 仍遵循原 `algorithm.ot_method`。

```bash
python reconstruct.py --config application.yaml
python reconstruct.py --config application.yaml --reference-config prepared/reference.yaml
```

各模式使用独立输出目录进行后续比较。常规发布事务保护已有结果，配置和代码身份参与 batch fingerprint。无共同基因、缺失候选、非法数值/权重或无运输质量明确失败；不隐式填均值。

## D1：confidence 与 reference 评分

公共数值函数置于中立工具模块。现有 Confidence 保留当前 annotation 分配的最大权重；max_median 对这些值取中位数，certainty/margin 保留。分阶段记录类别与来源，不用 assembly donor 权重覆盖 annotation Confidence。不同阶段的候选含义不同，不宣称数值可直接跨阶段比较。不新增 GA/LR 拆列、完整 posterior 持久化或 D2 方法修改。

## C1：明确输入来源，按已定线性契约交付

源表达声明是用户对处理历史的明确说明，不是由数值自动猜测。Application 可填写：

```yaml
inputs:
  st:
    path: raw_data/spatial.h5ad
    format: h5ad
    expression:
      identity: measured_expression
      scale: untransformed_nonnegative
  reference:
    path: raw_data/reference.h5ad
    format: h5ad
    expression:
      identity: measured_reference_expression
      scale: untransformed_nonnegative
```

这里的线性非负包括可信的非 log normalized 表达和小数，不等同于原始整数 counts。声明已知时还检查 X 有限且非负；log 声明拒绝。未声明则保持 unknown。源对象也可显式携带 `uns['revise_expression'] = {'identity': ..., 'scale': 'untransformed_nonnegative'}`；配置优先。

batch 的 ST 路径仍自动发现，允许 `inputs.st.expression`，不允许填写 ST 路径；直接 reference 支持同样 expression 字段。`inputs.reference.config` 仍与直接 reference 字段互斥。外部 reference 覆盖时舍弃旧 reference 声明，读取实际选中文件的显式 `uns` 记录；prepared paired 保留该记录，不伪造来源。

生产端按实际来源声明 Raw 与 SVC；已知线性数据不再携带旧 `scale: unknown` 阻塞消费者。sST 保留内部 normalization 与 parent-spot 校正说明，不宣称 raw counts；未知输入仍不能放行。分析库目前固定 finite/nonnegative/unlogged linear X，旧 log 声明已不适用；本批不修改独立消费者。

## R6：比较 Notebook

入口：[assembly_comparison.ipynb](../../../../reproduce/case/assembly_comparison.ipynb)。读取已生成产物，不自动重建。参数集中配置四种方法路径、原空间 subtype 基线、T/Macro/CAF 标签映射、seed 和 Leiden resolutions。对相同真实 ID/共同基因范围独立重建表达邻居图与 Leiden；原 cluster 不用于新分群。显示覆盖/缺失、空间图、列联表、ARI/NMI，保留各 resolution，不自动选择最高分或改默认方法。

## 验收与后续

工程验证：空间混合、两种权重方向/范围、数量不等、profile 均值、TACCO 真实小 fixture、资源上限、失败回滚、direct/batch、confidence 旧数值、真实分析 loader/表达步骤、Notebook fixture 执行与 R1 回归。最终 **964 passed、24 subtests passed、0 failed、0 skipped**，49 warnings，退出码 0。详情见[验证记录](../verification-ot-assembly-2026-09-19.json)。

D2 保持原有方法，不列本批修改；D3 out-of-spot/复杂映射未实现继续待做。真实规模、T/Macro/CAF 科学对照和最终默认值均后置，无用户必交新材料。

## 实际交付与验证复现

| 分类 | 结果 | 证据边界 |
|---|---|---|
| assembly | 两模式、配置、direct/batch、事务回滚通过 | 小 fixture 包含真实 TACCO；非真实大样本 |
| confidence | 共用数值与阶段来源通过 | 旧 Confidence 语义及 reference 排名不变 |
| Notebook | 副本执行、Leiden、图表和输入哈希通过 | 未跑真实 T/Macro/CAF 科学对照 |
| C1 | 正式消费者成功计算两侧 Moran 表达分析 | 线性浮点 fixture，不代表所有分析或所有真实输入 |
| 文档 | 本计划、索引、状态、历史修正与验证记录更新 | D2/D3 不误报为本批实现 |

执行中修复了缺失 subtype 的虚假声明、无效坐标 key、及新 source binding 对既有 orchestration mocks 的适配；未修改科学算法或真实数据来绕过失败。OpenMP 沙箱共享内存失败通过沙箱外同一 fixture 重跑解决，测试进程禁用可选 readline 避免本机导入崩溃。warnings 来自依赖弃用提示及既有边界 fixture，不是静默跳过。

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
MKL_NUM_THREADS=1 NUMBA_DISABLE_JIT=1 MPLCONFIGDIR=/tmp/revise-ot-mpl \
NUMBA_CACHE_DIR=/tmp/revise-ot-numba \
/Users/stephen/miniconda3/envs/python3.10/bin/python - <<'PYTEST'
import sys
sys.modules['readline'] = None
import pytest
raise SystemExit(pytest.main([
    'tests/application', 'tests/config', 'tests/backend', 'tests/batch',
    'tests/recon', 'tests/reference_preparation',
    'tests/analysis/test_assembly_comparison.py',
    'tests/integration/test_analysis_delivery.py',
    'tests/integration/batch/test_whole_sample_delivery.py', '-q',
]))
PYTEST
```

在 REVISE 仓库根执行；需要相邻分析库，其他路径用 `REVISE_ANALYSIS_ROOT` 指定。缺少消费者时会明确 skip，不能沿用本次联合验证结论。Notebook 需可用的 python3 Jupyter kernel。本批未 commit/push。
