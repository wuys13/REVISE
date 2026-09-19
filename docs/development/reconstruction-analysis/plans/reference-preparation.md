# Reference Preparation：独立准备与显式消费

[返回总览](../README.md) · [R7/R8 背景](../categories/04-future-extensions.md#r7)

用户已对齐本计划并授权整个接入批次。工程验证结果在本页末尾记录；真实规模和 reference 生物学质量不由工程测试证明。

## 两步之间的关系

第一步是独立的准备支线：已知配对时精确提取 reference；需要比较时对本地候选逐个运行 GA，生成 `reference.yaml` 和准备报告。第二步由调用者显式启动原有重建，只将该文件作为 reference 覆盖输入。没有外部配置时保持旧行为。

```text
paired 提取 / screen 候选 GA + 排名
                 ↓
        reference.yaml + report.json
                 ↓ 显式传给单样本或 batch
         原有 GA → LR → 样本发布
```

screen 对 N 个候选各运行一次 GA；正式重建对选中 reference 再运行一次。这是首版已接受的计算成本，不新增 GA 缓存、恢复或完整 posterior 交付协议。

## 第一步：生成 reference.yaml

从 REVISE 仓库根目录，使用已安装项目依赖的 Python：

```bash
python -m revise.reference_preparation --config /absolute/path/preparation.yaml
```

已知配对：

```yaml
schema_version: 1
mode: paired
source: ./reference-pool.h5ad
pair_column: pair_id
pair_key: "PAIR_001"
output_dir: ./prepared/PAIR_001
```

`pair_column` 和 `pair_key` 是明确的列和值，不推断 Patient/Sample 含义。精确匹配，不自动 trim 或转数字；输出保留所选 AnnData 行的其他内容。

候选筛选：

```yaml
schema_version: 1
mode: screen
reconstruction_config: ./application.yaml
candidates:
  directory: ./candidates
output_dir: ./prepared/screen_001
```

目录模式只看直接子级 H5AD。也可用 `candidates: {list: ./candidates.yaml}`；清单为 `schema_version: 1` 和 `candidates: [{id: candidate_a, path: ./candidate_a.h5ad}]`。准备配置、清单及 reference 配置中的相对路径分别相对于声明文件；Application 自身的 `root_dir` 规则保持不变。

`output_dir` 必须是新目录。paired 生成 `reference.h5ad`、`report.json`、`reference.yaml`；screen 生成报告和指向外部选中文件的配置，不复制候选。失败候选记录原因后继续；全部失败或无候选时不生成 reference 配置。再次运行使用新的输出目录。

## 第二步：原有重建显式使用结果

```bash
python reconstruct.py --config /absolute/path/application.yaml \
  --reference-config /absolute/path/prepared/screen_001/reference.yaml
```

Python 调用为 `run_application(application_path, reference_config=reference_yaml_path)`，原有返回类型保持不变。

batch 的样本或继承配置可写：

```yaml
inputs:
  reference:
    config: ./prepared/screen_001/reference.yaml
```

继续使用现有 batch 命令。该块与 `path`、`format`、filter 互斥；低层 reference 整块覆盖继承值，路径相对于声明它的 batch YAML。batch 不运行准备，不自动寻找 reference。

外部配置在原 reference 校验前替换 reference 路径、格式和 filter；不会重复执行旧筛选。预处理、注释列、route、solver、输出和 seed 仍由原 Application/batch 配置决定。无效外部配置直接报错，不回退旧 reference。

## 筛选分数与证据

三种分数均先逐 ST unit 计算，再取中位数：

| 分数 | 单元量 | 作用 |
|---|---|---|
| `max_median` | 最大类别概率 | 固定选择第一名 |
| `certainty_median` | 1 − entropy / log(类别数) | 报告分数及排名 |
| `margin_median` | 最大减次大概率 | 报告分数及排名 |

同分按候选 ID 升序。GA 概率必须有限、非负、逐行归一且至少两类。任何未分配行使该候选失败，不补概率、不删行、不换 solver。候选必须对应相同 ST 单位和顺序；类别覆盖、基因覆盖和 reference 细胞数允许不同并记录。高分只表示该次 GA 输出较集中，不表示准确率或生物学最优。

准备报告保存受控输入摘要、有效参数摘要/solver/seed、轴摘要和候选组成；不保存任意运行对象或完整 GA 矩阵。正式重建将外部配置、准备证据与最终输入身份纳入 provenance 和 batch 复用判断。消费工具产物时检查报告中的选中文件摘要，变化即拒绝；手写最小 reference 配置可使用，但不声称存在筛选证据。

当前 Confidence 列不变。筛选从 GA 概率直接计算指标，不读取可能被 LR 更新的最终 Confidence。

## 实施与验收

内部顺序为材料包纳入、GA-only 适配、单样本/batch 覆盖、联合验证。测试覆盖原材料包、真实小型 GA、候选失败与轴差异、摘要变化、旧路径兼容、batch 指纹和样本发布回归。不得只用模拟 callback 声称真实 GA 接线成功。

2026-09-19 工程验收完成：主回归 **917 passed、44 warnings、24 subtests passed**，无失败或跳过；独立真实 screen 命令入口检查 **1 passed**。合计 918 项测试加 24 个子用例。主回归包含配置、Application、backend、batch、recon、reference preparation，以及已有分析交接和完整样本交付集成测试。原材料包 33 项测试保留；命名空间迁移后隔离验证通过。

真实小型 POT 与 TACCO GA 均已执行，另以两个候选的真实 POT GA 验证准备命令生成可消费配置；paired → direct/batch → Raw/SVC/sample.yaml 的测试使用确定性重建计算替身，验证编排、产物与复用。没有运行真实全量重建，不据此证明规模性能或 reference 生物学质量。

初次整合回归暴露了新诊断回调未被测试替身调用、旧入口签名断言未更新，两项已修正后重跑通过；未更改科学算法。测试环境禁用可选 readline 和 pytest 插件自动加载，CLI/OpenMP 子进程在沙箱外运行。warnings 主要来自依赖弃用、POT 提示及既有边界 fixture。

可追溯证据见 [验证摘要](../verification-reference-preparation-2026-09-19.json)、[实现索引](../implementation-index.md) 与 [来源快照](../../../../revise/reference_preparation/SOURCE_PROVENANCE.md)。未提交、推送或修改 Atlas/分析库。

复现工程验收（在 REVISE 根目录、项目 Python 环境中）：

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
MKL_NUM_THREADS=1 NUMBA_DISABLE_JIT=1 MPLCONFIGDIR=/tmp/revise-ref-mpl \
NUMBA_CACHE_DIR=/tmp/revise-ref-numba python - <<'PYTEST'
import sys
sys.modules["readline"] = None
import pytest
raise SystemExit(pytest.main([
    "tests/application", "tests/config", "tests/backend", "tests/batch",
    "tests/recon", "tests/reference_preparation",
    "tests/integration/test_analysis_delivery.py",
    "tests/integration/batch/test_whole_sample_delivery.py", "-q",
]))
PYTEST
```

当前完整命令也包含最后单独执行的 screen 入口测试，预期总项数为 918。

自动公开数据检索、下载、reference 大小公平校正、新的 Confidence、cluster assembly 和大样本科学比较均继续后置。
