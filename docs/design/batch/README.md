# Batch reconstruction and analysis

这是批量重建与分析的统一入口。先确认下列输入输出结构和实现状态，再按问题读取协议或框架；代码开发顺序集中在计划中。

## Contract and status

**最终结构已经确定：标准 ST 目录作为输入，分层 YAML 管理配置，SC 通常按癌种分目录并通过路径共享；独立输出树保存重建、各方面分析和运行记录。** 本轮不增加样本清单、数据副本或另一套配置系统。

标准目录样式统一维护在[输入输出示例](input-output-example.md)：[输入与 SC 癌种分类](input-output-example.md#input-tree)、[iST 输出](input-output-example.md#ist-output-tree)、[hST/sST 输出](input-output-example.md#hst-and-sst-output-tree)。本页说明实现状态，不重复目录树。

| 约定 | 当前实现 | 当前边界 |
|---|---|---|
| ST、分层 YAML、reference 共享与输出布局 | 已实现 | 保留；不新增输入目录层级 |
| `paired / mean / random` 重建输出 | 已实现 | 保留文件语义；分析通过按需视图读取 |
| iST whole-sample `Raw + SVC + sample.yaml` delivery | 已实现工程交接 | GA/LR 全类型规则、Raw 保真和 sample-level 发布已接入；P2CRC 真实样本及科学/联合验收另行记录 |
| raw／paired／sST 分析输入 | `AnalysisInputs` 延迟加载 raw、重建、paired 和 sST 基线 | 辅助视图不做归一化或科学计算 |
| 单项与批量调用 | 两个批量命令和两个 Python 单任务接口已实现 | 单项调用只处理明确选择的样本、类型或方面 |
| 分析方面与返回记录 | 显式方面、资源、结构化产物和计算依据已实现 | 科学算法仍由方面适配器负责 |
| 恢复 | 重建与分析按相关指纹分别复用 | 输入、有效重建配置或代码变化会使依赖任务失效 |

外部调用方提供标准文件和项目配置，消费 `reconstruction.json` 与分析记录，不依赖临时文件或扫描文件名猜测状态。实际科学算法由方面模块提供；框架接入不等于算法已完成。

## Read by question

| 要确认什么 | 唯一定义位置 |
|---|---|
| 数据适配可直接参照的目录与配置 | [标准输入输出示例](input-output-example.md#input-tree)、[YAML 文件](input-output-example.md#configuration-files) |
| 输入目录、YAML 继承、reference 路径 | [protocol：配置](protocol.md#configuration-and-inheritance)、[标准数据](protocol.md#standard-st-and-reference-files) |
| 输出文件、分析视图与结果记录 | [protocol：输出](protocol.md#output-roles-and-current-result-semantics)、[分析输入](protocol.md#analysis-input-views)、[分析结果](protocol.md#analysis-results) |
| 单项／批量怎么调用，代码怎么组织 | [framework：目标结构](framework.md#target-architecture)、[单任务接口](framework.md#task-api) |
| 模块怎样接入，哪些变化需要重跑 | [framework：模块](framework.md#aspect-modules)、[复用](framework.md#reuse) |
| 先开发什么、为什么这样取舍 | [plan：实施顺序](plan.md#implementation-order)、[参考取舍](plan.md#reference-decisions) |
| 怎样验收 | [plan-validation：场景与测试](plan-validation.md#acceptance-scenarios) |

本目录由五篇说明文档和一份独立输入输出示例组成。原 `plan-inputs.md`、`plan-execution.md` 的内容已分别归并到协议与框架，不再并行维护。各章节明确标注实现状态；以下命令与恢复说明仍按**当前代码**书写。

## Set up and run

Copy the example tree under [configs/batch/](../../../configs/batch/) to your project. Supply real
`ST/CRC/S01/spatial.h5ad`, `ST/CRC/S02/spatial.h5ad` and reference files in
`SC/CRC/`. The template directories alone are not samples: discovery requires
`spatial.h5ad`. Remove the S02 override if it should share the CRC reference.
Review modality, coordinate units, reference label columns and protocol settings;
template QC thresholds are examples, not universal defaults.

目录概览见[标准示例](input-output-example.md#input-tree)，配置细节见[协议](protocol.md#configuration-and-inheritance)。

Install this checkout in the reconstruction environment with
`python -m pip install -e .`, then run:

```bash
revise-batch-reconstruct --config /data/project/batch.yaml
revise-batch-analyze --config /data/project/batch.yaml

# Equivalent source-checkout entrypoints:
python /path/to/REVISE/batch_reconstruct.py --config /data/project/batch.yaml
python /path/to/REVISE/batch_analyze.py --config /data/project/batch.yaml
```

The launch directory does not affect configured paths. See
[Batch input and output protocol](protocol.md#batch-input-and-output-protocol) for fields and inheritance, and
[Batch code framework](framework.md#batch-code-framework) for execution and adapter integration.

For one explicitly selected task, the Python API uses the same execution path:

```python
from revise.batch import run_analysis_task, run_reconstruction_task

run_reconstruction_task("/data/project/batch.yaml", "CRC/S01", cell_type="T")
run_analysis_task(
    "/data/project/batch.yaml", "CRC/S01", "carrier_inventory", cell_type="T"
)
```

The analysis aspect must be present and enabled in the effective YAML. Single-task
calls do not need a batch inventory; input views and structured result fields are
defined in the [protocol](protocol.md#analysis-input-views) and
[analysis-result](protocol.md#analysis-results) sections.

## Resume and inspect

Rerun the same command to resume. Successful tasks are reused only when their
relevant input identities, effective configuration, owning project and source
paths, runtime identity and artifacts still match. The complete configuration
chain is retained for audit and is not itself a reuse key; see the
[framework reuse rules](framework.md#reuse). Reconstruction and analysis have
separate reuse boundaries: an analysis-only setting or resource change reruns
the affected aspect without invalidating reconstruction, while an effective
reconstruction setting, input, or reconstruction-code change reruns
reconstruction and its dependent analysis. Changes to comments, whitespace, or
overrides that do not change effective reconstruction settings do not force
reconstruction. Failed or interrupted tasks rerun. Disabled samples and
removed iST types become inactive without deleting their scientific files.

Analysis adapter source changes are checked before publication. Because adapters
run in the current interpreter, restart that interpreter before retrying after a
source edit so an old imported module cannot be reused.

Consult `output/batch_status.json` and each task's `.revise/reconstruction.log`.
A failed rerun can leave older scientific files in place; require a current
successful `reconstruction.json` instead of checking filenames alone.
`output/analysis_status.json` and `analysis/analysis.json` report analysis
separately. An empty analysis directory never establishes analysis completion.

Reconstruction exits 0 on success, 1 for task failures and 2 for invalid invocation.
Analysis exits 0 only when requested tasks succeeded or were reused; missing
adapters and unmet prerequisites are incomplete outcomes. Execution is sequential;
one task failure does not stop subsequent tasks. Concurrent runs against the same
output root are rejected.

## Migrate the experimental input protocol

Schema version 2 replaces schema version 1, per-sample `sample.yaml` discovery
and `revise-prepare-sample`. Old batch configurations produce a migration error.
Single-run reconstruction configuration remains supported.

1. Choose separate input and output roots and set `schema_version: 2` in the
   project `batch.yaml`.
2. Put each already-standard ST file at `<input_root>/<sample>/spatial.h5ad`.
   Select existing standardized files deliberately; this release does not move
   or transform originals.
3. Move common reconstruction settings into the project or group `batch.yaml`.
   Move `sample.modality` to `modality` and coordinate metadata to
   `coordinates`. Reference paths can point to one shared file.
4. Keep sample overrides only where needed. Remove `inputs.st` and
   `preparation` settings from the new configuration. Convert nonstandard data
   before batch execution; the new entrypoint does not perform those conversions.
5. Run the new configuration. Old cached successes are not accepted directly as
   schema-2 successes. Existing data, preparation files and results are never
   automatically relocated or deleted.

The implemented scope is reconstruction orchestration, lazy analysis input views,
and structured analysis adapter handoff.
A separate data processing stage is planned but not implemented here.
Reconstruction-impact algorithms remain in their development branch; an aspect
without an adapter reports `not_implemented`.
