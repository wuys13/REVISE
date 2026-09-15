# Batch code framework

本文维护批量框架的代码职责、调用接口和恢复机制；数据字段与结果结构只在[协议](protocol.md)定义。[总入口](README.md#contract-and-status) · [实施依据](plan.md#implementation-order) · [验收](plan-validation.md#acceptance-scenarios)。

**本轮框架已实现并通过验收。** 具体测试范围与结果见[交付验证](plan-validation.md#verified-delivery)；科学分析算法仍由后续方面模块提供。

## Target architecture

```text
project batch.yaml + ST/<sample>/spatial.h5ad + shared reference
  → config: 分层配置解析；批量额外发现样本
  → sample: 只读协议检查
  → 重建执行器 → 既有单次重建引擎 → reconstruction.json
  → 分析执行器 → 按需输入视图 → 一个方面模块 → analysis/<aspect>/
```

单项和批量共用任务执行器。批量负责发现、展开和汇总，不另写一套重建或分析流程。hST／sST 每样本一个重建任务，iST 每样本 × 配置类型一个任务；一个分析任务对应其中一个重建任务的一个方面。

## Execution and module boundaries

| 职责 | 代码与关键符号 |
|---|---|
| 发现、继承及声明位置相对路径 | [config.py](../../../revise/batch/config.py)：`discover_samples`、`resolve_sample` |
| 标准 H5AD 只读检查与输入摘要 | [sample.py](../../../revise/batch/sample.py)：`read_sample`、`file_identity` |
| 重建执行、状态和交接 | [runner.py](../../../revise/batch/runner.py)：`run_batch`、`run_reconstruction_task`、`_run_task`、`verify_reconstruction_task` |
| 原有单次重建引擎 | [reconstruct.py](../../../reconstruct.py)：`run_application`；保留 QC、reference 筛选、基因对齐及 GA/LR |
| 分析输入辅助 | [inputs.py](../../../revise/batch/inputs.py)：`AnalysisInputs`、`InputView`、`PairedView`、`SSTBaseline` |
| 分析执行与发布 | [analysis.py](../../../revise/batch/analysis.py)：`run_analysis_batch`、`run_analysis_task`、`_run_aspect`、`AnalysisContext` |
| 公共 Python API | [__init__.py](../../../revise/batch/__init__.py)：两个单任务入口，延迟导入执行模块 |
| 两个批量命令 | [cli.py](../../../revise/batch/cli.py)：`main`、`analysis_main` |

重建执行器通过 `application_document` 转换到已有单次配置和路径契约，再调用子进程。分析模块在当前 Python 进程内运行，是可信项目代码，不是隔离沙箱。执行默认顺序进行，不包含跨类型 GA 缓存、并行或分布式调度。

## Task API

从 `revise.batch` 导入：

```python
run_reconstruction_task(config_path, sample_id, *, cell_type=None)
run_analysis_task(config_path, sample_id, aspect, *, cell_type=None)
```

| 参数 | 含义 |
|---|---|
| `config_path` | schema-2 项目 YAML；数据与输出路径由配置解析 |
| `sample_id` | 相对输入根的精确样本路径，例如 `CRC/S01` |
| `cell_type` | iST 必须指定一个已配置类型；hST／sST 必须省略 |
| `aspect` | 精确选择一个已配置分析方面，不执行其他方面 |

返回可序列化任务记录，包含状态、任务位置和实际产物信息；执行失败记录原因。错误的项目配置或非法任务选择抛出调用错误。单项调用不改写全批次汇总来声称其他任务已执行。

单项重建成功后可以直接单项分析，不要求先生成 `batch_status.json`。分析使用当前配置和实际重建状态；重建交接记录保留重建当时的来源事实，不作为当前分析配置来源。

公共单项入口和批量入口各获取一次输出根锁，内部执行器不重复加锁。所有路径选择都经共同的任务检查，防止单项绕过批量的目录和类型约束。

## Connect an analysis module

复制[示例模块](../../../configs/batch/analysis_example.py)，把[示例配置块](../../../configs/batch/analysis-example.yaml)合入项目或下级 `batch.yaml`，确保模块所在目录在 Python 导入路径中。该示例只生成载体清单，不是科学比较。

分析方面、参数、资源文件和返回结构的唯一详细定义见[分析配置](protocol.md#analysis-configuration)与[分析结果](protocol.md#analysis-results)。没有配置方面时不运行分析；明确请求而未提供实现时报告 `not_implemented`。

### Scientific adapter entrypoints

当前源码提供以下接入点。完整科学定义与待决选择仍在 [reconstruction-impact](../reconstruction-impact/README.md)，测试 fixture 的小规模参数不能直接作为真实数据推荐值。

| 方面 | 入口与必需配置 | 最小调用证据 |
| --- | --- | --- |
| `reconstruction_impact` | `revise.analysis.reconstruction_impact_adapter:run`；精确 `level1_column` / `parent_labels`、`partition_change`、`spatial_region`（需物理单位换算）、`raw_level2_mapping`；资源 `raw_level2_reference` | [真实计算与 batch handoff fixture](../../../tests/batch/test_scientific_aspects.py) |
| `pathway_activity` | `revise.analysis.pathway_activity_adapter:run`；显式资源 name/species/gene_id_type/version、coverage/基因支持阈值、seed、`provider_detected_gene_quantile` 与 `auc_threshold_quantile`；资源 `gene_sets` 为 pathway → gene ID 列表的 JSON | [AUCell 接入 fixture](../../../tests/analysis/test_pathway_activity_adapter.py)；provider 属可选 pathway 依赖 |

这两个方面设置 `requires_pairing: true`。impact 的 iST 输入为 paired spatial carrier；pathway 通过完整表达视图支持 paired 投影与原生 mean/random。Moran 当前只有接收外部权重图的数值比较核心，尚无可配置的 batch aspect；空间图和支持规则仍需确认。历史重建不能直接伪装成已验证 handoff。

### Render a published impact report

在源码 checkout 中调用 `revise.analysis.impact_report.render_impact_report(task_roots, output_dir)`，输出位置使用 task 或 sample 下的 `reports/reconstruction_impact/`。报告验证方面登记的 CSV/JSON 和内容摘要，在独立目录执行固定 [source notebook](../../../reproduce/case/reconstruction_impact/Reconstruction_Impact_Report.ipynb)，登记输入、执行 notebook 与图。它不重新读取 H5AD 或计算科学指标；最新任务失败时，旧有效产物明确标为 stale。

当前报告模板属于仓库，未随 wheel 分发；仅安装 wheel 的环境不能渲染此报告，需使用源码 checkout。`execute=False` 只准备可执行报告及输入，其 manifest 明确记录 execution skipped，不是执行验证。

## Aspect modules

`AnalysisContext` 提供：

| 字段 | 用途 |
|---|---|
| `reconstruction` | 已验证的重建交接记录，含输出角色、输入来源、坐标和原生配对说明 |
| `output_dir` | 本方面的 staging 目录，模块只向这里写产物 |
| `parameters` | 有效方面参数；模块返回实际使用的默认值和处理方式 |
| `resources` | 名称到已解析 `Path` 的映射；内容摘要由框架记录 |
| `inputs` | 按需建立的 `AnalysisInputs`；具体载体在调用对应方法时读取 |

`inputs.raw()` 读取原始 ST；`aligned_raw()` 提供 hST/iST 对齐视图；`paired()` 提供 iST donor 与 cluster 映射；`sst_baseline()` 构造父 spot 均分基线。`reconstructed()` 用于 hST/sST 或 iST mean/random，paired 下必须明确使用双载体视图。完整语义见[分析输入](protocol.md#analysis-input-views)。

一个模块负责本方面的两侧计算、必要的对象匹配、比较及产物写出。科学筛选、统计前提、指标含义和预处理顺序由模块负责。`requires_pairing: true` 要求原生观测配对；sST 均分基线不满足这一前提，支持基线的方法应明确调用 `sst_baseline()` 并记录比较依据。

## Identity and recovery

重建与分析使用同一个配置解析器，但相关依赖分别形成复用身份，见[复用规则](#reuse)。文件存在、非空目录或旧成功状态均不足以证明当前结果有效。

## Reuse

| 任务 | 复用依赖 |
|---|---|
| 重建 | 项目归属、样本／类型、重建有效配置、实际输入摘要、重建代码和运行依赖身份、发布产物校验 |
| 每个分析方面 | 相关方面配置、资源摘要、模块与公共分析代码、raw 和有效重建输入、发布产物校验 |

分析参数、资源或实现变化只影响相关分析；不触发重建。重建输入或有效配置变化使重建及依赖分析失效。框架按 adapter 与公共分析源文件识别实现；adapter 使用模块级 `CODE_DEPENDENCIES = ("package.module", ...)` 显式登记所调用的共享算法模块，其源文件摘要进入本方面指纹。框架不自动追踪完整的传递调用图；模块通过 `version` 明确标识未被自动覆盖的外部依赖变更。首次执行前，框架按当前源文件重新加载 adapter 及显式声明的依赖，避免先前导入的旧对象被登记为新代码结果；模块顶层代码须可重复导入。之后若 adapter、声明依赖或框架源文件变化，任务明确失败并要求重启解释器，再用新源码重跑；不会把缓存旧函数登记为新实现的结果。

每次重新解析沿途 YAML，包括新增和删除的文件。原文件摘要链用于审计，复用按相关有效配置判断；注释、空白或不改变有效值的覆盖不重算。旧指纹配方不会直接视为新成功缓存。

输入、配置及实现身份在结果接受前重新检查；分析资源也检查内容。发现运行中变化时拒绝把新产物登记为成功，稳定后重试。

## Recovery

任务日志和控制记录位于 `.revise/`；各方面状态与日志位于 `.revise/analysis/`。重建记录和分析记录分别表示完成情况，重建成功不等于分析已经运行。

模块结果先写 staging，检查路径归属与结果记录后整方面发布。发布失败恢复旧目录；失败或中断不成为可复用成功。模式切换沿用已有 iST 发布和受控旧文件清理，用户自有文件不清理。

批量中单项失败不阻断后续任务；不完整批次返回非零状态。禁用或移除的任务保留科学文件，但不得继续标为当前结果。单项分析汇总只说明本次选择的执行情况，不替未验证的其他方面背书。

## Verification and development boundaries

测试与场景对应统一见[验收文档](plan-validation.md#acceptance-scenarios)。其中单元测试检查输入与状态边界，真实小样本对照检查 hST/iST/sST 与原有单次引擎一致。

本框架不执行数据格式转换或标签 harmonization；科学算法由显式配置的方面 adapter 执行。示例 adapter 和调用回归证明框架可接入，不证明生物学结论有效。
