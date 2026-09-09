# Batch input and output protocol

本文统一定义批量输入输出协议，包括已实现的 schema 2、按需分析输入视图和结构化分析结果。协议中的文件、字段和状态语义以当前代码为准；科学算法仍由方面适配器负责。不再另维护一份计划输入协议。

[总入口与状态](README.md#contract-and-status) · [代码框架](framework.md) · [开发顺序](plan.md#implementation-order) · [验收](plan-validation.md#acceptance-scenarios)

| 协议部分 | 实现状态 |
|---|---|
| 分层配置、标准 ST、共享 reference、重建输出与控制目录 | 已实现，本轮保留 |
| 分析方面配置 | 已实现，方面、资源、启用状态和单项调用按显式配置执行 |
| raw／paired／sST 分析辅助视图 | 已实现，按需加载 |
| 分析结果返回结构与计算依据 | 已实现，适配器必须返回结构化产物和计算依据 |

## Configuration and inheritance

**已实现，本轮保留。**

The project `batch.yaml` declares `schema_version: 2`, `input_root` and
`output_root`. These keys are project-only. Roots must be separate and must
not contain one another. The input tree contains ST samples; shared SC files may
live beside that tree or elsewhere and are referenced by path.

SC 通常按癌种组织，例如 `SC/CRC/`、`SC/BRCA/`；这只是数据组织方式，不引入 SC 配置继承或自动 reference 匹配。标准样式与可复制 YAML 见[独立示例](input-output-example.md#configuration-files)。

Settings inherit in order: project, input root, each intermediate directory,
then sample. Each optional override file is named `batch.yaml`. Missing files
leave inherited settings unchanged. Directory names have no algorithmic meaning.

Mappings merge recursively; scalars and lists replace inherited values. An empty
mapping does not clear inherited fields. `inputs.reference` is the exception:
a lower-level reference replaces the entire block, including any filtering.
A replacement must supply its own path and any desired filter pair. This prevents
filters for one reference from leaking into another.

Every configured filesystem path resolves relative to the YAML declaring it,
before inheritance. Absolute paths are allowed. Adapter parameters are passed to
the adapter unchanged; they do not define additional framework filesystem fields.
The [group and sample YAML examples](input-output-example.md#configuration-files) demonstrate shared CRC reference paths and a sample-specific replacement. A replacement does not retain any upper-level reference filter, and no per-sample copy is created.

| Field | Meaning |
| --- | --- |
| `modality` | Required hST, iST or sST; maps to the existing single-run route. |
| `coordinates` | Required `unit: um \| pixel`; optional positive `microns_per_coordinate`. |
| `enabled` | Boolean, defaults to true; lower levels may override it. |
| `inputs.reference` | H5AD path and optional `filter_column` / `filter_value` pair. |
| `inputs.pm_on_cell` | Optional existing single-run mapping input, with a path relative to its YAML. |
| `algorithm`, `preprocessing` | Existing single-run solver and QC settings. |
| `global_anchoring`, `local_refinement` | Existing GA/LR settings; batch adds iST `cell_types`. |
| `output` | iST `ist_mapping` only; locations and filenames are batch-owned. |
| `execution` | Existing execution settings, including seed. |
| `analysis` | Per-aspect adapter specifications; see [Connect an analysis module](framework.md#connect-an-analysis-module). |

Required protocol fields cannot be null. `inputs.pm_on_cell` is optional but
cannot be null when declared; it must specify a path, as in the single-run
interface. Other single-run fields retain their own null rules, such as null QC
thresholds. Missing reference, modality or units
is an error; none is inferred from platform names, directory names or labels.
The full example configuration is under [configs/batch/](../../../configs/batch/). Refer to the existing
application reference for scientific parameter meanings.

## Standard ST and reference files

**已实现，本轮保留。**

A sample is a directory containing `spatial.h5ad` inside `input_root`. Its
identity is the relative directory path, such as `CRC/S01`. Two different
categories may contain S01 without collision. A sample cannot contain another
sample. Discovery does not follow directory symlinks and does not scan SC or the
output root. No sample list, sample ID field or per-sample YAML is required.

Both H5AD files must have nonempty unique observation and gene identifiers and
finite nonnegative expression in `X`. ST requires finite coordinates in
`obsm['spatial']`. Reference labels use configured GA/LR columns. Reference
filtering must select usable observations, and the inputs must share genes.
Choosing X does not prove that it contains counts; supply the representation
appropriate to the reconstruction protocol.

Checks are read-only. The runner does not write input H5AD, transform labels,
select layers, normalize coordinates or create preparation manifests. Prepare
nonstandard data separately. No platform converter or new preprocessing stage
is implemented in this release. Existing single-run QC, reference filtering and
gene alignment still run as configured.

Coordinate metadata describes existing values without scaling them. Pixel units
without calibration do not support physical-distance claims. Micrometre units
already imply `microns_per_coordinate: 1`; another scale is invalid. Label aliases are
not inferred: `Macro` and `Mono_Macro` remain distinct unless already
harmonized before batch execution.

## Output roles and current-result semantics

**已实现，本轮保留输出位置和重建载体。** 分析方面只按[显式方面配置](#analysis-configuration)执行；空的 `analysis` 不创建占位目录，目录非空也不代表完成。

Output directories mirror sample paths relative to `input_root`. hST uses
`sp-SVC`; iST uses `sc-SVC` cluster mode; sST uses `sc-SVC` sr mode.
For iST, each effective `local_refinement.cell_types` entry becomes an
independent task in a type subdirectory. If the setting is omitted, the current
runner default applies; this document does not prescribe a fixed biological
type list. Labels must be safe, unique path components. Missing types fail
explicitly; there is no name-based exclusion.

完整目录样式统一见[paired 输出示例](input-output-example.md#ist-output-tree)和[hST/sST 示例](input-output-example.md#hst-and-sst-output-tree)，此处只定义角色和状态语义。

hST/sST place `SVC.h5ad`, their handoff, controls and analysis at sample level.
Engine provenance remains available through references in the handoff.
iST `output.ist_mapping` selects assembly only:

* `paired` (default): spatial and expression H5AD carriers.
* `mean`: sparse cluster-mean expression on spatial observations in `SVC.h5ad`.
* `random`: seeded within-cluster donor assignment in `SVC.h5ad`, retaining
  donor identities and assembly provenance.

GA/LR are unchanged. Successful mode switches remove identified framework-owned
obsolete alternatives; failed publication restores prior files.

`reconstruction.json` records task identity, mode, actual carrier roles and
hashes, standard input sources, configuration provenance, coordinate meanings,
pairing checks and engine provenance. Require a successful current handoff before
analysis. iST expression observations are reference donors, not spatial
observations. sST generated units are not paired one-to-one with raw spots.
Unavailable pairing is explicit, not silently approximated.

Analysis is organized by aspect. Paired values on the same observations/windows
belong in one table with named baselines and deltas; independent axes belong in
separate files. Large carriers are referenced, not copied. There is no mandatory
raw/reconstruct/comparison directory layer. Specific scientific table schemas
remain the responsibility of the analysis modules.

## Analysis input views

**已实现。** 以下辅助视图按需在分析中使用，不改变前述落盘目录，不生成第二份标准输入。

### Carriers and alignment

辅助模块只读取并提供分析视图，不写回 H5AD，也不执行科学预处理。实现位于 [inputs.py](../../../revise/batch/inputs.py) 的 `AnalysisInputs`；对已验证的成功重建记录调用 `AnalysisInputs.from_reconstruction(record)`，分析上下文通过 `inputs` 属性延迟执行同一加载。不要求每个模块提前加载全部表达矩阵。

| 对象 | 含义与对应依据 | 当前代码支点 |
|---|---|---|
| raw | 标准 ST 原始文件的 `X`、基因和空间观测；不是 QC／基因对齐后的临时载体 | [sample.py](../../../revise/batch/sample.py) 的 `read_sample`、`BatchSample` |
| hST 重建 | `SVC.h5ad` 的实际输出观测；按 ID 对齐 raw 并核对坐标，输出顺序为分析顺序 | [runner.py](../../../revise/batch/runner.py) 的 `_handoff` |
| iST 空间载体 | 输出空间细胞及 `SVC_cluster`；paired 的 `spatial.h5ad` 不是完整 reference 表达载体 | [sc_svc_application.py](../../../revise/backend/runners/sc_svc_application.py) 的 `ScSVC.local_refinement` |
| iST donor | `expr.h5ad` 的 reference 侧观测，带 cluster 标签；不可视作空间观测 | 同上；[ist_assembly.py](../../../revise/application/ist_assembly.py) 的 `assemble_ist` |
| sST 父 spot | 生成 cell 输出中的 `spot_name`，关联 raw spot；生成 cell ID 与坐标来自实际输出 | [sc_svc_super_resolution_application.py](../../../revise/backend/runners/sc_svc_super_resolution_application.py)、[meta.py](../../../revise/backend/ops/meta.py) 的 `get_sc_obs` |

`raw()` 读取完整原始 ST；hST／iST 的 `aligned_raw()` 选择实际输出对应的空间观测，保留原始基因范围。检查 ID 唯一、映射完整、顺序和坐标一致；记录未进入输出的 raw 观测。无法建立对应关系时明确报不可用或输入错误，不用行号碰巧相同替代对应证据。

视图至少暴露表达矩阵及轴、空间 ID／坐标、基因范围、映射种类和来源。native ID 配对、cluster 映射和构造的父 spot 基线分别标识，不合并为一个含义不明的 `paired=true`。

### Paired view

读取原生双载体，验证空间与 donor 的 `SVC_cluster` 非空且集合一致。保留 donor ID，不添加虚构 donor 坐标。原始 raw 从标准 ST 读取：现有 paired 空间载体可能只保留对齐基因，不能据此缩小实测 panel。

提供两项最小辅助：

1. 对 donor 的原始 `X` 按 cluster 求稀疏均值，返回 cluster × gene 矩阵和轴。
2. 提供空间行到 cluster 行的映射，并允许把 cluster 结果映射回空间观测顺序。

均值定义与现有 `assemble_ist(..., mapping="mean", ...)` 对同一 donor 载体的行为一致；两者共用 `ist_assembly` 的内部均值计算，避免算法漂移。无需强制物化完整空间 cell × gene H5AD，也不隐式随机抽 donor。空间邻接使用实际空间细胞位置，不能用 cluster 中心代替。

模块可选择 donor 侧计算后汇总，或均值表达后计算；例如 EMT 两种顺序可能产生不同结果，框架不代选、不声称等价。模块记录采用的视图与处理顺序。若方法需要显式空间表达矩阵，由模块按需请求或构造，并承担内存需求。

`mean`／`random` 直接读取 `SVC.h5ad`，保留现有随机 seed 与 `revise_ist_donor_id`。不从 paired 自动切换模式。以后用不同输出根保存两种重建，再由明确的比较调用读取，首版不新增实验调度器。

具体入口为 `inputs.paired()`：`.donor` 与 `.spatial` 保留各自轴，`cluster_means()` 返回矩阵及 cluster 顺序，`map_to_spatial(values)` 将同一 cluster 顺序的一维分数或二维矩阵映射到空间行。`inputs.reconstructed()` 在 paired 下明确拒绝调用，避免把空间载体的 overlap panel 当作重建表达；其他模式直接返回 SVC 视图。

### SST baseline

sST 的原生 raw spot 与生成 cell 不是一对一。为在相同生成单位上比较，构造如下 raw 基线：

```text
n[s] = 实际发布结果中父 spot 为 s 的 cell 数
raw_baseline[c, g] = raw.X[parent(c), g] / n[parent(c)]
```

基线沿用生成 cell 的 ID、顺序和坐标；表达使用 raw 原始基因。以稀疏映射乘法和浮点除法计算，不做整数截断，不默认稠密化全矩阵。逐父 spot、逐基因检查其 cell 基线之和等于该 raw spot 表达（浮点容差内）。

坐标优先读取生成输出的 `obsm['spatial']`。当前 sc-SVC sr 引擎也可能把每个生成行的显式坐标写在输出 `obs[['x', 'y']]`；当 `obsm['spatial']` 缺失时，辅助视图只接受这组输出自身的字段，并在交接 provenance 中记录 `coordinate_source`。它不根据 raw 的 `spot_name` 父映射回填或猜测生成坐标。`get_sc_obs` 在没有 `revise_cell_locations` 时可能为生成行保留 parent-spot center；因此该来源表示 generated spatial-unit coordinates，可能是 parent-spot center 或显式 cell center，不构成真实 cell 几何证据。两种字段都缺失时，sST 基线明确不可用。

- 父 spot 必须来自实际输出映射；不采用预估 cell 数或最近邻推断。
- 生成 cell 缺少父映射、指向不存在的 raw spot、重复 cell ID，均拒绝构造基线并给出原因。
- 没有生成 cell 的 raw spot 记录为未覆盖；不除以零、不补造 cell。守恒声明仅覆盖已表示的父 spot。
- 明确标注这是 `parent_spot_equal_split` 构造基线，不是实测单细胞表达，也不是原生 cell 配对证据。
- 不修改重建算法、输出表达或 `rec_match_spot_sum`；基线守恒不要求重建总量与 raw 相等。

`inputs.sst_baseline()` 提供构造基线，`_handoff` 对 sST 的原生配对不可用判断保持不变。分析上下文提供这项构造基线能力，让支持这种比较的模块执行；需要真实 cell 配对的模块仍然不可用。

### Processing boundaries

辅助模块不负责归一化、log、EMT、Moran、CCI、对象匹配或显著性判断。raw 与重建可能使用不同表达尺度；模块必须声明实际处理，不能由辅助模块偷偷归一化成相同尺度。

iST raw 使用实测 panel，重建使用自己的基因范围。单基因差值只对实际共有且可比的基因计算；重建独有基因可报告新增可评估信息，不报告由补零制造的变化。其他分数的可比性、常量基因排除和最小样本条件由模块定义。

验证定义集中在[输入验收](plan-validation.md#input-scenarios)，计算视图与来源的落盘要求见[结果交接](protocol.md#analysis-results)。

## Analysis configuration

**已实现。** 分析方面采用显式、可追踪的适配器规格；当前可运行入口见[接入分析模块](framework.md#connect-an-analysis-module)。

配置示例：

```yaml
analysis:
  spatial_diversity:
    entrypoint: my_analysis.diversity:run
    version: "1"
    parameters:
      radius: 50
    resources:
      gene_sets: resources/gene_sets.csv
  cci:
    enabled: false
```

支持的方面字段为安全的方面名、`entrypoint`、非空 `version`、`parameters` 映射、`resources` 映射、`enabled` 和 `requires_pairing`。资源相对声明该字段的 YAML 所在目录解析后再继承，与现有路径原则一致；`resources` 是名称到文件路径的简单映射，内容摘要自动计算。资源缺失使该方面失败。

方面名为安全的单目录名。只执行明确配置且启用的方面；缺省或空 `analysis` 不自动运行或创建占位任务。递归继承中的空映射不会删除上层设置，禁用继承方面用 `enabled: false`。不提供方面间依赖调度。

## Analysis results

**已实现。** 每个启用的方面必须返回以下最小结构；只返回角色到文件路径的旧 adapter 需要迁移，框架不会把旧返回值解释成已提供计算依据：

```python
{
    "artifacts": {
        "comparison": {"path": "comparison.csv", "description": "实际检验结果"},
        "localization": {"path": "locations.csv", "description": "差异所在空间单位"},
    },
    "calculation": {
        "input_view": "cluster_mean",
        "parameters": {},  # 包含模块实际使用的默认参数、处理顺序及 seed（如适用）
        "comparison_basis": "模块说明比较单位、匹配依据、基因范围及处理尺度",
    },
}
```

路径均相对模块 staging 输出目录，继续使用越界／符号链接检查。示例角色不强制每个方面输出完全相同文件：模块声明其实际产物用途和必要比较依据，框架验证文件存在且可归属并记录摘要。来源引用由框架补充，包括输入、重建、资源和实现身份，避免模块重复抄写。`calculation.input_view`、`parameters` 和 `comparison_basis` 均为必填字段。

优先扩展现有 `analysis/analysis.json` 及方面控制记录，不新增平行 manifest。保留实际计算结果（包括不显著结果）与解释“哪里、多少”的必要空间证据；通过 ID 引用原有坐标或矩阵即可，不要求每次复制。按模块需要区分未测、筛除与未检验，不要求完整理论候选全集或所有中间图、表达矩阵。

分析状态记录在 `.revise/analysis/<aspect>.json`，批量分析汇总在
`analysis/analysis.json` 和输出根的 `analysis_status.json`。成功、复用、失败、
不可用和未实现均由记录区分；只有已发布的结构化产物才可标为成功。适配器在
当前解释器内导入执行；检测到源码变化时会拒绝继续，需重启解释器后再试。

科学筛选、匹配规则、统计前提、指标含义由模块声明。框架检查结构、来源与发布一致性，不判断科学结论有效性；“模块被调用”也不等于分析方法已验证。
