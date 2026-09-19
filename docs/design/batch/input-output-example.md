# Standard input and output example

这是输入输出目录样式的统一示例，供实现、数据整理和外部项目适配时直接对照。[返回入口](README.md#contract-and-status)。字段与继承规则以[协议](protocol.md#configuration-and-inheritance)为准，代码职责见[框架](framework.md#target-architecture)。

## Input tree

```text
project/
├── batch.yaml                         # 项目公共配置
├── ST/                                # input_root，仅在这里发现样本
│   └── CRC/
│       ├── batch.yaml                 # CRC 共用 reference 与重建设置
│       ├── S01/
│       │   └── spatial.h5ad
│       └── S02/
│           ├── spatial.h5ad
│           └── batch.yaml             # 可选：样本特有覆盖
├── SC/                                # reference 按癌种组织，按路径共享
│   ├── CRC/
│   │   ├── crc_reference.h5ad          # 所有适用的 CRC 样本共用
│   │   └── special_reference.h5ad      # 仅特殊样本使用
│   └── BRCA/
│       └── brca_reference.h5ad         # 供未来 BRCA ST 样本引用
└── output/                            # output_root，独立于输入
```

SC 通常按癌种分目录，也可以继续细分研究或保存多个 reference。SC 目录本身不参与 ST 样本发现或配置继承；使用哪个文件由 ST 沿途的 `batch.yaml` 显式指定，不因癌种同名自动匹配。允许引用外部路径，也允许不同分类共享同一 reference。

新增普通样本通常只需放入 `ST/CRC/S03/spatial.h5ad`，即可继承 CRC 配置。上图中的 H5AD 都是用户提供的真实数据，不随模板生成；BRCA 用于示范 SC 分类，不是必须建立的目录。

## Configuration files

可复制的 YAML 直接维护在以下文件中，文档不再重复配置代码：

| 文件 | 用途与解析结果 |
|---|---|
| [项目 batch.yaml](../../../configs/batch/batch.yaml) | 设置 `input_root: ST`、`output_root: output` 与公共重建参数 |
| [CRC batch.yaml](../../../configs/batch/ST/CRC/batch.yaml) | `../../SC/CRC/crc_reference.h5ad` → 项目下的 CRC 共用 reference |
| [S02 batch.yaml](../../../configs/batch/ST/CRC/S02/batch.yaml) | `../../../SC/CRC/special_reference.h5ad` → CRC 特殊 reference；完整替换上层 reference 与筛选 |
| [可选分析配置](../../../configs/batch/analysis-example.yaml) | 显式启用一个载体清单方面；由[示例适配器](../../../configs/batch/analysis_example.py)产生可重载 CSV |

不需要特殊 reference 时不放置 S02 覆盖文件。QC、标签列、单位等模板值需要与真实数据核对；模板不是任意数据均可直接运行的科学参数预设。分析方面只在配置显式启用；字段和结果记录见[分析配置](protocol.md#analysis-configuration)与[分析结果](protocol.md#analysis-results)。

示例中的 `cell_types` 与 `paired` 保留用于显式旧单类型兼容路径；whole-sample
交付应省略 `local_refinement.cell_types`，并使用 `random`（或明确选择
`mean`）。

## IST output tree

下图为显式配置 `local_refinement.cell_types` 或单类型兼容调用时的 iST
`paired` 布局示例。状态文件由对应调用产生，不要求用户预建；它不是
whole-sample 交付的目录。

```text
output/
├── batch_status.json                  # 批量重建汇总
├── analysis_status.json               # 批量分析汇总
└── CRC/S01/
    ├── .revise/sample.json
    └── T/                             # 示意：每个有效 cell_type 各一份
        ├── spatial.h5ad               # 空间载体
        ├── expr.h5ad                  # donor 表达载体，不是空间观测
        ├── reconstruction.json        # 重建交接入口
        ├── .revise/
        │   ├── application.yaml
        │   ├── reconstruction.log
        │   ├── task.json
        │   └── analysis/
        │       ├── <aspect>.json       # 方面状态
        │       └── <aspect>.log
        └── analysis/
            ├── analysis.json          # 分析汇总
            └── <aspect>/              # 实际分析产物，不固定文件名
```

`<aspect>` 是示意占位符，不是要建立的字面目录名。只有显式配置且启用的方面
才创建自己的目录；空配置不创建占位结果，目录存在也不表示分析已完成。适配器
返回的状态、产物和计算依据见[分析结果协议](protocol.md#analysis-results)。

`mean`／`random` 使用同一目录，将双载体替换为 `SVC.h5ad`；其余记录与分析位置不变。模式切换的受控清理和回滚规则见[输出协议](protocol.md#output-roles-and-current-result-semantics)。

省略 `local_refinement.cell_types` 时使用 sample-level whole-sample 路径：GA
只运行一次，实际 broad types 中符合 Level2 资格的类型合并到同一个 SVC，
默认使用 `random`。其目录还包含完整 Raw、sample 配置和同次运行记录：

```text
output/
└── CRC/S01/
    ├── SVC.h5ad                 # whole-sample random/mean assembly
    ├── raw.h5ad                 # 原始 X/轴/坐标/注释，含 QC 排除单位
    ├── sample.yaml              # 相对路径、标签、坐标和表达语义
    ├── reconstruction.json      # 同次运行 handoff、fingerprint 和 provenance
    ├── .revise/
    │   ├── sample.json
    │   ├── application.yaml
    │   ├── reconstruction.log
    │   └── task.json
    └── analysis/
        ├── analysis.json
        └── <aspect>/
```

Raw/SVC 的原始 annotation 保留；本次真实推断写入 `revise_Level1`、
`revise_Level2`，没有推断的位置保持缺失。`sample.yaml` 不把未知表达尺度
写成 raw counts 或 `log1p`。

## HST and SST output tree

hST／sST 无细胞类型目录，下面是样本级布局；批次汇总仍在输出根。

```text
output/CRC/S01/
├── SVC.h5ad
├── reconstruction.json
├── .revise/
│   ├── sample.json
│   ├── application.yaml
│   ├── reconstruction.log
│   ├── task.json
│   └── analysis/                      # 方面状态与日志
└── analysis/
    ├── analysis.json
    └── <aspect>/                      # 模块实际产物
```

这些是框架交接位置；引擎 provenance 等额外文件通过重建记录引用。不要用文件存在替代状态校验。sST 的生成 cell 与 raw spot 不是原生一一配对；支持按实际 `spot_name` 构造的[守恒基线](protocol.md#sst-baseline)，需要真实 cell 配对的方面仍应报告不可用。

## External data adaptation

1. 将已经满足协议的 ST 放到样本目录，固定命名为 `spatial.h5ad`；非标准数据先在外部处理。
2. SC 按癌种存放一次，在项目／分类配置填写路径；只为特殊样本增加覆盖。
3. 从项目配置启动调用，输出由框架生成，不将结果写回 ST 或 SC。
4. 通过 `reconstruction.json` 和分析记录消费真实产物。科学表格结构由模块声明，不要求创建统一的 raw／reconstruct／comparison 目录层。

本示例只展示文件组织；算法是否实现、结果是否有效，分别由模块能力和运行记录确定。
