# 输出、记录、报告与验收契约

完整 27 节点科学追溯树与状态见 [分析框架](analysis-framework.md#content-tree)。本页继续负责
字段、结论填写和展示契约，不复制科学树；节点内容通过本地静态映射与既有 question 关联。
本轮组织方案把主文 overview 固定为五行：分群、基因信息 + Moran、EMT、双 baseline、State；
基础定义、配对、完整基因和预处理通过固定目录的按需 audit 访问，不要求在主文中占行。
本轮实现与检查状态统一见[规范入口](README.md#report-display-status)。

本页负责“怎样保存和阅读已经计算出的事实”：结果身份、长表/JSON、状态、
结论 record、`node_narratives`、Notebook/HTML authoring contract、模板 provenance 和
R0–R4 验收。
科学定义分别见 [input-views.md](input-views.md)、[impact-analysis.md](impact-analysis.md)
和 [gene-and-function.md](gene-and-function.md)；未来批量边界见
[batch-integration.md](batch-integration.md)。本页不新增科学指标。

## 1. 结果身份和科学产物根

现有正式科学结果根只有：

```text
output/reconstruction_impact/P1CRC_VisiumHD/
output/reconstruction_impact/P2CRC_Xenium_ThreeParents/
```

schema 1 配置的 `output.dir` 是路径权威（允许已登记的环境覆盖）。每个 root 包含
`analysis/`、`figures/`、`notebook/` 和 `report.html`；旧 `Reconstruction_Impact_Report.ipynb`、
历史下载目录和开发 adapter 不是活动入口。正式运行的 source/executed identity、
stderr、checkpoint 状态、参数/实现/helper hash 由：

- `notebook/execution.json`：source/executed notebook digest、code-cell count、
  stderr、`checkpoint_used`；
- `notebook/validation.json`：source/executed cell identity、错误/警告检查和结果
  contract；
- `analysis/audit.json`、`parameters.json`：route、seed、cohort、expression、
  graph/resource/Region conditions 和输入文件 digest；
- `analysis/artifacts.csv`：analysis/figure/notebook/report 文件的相对路径、大小和
  hash。

阅读产物刷新不要求重新执行科学计算；在 source、参数、输入和科学表未改变时，应保留
原 `execution.json` 身份并只更新阅读记录/渲染证据。若其中任一科学依赖改变，必须另建
执行身份。证明当前阅读结果仍须读取这些记录；先前的执行次数或历史表不能代替对应的
执行证据。
VisiumHD 的 `All` scope 有独立 cohort/denominator；Xenium 没有 global scope。`All` 只消费
已有六类 global record：`foundation`、`complexity`、`matched_k`、`moran_all_valid`、
`moran_shared_valid`、`changed_units`；不为 EMT、局部邻域、State 或 Gain 补造 All record。
All 开头显示独立抽样、分群和 Moran 的各自分母。`changed_units` 专表按 Raw Level1 类型
口径保留 Wilson 区间；All 的 `provided_n` 归档缺失时保持 `NA` 并显式说明，另列已保存完整
载体计数 `18085 → 12926`，不把它填入 `provided_n`。

本节的目录是科学产物和证据的保存布局，不是 HTML 左侧阅读树。HTML 的固定左侧层级目录和
单页长页顺序见本页第 5 节；它按“总体 → 局部邻域 → 空间”组织，不能从结果目录、当前发现
或文件是否存在推导。

## 2. 记录粒度与状态

正式结论以 **科学问题 × `task_cell_type` × 必要的 `baseline`** 为粒度，
`scope` 是独立上下文。foundation 没有 baseline 时显式写 `baseline=null`；
Raw Leiden 与 Raw Level2 分别写 `raw_leiden` 和 `raw_level2`。一个记录至少包含：

```text
identity: record_id, sample_id, route_kind, scope, task_cell_type, question,
          baseline, layer, rule_id, spec_version
spec_version, rule_id, results
performance: {judgment, headline_eligible}
performance_judgment
evidence_condition: status, reason, interpretation_condition
conclusion: {zh, en}
limitations: {zh, en}
evidence: tables, figures, condition
run_identity/review: source and result identity, review status and digest
```

结果表和 record 不能把 status 压成一个 success flag：

| 层次 | 必须保留的状态 |
| --- | --- |
| 执行 | cell success/failure、source/executed mismatch、输入/资源错误分别记录 |
| 单侧科学结果 | `computed`、`unmeasured`、`low_coverage`、`insufficient_support`、`not_computable`，并写 reason |
| 比较 | `raw_status`、`reconstruction_status`；只有两侧有效才写 delta，否则保留单侧值和 `unavailable` |
| 特定规则 | `unmatched_cluster_complexity`、`no_stable_threshold`、`headline_eligible`；保存 threshold reason、valid/total bootstrap、CI width 和分母 |

预期数据不足可以作为带状态的科学结果发布；缺依赖、损坏资源、invalid axes 和
代码异常属于 `blocked`/失败。CSV 的空值重载为 `NA`，JSON 使用 `null`，不写非法
NaN。数值/规则判读可以自动产生；综合解释需要审阅，审阅状态和 record digest
绑定，内容变化后自动回到 `pending`。结论整理器沿用已批准的规则，从保存表整理阈值区间相对宽度、bootstrap 比例等事实并保存判读；不重新拟合阈值或运行科学分析。HTML renderer 只消费这份已保存记录，不重新判定 Region。

### 节点映射、同源叙述与兼容字段

`content_contract.py`仅保存节点顺序、双语短目的、方法链接及节点对应的既有question。
16类question身份不变；一个record可以服务多个节点，不能把内容树机械扩展为新的科学记录。
报告的阅读节点与Notebook的`impact_node`展示占位使用同一映射。旧`module-<question>`
锚点仍作为兼容入口保留。

本轮记录扩展保持旧字段类型和值：`foundation.results.cohort`仍为原字符串，
结构化队列审计放入`cohort_audit`；`selection`按本scope的Raw Level1对象登记覆盖、
纳入和排除，不能把全Raw其他cell type当作本parent的排除数。`preprocessing`读取
实际partition、Moran graph和EMT cutoff审计。配对核查范围保存在`pairing_audit`的
来源说明中，`n_checked`指抽样前完整空间载体，不是各分析的队列分母。
局部两baseline记录的`anatomy_summary`只转录既有分层CSV的对应中位数和窗口分母，
不重算统计；Kobs和evenness不得直接套用该record中基于Neff形成的方向判断。
overview 的“基因信息 + Moran”行读取 foundation/gene availability 中已有的 Raw/Recon
`provided_n`，并与 Moran 的 full-side/shared-valid 两个问题并列；`provided_n` 不由
Moran valid n 推算，也不新增一套基因分析记录。

记录包新增 `node_narratives`，作为 Notebook 与 HTML 共享的同源叙述列表。每项至少含
`node_id`、双语 `title`、`lead`、`interpretations`、`next`，以及 `record_ids` 和
`review`。不得把某一展示层临时写出的句子作为另一层的独立结论：

其中 `next` 继续保存在结果包并供 Notebook 使用；网页不沿用旧的 `next` 顺序组织 module 或
导航，网页顺序由固定左侧目录和页面阅读树决定。

```json
{
  "node_id": "2.1",
  "title": {"zh": "分群复杂度", "en": "Partition complexity"},
  "lead": {"zh": "本批先描述 K 的变化，再检查匹配条件。", "en": "This batch first describes K, then checks matching."},
  "interpretations": [
    {"zh": "Fibroblast 在同 resolution 下为 5 → 5。", "en": "Fibroblast is 5 → 5 at the same resolution."},
    {"zh": "该事实描述分群结构，不单独证明改善。", "en": "This describes partition structure and does not by itself prove improvement."}
  ],
  "next": {"zh": "转到 matched-K 条件。", "en": "Continue to matched-K conditions."},
  "record_ids": ["complexity:Fibroblast"],
  "review": {"status": "pending", "record_digest": "<digest>"}
}
```

`interpretations` 保持 1–3 条，按“事实 → 有条件的判断 → 限制”填写；没有对应事实时
保留空列表并在 `review` 标记 `pending`，不生成空 shell。旧 `section_summaries` 可继续
作为兼容投影，旧 16 类 `question` 身份、记录粒度和数值字段不变；新的
`node_narratives` 只增加节点级追溯与双语阅读内容。科学表变化或解释变化时，按既有
digest 机制重新审阅；仅刷新阅读产物不要求重新运行原始分析。

正式批次观察由同一批 records 填入对应节点；Notebook 即时计算表仍可展示，
但不得另拼一套正式科学判断。基础 audit 不提前执行 Moran、EMT、窗口或 Region
来填满主文 overview。

## 3. 科学产物布局（不等于 HTML 阅读树）

具体文件名以当前 writer 和 artifact index 为准；科学矩阵的数值身份不因 HTML
重绘改变。正式根的最小可重载布局为：

```text
<sample-root>/
├── analysis/
│   ├── audit.json, parameters.json, comparisons.csv, artifacts.csv
│   ├── conclusions/report_records.json, report_records.csv, review.json
│   ├── cross_parent_summary.csv, cross_parent_local_diversity.csv
│   ├── cross_parent_state_region_extent.csv
│   ├── inputs/{carriers.csv, gene_availability.csv, source_files.csv}
│   ├── integration/{source_inventory.csv, integration.json}
│   ├── reconstruction_impact/partition/
│   │   ├── summary.csv, assignments.csv.gz, mapping.csv, contingency.csv
│   │   ├── resolution_sweep.csv, complexity_sweep.csv, level1_summary.csv
│   │   └── <scope>_audit.json
│   ├── moran/{moran_by_cell_type.csv, moran_summary.csv,
│   │         moran_distribution_summary.csv, gene_availability.csv,
│   │         moran_graph_audit.csv, graphs/<scope>/...}
│   ├── pathway_activity/{scores.csv.gz, summary.csv, availability.csv,
│   │                    resource_genes.csv, resource.json, cutoff_audit.csv}
│   ├── anatomy/{anatomy_windows.csv, anatomy_context_summary.csv,
│   │           support_sensitivity.csv}
│   ├── changed_units/<parent>/matched_k_assignments.csv.gz
│   ├── spatial_fields/values.csv.gz
│   └── local_state/<parent>/
│       ├── raw_level2_labels.csv, raw_level2_posterior.csv.gz
│       ├── window_decision.json, unit_window_assignments.csv.gz, window_metrics.csv
│       ├── diversity_by_anatomy.csv, change_by_anatomy.csv
│       ├── state_threshold_bootstrap.csv, state_threshold.json
│       ├── gain_threshold_bootstrap_audit.csv, gain_threshold_audit.json
│       ├── state_region_extent_by_anatomy.csv, gain_region_extent_audit.csv
│       └── scale_sensitivity.csv
├── notebook/{<formal executed route>.ipynb, execution.json, validation.json}
├── figures/...
└── report.html
```

`analysis/<scope>/inputs/` 由实际 `save_table(..., scope=scope)` 生成时还登记
`observations.csv.gz`、Xenium `mapping.csv.gz`/`reference_clusters.csv.gz`、
`partition_observations.csv.gz` 和 `partition_features.csv`。大型表达矩阵引用原
H5AD/digest；不写 unit × full-gene 巨型 CSV。所有表显式写出 ID 列、单位、状态和
联合主键；排序不能代替 ID 对齐。

核心长表字段由方法文档定义，至少保留：

| 产物 | 行粒度 | 必须可重载的事实 |
| --- | --- | --- |
| partition | scope × edge × unit/cluster/resolution | 两组 K 及各自已有 ARI、unit_change、balanced（`balanced = 1 - macroF1`）、size、target/actual K、absolute/normalized contingency、mapping、assignment、changed、Wilson denominator/status；缺字段保留 NA/状态 |
| Moran | comparison × scope × gene | Raw/Recon `provided_n`（引用 gene availability/foundation）与 full-side/shared-valid 两视角、两侧值/delta/status、n_units/n_edges、graph ID、gene rule、expression view |
| EMT | scope × unit × pathway 或 scope × side × pathway | score/delta/status、resource_total、available/detected/coverage、provider、seed、rank cutoff、resource digest |
| anatomy/window | full-context unit/window 或 scope × scale × window | 坐标、Raw exact label、anatomy、unit/window counts、support curve、scale、denominator |
| local diversity | scope × window × baseline × metric × scale | Kobs、entropy、Neff、evenness、Raw/Recon/delta、draw count/seed、baseline/anatomy/support |
| Region | scope × scale × anatomy × region type | State/Gain、threshold status/reason、bootstrap、`valid_windows`/`region_windows`、`valid_units`/`region_units`、area、`area_fraction`、`unit_fraction`、分母 |

当前 writer 的关键长表名称和列语义保持如下；报告可以增加 named links，但不能以
展示需要改名或交换这些科学列：

| 表 | 必要列/语义 |
| --- | --- |
| `partition/contingency.csv` | absolute count、normalized value、对应 row/column denominator、scope/edge |
| `partition/resolution_sweep.csv` | target K、actual K、K difference、selection flag、status |
| `partition/complexity_sweep.csv` | Raw candidate resolution、Level1 ARI、K、selection reason |
| `local_state/<parent>/window_metrics.csv` | scope、scale、window ID、`valid_window`/support、Kobs、entropy、Neff、evenness、Raw/Recon/delta、baseline、draw count/seed、anatomy |
| `local_state/<parent>/support_sensitivity.csv` | `window_side_length`、`n_tissue_windows`、`n_valid_windows`、`valid_window_fraction`、`retained_parent_units`、`retained_parent_unit_fraction` |
| `local_state/<parent>/state_region_extent_by_anatomy.csv` | `level1_region`、`valid_windows`、`region_windows`、`region_area_um2`/`region_area_mm2`、`area_fraction`、`valid_units`、`region_units`、`unit_fraction`、`region_available`、`threshold_status`、scale |
| `local_state/<parent>/gain_region_extent_audit.csv` | 保持现有 valid/region window 与 unit counts、area/fraction、`region_available`；阈值状态/原因从同目录 `gain_threshold_audit.json` 读取，parent 与尺度关联该目录的窗口审计，不伪称 extent CSV 已内嵌这些列 |

公共比较身份保存在 `comparisons.csv`，固定至少包含
`comparison_id`、`sample_id`、`task_cell_type`、`scope`、`comparison_edge`、`raw_view`、
`reconstruction_view`、`label_source`、`observation_basis`、`gene_rule` 和
`normalization`。Global `All` 与 parent-internal comparison 使用不同的
`scope`/`comparison_edge`，不能靠展示 label 猜测。

`analysis/moran/moran_distribution_summary.csv` 的固定列为
`scope`、`comparison_id`、`gene_set`（`all_valid`/`shared_valid`）、`side`、
`n_valid`、`median`、`q1`、`q75`；`all_valid` 是各侧 full-side 分布，
`shared_valid` 是共同有效 gene 的两侧摘要。Moran graph 目录同时保存 raw weights、
row-normalized weights 和有顺序的 `observations.csv.gz`；graph audit 登记实际邻居
参数、坐标、components、isolates、权重变换和版本，后续 compare 只消费这套 shared
weights。overview 需要的 `provided_n` 仍从同一结果包的 gene availability/foundation
record 引用；它不改写 `moran_distribution_summary.csv` 的既有列，也不等同于 Moran 的
`n_valid`。

EMT 的 `availability.csv`/`resource.json` 必须保存当前批准 resource target
`resource_total=200`（若实际 digest 总数不同，写实际值和原因）、resource path/version/
digest、species/ID type、provider、seed、实际 per-side rank cutoff、score order、
available/detected hits 和 coverage。`cutoff_audit.csv` 保留动态 cutoff，不能把
`AUC_threshold=0.01` 写成固定前 1% genes。

partition writer 的 sweep 语义固定为：`resolution_sweep.csv` 是 matched-K target/
actual K、差距、选择标记和状态；`complexity_sweep.csv` 是 Raw candidate resolution、
Level1 ARI、K 和选择理由。历史文件名 `matched_k_resolution_sweep.csv` 与
`complexity_raw_resolution_sweep.csv` 不能交换含义。

State overview 只显示 `area_fraction` 一个 headline metric；area、
`region_windows/valid_windows`、`region_units/valid_units` 和分母进入正文条件列及
named data links。`no_stable_threshold` 时 Region extent/fraction 数值为 `NA`/`null`，
而 `valid_windows`、`valid_units` 等支持计数仍保留，不能强制填 0。

## 4. Notebook authoring contract（不要求同 HTML 排版）

Notebook 的核心参数、staged short calls、audit 输出、plot 调用和一段 brief batch
observation 必须可在 source 中审阅。按块的实际职责保留步骤：

```text
purpose + method link → parameters → short calls → necessary audit/table → relevant plot → brief observation
```

基础审计无需配图；`relevant plot`仅用于有科学解释价值的既有证据图，不要求每个节点新增图。
基础节点的审计结果按需引用，五行 overview 与正文优先用短段引用已有结果，不再生成
一套重复数值表。

Notebook 保留自己的执行/authoring 顺序，不复制 HTML 的固定左侧目录或单页长页；HTML
阅读树和 Notebook 排版相互独立。Notebook 仍须保留核心参数、短调用、必要检查、相关
图和同源批次观察，不能因为展示改造而删除已有分析。

长实现和绘图 helper 属于 `reproduce/case/reconstruction_impact/` 的本地模块；
Notebook 不复制长函数、不把 helper 默认值当作当前参数，也不把 HTML 结论写成新的
科学计算。Markdown/comments 使用英文，结论事实从保存 records 填充。

canonical source notebook 和 builder 为：

- [VisiumHD_sp_SVC_Reconstruction_Impact.ipynb](../../../reproduce/case/reconstruction_impact/VisiumHD_sp_SVC_Reconstruction_Impact.ipynb)
- [Xenium_sc_SVC_Fibroblast_Reconstruction_Impact.ipynb](../../../reproduce/case/reconstruction_impact/Xenium_sc_SVC_Fibroblast_Reconstruction_Impact.ipynb)
- [build_notebooks.py](../../../reproduce/case/reconstruction_impact/build_notebooks.py)

稳定英文 block 的唯一约定文件是 `templates/notebook-block.md`；它由 Notebook
formatter/generator 消费，source-generated notebook 的 intro/块输出 hash 进入
执行来源身份。不要为 foundation/overall/localization/region 各建一套模板。

## 5. HTML report contract

现有 renderer 是 [report_html.py](../../../reproduce/case/reconstruction_impact/report_html.py)，
事实和 evidence record 由 [report_records.py](../../../reproduce/case/reconstruction_impact/report_records.py)
提供。HTML 只消费保存的 CSV/JSON/图和 review status；不读取 H5AD、不重算指标、不
猜测图名、不补齐 NA。`report.html` 只属于同一 sample root。

Renderer-owned module anchors 使用 `module-<question>`，例如
`module-moran_all_valid`；HD All module 追加 `-all`；旧 scope-row 锚点保留或映射到对应新问题。record evidence 的字段名是 `evidence.tables`，每项
保存 `path`（以及存在性），由 renderer 作为具名 CSV link 消费；不要引入平行的
`evidence_tables` wire field。

### HTML 私有阅读树（本轮组织方案）

HTML 使用固定左侧层级目录和一张可连续滚动的长页。该阅读树独立于 27 节点科学注册表、
`content_contract.py` 的 question 映射和科学产物目录；它不因发现、record 数量、图是否
存在或状态动态重排。基础定义、范围/配对、完整基因和预处理放在固定目录的按需 audit，
缺少事实时在原位置显示 `pending`/`NA` 与原因。Notebook 保留执行/authoring 顺序，不要求
与 HTML 同排版。

```text
HTML private reading tree (fixed left navigation + one long page)
├─ 00 总览
├─ 01 总体
│  ├─ 分群
│  ├─ 基因信息 + Moran
│  └─ EMT
├─ 02 局部邻域
│  └─ 双 baseline
├─ 03 空间
│  ├─ anatomy
│  ├─ 分群图 / assignment
│  ├─ EMT spatial
│  ├─ local 分层
│  ├─ State 可靠性 → 连续场 → mask / coverage
│  └─ Gain audit（独立入口）
├─ 04 按需审计
└─ 05 HD All（仅已有六类 global record；仅 VisiumHD）
```

页面长页顺序固定为“总体 → 局部邻域 → 空间”；五行总览使用固定问题标题和固定 parent
顺序。HD `All` 只使用已有 `foundation`、`complexity`、`matched_k`、`moran_all_valid`、
`moran_shared_valid`、`changed_units` 六类 record，不为其他问题创造 global record。

### Overview

两份已生成 HTML 的 Overview 固定为以下五行：

| 行 | 固定问题标题 | 本轮组织方案中的可见 headline | 细节入口 |
| --- | --- | --- | --- |
| 1 | 分群 | 两组 K（`complexity`、`matched_k`）与成员改变预览（`unit_change`）；ARI 和 `balanced`（`1 - macroF1`）在正文分群表完整呈现 | 2.1–2.3、按需 pairing/partition audit |
| 2 | 基因信息 + Moran | Raw/Recon 已有 `provided_n`，并列 Moran full-side 与 shared-valid；只有 shared-valid 使用配对 delta | 1.4、2.4、完整 gene/graph audit |
| 3 | EMT | coverage 与 score 的主要摘要，分开显示 | 2.5、3.3、resource/cutoff audit |
| 4 | 双 baseline | Raw Leiden 与 Raw Level2 并列；Kobs、Neff、evenness 共用完整 2×3 图组，并显示窗口支持与尺度 | 3.4–3.8、window/anatomy audit |
| 5 | State | reliability → 连续场 → mask/coverage；Gain 只从独立 audit 入口进入 | 4.1–4.5、threshold/extent audit |

两组分群记录缺少某字段时保留 `NA`/状态，不从另一组推算；`provided_n` 取已有 gene
availability/foundation 事实，不用 Moran valid n 代替。parent 列固定
`Fibroblast`、`Mono/Macro`、`T`；每个 cell 只给一项主要事实和一句 brief conclusion，
并链接正文。HD `All` 作为独立 scope 放在适用行的独立 scope 区域；Xenium 不渲染 global。
Moran 的 full/shared、EMT 的 coverage/score 和局部两条 baseline 在对应五行内部保持
子问题分开，不能合并成总分。
Overview 的分群行只作两组 K 与成员改变的预览；完整 ARI、`unit_change` 和
`balanced = 1 - macroF1` 保留在正文分群表，不要求总览塞入全部诊断字段。

### 每个部分的说明段

每个固定部分先交代逻辑位置、采用理由和预期帮助形成的判断；完整内容在对应方法文档维护，
报告仅保留本模块所需说明，Notebook采用简洁英文目的与方法链接。预期判断不等于预设
本批改善。基础定义、范围/配对、完整基因和预处理只在按需 audit 中说明。常规组织按
既定阅读树固定；缺少事实保留状态，不以发现改变位置。共用要求见[节点填写与确认结构](analysis-framework.md#node-writing)。

### Body module

网页正文 module 的固定顺序是：

1. **固定问题标题（Question）**：问题、节点、scope 和 baseline 清楚可见；
2. **简短结果（Brief result）**：一句只引用本批保存事实的主线句，注明状态或解释条件；
3. **表/图（Table/Figure）**：只放该问题需要的数字、单位、分母、参数/状态和具名
   `evidence.tables`/`evidence.figures` 链接；图紧随表，能改变判断的图才嵌入；
4. **1–3 条 interpretations**：按“事实 → 有条件的判断 → 限制”填写，同一批 records
   同源双语。保存的 `node_narratives.next` 仍写入结果包并供 Notebook 使用，但网页不以
   `next` 排序 module 或导航。

双 baseline 部分必须把 Raw Leiden 与 Raw Level2 并列，保留 Kobs、Neff、evenness 三项
指标的完整共享 2×3 图组：每项均显示两个 baseline 行，以及 Raw、Recon、delta 三列。
同一 Recon 值可以在两条 baseline 中出现，但不能重复累计；窗口支持和尺度随图/表一并
给出。空间部分按 anatomy → 分群图/assignment → EMT → local 分层 → State 可靠性、
连续场、mask/coverage 展开，Gain 通过独立 audit 入口阅读。

Common figure 只在拥有它的 module/row 完整嵌入一次，其他 scope row 给 named link。
Parent-specific spatial images 必须放在对应 scope row；不能把一张 Fibroblast 图
复制到 Mono/Macro 或 T。caption 说明图回答的问题、scope 和比较方向；适用的单位、
坐标、分母及限制由图内标注和同一行的数值证明补齐，不在每张图注重复完整方法。
原图 link 使用 record evidence 中的已登记相对路径。

稳定 HTML shell 的唯一约定文件是 `templates/report.html`；它是 runtime renderer
读取的 shell，evidence tables 仍由 renderer/record 负责生成。模板路径、内容 hash
和 `report_html.py`/`report_records.py` source hash 一并进入 `report_source_digests`
或等价执行 provenance。不要再创建按层拆分的 HTML templates，也不要把 method/record
例子变成单独 template；示例只以内联 contract 形式存在。

<a id="report-examples"></a>
### 结果链接示例

以下链接指向同一结果根中的正式报告和已有证据。

| 要参考的情况 | 正式报告模块 | Notebook 的对应章节 |
| --- | --- | --- |
| 完整有效与共同有效 Moran 两个问题 | [HD Moran 分布](../../../output/reconstruction_impact/P1CRC_VisiumHD/report.html#module-moran_all_valid)、[共同有效基因](../../../output/reconstruction_impact/P1CRC_VisiumHD/report.html#module-moran_shared_valid) | HD 2.4 Gene spatial autocorrelation |
| matched-K 未成立，保留诊断和限制 | [HD matched-K](../../../output/reconstruction_impact/P1CRC_VisiumHD/report.html#module-matched_k)，查看 Mono/Macro、T | HD 2.2 Matched-complexity change |
| 无稳定阈值，连续场保留、覆盖 NA | [HD 阈值可靠性](../../../output/reconstruction_impact/P1CRC_VisiumHD/report.html#module-threshold_reliability)、[State Region](../../../output/reconstruction_impact/P1CRC_VisiumHD/report.html#module-region_extent) | HD 4.1–4.3 |
| 两条 Raw baseline 独立判读、共用 2×3 图 | [Xenium Raw Leiden](../../../output/reconstruction_impact/P2CRC_Xenium_ThreeParents/report.html#module-local_vs_raw_leiden)、[Raw Level2](../../../output/reconstruction_impact/P2CRC_Xenium_ThreeParents/report.html#module-local_vs_raw_level2) | Xenium 3.5–3.7 |

[HTML 模板](templates/report.html) 和 [Notebook 块模板](templates/notebook-block.md)
仍是两份 canonical shell；本页的文字模板负责内容顺序和填写职责。实际调用与出图格式参考
[HD source Notebook](../../../reproduce/case/reconstruction_impact/VisiumHD_sp_SVC_Reconstruction_Impact.ipynb)
和 [Xenium source Notebook](../../../reproduce/case/reconstruction_impact/Xenium_sc_SVC_Fibroblast_Reconstruction_Impact.ipynb)。

<a id="narrative-templates"></a>
### 文字模板与逐节点填写职责

文字模板只规定阅读内容，不增加科学计算。每个需要生成的节点 narrative 均按以下
顺序填写；foundation 节点默认收在按需 audit 中，不能因为节点存在而强制出现在主文：

```markdown
### <title.zh> / <title.en>
**简短结果 / Brief result.** <一句只引用本批保存事实的主线句。>

| scope | 事实、分母、状态和具名 evidence link |
| --- | --- |
| <parent or All> | <numbers, units, conditions, tables/figures> |

<已登记且能改变判断的图；没有图时保留数据入口，不补教学图。>

**Interpretations / 解读**
1. <事实及其双语表述。>
2. <在成立条件下的判断及其双语表述。>
3. <限制、反例或不可判断状态及其双语表述。>

**Next / 下一步（保存字段）.** <供 Notebook 和结果包追溯；网页不以此排序或导航。>
```

`node_narratives` 的 `title`、`lead`、每一项 `interpretations`、`next` 都必须同时
有 `zh` 和 `en`，并复用同一批 `record_ids`。`review.status=pending` 时只保存事实和
限制，不把未审阅句子当成结论。程序侧 `reading_placement` 的提示为：1.1–1.6=`audit`，
2.3=`reference`，其余节点=`main`；这是链接/归属提示，不是强制渲染 27 个节点的要求。

下面的职责表覆盖 27 个追溯节点；每行说明应填什么，不能用空标题代替：

| node_id | 填写职责 |
| --- | --- |
| 1.1 | 写明 Raw/Recon、空间/表达/reference 载体、观察单位、来源和解释边界；主文只给 audit 入口。 |
| 1.2 | 写明 parent exact label、覆盖/纳入/排除范围、抽样和各分析分母；主文只给 audit 入口。 |
| 1.3 | 写明 Raw/Recon ID 对齐、坐标单位/转换、覆盖和 paired 状态；主文只给 audit 入口。 |
| 1.4 | 写明两侧完整 gene space、表达来源及 unmeasured/zero/low coverage/not computable 语义；主文只给 audit 入口。 |
| 1.5 | 写明 Partition、Moran、EMT 各自输入、预处理、图/资源和状态；主文只给 audit 入口。 |
| 1.6 | 只汇总已保存的比较条件、限制传递和待检验项；没有证据保持 pending。 |
| 2.1 | 填 Raw/Recon K、resolution 或 fixed-final 条件、size/结构事实和诊断状态。 |
| 2.2 | 填 target/actual K、Hungarian mapping、change fraction、ARI、分母和 matched/unmatched 状态。 |
| 2.3 | 引用三个 parent 的已有变化和各自分母；作为 reference 归属，不重新计算一套 change。 |
| 2.4 | 分开填 Moran full-side 与 shared-valid 的有效数、Q75/分布、graph 条件和 paired delta。 |
| 2.5 | 分开填 EMT resource coverage、实际 cutoff、score 分布、paired n 和表达载体条件。 |
| 2.6 | 将分群、Moran、EMT 的事实并列到空间问题之前，不合成 improvement score。 |
| 3.1 | 填完整 Raw anatomy context、类别定义、覆盖和分母；图仅在能改变判断时使用。 |
| 3.2 | 复用已保存 mapping/change，填 changed units 的空间位置、scope 和 matched-K 限制。 |
| 3.3 | 复用 EMT score/delta 字段，填空间单位、坐标、carrier 和单侧/paired 状态。 |
| 3.4 | 填 Raw Leiden/Raw Level2/Recon、window support、selected scale、同 draws 和 baseline 角色。 |
| 3.5 | 分别填两条 baseline 下 Kobs 的 Raw/Recon/delta、window/scale 和状态。 |
| 3.6 | 分别填两条 baseline 下 Neff 的 Raw/Recon/delta、window/scale 和方向条件。 |
| 3.7 | 填 evenness 作为 richness/组成辅助事实，不套用 Neff 的方向判读。 |
| 3.8 | 填 anatomy 分层指标、delta、unit/window 分母和支持状态，区分数量与比例。 |
| 3.9 | 分开总结 assignment、EMT、局部 diversity 的空间事实，并指向 Region。 |
| 4.1 | 先填 threshold、bootstrap 有效率、CI width、window support 和可靠性状态。 |
| 4.2 | 填 Recon Neff 连续场、State threshold/mask 或 no-stable-threshold。 |
| 4.3 | 填 area、window/unit counts、area_fraction/unit_fraction 及各自分母；不可识别写 NA。 |
| 4.4 | 填同定义 rarefaction 指标/ delta 的尺度曲线和支持限制，不写成边界稳定性。 |
| 4.5 | 只在可靠性、状态和覆盖条件下说明 Region 能支持的解释；Gain 留作 audit。 |
| summary | 按五行结果和各自条件综合支持、相反证据或无法判断；不合成总分。 |

#### 既有已审阅解读示例（保留内容，仅改变阅读位置）

以下保留已有正式 HD 记录及已审阅中英解读，作为填写示例，数值不能复制到其他批次。
本节按科学节点编号引用，不改写这些解读；两份已生成 HTML 按固定的总体、局部邻域、空间
阅读位置挂载它们。小表展示数值，表外文字负责解释与承接；模板不是给所有字段再写一段定义。

**2.1 分群数量：Mono/Macro、T 增加，Fibroblast 群数不变。**
小表填三个 parent 的同 resolution Raw K → Recon K：`5 → 5`、`3 → 11`、
`3 → 19`。文字写：“群数增加主要见于 Mono/Macro、T；Fibroblast 群数相同，
仍不能推断成员没有重新分配。”接着引出 matched-K。旧的单侧 reconstructed-cluster
图描述后续局部分析使用的分群，不冒充同 resolution 数量比较图；保留在具名展开证据中。

**2.2 归属变化：仅 Fibroblast 可作 matched-K 比较。**
小表填匹配状态及有效结果：Fibroblast `5/5`，`13,705/23,151 = 59.2%`；
Mono/Macro、T 均为 `3/5`，不成立。文字写：“Fibroblast 的群数未变，成员分配却
发生变化；另外两类不能把群数差异与归属变化充分分开。”未匹配的诊断数值放展开证据，
不当作同等条件的 headline。承接：“分群结构变化之后，再看表达空间组织。”

**2.4 Moran：共同基因中仍可见空间相关性升高。**
分别填写各自有效基因分布及共同有效基因配对证据。Fibroblast Q75 为
`0.00495 → 0.18682`；共同有效 `12,728` 基因的配对 ΔI 中位数为 `+0.17038`。
文字写：“升高方向并不只出现在不同基因集的比较中；共同基因配对结果也支持它。
但 Moran 本身不能区分结构恢复与表达平滑。”另外两类按相同问题填写，随后转向 EMT。

**2.5 EMT：score 上移，同时覆盖与评分截断变化。**
Fibroblast 小表填 coverage `198/200 → 189/200`、score median
`0.01389 → 0.07676`、配对 Δ 中位数 `+0.06260`、有效排名长度 `20 → 193`。
文字写：“覆盖没有增加，评分仍上移；评分输入和截断条件同时不同。因此可描述现有
条件下的 score 变化，不能直接认定 EMT 生物学活动增强。”承接到空间位置，
不把总体 score 差异当作已经定位了 hotspot。

**3.1–3.3 空间位置：改变比例高与改变数量多要分开。**
Fibroblast assignment 小表列 Interface `222/314 = 70.7%`（8 有效窗口），
Normal `441/629 = 70.1%`（17），Tumor `2,076/3,348 = 62.0%`，
Other `10,966/18,860 = 58.1%`。文字写：“Interface、Normal 比例较高，但支持
范围小；Other 改变数量最多，同时具有最大的配对基数。”这些是 assignment 数据，
不能放进 State Region 覆盖表。EMT 图另写已审阅的定性位置观察，标注“图上观察”；
没有区域定量证据时，不命名 anatomy 富集或 EMT hotspot。

**3.4–3.8 局部状态：同一重建结果对两条 baseline 的含义不同。**
先给实际窗口支持，再分别展示 Kobs、Neff、evenness 的小表及各自 2×3 图。
Fibroblast 配对窗口 ΔNeff 中位数相对 Raw Leiden 为 `−0.9603`，相对 Raw Level2
为 `−0.2667`。文字写：“相对 Raw Leiden 支持局部有效状态更集中的方向，
却没有支持相对 Raw Level2 呈现更细状态。不能只选择其中一条判断改善。”
Kobs、evenness 填各自结果，不套用 Neff 判断；anatomy 分层继续核对区域方向和支持窗口。

**4.1–4.4 State：只有 Fibroblast 能给出稳定的高多样性范围。**
先填可靠性，再看连续场/mask 与覆盖。Fibroblast 阈值约 `1.780`，有效窗口面积
占比 `20.54%`；Mono/Macro、T 无稳定阈值，覆盖 `NA`。文字写：“可对 Fibroblast
展开高多样性状态范围，其余仍保留连续 Neff；20.54% 不是改善面积，也不是变化最大
区域。”最后用已有尺度结果说明局部指标方向如何依赖尺度，不把它改称边界稳定性。

**结尾：回到 impact，不重复所有数值。**
“本批重建表达呈现更强空间相关性；Fibroblast 在 K 不变时仍有成员重分配。
部分 parent 相对 Raw Leiden 更集中，但相对 Raw Level2 的更细状态方向未获支持。
EMT score 的变化受覆盖与截断条件限制；仅 Fibroblast 可识别稳定 State Region。”
各句链接对应节点；其他批次按实际支持、未支持及无法判断的结果改写，不沿用本段结论。

## 6. Authoring 与 provenance 细则

| 对象 | 作者 | 消费者 | 允许的内容 |
| --- | --- | --- | --- |
| 科学方法 | 本目录三份方法文档 | Notebook、record、HTML 链接 | Purpose、Computation、Meaning、Direction/conditions |
| 执行块 | source Notebook/builder | executed notebook、runner | 参数、短调用、audit、plot、brief observation |
| 长实现 | local helper/现有分析模块 | Notebook block | 可测试实现；不直接写 HTML 结论 |
| 结果事实 | notebook writer/record builder | HTML、下游 batch | 数值、分母、状态、证据路径 |
| 展示 | renderer + 两个 canonical templates | `report.html` | 固定左侧目录、单页长页、五行 overview、正文顺序、图、链接、transition |

本轮只调整 renderer/template 和正式 HTML 的阅读呈现，保留科学结果和已审阅解读。
不改 Notebook 科学步骤，不增加执行身份或科学重计算。

## 7. R0–R4 验收规范

以下保留既有科学执行与验收边界；本轮展示检查单独记于[规范入口](README.md#report-display-status)，
不把历史验收记录改标为本轮重新通过。

| 迭代 | 必须证明 | 证据 |
| --- | --- | --- |
| R0 | 固定左侧目录和单页长页；五行 overview、scope、旧图表与方法文档与框架只有一套规范；基础内容只在按需 audit | README/方法链接检查、无 competing active spec |
| R1 | 记录粒度正确；HD All 仅已有六类 global record 且独立；Xenium 无 global；Moran full/shared、EMT coverage/score、两 baseline 分开；NA/状态可见 | report_records JSON/CSV 重载、真实结果抽查 |
| R2 | Notebook 分块仍保留核心参数/短调用；cohort、seed、full-gene 和共同 draws 不变；长实现确在 helper | source/executed cell identity、对应结果与必要调用检查 |
| R3 | overview 恰为五行且每格一项事实/一句结论；网页正文按总体→局部邻域→空间及固定问题标题→简短结果→表/图→1–3 条解读；common 图只嵌一次；scope 图归属正确 | HTML parser/visual QA、named links |
| R4 | 阅读产物刷新保留原 execution identity；若科学依赖未变，不要求真实重跑；若依赖改变才需新身份；scientific CSV 数值身份可追溯 | execution/validation、artifact index、fresh reading evidence |

发布前还须检查：Moran 同组同图且最低 51 units；AUCell 同输入/参数且保留 coverage
与 cutoff；三 assignment 同窗口同 draws；anatomy、baseline、Region 分母独立；
threshold reliability 在 mask 之前；renderer/report 只读保存结果。纯呈现改造造成的
科学数值差异必须解释并修复，不能以 HTML 能打开代替验收。

<a id="acceptance-criteria"></a>
<a id="current-acceptance"></a>

### 本轮呈现修改

固定目录与正文采用同一阅读树；原有科学节点、问题编号、记录和审阅身份保持。
两份报告已完成本轮实现及展示检查，用户阅读效果尚未确认。检查范围与结果仅在
[规范入口](README.md#report-display-status)维护，避免重复状态台账。
