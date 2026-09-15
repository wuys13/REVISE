# Batch 接入与历史 handoff 边界

分析目的与共用节点见 [分析框架](analysis-framework.md#content-tree)。本页仅管理后续接入边界，
不把本轮文档刷新扩大为 batch 开发或真实重跑。

本页只规定 reconstruction-impact 进入既有 2.0 batch framework 的边界。它不是新的
科学层、不是当前 Notebook 入口，也不复制输入、Moran、EMT、局部状态或 Region
定义；这些定义分别归 [input-views.md](input-views.md)、
[impact-analysis.md](impact-analysis.md) 和 [gene-and-function.md](gene-and-function.md)。
当前输出/记录契约见 [outputs-and-test-plan.md](outputs-and-test-plan.md)。本地内容映射仅用于Notebook与报告展示，不是batch插件接口或新增调度配置。

## 1. Current boundary

当前正式执行仍由 schema 1 route source notebooks 驱动：

- [VisiumHD source notebook](../../../reproduce/case/reconstruction_impact/VisiumHD_sp_SVC_Reconstruction_Impact.ipynb)
- [Xenium source notebook](../../../reproduce/case/reconstruction_impact/Xenium_sc_SVC_Fibroblast_Reconstruction_Impact.ipynb)
- [builder](../../../reproduce/case/reconstruction_impact/build_notebooks.py)

每个 sample root 的 `report.html` 只读取保存的 CSV/JSON、conclusions、figures 和
review status；它不调用 batch adapter 或重算 H5AD。现有 2.0 的 `revise/batch/`
代码和 runner 可以作为未来接入参考，但不能把 adapter 的绿色测试写成当前 route
的科学证据。正式根仍由 schema 1 `output.dir` 决定：

```text
output/reconstruction_impact/<sample_id>/
```

不建立第二个开发输出树、不建立第二套 report records、不把未来 schema 2 文档写成
当前已完成能力。

## 2. Reusable 2.0 responsibilities

未来 batch 只复用已经通过 Notebook 和真实结果审阅的输入/结果契约：

| 2.0 能力 | 参考来源 | 接入边界 |
| --- | --- | --- |
| schema、分层 YAML、sample/path discovery | `revise/batch/config.py`、`sample.py` | 先映射已保存 route identity；不改变五行 overview 或 27 节点追溯 |
| reconstruction runner、success/failure handoff | `revise/batch/runner.py` | 只有 verified handoff 或明确 `imported_legacy` 才可消费 |
| aspect staging、恢复和 artifact ownership | `revise/batch/analysis.py` | 未来消费稳定 Notebook 表；不替代 Notebook 分块计算 |
| Raw/native/paired/sST read-only views | `revise/batch/inputs.py` | 先用 foundation audit 验证轴、cohort、projection |
| iST sparse cluster mean/mean/random publication | `revise/application/ist_assembly.py`、`publication.py` | 先保留 expression carrier、donor 和 seed provenance |
| root scripts/CLI | `batch_reconstruct.py`、`batch_analyze.py`、`revise/batch/cli.py` | 作为 future orchestration；不是当前 source notebook 总入口 |

2.0 的通用输入树、ST/SC protocol、runner 职责和 recovery 规则由现有 batch docs
负责；本页只说明 impact-specific 的 handoff 和 aspect 对齐。未来文档应链接本目录，
不要复制另一套 reconstruction-impact 计划。

## 3. Result-to-aspect mapping

稳定 Notebook 结果可以在未来映射为 aspect-owned artifacts，但 aspect 只是发布和
失效边界，不改变科学问题：

| Future aspect | 消费的已保存事实 | 明确不做的事 |
| --- | --- | --- |
| `reconstruction_impact` | partition、matched-K、anatomy/window/diversity、State/Gain | 不隐式启动 Moran/AUCell |
| `spatial_autocorrelation_by_cell_type` | Moran full/shared、availability、组内 graph 和逐 gene 结果 | 不生成额外的逐 gene spatial module |
| `pathway_activity` | EMT resource、coverage、unit score/summary 和 spatial fields | 不把 coverage/score 合成 biological validation |

每个 aspect 未来需要登记 input view、comparison basis、seed/parameters、axis/resource
identity、科学实现 digest、staging 文件和 artifact index。三个 aspect 仍各自保存
状态与 evidence；一个 HTML 页面把它们放在一起不代表它们变成一个分析类别。

未来公开调用可以继续采用 `run_reconstruction_task(config_path, sample_id, *,
cell_type=None)` 和 `run_analysis_task(config_path, sample_id, aspect, *,
cell_type=None)`，但调用名字不表示当前 Notebook 已经由它们驱动。安全 task ID
与 Raw exact label 必须继续显式映射：`Mono_Macro → Mono/Macro`。

## 4. Historical handoff

当前 impact 读取 `results/` 中的重建 H5AD。未来 batch 消费前必须得到标准任务的
`.revise/task.json`、`reconstruction.json`、输入/输出身份和 fingerprint；结构正确
的 H5AD 本身不等于 verified handoff。

| 选择 | 需要登记 | 代价/边界 |
| --- | --- | --- |
| 重新运行 reconstruction | 标准化输入、route/config、输入 digest、engine provenance、success record | 成本较高，但不会伪造历史 verifier；重建差异与 impact 差异分开核对 |
| 显式接管历史输出 | raw/reference/output hash、来源配置、route、轴、表达语义、`imported_legacy` 状态和缺失 provenance | 可复用昂贵输出，但不能把缺失 handoff 改写成 verified |

在该选择完成前，历史表只能作为迁移核对或 fixture。仅把 schema 1 路径改写为
schema 2，不足以声明所有旧结果可批量消费。

## 5. Failure、ownership 与 invalidation

未来 batch 必须把实际科学实现、Notebook helper、renderer/template、resource 和
输入 view 登记为依赖。依赖变更只使对应 aspect 失效；reconstruction fingerprint
和 analysis fingerprint 独立。资源/参数/科学实现变化不能继续沿用旧的 reviewed
结论。

低支持、未测量、不可计算和 `no_stable_threshold` 是 aspect 成功后可发布的科学
状态；缺依赖、资源损坏、invalid axes、代码错误和 staging 不完整是 failure/block。
失败尝试不能覆盖或冒充上一版成功产物；报告同时显示当前尝试状态和旧结果身份。
发布前必须验证声明 artifact 集合与实际文件集合一致，拒绝未经授权的额外文件，
不自动删除用户文件。

HTML shell 和 Notebook block template 仍只有两个 canonical paths：
`templates/report.html` 和 `templates/notebook-block.md`。template 内容 hash 属于
report/source provenance；未来 batch 不复制按层 template，也不在 batch runner
内嵌新的展示规则。

## 6. Future migration order

只有当前 R4 route 验收和独立 batch 授权都完成后，再按顺序推进：

1. 冻结 Notebook source/executed、输入、输出表、记录、template 和实现 digest，确定
   historical handoff 选择；
2. 用小 fixture 验证真实 ID/坐标、expression view、cluster mean、resource 和状态；
3. 将稳定 Notebook 表映射到 aspect-owned staging/publish tree，report 只消费已发布
   artifacts；
4. 分别接入 Moran 与 EMT aspect，保持 full/shared、coverage/score 和缺失状态；
5. 再迁移 schema 1 消费者，给旧入口增加迁移提示，并让 notebook/batch 同任务结果
   可逐表核对。

adapter/runner 测试不能替代 source/executed parity、真实 CSV/JSON reload、图/模板
   provenance 和 HTML/Notebook 验收。任何新的 batch 设计都必须回到本目录的唯一问题
   矩阵和三份方法文档，不另建活动规范。
