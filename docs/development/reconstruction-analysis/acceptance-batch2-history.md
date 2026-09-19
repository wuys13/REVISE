> 历史快照：第二批交付时的验收记录。旧表达支持结论不代表当前协议；当前状态见 [验收入口](acceptance.md)。

# 第二批验收清单：统一样本交付与分析衔接

> **协议更新（2026-09-19 OT assembly 批次）：** 当前消费者已固定线性非负、未取 log 的 X；旧 log 声明已拒绝。下文此前的表达 scale 支持范围属于历史快照，现行声明和本批处理以[新计划](plans/ot-assembly.md)为准。真实验收继续后置。
**状态：生产端实现及工程验收完成；真实样本与科学/联合验收待完成。** 这是 R2/R3/R4 的一页验收入口。最终工程套件已经通过；合成 fixture 的 loader/标签空间步骤也已验证。P2CRC_Xenium 命令保留为后续真实输入入口，本次文档更新不运行真实数据。

[返回总览](README.md) · [第二批总计划](plans/whole-sample-delivery.md) · [R2](plans/whole-sample-svc.md) · [R3](plans/raw-publication.md) · [R4/C1](plans/analysis-handoff.md)

## 使用边界

- 工程测试、真实输入消费、分析结果状态和用户科学审阅分开记录。
- 通过合成 fixture 或 loader 只能证明工程路径；不能证明真实样本的科学含义或生物学改善。
- 本批默认不执行昂贵的全量重建。下列 P2CRC_Xenium 命令是验收时使用的入口，当前未运行。
- `REVISE_Analysis_Agent` 是独立消费者；验收只调用其正式入口，不修改其并行工作树。
- sST 的 normalized nonlog 表达尚未具备完整消费方案；它必须在结果中标为 unavailable/partial，不得声明为 raw counts 或 log1p。

## 当前进展（2026-09-19）

- 最终生产端与交接套件通过 **759 passed, 43 warnings, 0 failed, 0 skipped**，退出码为 0；测试范围和源版本摘要见 [`verification-2026-09-19.json`](verification-2026-09-19.json)，完整日志为 `/tmp/revise-r2-final.log`。其中包含 16 项 sample-delivery 针对检查，以及两个交接/端到端集成检查。
- `tests/integration/test_analysis_delivery.py` 的正式消费者入口检查已通过：`load_sample` + `run_analysis(reconstruction_impact)` 完成标签/空间 diversity 和报告生成，输入哈希保持不变；由于 expression identity unknown，整体结果记录为 `partial`，没有把表达分析报成成功。全样本端到端测试还验证真实 application 预处理/pipeline/publication、GA 单次调用、多类型 skip/合并、Raw QC 排除单位保留、SVC 组成及 direct/batch 的 X/obs/var/坐标一致。
- 工程执行中处理的问题已留档：Python 的可选 readline 导入崩溃通过仅在测试进程中禁用该模块解决；CLI 沙箱内的 OpenMP 共享内存限制通过沙箱外重跑验证。测试另禁用插件自动加载及数值多线程，并修复 lazy kernel 缓存隔离，避免跨测试污染。43 个 warning 来自依赖弃用提示及刻意覆盖的 NaN/overflow、AnnData view/index、重复 ID 等边界 fixture，不构成失败。
- 上述证据使用合成 fixture，不替代 P2CRC_Xenium 的真实样本验收；当前没有运行真实 P2 全量重建，也没有完成用户科学审阅。坐标物理比例或表达 identity 缺少来源时，真实消费者应依据实际 `status`/`error` 判断可用性，不预先承诺所有空间或表达步骤只返回 `partial`。

## A. 工程验收

### R2：全类型重建与统一 SVC

- [x] 一个 whole-sample fixture 含至少两个 eligible broad types、一个只有单 Level2 的类型以及一个缺失/空白 Level2 的类型。
- [x] GA 只调用一次；所有 eligible 类型共享该结果，类型遍历顺序和 seed 可复现。
- [x] whole-sample 的单 Level2、缺失或空白类型产生可定位的 warning/skip 原因，不进入 SVC 合并；显式旧单类型调用对一个有效 Level2 保留成功，对零个有效 Level2 保持失败。
- [x] eligible 类型的 LR/assembly 失败使整个 sample failed；不发布空 SVC 或部分成功 SVC，历史完整结果保持不变。
- [x] whole-sample 默认 `random`；`mean` 仍可用；跨 broad type 的 cluster 标识不碰撞；合并前 gene axis 和 obs ID 校验有效。
- [x] 旧 `select_ct`、`select_cell_type`、batch `cell_type` 兼容调用仍保持单类型语义；旧单类型成功记录不会使 whole-sample 任务误判为已完成。

### R3：完整 Raw 与真实推断

- [x] fixture 含会被 QC 排除的原始单位；发布 Raw 仍包含完整原始 obs/var 轴、X、坐标和原始注释。
- [x] 发布前后 Raw 的 X、轴顺序、坐标和 observation ID 一致，证明没有使用 normalized 工作对象冒充 Raw。
- [x] 只有真实获得的 ID 回填 `revise_Level1`/`revise_Level2`；没有推断的位置继续缺失。
- [x] 原始同名标签不被覆盖；原始与推断冲突记录统计；输入已占用新增列名时明确失败。
- [x] Raw、SVC、`sample.yaml`、manifest 和运行记录来自同一成功运行；计算、写入或安装失败保护历史完整产物。

### R4/C1：交接与消费

- [x] `sample.yaml` 使用相对路径，层级 sample ID 映射确定、可逆且无碰撞；从其他工作目录调用仍能找到 Raw/SVC。
- [x] 配置指定推断列名，仅在存在时声明重建 cluster；未知单位或比例保持 unknown/缺失，不伪造标签。
- [x] 两侧 Raw/SVC 可独立加载，行数、顺序、ID 或基因集合不同不被强制伪配对；输入 H5AD 在分析前后字节或内容未被改写。
- [x] 在合成 fixture 上使用分析库正式入口完成至少一个已有前提的标签/空间分析；结果状态、错误或 unavailable 原因与文件一致。
- [x] sST normalized nonlog 的表达步骤单独标注未完成，不阻止与表达无关的标签/空间检查。

## B. P2CRC_Xenium 真实验收入口

当前仓库已存在以下真实输入（是否在验收环境可读仍需执行前核对）：

```text
raw_data/Real_application/P2CRC_Xenium.h5ad
raw_data/Real_application/adata_sc_all_reanno.h5ad
```

### 现有单类型控制命令

该命令核对当前旧 direct application 路径，只作为兼容基线，不能代替第二批 whole-sample 验收：

```bash
cd /Users/stephen/Documents/wuyushuai_project/wuyushuai_research_project/REVISE
/Users/stephen/miniconda3/envs/python3.10/bin/python reconstruct.py \
  --config configs/application/Xenium.yaml \
  --select-ct T
```

`configs/application/Xenium.yaml` 当前明确指向上述 P2CRC ST/reference，并显式设置
`output.ist_mapping: random`，因此这个单类型兼容命令输出一个 assembled
`SVC.h5ad`，不是 paired 双载体。若要核对旧的 paired 返回形式，应使用一份不声明
`output.ist_mapping` 的专用配置；不要修改维护中的模板。执行前确认输出目录
`results/sc_SVC_case/P2CRC_Xenium/T/` 是专用验收目录；不要用旧单类型产物冒充新
sample-level 结果。

### 第二批 whole-sample 命令

R2/R3 实施完成后，直接使用维护中的 whole-sample Application 配置。当前
`configs/application/Xenium.yaml` 已省略 `select_cell_type` 并将 mapping 设为
`random`，因此不需要先搭建 batch 目录：

```bash
cd /Users/stephen/Documents/wuyushuai_project/wuyushuai_research_project/REVISE
/Users/stephen/miniconda3/envs/python3.10/bin/python reconstruct.py \
  --config configs/application/Xenium.yaml
```

成功后预期的本批交付入口为：

```text
/Users/stephen/Documents/wuyushuai_project/wuyushuai_research_project/REVISE/results/sc_SVC_case/P2CRC_Xenium/raw.h5ad
/Users/stephen/Documents/wuyushuai_project/wuyushuai_research_project/REVISE/results/sc_SVC_case/P2CRC_Xenium/SVC.h5ad
/Users/stephen/Documents/wuyushuai_project/wuyushuai_research_project/REVISE/results/sc_SVC_case/P2CRC_Xenium/sample.yaml
```

这些路径只有在本次 whole-sample 运行成功后才可作为当前结果；命令不会在本轮
文档更新中执行。batch orchestration 仍由仓库已有的 `batch_reconstruct.py`
入口覆盖，但当前 P2CRC 的直接 Application 命令是本批代表性真实入口。

在独立消费者环境使用生成的 `sample.yaml`：

```bash
cd /Users/stephen/Documents/wuyushuai_project/wuyushuai_research_project/REVISE_Analysis_Agent
PYTHONPATH=/Users/stephen/Documents/wuyushuai_project/wuyushuai_research_project/REVISE_Analysis_Agent \
/Users/stephen/miniconda3/envs/python3.10/bin/python -m revise_analysis.cli run \
  --sample /Users/stephen/Documents/wuyushuai_project/wuyushuai_research_project/REVISE/results/sc_SVC_case/P2CRC_Xenium/sample.yaml \
  --analysis reconstruction_impact \
  --output /Users/stephen/Documents/wuyushuai_project/wuyushuai_research_project/REVISE_Analysis_Agent/output/P2CRC_Xenium_r2
```

这里的输出目录应使用专用验收位置。由于本生产输入当前没有进一步确认的
表达 identity/scale 和物理坐标校准，可检查已具备前提的标签能力；表达或需要物理尺度的步骤须依据实际
`status`、`unavailable` 和 `error` 判断，不预先承诺整条分析成功。不要从现有 consumer 迁移样本
的 `0.2125` 校准值倒推本次新交付的坐标含义。

## C. 已完成的工程验收与后续检查

| 验收层 | 本批结果 | 证据与边界 |
|---|---|---|
| R2 类型覆盖、GA 单次调用、SVC 合并 | 通过 | backend/application 测试及全样本 direct/batch 端到端 fixture；计算核使用确定性替代，未跑真实大样本 |
| R2 跳过、eligible 失败和兼容入口 | 通过 | 类型资格、全跳过、错误传播、显式单类型和恢复测试 |
| R3 Raw 完整性与标签来源 | 通过 | 原始 X/轴/坐标、QC 排除单位、冲突/缺失、源文件变化和内存别名保护 |
| R3 事务与恢复 | 通过 | 写入、安装、记录失败时旧产物不变；新指纹与旧入口分别处理 |
| R4 配置与真实消费者入口 | 通过 | 独立 Raw/SVC 轴、相对路径、正式 loader 和标签空间分析，输入哈希不变 |
| 表达分析 | 未完成 | fixture 整体 partial，表达 identity unknown；sST normalized nonlog 仍依赖 C1 |
| P2CRC_Xenium 真实运行与科学审阅 | 未执行 | 按上方命令运行后核对类型覆盖、Raw/SVC、缺失原因与科学解释 |

最终范围、版本与结果见[验证记录](verification-2026-09-19.json)。测试日志 `/tmp/revise-r2-final.log` 是本机临时证据，下面的命令可重新生成验证结果：

```bash
cd /Users/stephen/Documents/wuyushuai_project/wuyushuai_research_project/REVISE
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
MKL_NUM_THREADS=1 NUMBA_DISABLE_JIT=1 MPLCONFIGDIR=/tmp/revise-r2-mpl \
NUMBA_CACHE_DIR=/tmp/revise-r2-numba \
/Users/stephen/miniconda3/envs/python3.10/bin/python - <<'PYTEST'
import sys
sys.modules["readline"] = None
import pytest
raise SystemExit(pytest.main([
    "tests/application", "tests/config", "tests/backend", "tests/batch",
    "tests/integration/test_analysis_delivery.py",
    "tests/integration/batch/test_whole_sample_delivery.py", "-q",
]))
PYTEST
```

联测要求兄弟库 `REVISE_Analysis_Agent` 存在；另一个环境可通过 `REVISE_ANALYSIS_ROOT` 指定。若未提供消费者，联测会明确 skipped，此时不能沿用本次已验证消费的结论。真实样本未被本轮运行，分析库工作树未被修改，未提交、推送或创建外部 issues。
