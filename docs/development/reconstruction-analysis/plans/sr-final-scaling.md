# R1：sST 固定 spot 总量校正

**计划已对齐，R1 已获实施授权；实施与验证状态见[总览](../README.md#items)。**

[返回重建分析总览](../README.md) · [R1 小改动分类](../categories/01-small-changes.md#r1) · [共同约束](../decisions-and-sources.md#confirmed-decisions)

依据：[共同决策与暂停说明](../decisions-and-sources.md#execution-boundary)；旧本地草案仅为可选历史线索。

## 目标与范围

固定 sST 重建结束时的 parent spot、逐基因总量校正。

删除生成细胞各自缩放到 10,000 的最后分支。

保留 spot/reference 进入算法前的内部 `normalize_total`。

保留建图临时副本上的 normalization/log1p。

保留 posterior/reference allocation、quota 和 OT 图聚合。

不恢复原始 raw counts，不修改 Confidence，不修改 benchmark 算法。

此前暂停的 19 文件修改已纳入本批复核，验证结果在本文末尾记录。

## 已确认行为与代码证据

- [sST application runner](../../../../revise/backend/runners/sc_svc_super_resolution_application.py) 先建立共同基因轴和 normalized 工作矩阵。
- 同一 runner 的最终 `current_sum` / `ratio` 输出校正处；固定目标应是每个 parent spot 的共同基因逐列总和。
- [topology 建图](../../../../revise/backend/ops/topology.py) 的 normalize/log 只服务于图副本。
- [Application 配置](../../../../revise/application/config.py) 已有 unknown-key 拒绝入口。
- [engine loader](../../../../revise/config/loader.py) 已有 section unknown-key 校验。
- [Visium template](../../../../revise/application/templates/Visium.yaml) 不能继续暴露已经删除的二选一开关。

上述链接用于定位当前实现；它们不是本轮已经完成修改的证明。

## 具体模块与流程

1. 在 sST runner 中无条件执行现有 `current_sum`、逐基因 `ratio` 和 parent-spot 回写。
2. 从 `ApplicationScSuperResolutionConf` 和 adapter 构造中移除 `rec_match_spot_sum`。
3. 从 authority 默认、engine `sc` 允许键、application YAML 编译字段和 publication metadata 中移除开关。
4. 删除两个 Visium 模板中的开关；保留 graph、strength、PM 等已有输入。
5. 旧配置遇到该键时沿现有 unknown-key 路径显式失败，不静默忽略。

用户已确认：Application 和 Engine 旧键无论 true/false 均报错，明确提示删除该键；sST 固定采用 parent-spot 逐基因校正，不再支持最终逐生成细胞缩放到 10,000。Application 复用迁移消息表，Engine 在 sc 校验入口报告完整字段路径。

## 接口变化

拟移除内部 runner 字段 `rec_match_spot_sum`。

拟移除 application `local_refinement.match_spot_sum` 和 engine `sc.match_spot_sum`。

public output 继续使用现有 sST `svc` 角色；不新增 per-cell artifact。

不改变外部分析端的字段声明，也不把 normalized 输出声明成 raw count。

## 失败与恢复

逐基因校正出现空轴、非有限值或 ID 不一致时，当前样本失败，不安装新的结果。

发布沿用已有临时文件与 commit/rollback 回调；必要的输出校验在恢复实施时验证，不能把当前 publisher 描述为已经完成全部校验。历史产物受保护。

旧开关配置在编译阶段失败，避免不同调用者得到不同尺度。

不得通过回退到每细胞 10,000 或恢复 raw counts 来掩盖失败。

## 已对齐实施步骤

1. 复核开关的全库引用并建立最小 diff。
2. 先改 runner 与 runner config，再改 adapter、authority、loader、application config。
3. 更新模板、publication metadata 和工程文档。
4. 增加多生成细胞行为用例。
5. 用 `/Users/stephen/miniconda3/bin/python` 运行受影响的针对测试。
6. 核对 `git diff --check`、旧键拒绝和未改动的前置 normalization。

## 针对验证与验收

- source spot 的共同基因行和仍以 10,000 为算法内部目标。
- 对分配有非零支持的 fixture，多个生成细胞按 parent spot 汇总后逐基因在数值容差内等于 normalized spot X；零支持不应凭空造表达，也不能宣称此时精确守恒。
- 生成细胞不再全部各自等于 10,000。
- topology 的临时 normalization/log 行为未变。
- `match_spot_sum` 旧键被拒绝，模板不再包含该键。
- 不运行真实大样本，不以 benchmark 结果替代行为测试。

## 用户待决定、材料与必停边界

待决定：暂无。旧配置明确迁移文案已确认。

不需要用户提供科学材料；若发现输出尺度声明与分析端契约冲突，暂停并上提。

本批只授权 R1 代码和针对验证；不运行真实大样本，不启动 R2/R3。

不得扩大到 Confidence、Raw publication、benchmark 或 cluster assembly。

## 实施与验证记录（2026-09-19）

R1 已完成代码及针对行为验收，未提交。已复核此前 19 文件差异；算法仅固定原有 spot 逐基因校正，删除最终逐细胞 10,000 分支。Application/Engine 迁移提示覆盖 true/false，移除旧键后的配置正常编译；参数、模板和发布元数据已清理。内部 ST/reference normalization 与 topology normalize/log 源码未改。活动 application reference、batch protocol 说明及 runner 注释同步更正；历史归档不重写。

验证结果：**385 passed，14 warnings**，覆盖 `tests/application`、`tests/config/test_ot_config.py`、`tests/backend/test_application_column_contract.py`、`tests/backend/test_sc_local_ot.py`、`tests/backend/test_sc_sr_guidance.py`、`tests/backend/test_tacco_global_freshness.py`、`tests/batch/test_runner.py`。新 sST fixture 每个 spot 生成多个细胞，验证 normalized spot 逐基因总量及不再逐细胞固定 10,000；POT/TACCO 调用替身的现有 OT 用例增加了最终逐 spot、逐基因守恒与非逐细胞 10,000 断言，graph 相关用例继续通过。`git diff --check` 通过，活动代码中旧键只保留在拒绝与迁移提示路径。

环境限制：默认 Python 3.12 的 pytest 在导入 readline 时以 139 退出，关闭插件及限制线程仍复现；Python 3.10 也在 readline 导入处退出。最终使用现有 Python 3.10，在测试进程内将 readline 标为不可导入，保留 pytest capture/fixture 机制，关闭第三方自动插件、限制线程并禁用 Numba JIT 后完成测试。没有修改测试环境安装或仓库运行逻辑。首次实际执行为 382 passed/3 failed，失败来自新增测试缺少 route runtime 及 fixture 意外进入无关建图分支，修正夹具后全套通过。

可复现调用：

```sh
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_DISABLE_JIT=1 MPLCONFIGDIR=/tmp/revise-r1-mpl NUMBA_CACHE_DIR=/tmp/revise-r1-numba /Users/stephen/miniconda3/envs/python3.10/bin/python - <<'PYTEST'
import sys
sys.modules["readline"] = None
import pytest
raise SystemExit(pytest.main([
    "tests/application", "tests/config/test_ot_config.py",
    "tests/backend/test_application_column_contract.py",
    "tests/backend/test_sc_local_ot.py", "tests/backend/test_sc_sr_guidance.py",
    "tests/backend/test_tacco_global_freshness.py", "tests/batch/test_runner.py", "-q",
]))
PYTEST
```

警告为现有依赖弃用、故意溢出失败用例与 AnnData 索引转换。此结果只证明上述工程行为；未运行真实大样本、TACCO 集成 smoke 或 batch 真实样本 parity，不代表科学效果或正式发布验收。R2/R3 保持未实施。
