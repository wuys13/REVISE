# Batch validation and delivery

本文定义验收场景，并在[交付验证](#verified-delivery)记录本轮实际结果。[返回主计划](plan.md#implementation-order)。先读对应的[输入语义](protocol.md#carriers-and-alignment)或[执行约定](framework.md#task-api)，再按场景定位测试。以下是当前测试支点；场景清单与实际运行结果分开列出。

## Acceptance scenarios

### Input scenarios

| 场景与判定 | 当前测试支点 |
|---|---|
| 多个 ST 共用 reference；输入内容摘要和目录清单执行前后不变 | [test_sample.py](../../../tests/batch/test_sample.py)：`test_shared_inputs_are_read_only_without_new_files`；[test_runner.py](../../../tests/batch/test_runner.py)：`test_reconstruction_preserves_original_input` |
| 换启动目录，继承与资源路径仍按声明位置解析；新增／删除中间 YAML 被发现 | [test_config.py](../../../tests/batch/test_config.py)：`test_inheritance_paths_lists_and_reference_replacement`、`test_chain_detects_added_removed_and_changed_configuration`、`test_analysis_resources_resolve_relative_to_declaring_yaml` |
| paired 稀疏均值等于现有 mean；打乱 donor 顺序不破坏空间对应；donor 不能作为空间行 | [test_ist_publication.py](../../../tests/application/test_ist_publication.py)：`test_mean_is_sparse_and_uses_spatial_axis`；[test_inputs.py](../../../tests/batch/test_inputs.py)：`test_paired_view_keeps_donor_axis_and_maps_sparse_means_to_spatial_order`、`test_paired_view_rejects_mismatched_or_null_cluster_labels` |
| iST raw 保留实测 panel，reference 独有基因不在 raw 补零；输出 ID 子集正确对齐 | [test_inputs.py](../../../tests/batch/test_inputs.py)：`test_raw_alignment_uses_ids_and_coordinates_and_records_uncovered_rows`、`test_raw_alignment_rejects_missing_ids_or_changed_coordinates` |
| sST 用实际生成数均分：不等 cell 数、输出乱序、稀疏矩阵均逐 spot 逐基因守恒 | [test_inputs.py](../../../tests/batch/test_inputs.py)：`test_sst_baseline_equal_splits_raw_spots_in_output_order_and_preserves_coordinates` |
| sST 未覆盖 spot 明确记录；缺失／未知父映射失败；沿用输出 ID 与坐标且不改变重建表达 | [test_inputs.py](../../../tests/batch/test_inputs.py)：`test_sst_baseline_rejects_invalid_actual_output_mapping`、`test_sst_baseline_rejects_native_pairing_request` |
| mean/random 正常读取，随机可追踪；模式切换失败可回滚 | [test_ist_publication.py](../../../tests/application/test_ist_publication.py)：`test_random_donors_are_seeded_and_reference_order_independent`、`test_failed_switch_restores_old_pair`；[test_inputs.py](../../../tests/batch/test_inputs.py)：`test_mean_and_random_modes_use_svc_and_preserve_mapping_provenance` |

### Execution and reuse scenarios

| 场景与判定 | 当前测试支点 |
|---|---|
| 单项只运行一个样本／类型／方面；与批量调用同一执行器、同参数结果一致 | [test_reconstruction_api.py](../../../tests/batch/test_reconstruction_api.py)：`test_single_task_runs_only_the_selected_sample_and_type`、`test_single_task_reuses_without_batch_inventory`、`test_single_task_rejects_out_of_scope_selection`；[test_public_api.py](../../../tests/batch/test_public_api.py)：`test_public_task_signature` |
| 单项重建后可单项分析，无批次清单也能运行；锁不重入，并发仍拒绝 | [test_analysis_api.py](../../../tests/batch/test_analysis_api.py)：`test_single_analysis_does_not_require_batch_inventory`；[test_analysis.py](../../../tests/batch/test_analysis.py)：配对前提与归属检查 |
| 仅显式方面被调用；未实现、禁用、输入前提不满足分别记录 | [test_analysis_api.py](../../../tests/batch/test_analysis_api.py)：`test_empty_analysis_does_not_create_placeholder_task`、`test_disabled_analysis_is_inactive_without_running_adapter`、`test_removed_aspect_keeps_files_but_marks_control_inactive`；[test_analysis.py](../../../tests/batch/test_analysis.py)：`test_unimplemented_aspects_are_not_successes` |
| 适配器被调用并返回可重载文件；比较依据与实际参数保存 | [test_analysis_api.py](../../../tests/batch/test_analysis_api.py)：`test_single_summary_only_reports_requested_aspect`；[test_example_adapter.py](../../../tests/batch/test_example_adapter.py)：`test_example_adapter_records_actual_carriers`；[test_task_workflow.py](../../../tests/batch/test_task_workflow.py)：资源和重用记录 |
| sST 构造基线辅助可用；要求真实 cell 配对的模块仍不可用 | [test_inputs.py](../../../tests/batch/test_inputs.py)：`test_sst_baseline_equal_splits_raw_spots_in_output_order_and_preserves_coordinates`、`test_sst_baseline_rejects_native_pairing_request`；[test_analysis.py](../../../tests/batch/test_analysis.py)：`test_pairing_requirement_blocks_sst` |
| 单项失败后继续；部分输出不发布，发布失败恢复旧文件 | [test_analysis.py](../../../tests/batch/test_analysis.py)：`test_failed_analysis_does_not_prevent_other_samples`、`test_failure_continues_and_does_not_publish_partial_files`、`test_publication_failure_restores_previous_analysis` |
| 无变化复用；分析参数／资源／代码变化只重跑相应方面，不调用重建执行器 | [test_task_workflow.py](../../../tests/batch/test_task_workflow.py)：`test_resource_change_reruns_only_its_aspect`、`test_comment_only_config_edit_keeps_analysis_reusable`；[test_analysis_api.py](../../../tests/batch/test_analysis_api.py)：`test_resource_change_during_run_is_rejected_before_publish`、`test_code_change_during_run_is_rejected_before_publish`；[test_runner.py](../../../tests/batch/test_runner.py)：`test_analysis_configuration_does_not_invalidate_reconstruction` |
| 新增／删除／修改配置按有效相关字段失效；注释变化不重算 | [test_config.py](../../../tests/batch/test_config.py)：配置链变化；[test_runner.py](../../../tests/batch/test_runner.py)：`test_configuration_changes_follow_effective_settings`；[test_analysis.py](../../../tests/batch/test_analysis.py)：`test_analysis_rechecks_configuration_chain` |
| 输入或资源在运行中变化、产物被修改、旧指纹版本均不能复用；中断重跑 | [test_task_boundaries.py](../../../tests/batch/test_task_boundaries.py)：`test_output_root_change_during_execution_is_not_success`、`test_reference_change_during_handoff_is_not_success`；[test_runner.py](../../../tests/batch/test_runner.py)：`test_changed_input_during_reconstruction_is_not_published_as_success`；分析运行中变化由 `test_analysis_api.py` 覆盖 |

## Scientific regression boundary

保留 [test_real_sample_parity.py](../../../tests/integration/batch/test_real_sample_parity.py) 中 `test_real_batch_matches_direct_single_run` 的 hST／iST／sST 小样本对照；单任务 API、分析适配器交接和输入视图分别由上列 targeted tests 覆盖。使用同配置、受限样本和适当随机种子，不启动大规模计算；数值比较需注明容差与非确定性来源。

这些测试证明重建路由一致及辅助计算正确。测试 adapter 只证明框架真实调用、正确交接；EMT／CCI 等真实模块的科学有效性、统计前提和生物学解释需要其接入时另行验证，不能由本轮测试替代。

当前已有测试的旧预期若与新设计冲突，修改时须保留其原来防止的错误（如来源过期、项目串用），改用新边界验证；不简单删除测试来获得通过。

## Delivery checks

1. 两个单任务 API 和现有两个批量命令在安装环境可调用；同一小样本单项／批量示例中，任务选择不会扩大范围。
2. 示例 adapter 可执行、产物可重载，失败恢复示例可重复；配置示例中的资源路径按分层规则解析。
3. [protocol.md](protocol.md)、[framework.md](framework.md) 和 [README.md](README.md) 的字段、真实符号、复用边界和使用示例与当前实现一致；旧 adapter 返回形式及同进程源码缓存有明确迁移说明。
4. 尚未接入的科学模块保持未实现状态，重建完成与分析完成分别可判断；实际运行范围见下方交付验证，不将框架接入等同于科学验证。

## Verified delivery

以下安装和入口 smoke checks 已通过：

- clean `/tmp` source wheel `2.0.0rc1` 构建，并以 `--no-deps` 安装到目标环境；
- 两个已安装 console entrypoint 从不同工作目录执行 `--help`；
- 公共单任务 API 与 `AnalysisInputs` 相关类均可导入。

source-checkout wrapper 与已安装 wheel 入口保持分开。

2026-09-09 最终验收：

- `tests/batch tests/application`：**235 passed**，无失败或跳过。
- `tests/integration/batch/test_real_sample_parity.py`：**3 passed**，覆盖 hST、iST、sST 同配置单次／批量重建、单任务复用、实际 adapter 产物及输入视图。
- 真实输入来自 `raw_data/Real_application/` 的 P1CRC HD、P2CRC Xenium、共享 reference，以及 `raw_data/Sim2Real-ST-P5CRC/spot/part1/spot_50/xenium_spot.h5ad`；只读取受限切片，派生数据和结果写入测试临时目录。
- 运行环境为 Python 3.10，设置 `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1`、临时 `NUMBA_CACHE_DIR`／`MPLCONFIGDIR`，测试进程禁用本机有问题的 `readline` 扩展；在宿主环境运行科学子进程，避免沙箱 OpenMP 限制。两组分别有 17／15 条依赖及数值警告，不将其隐藏为无警告通过。
- `git diff --check`、`git diff --cached --check`、文档相对链接与章节引用检查通过。未运行整个仓库的无关测试，也未启动大规模重建或科学分析。

## Documentation checks

- 主计划独立覆盖目标、范围、关键决定、实施顺序和完成标准；[问题索引](plan.md#question-index)覆盖全部专题、现行规范及关键外部来源。
- 需要某项细节时，从主计划一次跳转可达；专题均返回主计划，必要前置章节直接链接。
- 规则详细定义只在对应专题，其他位置引用；统一规范按当前实现说明边界，不能把目标当作当前可用能力。
- 所有文档链接使用相对路径和稳定章节名；代码／测试链接指向真实文件，计划中的新增文件明确标为计划且不伪装成现有链接。
- 文档检查只证明文档交付；不将链接检查或计划中的测试清单报告成运行回归已经通过。

输入输出样式的唯一示例为 [input-output-example.md](input-output-example.md)。交付时检查其中的 SC 癌种路径与链接的 YAML 模板一致；目录样式不在各说明文档重复维护。
