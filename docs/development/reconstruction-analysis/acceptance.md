# 当前验收：2026-09-20 双仓 mini 真实联测

[协同入口](cross-repo-review.md) · [实施索引](implementation-index.md) · [比较协议](assembly-comparison-contract.md) · [机器记录](verification-mini-2026-09-20.json)

本轮完成代码修正、六次真实重建（三路线，其中 iST 四 assembly）、正式三文件交付及三个样本的 Notebook / batch / report。范围是原始样本约 1% 的连续空间正方形；完整 Raw 指该 mini 输入全部单位。全量与科学解释未验收。hST 表达来源未知、sST parent 校正缺口仍保留。

## 按本轮任务编号交付

编号对应 2026-09-20 用户工作包，不替代历史 R1–R8。

| 本轮编号 | 实际修改或复用 | 本轮证据与状态 |
|---|---|---|
| REVISE 1–2 | cell type `/ → _`，移除写法 collision 拒绝；类型参数去重；冲突计数先统一写法 | 原人工列、历史 cluster、NA 保留；真实 ID 唯一性、概率轴与发布保护仍在；回归通过 |
| REVISE 3–6 | 复用一次 GA、eligibility、失败传播、单一 SVC、完整 Raw 和发布事务；batch 默认 whole-sample random；显式表达/坐标声明 | 每个 iST run 一次 GA；四个 eligible 类型共 501 SVC；九个单 Level2 类型有原因地跳过。三路线 Raw 矩阵、轴、原标签、坐标与 mini 输入完全一致 |
| REVISE 7 | P2 random 主链；分析配置直接引用生产 sample.yaml | 3,409 Raw → 501 SVC；17 个 Notebook 代码单元执行；170 份科学表/JSON 与 batch 逐字节一致 |
| REVISE 8 / H2 | backed 提取三个历史 carrier 的 obs/spatial，构建零基因 baseline | 49,279 × 0；原 ID、各自 Level1、SVC_cluster、坐标逐值一致；来源 SHA256 保存 |
| REVISE 9–10 / H3–H4 | 四方法真实执行；显式 T、Mono_Macro、Fibroblast；逐 ID 严格坐标检查和共同 ID/gene 范围 | 四方法 upstream 配置、观测轴/标签/坐标一致，仅 assembly/output 不同。共同 ID 为 160/157/68，共同基因 13,088；坐标反例阻止对应类型评分 |
| REVISE 11 | comparison Notebook 补保存步骤并真实执行 | 36 个方法×类型×resolution 指标、36 份列联表、12 张图、coverage、参数、来源摘要及执行版；seed 42，resolution 0.6/0.7/0.8 |
| REVISE 12 | hST/sST 真实重建、交付、消费 | hST 5,077 Raw → 4,128 SVC；sST 43 spot → 567 SVC，parent/坐标通过；能力限制见下表 |
| Analysis 1–3 | 删除 label_aliases；类型与参数 `/ → _`；默认 All、Fibroblast、Mono_Macro、T | 主 SVC_cluster 身份保留；独立 Moran/pathway 沿用原生默认范围 |
| Analysis 4–5 | Raw Level2 使用 scope 内有效 ID及同 ID坐标，新增 coverage；K-control 缺主标签/窗口 unavailable | P2 All 有效 501/3409，排除 2908；目标类型分别 112/112、191/191、183/183；缺失/部分缺失/运行异常回归通过 |
| Analysis 6–8 | 六个 mini/full 项目配置；真实 Notebook/batch/report；报告补有效/排除数 | 三个 Notebook 各17个代码单元，无阶段异常；iST/hST/sST 分别170/11/44份科学表/JSON 与 batch一致；12份报告只读刷新，科学文件 hash 不变 |

## 真实结果与限制

| 路线 | 正式交付 | Impact | 独立 Moran / pathway | 解释边界 |
|---|---|---|---|---|
| iST | [P2 random](../../../output/mini-acceptance/20260920/delivery/P2CRC_Xenium_mini/random/sample.yaml)，Raw 3409×422，SVC 501×13088 | partial，无阶段异常 | partial / succeeded | 小 ROI 无稳定 State/Gain 阈值；QC 后部分 Raw broad 缺失使 Anatomy unavailable；常量基因 Moran 不可计算 |
| hST | [P1 HD](../../../output/mini-acceptance/20260920/delivery/P1CRC_HD_mini/default/sample.yaml)，Raw 5077×18085，SVC 4128×3451 | partial，无阶段异常 | skipped / skipped | 输入表达历史未确认，两侧 identity 保持 unknown，表达分支未验收；正常无 SVC_cluster |
| sST | [P2 Visium](../../../output/mini-acceptance/20260920/delivery/P2CRC_Visium_mini/default/sample.yaml)，Raw 43×18085，SVC 567×5021 | partial，无阶段异常 | partial / succeeded | 无 SVC_cluster 只限制相关分支；K-control 已开启验证；membership_same_units=false；0.73 µm/coordinate 暂定 |

sST parent 校正为 **partial**：139,457 个正目标 parent–gene 项中，821 项最终聚合为有效零（≤1e-12），涉及总目标表达量的 0.774%；其余138,636项在 rtol=1e-5、atol=1e-6 内一致，最大相对误差5.85e-7。单个 parent 最大总量差约1,953/10,000，不能只用全局比例淡化局部影响。现有 `X / (current_sum + 1e-10)` 乘法校正无法恢复无支持或极小支持的表达，本轮不改变算法。生成细胞总量约51.75–8625.29，中位数680.19，没有恢复最终逐细胞10,000缩放。[完整数值证据](../../../output/mini-acceptance/20260920/evidence/delivery-P2CRC_Visium-default.json)。

iST 坐标使用已确认0.2125 µm/coordinate；hST使用源 full-resolution metadata 的0.27380817798463214 µm/pixel，观测网格步长约8.00145/8.00017 µm，与8 µm bin一致。ROI中心、边界、原ID、来源摘要及数量见 [input_manifest](../../../output/mini-acceptance/20260920/input_manifest.json)。未平移坐标或改变 QC/reference 范围。

## 证据与测试

- [真实比较 Notebook](../../../output/mini-acceptance/20260920/assembly/comparison/assembly_comparison.executed.ipynb) · [HTML](../../../output/mini-acceptance/20260920/assembly/comparison/assembly_comparison.html) · [coverage](../../../output/mini-acceptance/20260920/assembly/comparison/coverage.csv) · [metrics](../../../output/mini-acceptance/20260920/assembly/comparison/metrics.csv)。ARI/NMI 对照指定历史重建标签，不是独立生物学真值，不自动选默认 assembly。
- 消费者 [iST Notebook](../../../../REVISE_Analysis_Agent/output/mini-acceptance/20260920/iST/notebook/01_reconstruction_impact.executed.ipynb) · [iST 报告](../../../../REVISE_Analysis_Agent/output/mini-acceptance/20260920/iST/notebook/report.html) · [hST 报告](../../../../REVISE_Analysis_Agent/output/mini-acceptance/20260920/hST/notebook/report.html) · [sST 报告](../../../../REVISE_Analysis_Agent/output/mini-acceptance/20260920/sST/notebook/report.html)。
- [Notebook/batch parity](../../../output/mini-acceptance/20260920/evidence/notebook-batch-parity.json) · [报告只读刷新](../../../output/mini-acceptance/20260920/evidence/report-refresh.json) · [baseline/四方法身份](../../../output/mini-acceptance/20260920/evidence/baseline-and-method-identities.json) · [四方法配置](../../../output/mini-acceptance/20260920/evidence/four-method-configurations.json)。
- 生产修改相关141 passed、额外调用链73 passed；消费者全套109 passed；比较与baseline套件最终13 passed（含fixture Notebook）；root联合交付/mini检查12 passed，随后sST校验补充3 passed。套件有重叠，不相加。历史759/964/75不是本轮通过数；当前命令与来源见机器记录。

`output/` 为本地运行证据，Git 忽略的大数据需保留或重建，缺文件不能沿用结论。首次 comparison kernel 被 socket 沙箱限制，失败记录已保留；最终执行版源码哈希不变、无cell error。普通启动出现过 OpenMP/numba 环境错误，隔离后重试成功，失败启动不作为科学通过证据。

## 服务器全量入口

两个仓库并列放置。在 REVISE 根目录激活包含 TACCO、POT、Scanpy、AnnData 的科学 Python 环境；Notebook python3 kernel 必须指向同一环境。输入根为 `raw_data/Real_application/`，需要 P2CRC_Xenium.h5ad、P1CRC_HD.h5ad、P2CRC_Visium.h5ad、adata_sc_all_reanno.h5ad。P2 reference 明确过滤 Patient=P2CRC；hST保持原reference范围。

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMBA_NUM_THREADS=1
export PYTHONHASHSEED=42 MPLBACKEND=Agg
python reconstruct.py --config configs/acceptance/full/P2CRC_Xenium-random.yaml
python reconstruct.py --config configs/acceptance/full/P2CRC_Xenium-mean.yaml
python reconstruct.py --config configs/acceptance/full/P2CRC_Xenium-within_cluster.yaml
python reconstruct.py --config configs/acceptance/full/P2CRC_Xenium-outside_cluster.yaml
python reconstruct.py --config configs/acceptance/full/P1CRC_HD-default.yaml
python reconstruct.py --config configs/acceptance/full/P2CRC_Visium-default.yaml
```

输出为 `results/acceptance-full/<sample>/<mapping>/`。四iST配置除assembly/output外相同，各自运行主流程，没有新增GA缓存。OT保留 `max_cost_entries=2_000_000`；超限应报告失败，不自动增大限制、抽样reference或降级assembly。完整iST稠密输出估计下界约23.6 GiB，hST相关全矩阵约34 GiB，尚不含副本、图和参考数据；建议独占至少64–128 GiB RAM节点并监视实际峰值。这是资源估算，未经全量验证。

在 Analysis Agent 根目录运行：

```bash
python -m revise_analysis.cli batch --config configs/p2_full_project.yaml
python -m revise_analysis.cli batch --config configs/p1_hd_full_project.yaml
python -m revise_analysis.cli batch --config configs/p2_visium_full_project.yaml
python ../REVISE/scripts/execute_acceptance_notebook.py --repo . --source notebooks/01_reconstruction_impact.ipynb --project configs/p2_full_project.yaml --output output/acceptance-full/iST/notebook
```

其余Notebook换相应项目配置和新输出目录。batch CLI在partial/skipped时返回1，需读batch_result.json与stage_errors区分前提不足和失败。全量comparison复制real JSON，显式替换四个method_paths、baseline_path、output_dir，使用执行脚本的 `--comparison-config` 指定，不覆盖mini证据。

mini使用 `configs/acceptance/mini/*.yaml` 和消费端 `configs/p2_project.yaml`、`p1_hd_project.yaml`、`p2_visium_project.yaml`。ROI准备脚本 `scripts/prepare_mini_acceptance.py` 要求全新的output/configs目录，拒绝覆盖已有输入和配置。hST来源确认、sST零支持校正策略、全量容量及科学解释为剩余事项。
