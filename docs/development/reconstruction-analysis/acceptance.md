# 当前验收：2026-09-20 双仓收口

[协同入口](cross-repo-review.md) · [实施索引](implementation-index.md) · [比较协议](assembly-comparison-contract.md) · [机器记录（当前见 closeout）](verification-mini-2026-09-20.json)

已完成 SVC-defined Anatomy、正式 Notebook 默认入口及 sST 正支持数值修正，并使用原固定 mini 输入完成真实联测。三路线共用分析 workflow，能力由字段与表达声明决定。实际重建完成的类型合并为单一 SVC；未重建类型不并入。工程收口完成；全量运行与科学解释未验收，sST 真零支持仍为 partial。

本轮基线 HEAD：REVISE `464ee26ac22978c3994e40af50325413e70cef13`，Analysis Agent `3e70ab364f89f9e24e981e3e9239e7fe563138a6`。本轮修改已提交并推送：REVISE `467af8b9a104328fba8db0e786471620eab4150f`，Analysis Agent `4c4278b3adde1907fd66474ff97957a9c7dab790`。用户审阅材料内容保留。旧 `output/mini-acceptance/20260920` 保留；新结果在两仓的 `output/mini-acceptance/20260920-closeout`。机器记录的原有字段保存上一阶段证据，`closeout` 是当前收口记录。

## A. 实际修改

| 仓库 | 文件／符号 | 修改内容 | 原因 |
|---|---|---|---|
| Analysis Agent | `ImpactWorkflow.stage_support` 及图／报告 | 完整 SVC broad 与坐标定义 Anatomy；沿用 Raw 原点；产物使用 `svc_anatomy_*` | 与重建后统一对象一致，Raw 推断标签缺失不再阻断 Anatomy |
| Analysis Agent | 连续 Impact Notebook、README | 默认正式 P2 project；保留显式 sample/project/override；按字段声明能力 | 日常入口与已验收正式交付一致 |
| REVISE | `_correct_sst_parent_gene_expression` | float64 先算 parent–gene 内份额，再乘目标；严格零支持不补值 | 去除固定 epsilon 缩放偏差并避免极小正数比例溢出 |
| REVISE | SVC `uns.sst_parent_gene_correction`、mini verifier | 保存校正前支持统计及稀疏零支持索引；分开核验正支持残差和 unresolved | 最终有效零不能证明校正前真零；旧文件无诊断为 unobserved |
| 双仓 | mini 配置、现有验收入口 | 新 sST 输出目录，三项目的新分析输出目录；科学参数不变 | 保留旧结果并明确本次证据对象 |

Anatomy 的 `Other` 仅表示已交付 SVC 窗口未观察到配置的 Tumor/Normal，不代表原组织没有这些类型。Raw/SVC 各自按坐标映射到此背景；没有 SVC Anatomy 网格的点为 `Unknown`。Raw 原始标签不改写，也不被这个背景重新注释。主 cluster 缺失仍只限制依赖它的分支。

## B. 删除或简化

| 原逻辑 | 当前处理 | 科学定义影响 |
|---|---|---|
| Anatomy 绑定完整 Raw 推断 broad | 直接使用完整交付 SVC 的 broad 与坐标；不新增来源选择／fallback | 按用户确认，背景改为 SVC-defined；不再称完整 Raw tissue |
| Notebook 默认为旧临时 carrier | 无显式路径时使用正式 P2 project；显式 sample 单独使用时不套用默认项目 | 不改变分析参数优先级或算法 |
| `X / (current_sum + 1e-10)` | 正支持采用 share-then-target，真零保持 unresolved | 恢复既定 parent–gene 目标；不定义新表达分配方法 |
| 以最终 `≤1e-12` 推断 zero support | 正式读取校正前诊断；旧文件只描述 effective-zero | 明确证据边界，不反推未观测状态 |

未新增 registry、planner、缓存、兼容框架、补标签或补表达策略。内部 normalization、OT、quota、seed、LR eligibility、发布事务、单一 SVC 和完整 Raw 保护继续复用。

## C. 实际验证

| 场景 | 实际检查 | 结果 | 能证明／不能证明 |
|---|---|---|---|
| sST 数值及 verifier | [28 项回归日志](../../../output/mini-acceptance/20260920-closeout/logs/sst-regression-tests.log)、[JUnit](../../../output/mini-acceptance/20260920-closeout/evidence/sst-regression.xml) | 28 passed | 普通／极小／真零支持、无最终逐细胞缩放、诊断保存、严格零输出及旧产物 unobserved 退出码；不证明真零补偿方法 |
| 消费者 Anatomy／交互／报告 | [聚焦日志](../../../output/mini-acceptance/20260920-closeout/logs/tests-consumer-focused-53passed.log) | 53 passed | 缺 Raw broad、无主 cluster、Unknown、网格不变及产物命名 |
| 正式发布器 → 真实消费者 | [联合检查](../../../output/mini-acceptance/20260920-closeout/evidence/joint-delivery-tests.xml) | 2 passed | 三文件直接加载，SVC Anatomy 与 Raw Unknown，输入保护 |
| Notebook 配置入口 | [5 个实际场景](../../../output/mini-acceptance/20260920-closeout/evidence/notebook-entry.json) | passed | 默认正式项目、单独 sample、显式 project、override、错误归属拒绝 |
| 固定 sST mini | [交付检查](../../../output/mini-acceptance/20260920-closeout/evidence/delivery-P2CRC_Visium-default.json)、[前后对照](../../../output/mini-acceptance/20260920-closeout/evidence/sst-before-after.json) | 发布成功，数值 partial | Raw 43×18085，SVC 567×5021；obs、基因轴、parent 坐标与旧结果一致；配置只变输出目录 |
| Notebook／batch | [逐文件对照](../../../output/mini-acceptance/20260920-closeout/evidence/notebook-batch-parity.json) | 三路各17个代码单元通过；0 stage errors | 参数与172／14／45份科学表JSON逐字节一致；iST 44份 State/Gain／窗口表与旧版一致 |
| 报告 | [只读刷新](../../../output/mini-acceptance/20260920-closeout/evidence/report-refresh.json) | 12份 | result、表、图哈希不变；只刷新读取已保存数据的HTML |
| 输入与旧证据保护 | [最终哈希](../../../output/mini-acceptance/20260920-closeout/evidence/final-input-hashes.json)、[旧证据保护](../../../output/mini-acceptance/20260920-closeout/evidence/previous-evidence-preserved.json) | 25个输入／旧交付、71份旧产物未变 | 旧证据与用户审阅保留；不代表全量容量通过 |

测试数分套件记录，聚焦子集不与完整套件相加。隔离线程和插件、禁用 readline 的科学环境用于有效验证；初期代理的 pytest capture 导入崩溃不计为通过。

### 本轮真实结果

| 路线 | 正式交付 | Anatomy | Raw 映射无覆盖 | Impact／独立 Moran／pathway |
|---|---|---|---:|---|
| iST | [复用 P2 random](../../../output/mini-acceptance/20260920/delivery/P2CRC_Xenium_mini/random/sample.yaml)，3409 → 501 | 96个 Other | 401 | partial／partial／succeeded |
| hST | [复用 P1 HD](../../../output/mini-acceptance/20260920/delivery/P1CRC_HD_mini/default/sample.yaml)，5077 → 4128 | 214窗：Other189、Tumor22、Normal2、Interface1 | 275 | partial／skipped／skipped |
| sST | [新 P2 Visium](../../../output/mini-acceptance/20260920-closeout/delivery/P2CRC_Visium_mini/default/sample.yaml)，43 → 567 | 1个 Interface | 0 | partial／partial／succeeded |

三路均无阶段异常，SVC Anatomy窗口/context均生成；sST的Anatomy stage因部分Raw parent scope缺失为partial。iST仍受小ROI稳定阈值支持限制；hST表达身份unknown；hST/sST正常无SVC_cluster。sST单窗结论依赖当前mini及暂定尺度，不推广到全量组织。

sST 139457 个正目标项中，139207 个正支持项全部在 rtol=1e-5、atol=1e-6 内守恒，最大相对误差约1.07e-7。旧821个有效零项中571个恢复；剩余250个是本轮直接观测的校正前严格零支持，涉及2个parent，未解决目标量1086.313，约占总量0.253%，最严重parent约缺9.90%。生成细胞总量范围约51.75–9903.48，中位数680.19；没有最终逐细胞10,000缩放。目标仍是共同基因轴上的内部 parent normalize_total(10000)，不是恢复原始counts。

新消费者 [iST Notebook](../../../../REVISE_Analysis_Agent/output/mini-acceptance/20260920-closeout/iST/notebook/01_reconstruction_impact.executed.ipynb) · [iST报告](../../../../REVISE_Analysis_Agent/output/mini-acceptance/20260920-closeout/iST/notebook/report.html) · [hST报告](../../../../REVISE_Analysis_Agent/output/mini-acceptance/20260920-closeout/hST/notebook/report.html) · [sST报告](../../../../REVISE_Analysis_Agent/output/mini-acceptance/20260920-closeout/sST/notebook/report.html)。完整路径与运行日志见机器记录。

## 上一阶段完成项与复用边界

名称统一与NA保护、collision拒绝删除、一次GA、eligibility与失败传播、单一SVC、完整Raw、正式sample.yaml、Raw Level2按scope有效ID、K-control局部unavailable均保留。本轮没有重做这些实现。

H2 baseline为49279×0，历史ID/Level1/SVC_cluster/坐标保留。H3四方法同ROI/reference/QC/GA/LR/seed；H4逐ID坐标及共同基因检查，T／Mono_Macro／Fibroblast共同ID为160／157／68，共同基因为13088。[比较Notebook](../../../output/mini-acceptance/20260920/assembly/comparison/assembly_comparison.executed.ipynb)保存36指标、36列联表、12图；相关源码和产物未变，本轮复用，不重新运行。历史cluster是指定比较基准，不是独立生物学真值。

## D. 剩余问题

| 问题 | 类型 | 当前处理 | 后续证据 |
|---|---|---|---|
| 250项校正前真零支持 | scientific／method | 保留unresolved，不新增分配策略 | 明确补偿语义及跨样本规模后单独讨论 |
| P1 HD表达来源 | input fact | Raw/SVC保持unknown | 可追溯的原始表达和预处理记录 |
| Visium 0.73 µm/coordinate | input fact | 保持暂定 | 对应样本的坐标标定来源 |
| 全量容量与科学解释 | scale／scientific | 本轮未运行、未验收 | 服务器实际峰值／失败记录和科学审阅 |
| 共坐标、Raw/SVC范围差异及四assembly选择 | scientific | 明示单位与覆盖，不自动给结论 | 更大范围真实数据与独立解释 |

正式三文件主链仍成立。本轮修正的是 Anatomy 来源与日常入口的不一致，以及确定的 sST 数值偏差；没有新建框架掩盖能力缺失或真零支持。

## 服务器全量入口

手动传输的六个输入、校验和及两仓库放置说明见[远程输入准备](remote-inputs.md)。四方法全量比较使用 `examples/assembly-comparison-full.json`，原 mini JSON 保持不变。

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

其余Notebook换相应项目配置和新输出目录。batch CLI在partial/skipped时返回1，需读batch_result.json与stage_errors区分前提不足和失败。全量comparison使用 `examples/assembly-comparison-full.json`，通过执行脚本的 `--comparison-config` 指定，不覆盖mini证据；完整命令见远程输入准备。

mini使用 `configs/acceptance/mini/*.yaml` 和消费端 `configs/p2_project.yaml`、`p1_hd_project.yaml`、`p2_visium_project.yaml`。ROI准备脚本 `scripts/prepare_mini_acceptance.py` 要求全新的output/configs目录，拒绝覆盖已有输入和配置。hST来源确认、sST零支持校正策略、全量容量及科学解释为剩余事项。
