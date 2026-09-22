# 首轮真实交付：主代理与独立审视记录

状态：生产交付 `technical_verified`；下游分析与最终科学验收仍待完成。依据为本地已取回、与远程哈希一致的 [交付审计](../../../../results/p2-xenium-impact-20260922/iter-001/audit/delivery-audit.json)、[运行回执](../../../../results/p2-xenium-impact-20260922/iter-001/receipt.json) 和 [容量证据](../../../../results/p2-xenium-impact-20260922/iter-001/capacity/probe-002.json)。

正式代码 `revise-2.0@0abf8515b2d80340e33b530719878a95182b2975`，QZ job `20260923-012058-l2magz`，cpu；解释器为该仓 `.venv/bin/python`。正式运行回执时间为 2026-09-22 17:21:12 至 17:41:49 UTC，重建及审计均退出 0。输入、命令、依赖记录和完整哈希见回执。进程观测峰值 RSS 3,569,844,224 B，cgroup working set 峰值 7,394,627,584 B；90% 保护未触发，memory.events max 保持 80，OOM/kill 均为 0。cgroup 包含同容器其他基础进程，不等同于本任务 RSS。

## 已核实事实及适用范围

| 观察 | 证据与分母 | 当前判断与限制 |
|---|---|---|
| Raw 原对象保护 | 340,837 × 422；原 X、轴、已有标签、坐标逐项核验 | 原始输入未被 SVC 子集覆盖；新增 `revise_*` 标签须另追溯 |
| 预处理保留 | 242,294 个空间细胞；30,774 个参考细胞全部为 P2CRC | 98,543 个 Raw 单位未进入推断流程；不使用原标签静默填补新增推断标签缺失 |
| 最终 SVC 覆盖 | 50,626 × 13,088；DC 1,347，Fibroblast 17,455，Mono_Macro 17,605，T 14,219 | 约占 Raw 14.85%、预处理后单位 20.90%。full 输入执行不等于全组织输出覆盖 |
| 10 类跳过 | B、Gliacyte、Intestinal Epithelial、Lymphatic EC、Mast、Pericyte、Plasma、SMC、Tumor、Vascular EC；每类仅一个有效参考 subtype | 明确的 full-route 方法边界，实际没有进入 random 载体。不能代理自行增加 singleton passthrough 或伪造 subtype |
| 三个 named parent | 两侧推断 broad 计数一致；候选与最终 SVC 的各 parent ID 集合哈希一致，顺序不保证一致 | 可以继续 parent 内描述与必要敏感性分析；parent 来自同一 Global Anchoring，不能称独立标签验证 |
| Raw Level2 来源 | `delivery._inferred_labels` 从同一 Local Refinement 的 `sc_svc_spatial` 回填 | 只能作为来源明确的同流程 comparator。Raw Leiden 使用 Raw X，但 parent 选择仍共享 GA |
| All 比较 | Raw 340,837 与 SVC 50,626，覆盖和类型组成明显不同 | 保持既定独立对象，不自动过滤 Raw。All 差值不能单独归因为重建改善 |
| Anatomy 来源缺失 | SVC 无 Tumor、无 Intestinal Epithelial | 预期 SVC-defined Anatomy 缺乏 Tumor/Normal/Interface 区分力；实际点表仍待 Impact 核验。不得把 Other 解释为原组织不存在这些类型 |
| 基因身份 | Raw/SVC 共同 343；SVC 13,088 全部在参考；Raw/参考共同 407 | 原生基因全集与深度不同；跨侧 Moran 仅按共同基因证据描述。参考扩展不等于空间真值 |
| random 信息转移 | 50,626 次分配全部逐行与来源 donor 表达精确一致；10,021 个独立 donor，复用最大 34 | 80.21% 分配为超过独立 donor 数的重复分配；96.99% 单位落在复用 donor 上。表达变化不能称新增独立信息或 Xenium 真值验证 |

## 后续验收

独立审视与主代理均同意交付可进入既定 Impact，不能据此授予科学接受。完整覆盖链及 ID 集合补证已完成，见 `audit/coverage-identity-supplement.json`、`coverage-by-broad.csv` 和 `broad-label-transitions.csv`。原标签与 GA 的 1,805 处冲突中，996 个最终四类单位的原标签不同，其中包含四类内部改判；另有 500 个原四类单位被改判至这四类之外，360 个原其他类单位迁入四类；不能只用整体冲突比例掩盖局部敏感性。若关键结论异常或接近判断边界，再检查原标签口径敏感性。下游须检查完整底表、真实分母、原生评分背景与 donor 复用；需要时执行预先说明问题的独立诊断。Anatomy 若实际全 Other，报告应明确其输入来源限制，而不是以计算成功授予组织区分结论。

报告中“数量差不表示单位丢失”的笼统说法已改为结合输入筛选与交付覆盖解释，并补充 Raw Level2 与 parent 来源的独立性限制；这属于展示层纠错，不改变科学计算。
