# P2 Xenium：真实重建—分析—审视报告

Task ID：`p2-xenium-impact-20260922`。状态：`scientific_review_complete_with_limits`。本状态不等于方法生物学有效性已获验证。

## 主要判断

全量真实输入已完成 random 重建、三文件交付、原生 Impact、K-control 全量 Impact、共同有序基因背景诊断和 donor-aware 审计。输入、运行、抽样对象和产物均可追溯。所有正式进程正常退出，无新增 OOM；原生 Impact 的 `partial` 仅指部分 SVC 常量基因无法计算 Moran，未被补零，也没有阶段错误。

**本轮没有足够证据支持“重建整体更好”。** 已验证的变化是：标签多样性差对 Raw 分群粒度敏感；共同基因空间自相关大多降低，统一基因背景后方向保留；原生及共同背景 EMT 分数都不足以支持跨侧定量活性改善。这些是有效的审视结果，不通过调整参数追求正面结论。

从方法机制看，random 在同群内选用参考细胞表达载体，所有 SVC 表达行均与来源 donor 精确一致。它扩展了表达基因范围，但没有增加独立原位测量。Moran 下降与这种信息转移及较低共同轴检测支持相容；本轮未分离随机赋值、参考选择、跨 assay 稀疏性等因素的各自作用，不能把它直接称为生物学恶化或算法造成的空间损伤。

范围必须明确：Raw 保留 340,837 个单位，SVC 为 50,626 个（14.85%），包含 DC、Fibroblast、Mono_Macro、T；SVC 不含 Tumor/Intestinal epithelial。三个 parent 的完整单位 ID 集合一致，All 的范围与抽样口径不同。因此 parent 是较清楚的局部比较，All 不能代表全组织重建改善。

## 统一问题与判断表

| 问题／观察 | 前因与机制 | 证据与比较口径 | 已采取措施 | 当前判断 | 对主结论的影响 | 需要用户判断或补充 |
|---|---|---|---|---|---|---|
| 输入、坐标及表达来源是否正确 | 原位 Raw 与参考表达分别保留；random 按同群分配 donor | 输入双端哈希；Raw X/轴/坐标/原标签保护；50,626 assignment 与原始参考逐行精确核验；两源 stored X 全为有限、非负整数 | 修复交付审计中错误的全基因覆盖要求、补充原始 donor 来源审计及运行回执 | 身份及技术链路可信 | 证明来源正确，不证明恢复了空间真值 | 无当前决策 |
| 重建覆盖只占 Raw 的一部分 | 98,543 单位未过预处理；191,668 候选因参考只有一个有效子类型而跳过 | 242,294 个候选最终交付 50,626；逐大类排除有底表 | 保留完整 Raw；不虚造子类型或悄悄补入被跳过对象 | 覆盖受限、原因可追溯 | 不可将 All 差异外推全组织；肿瘤/上皮机制不在本轮 SVC 证据范围 | 若后续目标必须包含这些对象，需要决定参考或方法语义；本轮不自动更改 |
| 正 Gain 是否意味着改善 | Raw Leiden 与 SVC 标签的粒度和群大小不同 | 40µm、同一共同有效窗口场；预设 grid 的 K-control 结果见下表 | 完整重跑 iter-002，仅开启 raw_k_control；保存全部 grid/selection/正负 delta | 对粒度敏感；只有 Mono 精确匹配 K，其他是 nearest | 不支持稳健的整体多样性改善；不能把 nearest 当等 K 因果控制 | 精确匹配其他 scope 是否值得追加属后续问题，本轮不调 grid 追结果 |
| Membership 变化比例能否作改善指标 | All 比较是 Raw 30k 抽样与 SVC 覆盖交集；两种标签系统不等价 | All 仅 4,506 单位（Raw样本15.02%、Raw全体1.322%、SVC全体8.901%），交集中 Raw6群/SVC31群；85.8%变化 | 前置共同 cohort 与实际群数；核对 K-control 不生成受控 membership | 仅原生标签对应关系描述 | 不是全组织重建改变率，更不是改善率 | 无需将其升级为核心改善证据 |
| Moran 大多下降是否只是基因背景造成 | 原生422/13088基因轴不同；同时存在参考载体和检测支持差异 | 同343基因轴、保留原 cohort/k6/方法后，parent负差307/333、315/338、305/339 | 独立共同有序轴诊断；逐基因及排除原因完整保留 | 下降方向不是仅由基因背景差异造成 | 支持描述方向，不授予生物优劣或单一原因；父群同ID仍非同一表达测量 | 分离随机赋值与跨assay稀疏性的贡献需新设计；独立空间真值尚缺 |
| EMT 分数变化是否可信地反映活性变化 | 基因覆盖、排名背景、cutoff和评分单位均有差异 | 原生31/200 vs187/200、cutoff3/250；共同轴29/200 vs29/200，cutoff仍2/6；Raw共同轴仅2种分值 | 保存两种背景完整单位评分、detection、provider/adapter/seed身份；不改cutoff“救活”比较 | 跨侧定量活性结论不可授予 | 不能把更高分称EMT增强或改善；本轮也没有SVC上皮/Tumor支持转化判断 | 若需要此类主张，需明确对象与可比评分/独立验证设计 |
| donor 复用是否抬高平均 EMT 分 | 50,626 assignment由10,021参考细胞产生；不同大类复用权重不同 | 已评分30,000单位/8,996参考细胞；均值assignment .0556223，保留大类权重的donor等权 .0575542，全donor等权 .0439714 | 完整独立标签/ID/cluster/score联接与两层权重分解；Kish仅作权重描述 | 本次保留大类占比口径下未见复用抬高All均值 | 不等于reference选择无偏，也不把参考细胞当患者或独立生物重复 | 无独立真值时保留信息转移限制 |
| Other/Unknown 是否具有生物组织意义 | Anatomy由交付SVC派生；没有肿瘤/上皮支持 | SVC50,626全Other；Raw231,604 Other/109,233 Unknown；Unknown全部等于未被派生背景覆盖 | 核完整点表、覆盖标记和条件表 | 组织分层不可用；Unknown是覆盖不足 | 不支持肿瘤/正常/界面定位或独立组织真实性证明 | 需要独立组织信息才能检验这些主张 |
| State/Gain 二元区域是否稳健 | min4稀释引入离散Neff原子点；阈值区间可跨越大量同值窗口 | T State窗口在保存阈值区间端点为296/168、Gain253/49；Mono/T阈值附近窗口均为4单位 | 只读完整连续场与保存bootstrap阈值；不改阈值、不补零 | 精确边界敏感；阈值ok不等于掩膜精确 | 区域保留描述性使用；端点计数不是区域数量总体置信区间 | 不以这些掩膜承担稳健机制结论；若未来需要再设计尺度/空间稳健性验证 |
| 环境、回执及报告是否可靠 | 8GiB实际容器配额；原AUCell多rank副本；服务器中文字体缺失 | 真实adapter严格数值等价；各运行guard/receipt/输入及产物哈希；科学底表与展示副本一致 | 两仓独立环境、已验证低内存rank缓冲、受控运行；从保存表刷新中文图并前置限制 | 已修复影响本次执行及阅读的问题 | 技术成功与科学结论分开；未重复展示Notebook全量计算 | GitHub推送未进行，见交付说明；本地与启智执行不受影响 |

## K-control：灵敏度，非一概等群数

共同窗口集合与原生比较完全一致。预设 resolution grid 为0.2/0.4/0.6/0.8/1.0，未为匹配或结果方向追加参数。下列均值/中位数为完整正负 delta 场，不是仅正Gain候选。

| Scope | 目标K / 实际K | 状态 | 共同窗口数 | 原生均值 → 诊断均值 | 原生中位数 → 诊断中位数 |
|---|---:|---|---:|---:|---:|
| All | 31 / 15 | nearest | 1417 | 0.6052 → 0.0201 | 0.6518 → 0.0161 |
| Fibroblast | 8 / 6 | nearest | 1547 | 0.1055 → −0.1548 | 0.0159 → 0 |
| Mono_Macro | 8 / 8 | exact | 1635 | −0.2683 → −0.3466 | −0.1690 → −0.2888 |
| T | 8 / 9 | nearest | 1265 | 0.0684 → −0.2619 | 0 → −0.2452 |

这足以表明原生正向结论依赖标签粒度。只有Mono可称精确K匹配，其他不能据此定量分解“多少差异纯由K造成”。All还保留范围/抽样混杂，不能升级为全组织结论。现有实现不产生K-controlled membership，也未保存control逐单位标签及群大小，不能追溯具体群大小如何贡献差值；该限制不阻断对预设分区分辨率敏感的保守结论。即使K相等，也不意味着生物学粒度相同。

## Moran：控制基因背景后方向保留

| Scope | 共同可计算基因 | 共同轴负差数 | 原生差中位数 | 共同轴差中位数 |
|---|---:|---:|---:|---:|
| All | 343 | 298 | −0.056417 | −0.056678 |
| Fibroblast | 333 | 307 | −0.032709 | −0.032896 |
| Mono_Macro | 338 | 315 | −0.038609 | −0.038049 |
| T | 339 | 305 | −0.023295 | −0.023660 |

三个parent使用相同完整单位ID，All仍是各侧独立30k。常量基因不进入可计算分母。共同轴的count-like总量/检测基因数中位数仍不同：Fib Raw112/46 vs SVC35/18，Mono96/50 vs56/23，T84/51 vs37/21。它们来自不同assay及表达载体，不能称可直接比较的真实测序深度，也不能认为统一基因轴已消除所有混杂。

## 验收与保留未知

运行完成不等于证明方法更好。主代理与独立审视已核对最终代码和完整底表：当前没有仍影响上述保守结论、且必须继续解决的重大问题。本次必须完成的执行、必要诊断和审视已达到；生物学改善主张不获接受。保留未知包括：未覆盖肿瘤/上皮的后续目标、精确K匹配其他scope的意义、reference/random与稀疏性的因果分解、独立空间真值及生物重复。这些未知限制改善主张；不通过参数调整或换分母掩盖它们。

## 运行与证据入口

| 运行 | 代码提交 | 实例 / QZ job | 结果 |
|---|---|---|---|
| 全量random重建 | `0abf8515b2d80340e33b530719878a95182b2975` | cpu / `20260923-012058-l2magz` | 三文件及来源审计通过 |
| 原生Impact iter-001 | `c9e4b7454d6f42708990690eaf1ab34867f4515a` | cpu-4 / `20260923-023200-biy79t` | exit0，322登记产物，partial仅常量基因 |
| K-control完整Impact iter-002 | `bf3fe5f065d291b9f03259185e7ab40db79d5b6a` | cpu-4 / `20260923-041123-s9hvwo` | exit0，342登记产物，原生189表仅配置开关记录变化 |
| 共同ordered-gene背景 | `bf3fe5f065d291b9f03259185e7ab40db79d5b6a` | cpu / `20260923-041743-01rkwg` | exit0，29登记产物，cohort身份保持 |

两仓分别使用其远端`.venv/bin/python`；启智代码根为`/inspire/hdd/project/sais-phenotype-fm/public/wuyushuai/zz/0_REVISE_developing/REVISE_code`。分支REVISE=`revise-2.0`，Analysis=`sl`。实际解释器版本、依赖、命令、输入完整SHA、耗时和资源以以下不可变receipt为准。正式输出未被展示修复覆盖。

- [生产运行与完整donor审计](../../../../results/p2-xenium-impact-20260922/iter-001/receipt.json)
- [原生运行回执](../../../../../REVISE_Analysis_Agent/output/p2-xenium-impact-20260922/iter-001/_retrieved_final/20260922T194202Z/extracted/iter-001/_run/run_receipt.json)
- [K-control运行回执](../../../../../REVISE_Analysis_Agent/output/p2-xenium-impact-20260922/iter-002/_retrieved_final/20260922T203950Z/extracted/iter-002/_run/run_receipt.json)
- [共同背景运行回执](../../../../../REVISE_Analysis_Agent/output/p2-xenium-impact-20260922/iter-001/diagnostics/_retrieved/common-gene-background-001-full-evidence/common-gene-background-001/run_receipt.json)
- [原生中文审阅报告](../../../../../REVISE_Analysis_Agent/output/p2-xenium-impact-20260922/iter-001/review/native-display-001/report.html)
- [K-control中文审阅报告](../../../../../REVISE_Analysis_Agent/output/p2-xenium-impact-20260922/iter-002/review/k-control-display-001/report.html)
- [逐项发现和修复记录](findings.md)
- [生产交付审视](producer-review.md)
- [原生与诊断详细解读](early-impact-review.md)
- [K灵敏度完整摘要](../../../../../REVISE_Analysis_Agent/output/p2-xenium-impact-20260922/iter-002/review/main-k-control-001/k-sensitivity-summary.csv)
- [共同背景完整底表](../../../../../REVISE_Analysis_Agent/output/p2-xenium-impact-20260922/iter-001/diagnostics/_retrieved/common-gene-background-001-full-evidence/common-gene-background-001/artifact_hashes.json)
- [donor-aware完整审计回执](../../../../../REVISE_Analysis_Agent/output/p2-xenium-impact-20260922/iter-001/review/donor-aware-emt-native-002/receipt.json)
- [逐单位表达支持审计](../../../../../REVISE_Analysis_Agent/output/p2-xenium-impact-20260922/iter-001/review/depth-audit-001/receipt.json)

GitHub推送重试被自动审批拒绝，原因是提交含内部服务器路径及运行/资源元数据，未建立向该外部目的地发送这些内容的明确授权；未重试。已用本地提交与定向Git bundle同步启智，不影响本次计算或证据交付。

最终接手入口：[主代理关闭回执](../../../../../REVISE_Analysis_Agent/output/p2-xenium-impact-20260922/review-closeout/receipt.json)、[独立审视记录](../../../../../REVISE_Analysis_Agent/output/p2-xenium-impact-20260922/review-closeout/independent-review.md)。最终代码、完整底表与本报告/两份审阅HTML已经过独立复核，无未关闭重大项。
