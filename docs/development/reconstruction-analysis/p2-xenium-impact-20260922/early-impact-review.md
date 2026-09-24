# 首轮 Impact 底表审视

状态：`scientific_review_complete_with_limits`。以下不构成整体科学验收。来源为 Analysis Agent 的两个独立只读快照：`iter-001/_running_snapshot/20260922T190702Z`（35个早期摘要/manifest）与 `20260922T191731Z`（24个完整Moran/区域/阈值CSV）。快照保存了staging来源、QZ job及归档哈希；最终发布后主代理已核对59项同名文件：58项一致，仅cohort_manifest补齐pathway记录；本文数值对应的科学底表均一致。

主代理直接读取完整小底表与实现；独立审视另行核对。下列诊断保留首轮结果，不以得到正面方向为目的。

| Scope | Raw/SVC K | 40µm共同有效窗口 / SVC有效窗口 | 正 / 零 / 负 Gain窗口 | Gain中位数 | Moran负差基因 / 共同可计算基因 | Moran差值中位数 |
|---|---:|---:|---:|---:|---:|---:|
| All | 10 / 31 | 1417 / 4697 | 1069 / 27 / 321 | 0.651788 | 295 / 343 | -0.056417 |
| Fibroblast | 5 / 8 | 1547 / 1547 | 806 / 209 / 532 | 0.015940 | 307 / 333 | -0.032709 |
| Mono_Macro | 7 / 8 | 1635 / 1635 | 432 / 294 / 909 | -0.169012 | 315 / 338 | -0.038609 |
| T | 5 / 8 | 1265 / 1265 | 574 / 175 / 516 | 0 | 307 / 339 | -0.023295 |

Gain是SVC重建标签与Raw Leiden标签的稀释后有效多样性差。正Gain阈值只拟合正差窗口，不应遮蔽完整正负分布。All还受观察人群和覆盖差异影响：共同窗口内Raw7272个单位、SVC14614个单位；不能将All差值提升为全组织改善。三个parent的原生单位cohort与窗口支持相同，可更直接检查标签系统差异，但仍不是独立真值。

Moran表目前仅在各自完整原生轴归一化后取共同343基因；Parent中常量/不可计算基因不进入有效分母，不能一律以343为分母。多数下降是实际原生结果方向，归一化背景控制前不授予生物学因果解释。

## 已触发的必要诊断

1. **Raw K-control**：三个parent特别是Fib/T的群数差可改变Gain解释。新增明确只允许`raw_k_control=false→true`的入口；正式首轮冻结门保持不变，`iter-002`仍消费producer `iter-001`并完整运行Impact。使用已有resolution grid，报告目标K、实际K和matched/nearest状态，不调网格追求结果。All控制后仍受覆盖混杂。
2. **共同基因背景 Moran 与 EMT**：保持首轮每侧每scope实际单位cohort、seed和kNN参数，先统一343共同基因，再执行已有归一化和Moran；EMT保留首轮每侧一次的30000单位union cohort、原provider和评分参数，同时记录共同背景下的coverage、detection与实际cutoff。cutoff若不可用则显式保留，不更改阈值补救。独立局部诊断，不改主流程、不统一All单位、不新增全局配对；记录逐cohort哈希和输入/方法身份。
3. **State边界只读敏感性**：用已保存的bootstrap阈值区间端点对同一连续场重新计数，只描述边界敏感性，不更换正式阈值、不重算原始数据。

## State边界检查

| Scope | 正式State窗口数 | 阈值区间下端 / 上端对应窗口数 | 对应单位数下端 / 上端 | 解释 |
|---|---:|---:|---:|---|
| All | 2286 | 2353 / 2213 | 23272 / 21546 | 区域计数随已估计阈值范围变化 |
| Fibroblast | 540 | 587 / 474 | 5245 / 4045 | 不把单一掩膜视为精确边界 |
| Mono_Macro | 164 | 276 / 164 | 1527 / 770 | 102个窗口的Neff位于阈值附近，均只有4个单位 |
| T | 296 | 296 / 168 | 1649 / 1118 | 123个窗口恰处Neff≈2.828427，均只有4个单位，存在离散边界敏感性 |

这些区间来自现有窗口bootstrap；这里只检查保存的阈值变化如何影响分类，不额外声称其已验证为独立空间观测的总体置信区间。稳定阈值门通过不等于区域边界稳定，也不通过放松门或挑阈值修饰发现。

## Gain边界检查

| Scope | 正式Gain窗口数 | 阈值区间下端 / 上端对应窗口数 | 对应单位数下端 / 上端 |
|---|---:|---:|---|
| All | 357 | 441 / 225 | Raw 2230 / 1111；SVC 5083 / 2586 |
| Fibroblast | 253 | 427 / 117 | 3704 / 925 |
| Mono_Macro | 73 | 166 / 73 | 811 / 333 |
| T | 163 | 253 / 49 | 1379 / 248 |

与State检查相同，这些是固定观测窗口场在保存阈值区间端点下的分类数量，不是区域数量的bootstrap置信区间。尤其T的正式163个Gain窗口对阈值边界敏感；不能用阈值状态`ok`代替区域精确性声明。

正式结果补充：Raw EMT覆盖31/200、cutoff3；SVC187/200、cutoff250，原生绝对分数不可直接比较。Raw三个parent评分单位1558/1526/1300，SVC为10249/10503/8474；各侧All均评分30000，其余单位保持未评分。Anatomy中SVC50626全Other；Raw为231604 Other、109233 Unknown，不能提供独立组织分层验证。All membership仅4506个共同ID（Raw All抽样与SVC覆盖交集），parent分别17455/17605/14219，不能把All membership变化比例称为全组织改变率。

完整322产物已由主代理独立核hash；无阶段错误。Moran部分SVC常量基因明确不可计算，故整体result为partial。服务器PNG中文字体缺失，已从原底表建立独立中文展示刷新副本，不重算科学结果。

## 表达支持补证

原始Raw与reference全量stored X均有限、非负、整数，源文件哈希与正式输入一致，支持count-like声明，但不证明跨assay可比。用已核对的donor行与SVC交付基因轴重建每单位描述统计，在同一343基因轴和相同parent ID集合中：

| Parent | Raw total / detected中位数 | SVC total / detected中位数 |
|---|---:|---:|
| Fibroblast | 112 / 46 | 35 / 18 |
| Mono_Macro | 96 / 50 | 56 / 23 |
| T | 84 / 51 | 37 / 21 |

这些是共同轴count-like表达总量与检测支持，不能当作可直接比较的真实测序深度。即使共同基因背景诊断后Moran方向保留，也不能把剩余差异全部窄归因空间重建；跨assay载体与稀疏性仍不同。完整per-unit底表和输入/输出hash见Analysis该轮`review/depth-audit-001`，独立审视已复核。

## K-control的适用范围

现有实现保存resolution sweep、selection、窗口多样性与相对SVC的delta，只诊断Gain标签粒度。它不输出K-controlled membership。原生All membership的4506单位交集中实际K为6 vs31（全Raw partition才是10），占Raw 30k样本15.02%；其85.8%变化不能在K诊断后改称已控制群数的结果。

## 共同有序基因背景诊断结果

诊断commit `bf3fe5f065d291b9f03259185e7ab40db79d5b6a`；cpu job `20260923-041743-01rkwg`，exit0，无新增OOM或memory max事件。完整29制品经执行代理、独立审视及主代理核验。仅改变ordered gene axis，其余单位cohort/seed/k6/原方法不变。

| Scope | 共同轴负Moran差 / 可计算 | 原生差中位数 | 共同轴差中位数 |
|---|---:|---:|---:|
| All | 298 / 343 | −0.056417 | −0.056678 |
| Fibroblast | 307 / 333 | −0.032709 | −0.032896 |
| Mono_Macro | 315 / 338 | −0.038609 | −0.038049 |
| T | 305 / 339 | −0.023295 | −0.023660 |

可以排除“下降方向仅由完整基因背景差异造成”这一窄解释。不能据此分离random/reference信息转移、跨assay检测支持与真实空间结构的各自作用；没有独立空间真值，亦不能将较低Moran等同于生物学恶化。

EMT共同轴两侧均29/200基因、同provider/seed/tie permutation，但cutoff为2与6；Raw仅2个分值（24898/30000为0），SVC16个分值（13585/30000为0）。均值0.002932与0.019707仍不是可直接比较的生物活性差。此诊断已足够判断本次设计不支持跨侧EMT定量改善，不为得到可比正面结果另调cutoff。

完整入口：Analysis `output/p2-xenium-impact-20260922/iter-001/diagnostics/_retrieved/common-gene-background-001-full-evidence`。

## 最终闭环

K-control完整重跑和独立复核已完成；仅Mono达到exact K，其他nearest，原生Gain对预设分区分辨率敏感。所有最终数值、限制、未知与停止依据统一见[最终审视报告](final-review.md)。主代理确认当前没有影响该保守结论且必须继续解决的重大问题；无需为获得阳性或exact K扩grid、改cutoff或追加其他平台。
