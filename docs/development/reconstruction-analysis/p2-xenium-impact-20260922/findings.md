# 问题、诊断与判断

当前状态：`scientific_review_complete_with_limits`；最终结论和按重要性排序的统一表见[最终审视报告](final-review.md)。本表保留实施细节，区分代码错误、比较混杂和科学未知。`pending` 表示尚未检查，不等于通过。每轮审视需引用实际轮次与证据文件；没有证据的解释保持假设。

| ID | 问题／观察 | 前因与机制 | 证据与比较口径 | 已采取措施 | 当前判断 | 对主结论的影响 | 需要用户判断或补充 |
|---|---|---|---|---|---|---|---|
| PRE-001 | 新增交付审计误要求 Raw 覆盖 SVC 全基因轴 | 将 panel overlap 错当输出合法性条件，会拒绝合法 reference 扩展 | 独立审视与主代理核对 `verify_full_ist_delivery.py`；历史mini为422 vs13088基因 | 执行代理修正为参考覆盖断言与Raw交集记录，并补回归 | corrected / independently_reviewed | 若不修会必然误判真实交付 | 无 |
| PRE-002 | 新增 donor 审计仅比较交付内两列 | 同一产物的自洽不能替代来源核验 | 独立代码审视 | 补全部assignment原reference成员/P2身份、Level1、逐行X及provenance hash核验 | corrected / independently_reviewed | 影响donor来源解释 | 精确一致只证明参考信息转移身份，不证明恢复了真实空间表达 |
| PRE-003 | 新运行器代码/参数/缺失产物启动与完成门不足 | HEAD+dirty路径不能还原未提交实现，登记文件缺失不能报完成 | 独立代码审视，主代理已核对 | clean-worktree门、固定科学配置、missing-artifact失败语义 | corrected / independently_reviewed | 影响复现与验收完整性 | 无 |
| PRE-004 | 审视摘要未登记cohort manifest | 新抽样证据若未进入摘要会漏审 | 主代理补充显式read；历史缺失与注册manifest两情形已检查 | 摘要新增cohorts和源文件hash；缺失不补空列表 | corrected / independently_reviewed | 不改变科学计算 | 无 |
| RESOURCE | 两端容器8GiB | 稀疏random与AUCell排名的峰值来源不同 | 重建/原生/K-control/共同轴峰值WS约7.395/5.193/5.029/4.120GB（十进制），均无新增OOM；原rank三副本约4.712GB | 重建保留90%保护；分析启用已审视且严格等价的单rank缓冲适配 | runtime_verified / scientific_review_pending | 环境技术验证不等于科学接受 | 暂无 |
| IDENTITY | 输入、坐标及donor身份 | Raw保护及同群reference信息转移 | 全50,626 assignment逐行精确核对、三文件与补证哈希 | 主代理与独立审视核对完整audit | technical_verified | 支持来源身份，不证明空间表达恢复 | 无 |
| COVERAGE | Raw340837与SVC50626范围不同；三个parent输入ID集合一致 | 98543单位未过eligibility；191668候选因reference子类型不足跳过；All另有独立抽样 | producer coverage supplement与逐broad底表；consumer完整底表已核：All共同有效窗口1417/4697；parents1547/1635/1265 | 完整记录输入、候选、交付与跳过，Raw保持完整 | coverage_verified / interpretation_limited | All差异不能单独归因重建；parent各阶段实际抽样及哈希已核 | 若研究目标必须包含Tumor/上皮，需另行决定方法或参考输入，不能自动补入 |
| GRANULARITY | State/Gain 与标签粒度 | Raw Leiden K为10/5/7/5，SVC为31/8/8/8；标签数量及大小分布不同 | 原生Gain中位数All/Fib/Mono/T为0.652/0.016/−0.169/0；完整正负场已核 | 必要Raw K-control已完成iter-002并独立复核，仅改raw_k_control=true；All/Fib/Mono/T实际K=15/6/8/9，只有Mono精确达到目标 | diagnostic_complete / sensitivity_confirmed | 原生Gain不等于生物改善；All还受覆盖影响 | 不为exact或阳性扩展既定grid |
| ANATOMY | 无可用组织分层 | SVC未覆盖Tumor/上皮，派生Anatomy只形成Other；Raw支持外为Unknown | 完整点表SVC50626 Other，Raw231604 Other/109233 Unknown；三种目标组织区域均无支持 | 核对点表与派生条件表，保留Unknown | limited / no_independent_anatomy_validation | 不能支持组织区域恢复或区域机制 | 若需要该科学主张，须独立组织信息及覆盖对应对象 |
| GENE_AXIS | 多数共同可计算基因原生Moran下降 | 各侧先按422/13088全轴归一化再取343共同基因，背景不同 | 负差295/343、307/333、315/338、307/339；不把常量基因计入可计算分母 | 共同ordered 343轴诊断完成并独立复核：负差298/343、307/333、315/338、305/339；中位数−0.056678/−0.032896/−0.038049/−0.023660 | diagnostic_complete / direction_persists | 已排除仅由基因背景造成负方向；仍受参考载体、共同轴检测支持与稀疏性影响，不能称生物优劣或归因于算法 | 因果分解与独立真值属于后续研究选择 |
| EMT | 原生AUCell分数范围不可直接跨侧比较 | gene coverage、排名背景、实际cutoff及评分单位不同；All不包含重建后的上皮/Tumor | 本轮Raw EMT31/200、cutoff3；SVC187/200、cutoff250；parents实际评分Raw1558/1526/1300、SVC10249/10503/8474 | 共同ordered轴诊断完成：两侧29/200，但cutoff仍2/6；Raw仅2种分值，24898/30000为0；保留所有实际结果 | diagnostic_complete / cross_side_activity_not_supported | 分数升高不能称EMT增强或改善；此任务没有重建上皮/Tumor状态以支撑转化判断 | 若未来需要肿瘤EMT主张，须重新明确对象与独立证据 |
| DONOR | 50626 assignment来自10021参考细胞，最大复用34次 | random同cluster载体抽样，SVC表达逐行等于reference | 实际30000评分单位对应8996参考细胞；完整ID/cluster/score严格一致 | donor-aware保存表诊断完成：均值assignment0.0556223；保留broad权重的donor等权0.0575542；全donor等权0.0439714 | diagnostic_complete / independently_reviewed | 后者改变大类组成，不能把全部下降归因复用；所有参考细胞仍非生物独立重复 | 无独立真值时不把参考表达结构当独立恢复验证 |
| STABILITY | 阈值通过不代表边界稳定 | 最小4单位稀释产生离散Neff，窗口阈值区间可跨越大量相同值 | Mono102/T123窗口处于Neff≈2.828427原子点；T State计数在保存阈值区间端点为296/168，Gain253/49 | 从保存bootstrap阈值与完整连续场只读检查，无改阈值 | limited / boundary_sensitive | 区域掩膜只作描述；端点计数不是区域数量总体置信区间 | 不需为了正面结论调参；独立空间稳健性仍未知 |
| PRESENTATION | 服务器PNG中文字体缺失 | Linux环境未找到绘图配置中的中文字体 | 主代理视觉检查标题/轴标签为方框 | 独立native-display-001重绘129图，k-control-display-001重绘133图；各自189/205张表及result哈希不变，中文及解释限制已复核 | corrected / technical_visual_verified | 不改变科学计算及原始发布证据 | 无 |

生产端已核事实及科学边界见 [producer-review](producer-review.md)；consumer正式结果已取得并独立核验322项哈希；共同背景、donor、表达支持、K-control及最终整体复核均已完成，状态scientific_review_complete_with_limits；方法生物改善未获接受。原Level1→GA转移的四类边界经主代理复算：内部636，外部迁入360，迁出500。

分析证据统一入口为Analysis `output/p2-xenium-impact-20260922/iter-001/`：正式回收`_retrieved_final/20260922T194202Z`、主代理核验`review/main-native-001`、donor审计`review/donor-aware-emt-native-002`、展示修复`review/native-display-001`。失败的donor首次schema尝试保存在`review/donor-helper-source-001/actual-run.log`，未生成伪汇总；第二版增加已注册full reconstruction labels作为独立标签来源。
