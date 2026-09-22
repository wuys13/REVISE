# 问题、诊断与判断

此表区分代码错误、比较混杂和科学未知。`pending` 表示尚未检查，不等于通过。每轮审视需引用实际轮次与证据文件；没有证据的解释保持假设。

| ID | 问题／观察 | 前因与机制 | 证据与比较口径 | 已采取措施 | 当前判断 | 对主结论的影响 | 需要用户判断或补充 |
|---|---|---|---|---|---|---|---|
| PRE-001 | 新增交付审计误要求 Raw 覆盖 SVC 全基因轴 | 将 panel overlap 错当输出合法性条件，会拒绝合法 reference 扩展 | 独立审视与主代理核对 `verify_full_ist_delivery.py`；历史mini为422 vs13088基因 | 执行代理修正为参考覆盖断言与Raw交集记录，并补回归 | corrected / independently_reviewed | 若不修会必然误判真实交付 | 无 |
| PRE-002 | 新增 donor 审计仅比较交付内两列 | 同一产物的自洽不能替代来源核验 | 独立代码审视 | 补全部assignment原reference成员/P2身份、Level1、逐行X及provenance hash核验 | corrected / independently_reviewed | 影响donor来源解释 | 精确一致只证明参考信息转移身份，不证明恢复了真实空间表达 |
| PRE-003 | 新运行器代码/参数/缺失产物启动与完成门不足 | HEAD+dirty路径不能还原未提交实现，登记文件缺失不能报完成 | 独立代码审视，主代理已核对 | clean-worktree门、固定科学配置、missing-artifact失败语义 | corrected / independently_reviewed | 影响复现与验收完整性 | 无 |
| PRE-004 | 审视摘要未登记cohort manifest | 新抽样证据若未进入摘要会漏审 | 主代理补充显式read；历史缺失与注册manifest两情形已检查 | 摘要新增cohorts和源文件hash；缺失不补空列表 | corrected / independently_reviewed | 不改变科学计算 | 无 |
| RESOURCE | cpu 容器只有 8 GiB | 宿主机内存不等于容器配额；旧稠密估算又不等于 random 稀疏峰值 | 重建代理读取 cgroup memory.max=8589934592、cpu.max=4核，baseline约2.06GiB | 保持full gate关闭；补查真实矩阵编码/nnz与复制热点 | investigating | 尚不证明可安全执行full，也不能仅凭旧估算宣布不可行 | 若等价实现优化仍不足，再提供具体资源方案 |
| IDENTITY | 输入、坐标及 donor 身份 | full random 需保持 Raw 并追溯同群 donor | 待 full 输入/交付审计 | 准备审计 | pending | 错误则阻断 | 无 |
| COVERAGE | Raw/SVC 覆盖与分母 | eligibility、过滤、独立抽样改变观察范围 | 待各 scope 单位与窗口支持表 | 保持真实缺失与跳过 | pending | 决定比较适用范围 | 科学覆盖要求待结果后说明 |
| GRANULARITY | State/Gain 与标签粒度 | Raw Leiden 与 SVC 主标签含义及 K 不同 | 待分群及多样性表 | 如影响解释，触发单独 K-control | pending | 不可直接称生物改善 | 科学意义待结果审阅 |
| ANATOMY | Other/Unknown 的来源 | Anatomy 由交付 SVC 定义，不是独立组织真值 | 待点与窗口组成表 | 不补标签、不混用分母 | pending | 约束空间解释 | 必要时需独立组织信息 |
| GENE_AXIS | 共同观测与参考扩展基因 | Xenium panel 与 SVC 全表达基因范围不同 | 待基因交集、表达身份及 shared-gene Moran | 原生与跨侧结论分开 | pending | 不得声称未测基因已独立验证 | 额外验证证据若需要再列 |
| EMT | AUCell 背景与 cutoff | 原生基因轴和覆盖可能不同；provider参数是检测数分位而非固定排名比例 | 历史mini底表：Raw 422基因/EMT31/200/cutoff5，SVC13088基因/EMT187/200/cutoff255；本轮full仍待核验 | 必要时共同背景诊断，仍单独核验实际cutoff | historical_risk_confirmed / full_pending | 原生绝对分数差不能直接代表改善 | EMT 生物意义待审阅 |
| DONOR | reference 信息转移及复用 | random 从同 cluster expression carrier 抽样 | 待 donor 唯一数与复用分布 | 必要时 donor-aware 诊断 | pending | 限制表达变化的因果解释 | 无独立真值时保留限制 |
| STABILITY | 参数、抽样与阈值稳定性 | 单个 seed/窗口不证明稳健 | 待支持曲线、bootstrap 与 cohort 身份 | 对具体疑点做单变量诊断 | pending | 翻转则报告不稳健 | 必要时选择后续研究方向 |
| PRESENTATION | 图文与底表一致性 | 报告摘要可能遗漏限制 | 待 result 与完整登记底表 | 独立复核 | pending | 错误解释必须修复 | 无 |
