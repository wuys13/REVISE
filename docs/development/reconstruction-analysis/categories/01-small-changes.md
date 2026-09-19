# 一、范围已定的小改动

[返回总览与状态](../README.md) · [落实文件索引](../implementation-index.md) · [共同决策](../decisions-and-sources.md#confirmed-decisions)

这一类放能够明确限定行为、可独立验证的变化。“改动小”不等于可跳过执行计划与用户对齐。目前只有 R1 满足这一分类；random 默认与全类型输出放在第二类一起交付。

<a id="r1"></a>
## R1：删除 sST 最终逐生成细胞缩放到 10,000

**为什么做。** 当前 sST 算法内部尺度处理与最终输出处理有不同目的。用户明确要求移除最后对每个生成细胞独立缩放到 10,000 的行为，避免把它与 spot 内表达分配混为一谈。

**已知事实。** sST 先将共同基因上的 spot 表达、reference 表达归一化，再进行分配与邻域调整；建图 normalization/log1p 作用于副本。原来最后有两个分支：按 parent spot、按基因校正回目标表达，或将每个生成细胞单独缩放到 10,000。目标 spot 本身已经归一化，保留前一分支并不恢复 raw counts。

**用户确认。** 仅删除最终逐细胞分支，保留前面的算法内部归一化；按照讨论后的方案保留 spot 校正，移除二选一开关。不改 confidence，不在此次重新定义表达分配算法。

**已确认兼容方案。** Application/Engine 的旧键无论 true/false 均明确报错，提示删除该键并说明固定 spot 逐基因校正；模板、参数传递、元数据及相关测试一起清理。

**依赖与节奏。** 不依赖 R2、cluster 或 confidence 新方案。执行计划已获用户确认；复核暂停修改后完成针对性验证。

**需要用户参与。** 新增科学讨论及必需材料：暂无。工程兼容方案及验证范围已对齐。现有数据如用于真实案例验证，仅作补充，不先要求用户准备大型重建输入。

**下一步。** 查看[执行计划与验证记录](../plans/sr-final-scaling.md)。本批实施状态在[总览](../README.md#items)维护；后续只讨论第二批，不自动执行。

**实现依据。** [sST application runner](../../../../revise/backend/runners/sc_svc_super_resolution_application.py)、[建图副本](../../../../revise/backend/ops/topology.py)、[配置入口](../../../../revise/application/config.py)。相关差异及验证边界见执行计划。
