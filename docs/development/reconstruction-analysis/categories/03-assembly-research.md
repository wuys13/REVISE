# 三、assembly 方法与科学比较

[总览与状态](../README.md#items) · [落实索引](../implementation-index.md) · [本批执行计划](../plans/ot-assembly.md)

<a id="r5"></a>
## R5：两种 OT assembly

**为什么做。** 在统一空间轴上利用局部空间信息与 SC 表达加权，提供 random 和固定 mean 之外的分配方式。

**已确认。** 外部提供 SC Level2；现有 LR 与 SVC_cluster 保留。ST 匹配副本按同 cluster 邻居 0.2/自身 0.8 混合。within_cluster 对同群 SC 单细胞做 OT；outside_cluster 对同 broad type 全部 SC cluster mean profiles 做 OT。TACCO 权重按行归一化后乘全基因表达。默认仍 random。

**修正记录。** 早先两图结构匹配、代表状态插值是讨论候选；用户随后明确采用 Benchmark 风格“空间表达加权→OT→表达加权”。本批不叠加显式图结构迭代。outside_cluster 也重新计算 ST×profile OT，不复用方向相反的 SC×cluster 注释矩阵。

**依赖。** 复用已完成的 R2/R3 配对载体和发布路径；不依赖新 posterior、D2 或真实验收。

**需要用户参与。** 本批暂无必需材料或待选算法。科学效果后续审阅。

**下一步与证据。** 配置、数值实现、边界测试与文档已完成；下一步是真实比较，见[比较协议](../assembly-comparison-contract.md)，工程证据见落实索引。

<a id="r6"></a>
## R6：比较 Notebook 与最终默认值

**为什么做。** 工程上输出表达不证明空间亚群结构保持。用户希望观察 T、Macro、CAF 重建后能否找到与原空间 sc-SVC subtype 对应的亚群。

**已确认。** 独立 Notebook 读取 mean/random/within/outside；每种重建表达重新做 Leiden，与旧空间 subtype 比较空间图、列联表和 ARI/NMI。统一参数，显示缺失与共同范围，保留多 resolution；不把旧 cluster 直接作为新分群，不自动挑最高分或修改默认值。

**依赖。** Notebook 开发依赖统一 SVC 文件协议，可用小 fixture；真实运行依赖各方法产物，暂不启动。

**需要用户参与。** 开发暂无；真实样本、规模、科学解释随后由用户逐项审阅并安排。

**下一步与证据。** [Notebook](../../../../reproduce/case/assembly_comparison.ipynb)与[计划](../plans/ot-assembly.md)。工程执行不替代真实生物学比较。
