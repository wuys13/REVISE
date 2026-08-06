REVISE 1.0 → 2.0 迁移与完善任务清单
1. 1.0 优化同步与整体设计保持
1.1 1.0 已有优化迁移
梳理 1.0 中已完成的优化内容，确保 2.0 架构同步更新。
包括但不限于：
参数优化；
流程优化；
模块接口调整；
已验证的稳定性改进。
1.2 2.0 三大类 SVC 设计整合
当前 2.0 已实现 SVC 归入三大类别的最终设计方向。
后续需要：
将已有 1.0 优化内容迁移并 apply 到 2.0 三类 SVC 架构中；
保持已有模块设计一致性；
避免迁移过程中破坏 2.0 的统一接口和扩展性。

2. 规范化重建评估体系
2.1 重建质量监控
建立统一的 reconstruction evaluation framework：
增加标准化 warning / logging 机制；
对每次 reconstruction 自动评估：
是否存在异常情况；
是否满足预设质量标准；
是否需要用户进一步检查。
2.2 评估指标体系完善
后续进一步明确：
必要 warning 类型；
阈值设置；
不同数据类型（sp-SVC / sc-SVC / sc-SVC-sr）的评价标准；
失败模式分类与反馈。

3. 指标体系与接口重新对接
3.1 通用评价指标接口
重新整理并接入已有 benchmark / evaluation 指标：
聚类与空间结构：
ARI 等 clustering metrics；
spatial domain consistency metrics。
3.2 Cell type / marker 评价接口
重新接入细胞类型相关指标：
cell type marker specificity；
TMP、MER 等 marker evaluation metrics；
相关 annotation quality assessment。
3.3 下游应用分析接口
恢复并统一应用层分析接口，包括：
CCI 分析
ligand-receptor interaction；
cell-cell communication analysis。
Pathway 分析
pathway activity；
functional program evaluation。
Cell type neighborhood 分析
spatial neighborhood；
niche / microenvironment analysis。

总体目标
完成从 REVISE 1.0 → REVISE 2.0 的工程与分析体系迁移：
保持 1.0 已验证优化，同时适配 2.0 统一 SVC 架构；建立标准化 reconstruction quality control；恢复完整 evaluation 与 downstream analysis ecosystem。