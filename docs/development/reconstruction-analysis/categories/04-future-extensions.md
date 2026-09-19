# 四、等待方案或材料的扩展

[返回总览与状态](../README.md) · [落实文件索引](../implementation-index.md) · [共识与来源](../decisions-and-sources.md)

本类是后续讨论的固定入口，不是从近期范围中消失的事项。每一条保留当前认识、待讨论内容、材料和再次启动条件；没有具体设计时不先创建空接口或 TODO 代码。各项是否互相依赖取决于最终方法，不强行排成一条串行链。

<a id="d1"></a>
## D1：confidence/posterior 保存与消费

**为什么做。** 重建 annotation 与 reference 排名需要一致的基础数值定义，避免同名 confidence 漂移。

**已确认。** D1 与 R7 参考选择合并。保留当前列最大分配权重；max_median 为单位 confidence 中位数，certainty/margin 保留。共享中立数值函数，记录阶段和类别含义，不以 assembly 权重覆盖已有 confidence，不新增拆列或完整 posterior 保存协议。

**依赖与参与。** 已有 reference 材料和实现足够，本批不再等用户另外提供方案；不依赖 D2。

**下一步。** 已按[本批计划](../plans/ot-assembly.md)完成实现和数值回归；不再等待用户补充 confidence 方案，证据见落实索引。

<a id="r7"></a>
## R7：Reference evaluation、ranking 与 selection

**为什么做。** 将 reference 选择作为独立可解释的准备步骤，保存候选评分及输入证据，再显式交给既有重建流程。

**已确认。** 用户提供 reference-preparation 材料包后，已对齐[接入计划](../plans/reference-preparation.md)：本地 whole-file 候选逐个 GA，保留三种分数，max_median top-1，同分按 ID；准备和重建独立，首版接受重复 GA。未分配或非法概率使候选失败，不静默修复。所有排名候选 ST 轴相同，类别/基因/细胞覆盖差异记录。

**工程与科学边界。** 接入不改变 Confidence；筛选直接使用 GA 概率。当前评分表示概率集中程度，不证明 reference 生物学质量，也不校正大小或类别覆盖偏差。完整 GA 缓存、分层公平评分、top-k/ensemble 未获本批实施授权。

**依赖与参与。** 实施使用已交付材料，无用户必交补充材料；真实科学比较后续需要代表样本及判断标准。Agent 负责 GA-only 适配、配置覆盖、batch 指纹和验证，不修改 Atlas 工作树。

**下一步与证据。** 本批工程状态见[总览](../README.md#items)，运行方式和验收见[执行计划](../plans/reference-preparation.md)。历史 A16/A19/A20 与 B2–B4 仍由[来源映射](../decisions-and-sources.md#source-map)追溯。

<a id="r8a"></a>
## R8a：显式配对来源（含 Patient ID 场景）

**为什么做。** 对已有匹配关系的 reference pool 精确取出候选，避免重复手工筛选或在重建时误用旧 filter。

**已确认。** paired 使用用户指定的 pair_column 和 pair_key，不根据 Patient/Sample 名称推断含义。提取保留 AnnData 内容，输出 reference.h5ad 和最小 reference.yaml；主流程消费时清空旧 filter。

**依赖与参与。** 通用实现无需额外材料；具体运行需已有 pool、明确列和值。用户无需现在整理新的 pool，Agent 在实际运行前核对现有输入。

**下一步与证据。** 已与 R7 同批完成接入、验证，见[两步操作与计划](../plans/reference-preparation.md)。自动推断患者匹配仍不在范围内；历史 A17 保留。

<a id="r8b"></a>
## R8b：本地候选与 CELLxGENE 来源

**为什么做。** 允许已准备好的公开 reference 参与统一筛选，而不为不同来源重建另一套后半段。

**本批已定。** 支持本地目录直接 H5AD 和显式候选清单。文件可来自 CELLxGENE，但不因此包含在线搜索、下载、metadata eligibility 或 catalog 能力。

**仍待讨论。** 若要自动寻找公开数据，再核对既有数据仓库、Census/下载接口、筛选条件及资源限制；这些没有被本地筛选实施替代。

**依赖与参与。** 本地接口实施暂无额外用户材料；在线来源需要启动时再讨论范围和现有入口，不要求现在准备。

**下一步与证据。** [本地筛选接入](../plans/reference-preparation.md)已完成；在线来源扩展继续保留为后续事项。历史 A18 仍映射到此，不把本地文件枚举称为完整 CELLxGENE discovery。

<a id="d2"></a>
## D2：gene-wise uncertainty

**最新决定：保持之前的方法，本批不修改，也不列入近期开发；以下定义问题仅留作未来扩展记录。**

**为什么做。** 若重建真正计算出昂贵的基因级 uncertainty/posterior，下游需要时应保存，不能要求分析重新运行重建。

**已有原则与未决问题。** 普通 gene overlap 可现场推导，与本项不同。具体方法、数值含义、gene/observation 轴、保存哪些量以及谁消费均未确定。仓库已有 [gene_uncertainty 实验实现](../../../../revise/backend/kernels/gene_uncertainty.py)，不等于本轮目标方案已经被确认。

**依赖与参与。** 需要用户的方法/公式或已有实验思路，并核对下游真实需求；若来自已有代码，Agent 先阅读和提取差异。是否依赖 D1/R5 由选定方法决定，不预设合并。

**下一步。** 讨论估计对象和可验证性后再确定保存边界与执行计划，不提前设计大型 uncertainty 矩阵。[历史映射](../decisions-and-sources.md#source-map)保留 A15。

<a id="d3a"></a>
## D3a：out-of-spot 重建与 provenance

**最新决定：当前 runner 仍只有 TODO，尚未实现；继续待做，不进入本批。现有位置的表达插值不等于生成 out-of-spot 单位。**

**为什么做。** 未来若生成没有真实 parent spot 的细胞，需说明它们如何产生及证据来自哪里；不能为了兼容分析而制造一个假的 spot_name。

**共识与未知。** 不伪造 parent spot 的原则已明确；生成方法、空间范围、表达来源、标识和验收尚未确定。当前 [sST runner](../../../../revise/backend/runners/sc_svc_super_resolution_application.py) 的 TODO 不是实现证据。

**依赖与参与。** 用户需要参与定义重建目标；已有代码、位置/表达示例或预期结果可供讨论。来源协议与生成算法一起设计，不独立发明字段。现有 sST 真实 parent spot 的保留不等待本项。

**下一步。** 有具体生成方案或案例后再展开方法与 provenance 计划。[历史映射](../decisions-and-sources.md#source-map)保留 B1。

<a id="d3b"></a>
## D3b：复杂 Raw↔SVC correspondence

**为什么做。** 真实存在多对多或其他非平凡关系时，不能丢掉重建已知映射；也不应仅为了统一协议强迫所有技术使用同一 source_unit_id。

**当前边界。** 先由 R3 保留真实单位 ID 或 spot_name/cell_id。是否需要独立 sidecar 取决于具体关系和分析需求；目前没有已确认 schema，不新增通用 pairing 前置条件。

**依赖与参与。** 当前暂无用户必交材料。出现复杂关系时提供实际关系例子及消费问题，再判断是否需要新表示；它可能与 D3a 有关，但二者不是同一个任务。

**下一步。** 有例子后比较现有字段是否足够，再决定是否立独立实现计划。当前关系视图参考 [batch inputs](../../../../revise/batch/inputs.py)，不能将其直接变成分析库必须兼容的入口。[历史映射](../decisions-and-sources.md#source-map)保留 A12 中后置部分。
