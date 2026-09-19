# Reconstruction 与 Analysis 对齐：决策、来源与执行边界

本页保存 2026-09-19 的决策与来源；[总览](README.md)是文档入口和唯一状态表。它把历史事项压缩成可引用的决策地图；下面的编号、原始含义和当前落点足以独立阅读，不以未跟踪的历史原文存在为前提。

状态词的含义是：**已确认**只记录本轮用户明确给出的行为；**提案**记录实现建议或历史材料中的方向，不能当作用户同意；**待方案**表示必须等待用户的科学或协议方案。`other_reviews/GPT.md` 是可选的未跟踪历史线索，不是本页的必要输入。

<a id="confirmed-decisions"></a>
## 已确认行为与边界

- **Reference Preparation 新批次**：用户已授权[独立准备与显式消费](plans/reference-preparation.md)。paired/screen 均纳入；准备与重建分开；单样本和 batch 消费 reference.yaml，校验前覆盖旧 reference/filter；首版接受重复 GA。保留 max_median top-1，未分配候选失败，跨候选 ST 轴一致，类别/基因/细胞差异记录而不新造公平性校正。当时 Confidence、在线发现和完整 GA 缓存后置；最新 D1 已与 R7 合并共用数值定义，其他仍后置。

- **iST assembly**：默认仍为 `random`，保留 `mean`；新增已对齐的 `within_cluster` / `outside_cluster`，详见[本批计划](plans/ot-assembly.md)，最终默认值留待科学比较。
- **类型遍历与纳入**：whole-sample 在 GA 后遍历样本实际出现的全部 broad cell types；对当前有效 reference 的对应 broad type 检查 Level2。只有超过一个不同的、非空有效 Level2 类别才进行该类型的 LR；不满足时提示并记录原因，提示应提前完善 reference Level2 注释；最终 SVC 只合并符合条件的类型。显式旧单类型调用保留既有兼容规则：一个有效 Level2 可继续运行，零个仍失败。
- **失败语义**：已满足纳入条件的类型在计算中报错时，整个 sample 失败，不发布新结果，也不能把运行故障改写成类型跳过；历史结果受保护。
- **sST 边界**：只删除最后把逐生成细胞缩放到 10,000 的分支；保留内部 ST/reference normalization、建图副本 normalization/log1p 和按 parent spot/基因的总量校正。固定 spot 校正开关的移除沿用本轮已确认的工程方向。
- **Raw 与缺失**：Raw 保留完整原始观测轴；标签只回填本次真实获得的推断，没有推断的单位保持缺失，不为覆盖率补造标签。
- **Confidence**：当前复用已有 `Confidence` 列；GA/LR 拆分、posterior 和 uncertainty 的具体协议等用户方案，不能由执行者自行补齐，也不阻塞其他已定工作。
- **消费者边界**：Analysis 作为独立消费者协同；生产端交付可读的 Raw/SVC 与配置，不要求消费者导入 reconstruction backend、cell-type 目录或 paired carrier。消费者协同不等于本页代改兄弟库。
- **范围演进**：首次交付只落分层文档；后续用户已授权 R1 实施及针对测试，并进一步批准将 R2/R3/R4 合并为第二批实施；`categories/01-small-changes.md` 与 `categories/02-sample-delivery.md` 负责前两类详细规划，其余第三、四类保留问题和依赖入口，R7/R8 按上述新增授权实施。
- **交付方式**：纳入 Git 工作树，不 commit、不 push、不发布外部 issues。R1 已完成针对测试；第二批生产端实现和工程套件已完成，真实 P2CRC_Xenium 验收、分析库联合验收和科学解释分开记录，默认不运行昂贵的全量重建。
- **历史边界**：旧 whole-SVC 计划中的 `mean` 默认、route-qualified 输出等不自动成为本轮决定；Raw coverage 字段只是可选提案，尚未锁定为公共协议。

<a id="evidence"></a>
## 实现依据与核对边界

| 来源 | 本页用途 | 当前状态 |
|---|---|---|
| [`revise/backend/runners/sc_svc_super_resolution_application.py`](../../../revise/backend/runners/sc_svc_super_resolution_application.py) | sST 生成、内部 normalization、spot 校正与 out-of-spot TODO 的代码事实 | producer 当前工作树来源；先前暂停的 19 个文件已获 R1 恢复授权；验证见 R1 计划 |
| [`revise/application/ist_assembly.py`](../../../revise/application/ist_assembly.py) | 现有 `mean`/`random` assembly 和 donor/cluster 轴语义 | 复用现有实现的来源；不代表 whole-sample 改造已完成 |
| [`revise/application/publication.py`](../../../revise/application/publication.py) | 当前 H5AD publication、paired carrier 与事务性发布边界 | 现状依据；R2/R3 需要在此边界上扩展 |
| [`revise/batch/runner.py`](../../../revise/batch/runner.py) | 当前 sample、cell-type task、fingerprint、handoff 和恢复状态 | 现状依据；R2 将改变 task ownership，不能把目标写成已实现 |
| [`revise/application/preprocess.py`](../../../revise/application/preprocess.py) | broad/Level2 label normalization 与 reference 过滤的现有入口 | 只作为既有处理依据；有效 Level2 规则不在本页发明新阈值 |
| [`tests/application/test_ist_publication.py`](../../../tests/application/test_ist_publication.py) | 现有 paired/mean/random、seed、发布回滚和碰撞检查的验证入口 | 既有验证入口；R1 针对套件包含 application 测试，详见 R1 记录 |
| [`configs/batch/batch.yaml`](../../../configs/batch/batch.yaml) | 当前 batch 的 `cell_types`、`Level2` 和 `paired/mean/random` 配置事实 | 旧输入事实；不等于本轮 sample-level 目标已落地 |
| `other_reviews/GPT.md` | A1–A22/B1–B4 的历史事项来源 | 可选、未跟踪；本页已完整转述其事项，不依赖该文件 |
| `docs/plans/2026-08-03-002-feat-v2-unified-svc-output-plan.md` | 旧 whole-SVC 方案的可选本地线索 | 不作为本轮依据；其中 mean 默认和 route-qualified 方向不继承 |
| Analysis consumer working tree | 当前消费者尺度/入口状态的外部事实 | 2026-09-19：已取消 global gate、改为 per-side；normalized non-log 仍欠缺；本页不再写“全局门槛仍在” |

分析库可选本地线索位于兄弟目录 `REVISE_Analysis_Agent`：`pyproject.toml`、`docs/migration.md` 与包源码用于核对 backend 解耦；`revise_analysis/io.py`、`methods/partition.py`、`methods/spatial.py` 和 `docs/input-output.md` 用于核对分侧矩阵声明及计算支持。关键观察已在本页和 R4 计划转述，即使该目录不存在也能理解待办；联合验收时仍需取得实际消费者版本。

本轮新增但不在历史编号中的决定：sST 最终逐细胞 10,000 分支删除归入 [R1](categories/01-small-changes.md#r1)；逐类型 Level2 资格和全样本失败归入 [R2](categories/02-sample-delivery.md#r2)；文档先行、前两类详细规划、第三四类保留、消费者独立协同及版本管理方式均在本页保留。

<a id="chronology"></a>
## 重要讨论与更正

1. 较早的 2026-08-03 计划提出 mean 默认与 route-qualified 输出；本次提供的 GPT.md 倾向 random、sample 级输出。本轮用户明确选定 random 默认、保留 mean，不把两份历史材料混为同一共识。
2. 先前把消费者表达尺度问题概括成全局门槛的说法已更正：当时 consumer 是 per-side 条件；最新已固定线性非负非 log X，C1 按源声明如实报告可用性，旧 log1p 支持快照已失效。
3. 主代理在整体拆分尚未对齐时提前启动 sST 实施。用户纠正后已停止：19 个文件修改未测试、未提交。文档整理阶段仅保护现场；后续用户确认 R1 后才恢复实施。
4. 当前文档将“科学选择待用户方案”和“工程衔接可先做”分开：R1–R4 属于前两类详细规划；R2/R3/R4 已合并为第二批并进入实施，当时 R5–R8、D1–D3 只保留未来入口；后续 R7/R8 已单独对齐并授权，见本页新增决定。
5. 先前曾把 Raw coverage 字段写成 R3 的必备验收项；本页已更正为可选提案，不把它写入本轮公共协议。
6. 2026-09-19 用户批准第二批 R2/R3/R4 一次实施：旧单类型入口保留为兼容调用；whole-sample 中部分 Level2 缺失的 reference 行只在 LR 资格判断中排除且要求超过一个有效类别，显式旧单类型允许一个有效类别、零个仍失败；原始标签保留、推断写入独立列。随后生产端实现和工程套件完成（验证摘要见 [`verification-2026-09-19.json`](verification-2026-09-19.json)）；真实样本、科学解释、sST 表达和分析库联合验收仍分别记录，不能由工程通过代替。

## 如何解释历史材料

- `other_reviews/GPT.md` 中的“基本想清楚”“事情确定”是整理文档的状态词，不是逐轮用户确认记录。
- A1–A15、A20–A22 的当前落点优先依据本轮明确行为；原材料中的实现顺序、字段名和默认建议只有在本页再次标为已确认时才有效。
- A16–A19、B2–B4 对应 R7/R8；后续已对齐本地准备与 max_median top-1 接入，在线 discovery 和其他选择策略继续后置，也不是前两类实施的前置条件。
- D1 保留已有 `Confidence`，最新已与 R7 合并定义；不再等待独立方案，也不设计 GA/LR 分列。
- D2、D3a、D3b 不为本轮 R2/R3 增加 gene-wise uncertainty、out-of-spot 或复杂 sidecar 的伪接口。
- “当前落点”列表示文档归属，不表示代码、数据、测试或消费者接线已验收。
- 代码链接用于定位现有行为和边界；它们不能被解释为对修改后行为的证据。
- 本页不替代用户逐项确认；若后续方案改变 Raw/SVC 公共 contract，应记录决策变化、更新类别内容，并同步总览状态。

<a id="source-map"></a>
## 历史事项对应表

下表完整映射历史材料的 A1–A22/B1–B4。每个编号一行；“当前落点”是本页的文档归属，不表示对应代码已经完成。

| 编号 | 原始含义简述 | 当前落点 | 确认/建议边界 |
|---|---|---|---|
| A1 | 所有 route 向 downstream 提供完整 `SVC.h5ad` | [R2 全类型 SVC](categories/02-sample-delivery.md#r2)、[R3 完整 Raw](categories/02-sample-delivery.md#r3)、[R4 分析交接](categories/02-sample-delivery.md#r4) | whole-sample 交付方向已确认；具体字段和 route 兼容性仍是工程工作 |
| A2 | 取消 sc-SVC 的 per-cell-type 最终输出结构 | [R2 全类型 SVC](categories/02-sample-delivery.md#r2) | sample-level task 与合并目标已确认；内部 carriers 可保留 |
| A3 | sc-SVC 遍历所有可重建 broad cell types | [R2 全类型 SVC](categories/02-sample-delivery.md#r2) | “GA 后所有实际出现类型”已确认；有效 Level2 判定按本页规则，不增设细胞数阈值 |
| A4 | 保留 mean/random，增加 cluster assembly | [已有 assembly](categories/02-sample-delivery.md#existing-assembly)、[R2 全类型 SVC](categories/02-sample-delivery.md#r2)、[R5 cluster assembly](categories/03-assembly-research.md#r5) | random 默认、mean 保留；最新已对齐 within/outside 两种 OT assembly |
| A5 | cluster 完成前可把 random 设为默认 | [R2 全类型 SVC](categories/02-sample-delivery.md#r2) | 本轮已确认 random 默认；旧稿中 mean 默认不继承 |
| A6 | 后续比较 mean/random/cluster 并选择最终默认 | [R6比较](categories/03-assembly-research.md#r6) | Notebook 目标已定：T/Macro/CAF 重建后 Leiden 对照原空间 subtype；真实比较与默认值后置 |
| A7 | 生成 analysis-ready Raw H5AD | [R3 完整 Raw](categories/02-sample-delivery.md#r3) | Raw 交付方向已确认；具体公共字段不扩展为未确认协议 |
| A8 | Raw 写入 Level1/Level2，避免 Analyst 重算 | [R3 完整 Raw](categories/02-sample-delivery.md#r3) | 只回填真实获得的推断；缺失保持缺失 |
| A9 | 保存昂贵 confidence/posterior | [D1confidence](categories/04-future-extensions.md#d1) | 与 R7 共用数值定义已确认；本批不新增 GA/LR 分列或完整 posterior 保存 |
| A10 | Raw 保留未进入 reconstruction 的全部 ST units | [R3 完整 Raw](categories/02-sample-delivery.md#r3) | 已确认；不能为 pairwise 裁剪 Raw |
| A11 | 标记 Raw unit 是否进入最终 reconstruction | [R3 完整 Raw](categories/02-sample-delivery.md#r3) | coverage 字段只是可选提案，未锁为本轮公共协议 |
| A12 | 保留真实 Raw↔SVC correspondence，不强制 universal ID | [R3 完整 Raw](categories/02-sample-delivery.md#r3)、[D3b 复杂映射](categories/04-future-extensions.md#d3b) | 真实关系保留方向已确认；复杂 sidecar 仍待方案 |
| A13 | sST 保留 parent spot 信息 | [R3 完整 Raw](categories/02-sample-delivery.md#r3)、[D3a out-of-spot](categories/04-future-extensions.md#d3a) | `spot_name`/`cell_id` 关系保留；不制造一对一或虚假 parent |
| B1 | 为 out-of-spot SVC 建 provenance | [D3a out-of-spot](categories/04-future-extensions.md#d3a) | 不制造 parent spot 已确认；生成方式和字段待设计 |
| A14 | 不持久化普通 gene-overlap provenance | [no gene overlap](categories/02-sample-delivery.md#no-gene-overlap) | 便宜推导不新增字段已确认 |
| A15 | 保存未来真正昂贵的 gene-wise uncertainty | [D2geneuncertainty](categories/04-future-extensions.md#d2) | 原则保留；具体方法和字段等待用户方案 |
| A16 | Reference Selection 成为 reconstruction 前置阶段 | [R7reference评估](categories/04-future-extensions.md#r7) | 后续已确认独立 preparation 支线，显式消费 reference.yaml |
| A17 | 支持按 Patient ID 找 reference candidates | [R8aPatient](categories/04-future-extensions.md#r8a) | 后续来源；pool 组织和选择策略待方案 |
| A18 | 支持按 CELLxGENE 找 reference candidates | [R8bCELLxGENE](categories/04-future-extensions.md#r8b) | 后续来源；retrieval/filtering/sampling 待方案 |
| A19 | 两类 reference 来源汇合到统一后半段接口 | [R7reference评估](categories/04-future-extensions.md#r7)、[R8aPatient](categories/04-future-extensions.md#r8a)、[R8bCELLxGENE](categories/04-future-extensions.md#r8b) | 后续已确认 paired/本地 screen 汇合 reference.yaml；在线 discovery 后置 |
| B2 | 用 ST Unit Confidence 对候选 reference 排序 | [R7reference评估](categories/04-future-extensions.md#r7) | 后续接受材料包三种中位数分数，max_median 为本批选择规则 |
| B3 | 持久化 candidate ranking 中间文件 | [R7reference评估](categories/04-future-extensions.md#r7) | 本批沿用 report.json 并补受控输入/GA 证据 |
| B4 | 决定 top-1、top-k 或 ensemble | [R7reference评估](categories/04-future-extensions.md#r7) | 后续已确认本批 max_median top-1；top-k/ensemble 后置 |
| A20 | 将 selected reference 纳入 provenance/fingerprint | [R7reference评估](categories/04-future-extensions.md#r7) | 本批已确认实际 reference 与外部配置/报告证据进入 provenance/fingerprint |
| A21 | 简化 Analyst-facing output 为 Raw + SVC | [R4 分析交接](categories/02-sample-delivery.md#r4)、[C1分析协同](categories/02-sample-delivery.md#c1) | Analysis 入口边界已确认；消费者独立协同，不代替其内部方案 |
| A22 | 移除 Analyst 对 REVISE backend 的运行依赖 | [backend independence](categories/02-sample-delivery.md#backend-independence) | 目标边界已确认；包依赖与源码已静态核对；不等于新交接已联合验收 |

<a id="execution-boundary"></a>
## 当前实施边界

本页只承担决策和来源索引，不承担代码实现、测试、数据运行或 external issue 发布。第二批的执行细节和验收记录分别落在[总批次计划](plans/whole-sample-delivery.md)和[验收清单](acceptance.md)；允许的下一层落点如下：

- `categories/01-small-changes.md#r1` 详细记录 R1 的 sST 分支删除；保留内部 normalization 和 spot correction，不恢复 raw counts，不修改 confidence。
- `plans/whole-sample-delivery.md` 记录第二批 R2/R3/R4 的统一执行顺序、发布边界和验收入口；生产端及工程验收已完成，剩余真实样本和科学/表达联合验收仍在该页与[验收清单](acceptance.md)跟踪；`categories/02-sample-delivery.md#r2` 详细记录 R2 的 sample-level iST 合并、类型资格、全样本失败与历史结果保护。
- `categories/02-sample-delivery.md#r3` 详细记录 R3 的完整 Raw、真实 annotation、缺失保留和 correspondence；coverage 只作为可选提案。
- `categories/02-sample-delivery.md#r4` 与 `#c1` 记录已验证的生产端 handoff 和消费者独立协同；不要再写全局尺度门槛，按最新线性契约区分每侧实际来源与可用性。
- `categories/02-sample-delivery.md#existing-assembly`、`#backend-independence`、`#no-gene-overlap` 只承接已存在的工程边界，不把旧 whole-SVC plan 写成当前决定。
- `categories/03-assembly-research.md#r5`、`#r6` 承接已批准的两种 OT assembly 和比较 Notebook，科学解释另验收。
- `categories/04-future-extensions.md#r7`、`#r8a`、`#r8b` 由已授权的 [reference preparation 计划](plans/reference-preparation.md) 承接本地实现；`#d1` 与 R7 合并实施；`#d2` 保持现状；`#d3a`、`#d3b` 保留后续需求，不统一阻塞前两类；具体分析或评分若依赖新语义，再记录局部依赖。

当前不在本页决定：复杂 sidecar schema、coverage 字段最终命名、GA/LR confidence 拆列、posterior/uncertainty 指标、新的 reference 公平性校正、top-k/ensemble、额外两图结构优化，以及任何 CLI 全面改名或发布版本决策。第二批已确定 `sample.yaml` 使用相对路径、层级 sample ID 使用可逆单段映射，工程验证已完成；真实样本和科学/表达联合验收仍待进行。

## 各层的内容责任

| category | 本页交给下一层的最小内容 | 不得越界的内容 |
|---|---|---|
| 01-small-changes | R1 的删除范围、保留的 normalization/spot correction、旧产物保护 | 不恢复 raw counts，不把旧配置兼容策略写成科学共识 |
| 02-sample-delivery | R2 的资格/跳过/失败语义；R3 的 Raw 保留与真实 annotation；R4/C1 的生产到消费者路径 | 不加入 coverage 必选字段、不恢复 global gate、不代改 consumer |
| existing-assembly | 现有 mean/random carrier 及 seed/rollback 事实 | 不把旧 mean-default plan 当当前默认，不实现 cluster |
| backend-independence | Analysis 只读交付产物的边界及已核对依赖 | 不把 backend 解耦当作新交接的联合验收证据 |
| 03-assembly-research | cluster 的问题定义和后续比较入口 | 不猜空间连续性、donor 分配或最终默认 |
| 04-future-extensions | Reference、confidence、uncertainty、out-of-spot、mapping 的依赖清单 | 仅按已对齐计划接入本地 paired/screen；不新增在线 discovery、posterior 或复杂 sidecar |

首次交付包括总览、四类事项页、四份执行计划和本页。后续用户已确认并授权第一批 R1：补齐旧键迁移提示，复核删除范围并完成针对测试；随后授权第二批 R2/R3/R4 合并实施，生产端和工程套件已完成。第二批不修改分析库并行工作树、不发布外部 issues、不 commit、不 push；真实命令保存在[验收清单](acceptance.md)，P2CRC_Xenium、科学解释和 sST 表达结果继续回填。

### R1 后续确认（2026-09-19）

用户选择旧键报错并提供明确迁移说明，随后批准 R1 执行计划。Application 与 Engine 均拒绝 true/false，固定现有 spot 逐基因校正，保留内部归一化。实现与验证记录见 [R1](plans/sr-final-scaling.md)。

## 2026-09-19 后续修正：OT assembly 批次

用户批准[两种 OT assembly 与配套计划](plans/ot-assembly.md)。R5 固定 SVC_cluster；空间表达混合限定同群；within 对 SC 单细胞，outside 对同 broad type 全部 cluster mean profiles；均重新计算 TACCO OT 并用全部归一化权重重建全基因表达。默认 random 不改。

D1 与 R7 共享 confidence 定义，不再等待独立方案。R6 先交付 T/Macro/CAF 重新 Leiden 与原空间 subtype 对照 Notebook，真实运行后置。D2 保持现状，不列修改；D3 待做。C1 按当前消费者已定 finite/nonnegative/unlogged linear X 契约更新 producer，之前“接受 log1p、normalized nonlog 待定”的快照已被取代。历史快照不再作为现行消费规则。


## 双库审阅前收尾（2026-09-19）

- 仅做小型工程浅验收，保留产物与版本证据；真实全量与科学解释后置。
- 用户指定 Real_application 既有重建文件的 SVC_cluster 为横向比较 Ground Truth。旧 sc-SVC 是标签基线，不是第五种表达方法；四种 assembly 独立重跑 Leiden。
- 当前协议与审阅对象见 [双库审阅入口](cross-repo-review.md)，Notebook 输入、角色及元数据核对见 [比较约定](assembly-comparison-contract.md)。
- 历史第二批验收归档，当前入口不再沿用旧线性表达不支持的结论；D1 已完成，D2 保持现状，D3 与在线来源暂缓。

- 只读对照发现 sST 已有 x/y 未写入 obsm.spatial、已有 cell_type 未作 broad 别名来源；本轮仅修复这两处生产端交接，不修改生成位置、表达或标签值。

- 用户再次澄清 Ground Truth 指重建后已有数据。已定位 T/Macro/CAF 三份 spatial.h5ad 且均含 SVC_cluster；不是 Real_application 原始文件。材料已存在，汇集基线属于后续输入准备，不再列为用户待补材料。
