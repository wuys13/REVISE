# Assignment Guidance 行为标准

## 行为画像

所有 SVC 路线在 Global Anchoring（GA）结束后共享同一个最小 assignment 契约：完整 posterior 矩阵 `Q`，以及由 `argmax(Q)` 得到的硬标签。Assignment Guidance 只负责验证这两个输入，并按照 SVC 路线决定 `Q` 是否以及如何影响后续局部步骤；它不再承担通用状态管理、provenance、fallback 或任意 assignment 变换。

三条路线的目标行为如下：

- **sp-SVC**：先使用 GA 硬标签完成分组；在组内局部邻居图的支持上，用 `Q` 的 posterior compatibility 调整局部 OT cost，最终影响 transport coupling 和表达平滑。
- **sc-SVC-sr**：保留 SR 路线自身的 virtual-cell allocation；将 spot-level `Q` 投影到对应 virtual cells，并在 virtual-cell 局部表达平滑阶段调整局部 OT cost，最终影响 transport coupling。virtual-cell 的分组标签仍由 SR allocation 逻辑决定，不强制等同于投影后 `Q` 的 argmax。
- **sc-SVC**：GA 的 `Q` 当前只通过 `argmax(Q)` 参与 Level1 分流；完整 `Q` 保留在 GA 输出中，但不参与后续 GraphCluster 边重加权。每个 Level1 分组继续独立执行既有 local refinement、聚类和 annotation。

当前固定采用的 posterior compatibility 为：

\[
A_{ij}=Q_iQ_j^\top
\]

局部 OT 使用当前 cost-conditioning 语义：

\[
C'_{ij}=C_{ij}+\lambda\left[-\log\left(\max(A_{ij},\epsilon)\right)\right]
\]

因此，posterior 并非直接改写原始 adjacency，而是在已有局部图支持上调整 OT cost，继而改变 transport coupling；用户层面可将其理解为影响局部平滑所使用的有效边权。

## 标准修订点

当前实现隐含的标准是：不同路线都通过一个通用 `AssignmentState` 和统一 Assignment Guidance 策略消费 assignment，并围绕该过程维护 availability policy、fallback、运行阶段状态和证据字段。

修订后的标准是：

1. 只统一 GA 的 assignment 产物，不统一各路线对 posterior 的下游使用方式。
2. Assignment Guidance 是一个小而明确的路线分派入口，允许直接使用 `if / elif / else` 表达三种封闭行为。
3. 核心 assignment 数据只包含硬标签与完整 posterior `Q`。索引和 category 名称由 `Series`/`DataFrame` 自身承载。
4. `source`、`level`、`lineage`、hash 和运行事件属于 provenance/diagnostics，不属于算法接口。
5. GA 已承诺产生 `Q`，因此缺少或错位的 `Q` 是契约错误，不允许从硬标签静默重建 one-hot posterior。
6. 第一版不保留通用的 `off | prefer | require` 和 `cost | reference` 组合；sp-SVC 与 sc-SVC-sr 固定使用 local-OT conditioning，sc-SVC 固定为 store-only/no-op。可调参数只保留局部 OT guidance strength `lambda`。

## 最小接口

概念上的最小数据结构为：

```python
@dataclass(frozen=True)
class AssignmentState:
    labels: pd.Series
    posterior: pd.DataFrame
```

统一入口的概念行为为：

```python
def apply_assignment_guidance(route, assignment, local_ot=None, *, strength):
    validate_assignment(assignment)

    if route == "sp_svc":
        return condition_local_ot(local_ot, assignment.posterior, strength)
    elif route == "sc_svc_sr":
        return condition_local_ot(local_ot, assignment.posterior, strength)
    elif route == "sc_svc":
        return local_ot  # explicit no-op; Q remains in the GA output
    else:
        raise ValueError(f"Unsupported SVC route: {route}")
```

这段代码只定义接口方向，不要求照抄具体类型或函数名。实现时应复用 sp-SVC 和 sc-SVC-sr 已有的局部 OT 数据结构，避免为统一入口新增一组宽泛的可选字段。

## 验收标准

### GA 输出契约

- 每条 SVC 路线完成 GA 后，都能取得与 observation 顺序严格对齐的 `Q`。
- `Q.index` 与目标 observation identifiers 一致，`Q.columns` 是明确且无歧义的 category labels。
- `Q` 的每个元素有限且非负，每行在允许误差内归一化。
- GA 硬标签等于对应行 `Q.idxmax(axis=1)`。
- `Q` 缺失、轴错位、category 冲突或硬标签不一致时明确失败，不进行 one-hot fallback。

### sp-SVC

- 局部 refinement 仍按 GA 硬标签分组。
- posterior compatibility 只在当前局部邻居图支持上计算。
- 当两个有效 posterior `Q` 不同时，启用非零 `strength` 可以改变局部 OT cost/coupling。
- `strength=0` 时结果与未应用 Assignment Guidance 的局部 OT 一致。

### sc-SVC-sr

- spot-level `Q` 按 virtual-cell-to-spot 映射正确投影，不产生新的独立 posterior 来源。
- virtual-cell allocation 和分组语义保持由 SR 路线负责。
- posterior 只在 virtual-cell 局部表达平滑的 OT 阶段生效。
- `strength=0` 时结果与未应用 Assignment Guidance 的局部 OT 一致。

### sc-SVC

- GA 硬标签继续决定 Level1 分流。
- 完整 `Q` 保留并可追溯，但不改变 GraphCluster edge weights、Leiden cluster 或 Level2 local anchoring。
- 移除 Assignment Guidance 后新增的 posterior graph reweighting，不改变既有的组内 local refinement 流程。

### 配置与证据

- 核心运行配置不再暴露 `off | prefer | require`、`cost | reference` 的组合状态空间。
- local-OT 路线最多暴露一个 `strength` 参数；其默认值必须由现有实验或兼容性要求明确决定。
- 运行证据只需说明路线、是否实际应用以及 `strength`；详细 hashes、lineage 和 fallback 诊断如仍需保留，应位于独立 provenance/diagnostics 输出中。

## 约束域

实现该标准时，预计需要检查和外科手术式修改以下范围：

- `revise/backend/ops/assignment.py`
- `revise/backend/ops/assignment_guidance.py`
- `revise/backend/ops/posterior_conditioning.py`
- `revise/backend/ops/sr_allocation.py`
- `revise/backend/runners/sp_svc_assignment_guidance.py`
- `revise/backend/runners/sp_svc_application.py`
- `revise/backend/runners/sc_svc_sr_application.py`
- `revise/backend/runners/sc_svc_application.py`
- `revise/backend/kernels/graph_cluster.py`
- `revise/backend/adapters.py`
- `revise/config/loader.py`
- 与 assignment guidance、三条 SVC 路线及配置兼容相关的测试和文档

该范围是重构检查边界，不表示所有文件都必须产生修改。每一处实际改动都必须能追溯到上述行为标准。

## 关键决策

- 统一点放在 GA 输出 `Q + argmax label`，而不是下游使用策略，因为三条 SVC 的 local refinement 语义本来不同。
- 保留 `Assignment Guidance` 名称，但将它收缩为显式路线分派和 posterior 应用，不再作为通用 assignment 状态平台。
- 使用直接的 `if / elif / else`，因为当前只有三个稳定且封闭的算法分支；暂不引入 adapter registry 或策略类层级。
- sp-SVC 和 sc-SVC-sr 共享局部 OT conditioning 的数学实现，但各自保留数据准备和分组责任。
- sc-SVC 暂不传播软 posterior，避免 GA posterior 或 Level2 posterior 被重复用于 Graph clustering，改变原始算法因果链。
- SR allocation 是 sc-SVC-sr 自身的必要步骤，不属于可选 Assignment Guidance。
- 不使用 one-hot fallback，因为它会丢失 posterior uncertainty，并可能让同一硬分组内的 compatibility 退化为常数 1。
