# REVISE × REVISE_Analysis_Agent 当前收口任务

## 0. 任务定位

你现在需要继续落实两个仓库的协同收口：

- `REVISE`
  - 目标分支：`revise-2.0`
  - 本次审阅基线 commit：`464ee26ac22978c3994e40af50325413e70cef13`
- `REVISE_Analysis_Agent`
  - 目标分支：`sl`
  - 本次审阅基线 commit：`3e70ab364f89f9e24e981e3e9239e7fe563138a6`

开始前必须先检查当前 HEAD 和工作树；如果已经有更新，以当前工作树为准，不机械重复下面已经完成的修改，也不要覆盖用户已有未提交内容。

这轮不是重新设计两个项目，也不是继续扩大分析能力。

当前共同主线已经基本建立：

**REVISE 重建 → 正式发布 `sample.yaml + raw.h5ad + SVC.h5ad` → Analysis Agent 独立加载 → Notebook / batch / report 使用同一分析逻辑。**

现在的目标是：

> **检查并修正仍然存在的少数语义或数值问题，让已经建立的主链更加自洽、简单、可直接使用。**

------

# 1. 重要设计原则

## 1.1 做减法，不做加法

如果现有逻辑已经能满足目的，优先复用。

不要主动增加：

- 通用兼容层
- fallback / 自动猜测
- 自动补标签
- 自动修复输入
- 新的 registry / planner
- 新的缓存或任务调度框架
- 新的通用 correspondence 层
- 另一套 Notebook / report 体系
- 为了让测试“全绿”而改变科学定义

如果某个条件真的不满足：

**把问题暴露出来。**

不要通过改变输入、换算法、自动补值等方式让它看起来成功。

------

## 1.2 大规模数据才能回答的问题，本轮后置

这轮不要求通过 mini 数据获得完整生物学结论。

例如以下内容不阻塞当前收口：

- State / Gain 在全量数据上是否形成稳定区域
- 四种 assembly 哪个最终最好
- 大规模生物学解释
- 跨平台 / 跨器官 / 跨疾病结果
- 全量运行的资源和容量表现
- sST 共坐标对空间统计的真实影响
- 更完整的 uncertainty 理论
- CCI / niche / recoverability 等后续分析

这些只需要在已有文档中的“后续检查”里明确保留。

**不要为了这些问题扩大 mini ROI、降低阈值、改科学参数或新增方法。**

------

# 2. 当前已经完成的部分

先核对当前代码；如果以下内容仍成立，不要再次修改。

## REVISE 已基本完成

- cell type 主流程统一 `/ → _`
- 缺失值 NA 保留
- 不再通过 `label_aliases`
- 已删除 slash normalization 后的 collision 拒绝
- whole-sample iST 重建
- 单一 SVC 输出
- 完整 Raw 保留
- `sample.yaml + raw.h5ad + SVC.h5ad` 正式交付
- 事务发布与失败回滚
- mean / random / within_cluster / outside_cluster 四种 assembly
- H2 历史 baseline 汇集
- H3 四 assembly mini 真实执行
- H4 exact ID / gene / coordinate 检查
- hST / sST mini 正式重建与交付
- full configs 和服务器执行入口

## Analysis Agent 已基本完成

- 删除主流程 `label_aliases`
- cell type `/ → _`
- 默认 scope：
  - `All`
  - `Fibroblast`
  - `Mono_Macro`
  - `T`
- Raw Level2 按 scope 内有效 ID 使用
- Raw Level2 coverage 输出
- 无 `SVC_cluster` 时 K-control 等相关能力按 unavailable 处理
- 意外运行错误仍记录为真实 error
- Notebook / batch / report 使用同一 workflow 和有效参数
- 项目配置直接引用 REVISE 正式 `sample.yaml`
- mini iST / hST / sST 均已实际消费
- report 只读取已保存结果，不重新计算

如果这些已经在当前工作树中成立，视为完成项，不继续“优化”。

------

# 3. 当前真正需要继续处理的事项

本轮重点只有少数几项。

------

# A. Analysis Agent：Anatomy 标签来源

## 背景

当前生产端完整保留 Raw。

REVISE 的 `revise_Level1` 是推断结果，只对实际经过对应推断的 Raw 单位回填，因此合法情况下会存在 NA。

Analysis Agent 当前 Anatomy 构建路径仍大致是：

```python
raw_broad = sample.labels("raw", sample.broad_key)
```

然后用这组标签构建完整 Raw Anatomy。

这可能存在一个语义问题：

**Anatomy 的目的本来是描述完整 Raw tissue 的组织背景，而不是描述“哪些单位成功获得了重建 broad inference”。**

因此当前应先确认：

> Anatomy 应该读取哪一个 Raw 标签来源？

尤其检查完整 Raw 中是否已经保留原始人工 `Level1`，以及当前 `sample.broad_key` 是否实际指向 `revise_Level1`。

## 要做的事情

1. 沿正式 REVISE delivery → sample.yaml → Analysis Agent `Sample` → `stage_support()` 检查实际 Anatomy 使用的列。
2. 根据现有设计目的，判断 Anatomy 是否应该明确使用完整 Raw 的原始 Level1。
3. 如果是：
   - 只在 Anatomy 这个局部调用点修正来源。
   - 不改变 Reconstruction Impact 中其他 scope / broad routing 的列。
   - 不把原始 Level1 和 `revise_Level1` 自动混合填充。
4. 使用一个很小的 fixture 验证：
   - 原始 Level1 完整；
   - `revise_Level1` 存在部分 NA；
   - Anatomy 仍按照设计指定的来源工作。

## 不要做

不要先实现：

- 通用 missing-label Anatomy 系统
- 自动用原始标签补 `revise_Level1`
- 多标签来源 fallback
- 根据“哪列更完整”自动选择
- 全局放宽所有空间窗口函数对 NA 的约束

如果正确标签来源修正后问题自然消失，就到此为止。

如果正确来源本身仍有缺失，只把真实缺失暴露出来，再另行讨论。

------

# B. Analysis Agent：Notebook 正式默认入口收口

## 背景

当前正式 mini/full 分析已经通过 project YAML 直接引用 REVISE 发布的 `sample.yaml`。

但 `notebooks/01_reconstruction_impact.ipynb` 仍保留历史默认：

```python
data/P2CRC_Xenium/sample.yaml
```

该文件还是旧的临时 P2 carrier 入口。

因此现在存在：

- 验收脚本使用新正式主链
- 用户直接打开 Notebook 时仍可能进入旧路径

这不符合项目已经迁移到正式交付的状态。

Notebook 中也可能仍有类似：

> 样本必须包含完整 SVC reconstruction label

这样的旧表述，而当前 hST / sST 已允许无 `SVC_cluster` 时按能力 partial 运行。

## 要做的事情

1. 收敛 Notebook 的默认使用方式到当前正式主链。
2. 最简单的方案优先：
   - 默认明确使用一个正式 project config；
   - 或明确要求用户设置 project，而不是默默落到旧临时 sample。
3. 保留：
   - `PROJECT_YAML`
   - `SAMPLE_YAML`
   - override
     这些现有能力。
4. 删除或修改与当前实际行为矛盾的旧说明。
5. README 与 Notebook 首部保持一致。

## 不要做

不要增加：

- 自动扫描项目
- 自动发现最近的 sample
- 缺文件时 fallback 到 example
- 多套默认策略
- 第二本“正式 Notebook”

目标只是：

> **日常入口和已经验收的正式路径一致。**

------

# C. REVISE：sST parent–gene 校正式

## 背景

当前 sST 最终校正代码大致是：

```python
current_sum = ...
ratio = X / (current_sum + 1e-10)
SVC_X = SVC_X * ratio[spot_indices]
```

当前 mini 验收已经发现：

- 大多数正支持 parent–gene 项守恒良好
- 仍有一部分正目标最终聚合接近零
- 当前状态被诚实记录为 partial

这里要把两类问题分开。

### 1. 正支持但非常小

如果 `current_sum > 0`，加入固定 `1e-10` 会使校正结果偏离严格的目标比例。

这属于可以直接从公式层面判断的问题，不需要大规模数据。

### 2. 真正零支持

如果某个 parent–gene：

```text
target > 0
current_sum == 0
```

纯乘法校正无论如何都不能恢复该表达。

这属于算法语义问题，不能由 Codex 自行选择一种补偿方式。

## 要做的事情

1. 检查校正前的 `current_sum`，不要仅根据最终输出 ≤1e-12 就把所有情况叫 zero support。
2. 区分：
   - strictly positive support
   - true zero support
3. 对 **positive support**：
   - 明确当前 epsilon 是否只是为了避免除零。
   - 如果是，只对正支持项使用不会引入系统缩放误差的比例计算。
4. 对 **true zero support**：
   - 不新增自动分配策略。
   - 保留并明确报告。
5. 更新已有验证代码，使报告明确区分：
   - positive-support residual
   - zero-support unresolved target
6. 用极小 synthetic case 验证：
   - 普通正值
   - 极小正值
   - 真零

## 不要做

不要自行实现：

- 给零支持 gene 加 pseudocount
- 均分给所有虚拟细胞
- nearest-cell 补值
- reference 表达补偿
- 修改 OT
- 新的 imputation 策略

这些都属于后续科学决策。

本轮目标只是：

> **把确定的数值问题修正，把真正没有定义好的问题暴露出来。**

------

# D. 来源声明：只核对，不开发新方法

目前仍有两个明确未完全确认的输入事实。

## hST expression identity

当前 P1 HD 的表达历史仍为 unknown。

不要根据：

- 非负
- 整数
- 数值范围
- 看起来像 counts

自动把它升级成 raw_counts。

要做的只有：

1. 检查已有数据来源或 preprocessing 记录；
2. 能确认则更新正式配置；
3. 不能确认则继续保持 unknown。

unknown 是有效结果，不是错误。

------

## sST physical scale

当前 P2 Visium 配置使用：

```yaml
microns_per_coordinate: 0.73
```

但现有注释说明该值仍 provisional。

这里只需要核对已有来源。

- 能确认 → 保留并在文档中说明来源。
- 不能确认 → 不要把依赖这一值的物理尺度解释写成已确认。

不要新开发尺度估计算法。

------

# 4. 当前不要求解决的内容

以下内容只更新已有“remaining / future checks”文档，不开发。

## 后续真实数据检查

- State/Gain 在大规模数据上能否形成稳定区域
- region threshold 的真实稳定性
- State × Anatomy 的生物学解释
- EMT / Moran / pathway 的真实结论
- 四种 assembly 的最终科学选择
- K-control 对真实结论的影响
- 全量 P2 Xenium / HD / Visium 的最终结果
- 全量资源和峰值内存
- sST 共 parent 坐标对 Moran 等空间统计的科学影响
- Raw/SVC `All` 覆盖范围不同对解释的影响
- sST zero-support 的跨样本规模和科学后果

## 后续方法方向

- unified uncertainty
- gene-wise uncertainty
- online CELLxGENE
- out-of-spot
- 通用 Raw↔SVC correspondence
- CCI / niche / cross-platform正式流程

本轮不要因为看到这些 TODO 就开始实现。

------

# 5. 关于 mini 验收应该如何理解

当前 mini 的作用主要是：

### 可以证明

- 真实数据能够经过当前 REVISE 路线运行
- 正式三文件可以发布
- Analysis Agent 可以直接消费
- 不同 route 的能力限制不会被错误兜底
- Notebook/batch/report 使用同一套逻辑
- 四种 assembly 比较工具链实际能执行

### 不要求证明

- 小 ROI 一定能找到稳定 State/Gain region
- mini 数据能支持生物学结论
- 四方法 mini ARI/NMI 能决定最终方法
- hST unknown 表达能做分子分析
- 全量数据一定能在当前资源限制下运行

因此，不要为了使 mini 结果更“完整”而：

- 降低 200-window threshold
- 改 window size
- 改 min_window_units
- 选择更有利的 ROI
- 改科学参数
- 自动填补标签

mini 如果因为支持不足而 unavailable，这是合理结果。

------

# 6. 测试要求

不要新增大规模防御测试。

优先复用现有测试。

只针对实际修改补最小测试。

## Analysis Agent

### Anatomy

至少覆盖：

```text
raw original Level1 完整
revise_Level1 局部缺失
→ Anatomy 使用确定的正确来源
```

### Notebook

检查：

```text
正式 project
→ 正式 REVISE sample.yaml
→ 与 batch 使用相同 effective parameters
```

不需要重新运行全部真实重建来测试 Notebook 路径修改。

------

## REVISE

### sST correction

最小数值测试：

```text
case 1: positive ordinary support
case 2: extremely small positive support
case 3: true zero support
```

需要明确检查：

- 正支持分支是否满足预期守恒
- 真零是否被明确识别为 unresolved
- 不产生 NaN / inf
- 不重新引入 final per-cell 10000 scaling

除此之外，不扩展测试矩阵。

------

# 7. 文档收口

本轮完成代码修改后，只更新现有文档。

不要创建新的 review / handoff 文档体系。

建议只改现有：

## REVISE

```text
docs/development/reconstruction-analysis/
    acceptance.md
    cross-repo-review.md
```

若 sST 校正状态发生变化，同时更新机器可读验证记录或相应 verifier 输出定义。

## Analysis Agent

```text
docs/cross-repo-review/evidence-and-gaps.md
docs/input-output.md
docs/analyses/reconstruction-impact.md
notebooks/README.md
```

仅修改受本轮实际变化影响的部分。

------

# 8. 实施顺序

建议按这个顺序。

## Step 1：先核对当前工作树

输出简短表格：

| 项目                           | 当前是否仍存在 | 是否需要修改 |
| ------------------------------ | -------------- | ------------ |
| Anatomy source                 |                |              |
| Notebook old default           |                |              |
| Notebook outdated cluster text |                |              |
| sST positive-support epsilon   |                |              |
| hST expression source          |                |              |
| sST scale source               |                |              |

已经解决的不要重做。

------

## Step 2：完成 Analysis Agent 的两个收口

1. Anatomy label source
2. Notebook 正式默认入口和文案

跑聚焦测试。

------

## Step 3：完成 REVISE 的 sST 数值核对

1. 区分 positive / zero support
2. 修正能明确判断的 positive-support 计算
3. zero-support 保留 unresolved
4. 更新 verifier
5. 跑聚焦测试

------

## Step 4：核对已有 provenance

检查 hST expression、sST scale。

没有新证据就不要改声明。

------

## Step 5：更新现有文档

最终状态分成：

```text
Engineering completed
Current local fixes completed
Known unresolved input/method facts
Deferred large-data/scientific checks
```

------

# 9. 明确禁止的实现方式

以下行为都不要做：

- 自动补 Raw 标签
- 自动从原始 Level1 填 `revise_Level1`
- 自动根据列完整程度选 Anatomy source
- 自动猜 expression identity
- 自动猜 physical scale
- 自动降低 region threshold
- 自动调整 scientific parameters 让 mini 成功
- 为 zero-support expression 自动补值
- 自动扩大 OT memory limit
- 自动缩 reference
- 自动 fallback 到其他 assembly
- 新增 global compatibility layer
- 新建第二套分析入口
- 新建第二套 report
- 把四方法比较复制到 Analysis Agent
- 为后续 uncertainty / CCI / niche 提前写框架

------

# 10. 最后汇报格式

完成后不要只说“all tests passed”。

请给出以下四张表。

## A. 实际修改

| 仓库 | 文件/符号 | 修改内容 | 为什么 |
| ---- | --------- | -------- | ------ |
|      |           |          |        |

## B. 删除或简化

| 原逻辑 | 当前如何简化 | 是否改变科学定义 |
| ------ | ------------ | ---------------- |
|        |              |                  |

## C. 实际验证

| 场景 | 命令/测试 | 结果 | 能证明什么 | 不能证明什么 |
| ---- | --------- | ---- | ---------- | ------------ |
|      |           |      |            |              |

## D. 剩余问题

| 问题 | 类型                                             | 当前为何不解决 | 后续需要什么证据 |
| ---- | ------------------------------------------------ | -------------- | ---------------- |
|      | implementation / input fact / scientific / scale |                |                  |

最后明确回答：

1. 当前两个仓库的正式主链是否仍然成立？
2. 本轮有没有发现新的真实接口错位？
3. 哪些问题已经通过做减法解决？
4. 哪些问题被明确保留，而没有用防御性实现隐藏？
5. 是否新增了任何不必要的框架；如果有，删除它。
