# REVISE 2.0：面向 REVISE_Analysis_Agent 的 Reconstruction 开发事项

## 1. 总览：`revise-2.0` 接下来具体要做什么

| 编号      | 开发事项                                                             | 目标                                                                         | 当前确定程度                          | 下一步                                                  |
| ------- | ---------------------------------------------------------------- | -------------------------------------------------------------------------- | ------------------------------- | ---------------------------------------------------- |
| **A1**  | **统一 reconstruction 公共输出**                                       | 所有 route 最终都向 downstream 提供一个完整 `SVC.h5ad`                                 | **基本想清楚**                       | 修改 publication / batch output contract               |
| **A2**  | **取消 sc-SVC 的 per-cell-type 输出结构**                               | reconstruction 内部可按 cell type 计算，但最终一个 sample 只输出一个 whole-tissue SVC       | **基本想清楚**                       | 把 batch task 从 `sample × cell_type` 改成 sample-level  |
| **A3**  | **sc-SVC 遍历全部 broad cell types**                                 | 不再只 `select_cell_type=T/Fibroblast/...`，而是把所有可重建 Level1 都跑完                | **事情确定，具体遍历策略待定**               | 等你提供以前遍历 cell type 的参考代码，再落实                         |
| **A4**  | **sc-SVC 统一 assembly：保留 mean/random，新增 cluster mode**            | 把 spatial carrier 与 reconstructed expression 合成一个真正可分析的 spatial `SVC.h5ad` | **事情确定；cluster 算法未定**           | 第一阶段先打通 random；再单独设计 cluster                         |
| **A5**  | **当前默认 assembly 可先改为 random**                                    | 在 cluster mode 尚未完成前，先有一个统一可运行的 SVC 输出                                     | **基本想清楚**                       | 修改默认配置及 publication                                  |
| **A6**  | **mean/random/cluster 后续做系统比较**                                  | 决定最终哪种 assembly 更合理、是否把 cluster 设为默认                                       | **事情确定；benchmark 方案待定**         | 等三种模式都可运行后再设计比较                                      |
| **A7**  | **生成 analysis-ready Raw H5AD**                                   | Raw 保留原始实测 expression，同时增加 reconstruction 已经算出的重要 annotation               | **基本想清楚**                       | 增加 Raw publication/finalization                      |
| **A8**  | **Raw 写入 Level1 / Level2 等标签**                                   | Analyst 不再重新调用 OT/TACCO 获得 Raw Level2                                      | **基本想清楚**                       | 把 GA/LR annotation 映射回完整 Raw obs                     |
| **A9**  | **保存高成本 confidence/posterior**                                   | reconstruction 已经计算出的昂贵 inference 不要让 Analyst 重算                           | **原则明确；具体保留字段待定**               | 结合 uncertainty 重构方案确定最终字段                            |
| **A10** | **Raw 保留没有进入 reconstruction 的 ST units**                         | 不为了方便 pairwise 而裁剪 Raw                                                     | **基本想清楚**                       | Raw publication 必须基于完整原始 observation axis            |
| **A11** | **标记 Raw unit 是否进入最终 reconstruction**                            | 某些需要 comparable subset 的分析可以按需筛选                                           | **基本想清楚，但字段名可后定**               | 可加入简单 boolean 标记                                     |
| **A12** | **保留真实 Raw↔SVC correspondence**                                  | 后续需要 pairwise 时有真实依据，但不强制全局 pairing                                        | **原则想清楚；具体载体待定**                | hST/iST 尽量保留 ID；sST 保留 `spot_name`；复杂关系以后再定          |
| **A13** | **sST 保留 parent spot 信息**                                        | 必要时 Analysis 可在单细胞轴构造 parent-spot equal-split baseline                     | **基本想清楚**                       | 延续 `spot_name` / `cell_id`                           |
| **B1**  | **sST out-of-spot SVC provenance**                               | out-of-spot cells 不制造假的 parent spot                                        | **事情确定；生成方式和字段待定**              | 等 out-of-spot reconstruction 本身开发时一起设计               |
| **A14** | **不额外保存普通 gene overlap provenance**                              | `Raw ∩ SVC genes` 等便宜信息 Analysis 现场计算                                      | **明确不做**                        | 不增加 `measured_in_raw` 等冗余字段                          |
| **A15** | **未来真正昂贵的 gene-wise uncertainty 要保存**                            | 避免 Analyst 重算 reconstruction uncertainty                                   | **原则明确；内容依赖后续 uncertainty 方法**  | 和新的 gene-wise reconstruction 一起设计                    |
| **A16** | **Reference Selection 作为 reconstruction 前置阶段**                   | reconstruction 不再默认由用户直接指定唯一 reference                                     | **整体架构基本想清楚**                   | 新增 candidate → score → selection → reconstruction 流程 |
| **A17** | **Reference 来源支持 By Patient ID**                                 | 从患者匹配 reference pool 中找到候选                                                 | **事情确定；数据组织细节待定**               | 后续结合实际 reference 目录设计                                |
| **A18** | **Reference 来源支持 By CELLxGENE**                                  | 自动寻找与当前 ST 匹配的公开 scRNA reference candidates                                | **事情确定；candidate discovery 未定** | 后续单独开发 CELLxGENE retrieval/filtering                 |
| **A19** | **统一两种 Reference 来源的后半段接口**                                      | 不做 Patient pipeline / CELLxGENE pipeline 两套 reconstruction                 | **基本想清楚**                       | 两者都先产出统一 selection/config artifact                   |
| **B2**  | **用 ST Unit Confidence 排序候选 Reference**                          | 判断哪个 sc reference 最适合当前 ST                                                 | **核心思想确定；score 定义未定**           | 单独研究如何从 unit confidence 聚合成 reference-level score    |
| **B3**  | **生成 reference ranking 中间文件**                                    | 保存昂贵 candidate evaluation，供后续选择和审计                                         | **基本方向明确；格式细节待定**               | 可先考虑 ranking table + selected reference config       |
| **B4**  | **Reference selection 最终选 top-1 / top-k / ensemble**             | 决定 reference 如何进入正式 reconstruction                                         | **暂未决定**                        | 第一版可先 top-1，后续根据实验决定                                 |
| **A20** | **Reference selection 写入 reconstruction provenance/fingerprint** | reference 改变必须视为 reconstruction 改变                                         | **基本想清楚**                       | 接入 2.0 fingerprint/provenance                        |
| **A21** | **简化 Analyst-facing output**                                     | Analyst 主要只看 Raw + SVC，不理解 paired / route / cell-type folders              | **基本想清楚**                       | reconstruction.json 可继续内部保留，但不成为复杂分析入口               |
| **A22** | **移除 Analyst 对 REVISE backend 的依赖**                              | Analysis repo 不需要 `revise.backend.kernels.ot` 等 reconstruction 实现          | **基本想清楚**                       | 把必要 inference 全部在 reconstruction 阶段做完并保存             |

---

# 2. 其他分支只需要这样理解

后续开发主线：

```text
revise-2.0
```

负责真正完成上述 reconstruction 改造。

`reconstruction-impact`：

* 暂时冻结；
* 不继续添加 reconstruction 功能；
* 其 analysis / Notebook / Web 能力迁入 `REVISE_Analysis_Agent`；
* 同时可作为历史结果与 parity reference。

`main`：

* 只在已有 reconstruction backend / unified pipeline 实现确实更成熟时作为参考；
* 不在这一轮和 `revise-2.0` 双线同时开发同一套 Analyst 对接功能。

---

# 3. 第一条主线：把 reconstruction 输出真正统一

## 3.1 最终公共接口

目标是每个 sample 最终至少形成：

```text
output/<sample>/
├── raw.h5ad
├── SVC.h5ad
├── reconstruction.json
└── .revise/
```

其中给 `REVISE_Analysis_Agent` 的真正核心数据产品只有：

```text
raw.h5ad
SVC.h5ad
```

Analysis 不再需要理解：

```text
sp-SVC / sc-SVC / sc-SVC-sr
paired
expr.h5ad
spatial.h5ad
T/
Fibroblast/
Mono/
```

这些都应该终止在 reconstruction 内部。

---

# 4. sc-SVC：当前最核心的改造

## 4.1 从单一 cell type reconstruction 改为 whole-sample reconstruction

现在的逻辑本质是：

```text
select_cell_type
        ↓
只重建一个 broad cell type
        ↓
输出该 cell type 的结果
```

以后应该：

```text
ST sample
   ↓
确定需要重建的 broad cell types
   ↓
逐个执行 sc-SVC LR
   ↓
合并全部 reconstructed spatial cells
   ↓
whole-sample SVC.h5ad
```

也就是说：

> cell type 可以继续作为内部计算分块，但是不能继续作为目录和最终文件分块。

这也是未来做 whole-tissue：

```text
CCI
TLS
Niche
Treg ↔ Mreg
CAF ↔ Tumor
```

的前提。

---

## 4.2 “遍历哪些 cell types”目前暂不锁死

这件事确定必须做，但具体策略还需要实际开发时讨论。

可能要考虑：

```text
GA 中存在的 Level1
∩
Reference 中有足够细胞的 Level1
```

以及：

```text
最少 spatial cells
最少 reference cells
Unknown
特殊 broad classes
```

你之前应该已经有遍历所有 cell type 的参考代码。

因此这里的具体行动是：

> **等实际开始 A3 时，你提供之前相关实现，然后在 `revise-2.0` 当前 pipeline 上重构，而不是现在重新猜一套。**

---

# 5. sc-SVC：统一 spatial expression assembly

当前 sc-SVC 的核心历史问题是：

```text
spatial cells
```

和：

```text
完整 reconstructed expression
```

不在同一个 observation axis。

最终必须得到：

```text
SVC.h5ad

obs = spatial SVC cells
X = 对应的 reconstructed expression
obsm["spatial"] = cell locations
```

这是 Analyst 能真正只读一个 SVC 文件的关键。

---

# 6. Assembly 第一阶段保留三种模式

明确保留：

```text
mean
random
cluster
```

其中：

### mean

cluster/reference expression 的平均表示映射回 spatial cells。

### random

从对应 cluster 的 reference/reconstructed donors 中随机选择表达映射给 spatial cell。

### cluster

新增模式。

核心思想是：

> 利用 reconstructed cluster 以及空间连续性，把 cluster-level reconstructed expression 合理分配到空间中的实际 SVC cells。

---

# 7. 当前工程推进顺序：先 random，再 cluster

不要因为 cluster mode 尚未设计完成而阻塞整个统一输出。

可以先：

```text
遍历所有 cell types
        ↓
random assembly
        ↓
merge
        ↓
SVC.h5ad
```

把完整链路跑通。

因此当前：

```text
default = random
```

可以作为阶段性默认。

之后再：

```text
mean
vs
random
vs
cluster
```

系统 benchmark。

---

# 8. Cluster mode 是确定要做，但当前最明显的“待设计事项”

现在只固定目标，不固定算法。

目标：

```text
cluster/state reconstruction
        +
spatial continuity
        ↓
cluster → spatial-cell assignment
        ↓
cell-level reconstructed expression
```

但以下问题都留待开发时讨论：

* spatial continuity 怎么定义；
* neighbor graph 是否参与；
* cluster 是 hard assignment 还是 soft assignment；
* 一个 spatial state 中多个 donor 如何分配；
* 是否继续用 cluster mean；
* 是否利用 posterior；
* 是否最终直接和新的 gene-wise uncertainty reconstruction 统一。

这部分应该单独作为一个开发专题。

---

# 9. Raw H5AD：真正要保存的是“昂贵 inference”

Raw 的目标不是复制一个干净输入。

应该是：

```text
原始实测 ST expression
+
REVISE 已经计算出的 downstream 重要 annotation
```

其中 `X`：

> 保持真实 measured expression。

不要变成：

* QC 后矩阵；
* reference-overlap genes；
* normalized；
* log-transformed；
* reconstruction 临时 carrier。

---

# 10. Raw 第一阶段明确需要增加什么

至少：

```text
Level1
Level2
```

如果 reconstruction 已经计算：

```text
ST Unit Confidence
posterior probabilities
```

而后续 analysis / uncertainty 需要，则应该一起保存。

原因是：

> 这些信息重算意味着再次运行昂贵的 reference mapping / OT / annotation，而且容易因参数变化造成 Raw baseline 漂移。

因此 Raw Level2 以后不能再由 Analyst 重算。

---

# 11. Raw 仍然保留全部原始 ST units

即使：

```text
100,000 raw units
```

最后只：

```text
92,000 reconstructed
```

Raw 仍然应该是 100,000。

不要为了方便 comparison 裁掉低质量、未重建单位。

因为很多分析只比较：

```text
distribution
state structure
pathway
spatial pattern
```

不需要 observation-level pairing。

---

# 12. 可以增加一个很轻的 reconstruction coverage 标志

例如概念上：

```text
raw.obs["revise_reconstructed"]
```

表示该 Raw unit 最终是否进入 reconstruction。

这样真正需要 matched subset 时：

```python
raw[raw.obs["revise_reconstructed"]]
```

即可。

具体字段名以后可以再定。

如果 filtering reason 很容易获得，可以保留；否则第一阶段没必要增加复杂 reason ontology。

---

# 13. Raw↔SVC correspondence：保存真实关系，但不要统一强制成一个字段

总体原则已经确定：

> 没必要把所有技术强制成 universal `source_unit_id`。

但是 reconstruction 本身已经计算出来的昂贵真实映射不能丢。

---

## hST / iST

如果可以：

```text
SVC.obs_names
```

继续保留原 spatial unit ID。

那这是最简单、最好的 correspondence。

Raw 中没有进入 reconstruction 的单位自然没有 SVC counterpart。

---

## sST

继续保留：

```text
spot_name
cell_id
```

让 Analysis 在真正需要的时候知道：

```text
一个 raw spot
→ 哪些 reconstructed cells
```

---

## 更复杂的 correspondence

如果未来确实不能用 obs columns 表达，再考虑：

```text
mapping table
```

但目前不要提前设计文件格式。

原则只有：

> **已经计算出的真实 correspondence 保存；便宜且可推导的 pairing 不提前形式化。**

---

# 14. sST 的 out-of-spot 是一个明确待办，但尚未完全设计

未来可能会 reconstruct：

```text
out-of-spot single cells
```

这些天然没有 parent raw spot。

因此不能：

```text
out-of-spot
→ nearest spot
```

人为制造 relation。

概念上应该类似：

```text
spot_name = missing
is_out_of_spot = true
```

但具体怎么生成 out-of-spot、怎么保存 provenance，等 reconstruction 本身实现这部分时再一起决定。

---

# 15. Gene provenance：这一轮明确不要过度保存

例如：

```text
Raw 和 SVC shared genes
measured genes
reconstructed-only genes
```

Analysis 端直接：

```python
raw.var_names.intersection(svc.var_names)
```

即可。

没有理由为了它增加新的 persistence contract。

因此第一阶段：

> 不做 `measured_in_raw` 等普通 gene provenance。

---

# 16. 但 gene-wise uncertainty 是另外一回事

如果后续 REVISE 真正计算得到：

```text
gene-specific reconstruction uncertainty
gene posterior
```

这是昂贵的 reconstruction inference。

必须保存。

所以原则是：

```text
便宜的 derivation
→ 不保存

昂贵的 inference
→ 保存
```

---

# 17. 第二条主线：Reference Selection

这是 reconstruction 端另一块独立的大开发。

目标从：

```text
用户指定一个 reference
→ reconstruction
```

转成：

```text
ST
→ candidate references
→ reference evaluation
→ reference selection
→ reconstruction
```

---

# 18. 至少支持两种 candidate 来源

## 18.1 By Patient ID

输入 patient/sample identity 后，从已有 reference pool 中找到：

```text
patient-matched reference candidate(s)
```

如果只有一个，也仍然建议走统一 selection artifact。

如果多个，则可以进一步比较。

---

## 18.2 By CELLxGENE

从 CELLxGENE 中按照：

```text
species
tissue
disease
...
```

获得一批候选 reference。

然后再评价它们和当前 ST 的匹配程度。

具体 CELLxGENE retrieval / metadata filtering / sampling 后续再设计。

---

# 19. 两种 Reference 来源最终必须汇合

后半段不要维护两套 reconstruction。

应该：

```text
By Patient ID ──┐
                ├→ candidate references
By CELLxGENE ───┘
                       ↓
              Reference Evaluation
                       ↓
              Reference Selection
                       ↓
               unified config
                       ↓
                Reconstruction
```

所以 reconstruction 主体最终只消费：

```text
selected reference
```

而不知道它是怎么找到的。

---

# 20. Reference Evaluation 的核心信号：ST Unit Confidence

已经确定的核心思想：

> 对每一个 candidate reference，运行第一阶段必要的 ST-unit mapping / confidence inference，用得到的 ST Unit Confidence 衡量它和当前 ST 的匹配程度。

概念上：

```text
reference A
↓
ST Unit Confidence distribution
↓
reference score

reference B
↓
ST Unit Confidence distribution
↓
reference score
```

然后排序。

---

# 21. Candidate evaluation 的结果应该持久化

因为：

```text
N candidate refs
×
ST Unit Confidence
```

计算代价不低。

因此值得形成一个中间 artifact。

例如：

```text
reference_ranking.csv
```

具体 schema 后面再定。

然后形成：

```text
reference_selection.yaml/json
```

供正式 reconstruction 消费。

---

# 22. Reference selection 当前最主要的待设计问题

这几个要记下来，之后单独讨论：

### 如何从 ST Unit Confidence 得到一个 reference-level score？

可能涉及：

```text
mean
median
low-confidence fraction
entropy
coverage
unknown fraction
```

现在还没必要选。

### reference size 是否会导致不公平？

大 reference 是否天然更容易得到高 confidence，需要实测。

### 是否需要按 cell type 分开评价？

可能出现：

```text
overall reference A 最好
Fibroblast reference B 更好
```

目前是否允许这种情况尚未决定。

### CELLxGENE candidate 数量如何控制？

需要考虑：

```text
dataset filtering
cell sampling
batch effect
download/cache
```

后续再定。

### top-1 还是 top-k？

第一版最简单可以：

```text
select top-1
```

以后是否 ensemble 暂时不定。

---

# 23. Reference selection 必须进入 provenance/fingerprint

这一点基本明确。

如果：

```text
reference A → reference B
```

即使 ST 和其他参数不变，

也必须视为：

```text
different reconstruction
```

因此至少：

```text
selected reference identity
selection artifact
```

需要纳入 reconstruction provenance/fingerprint。

---

# 24. `revise-2.0` batch 结构也要跟着变化

当前 iST 类似：

```text
sample
├── T/
├── Macro/
└── Fibroblast/
```

以后改成：

```text
sample/
├── raw.h5ad
├── SVC.h5ad
└── reconstruction.json
```

内部仍可：

```text
T LR
Macro LR
Fibroblast LR
...
```

但 batch task 是：

```text
sample
```

而不是：

```text
sample × cell type
```

---

# 25. Analyst-facing 和 reconstruction-internal 文件要区分

## Analyst-facing

核心：

```text
raw.h5ad
SVC.h5ad
```

未来复制到：

```text
REVISE_Analysis_Agent/data/<sample>/
```

即可。

---

## REVISE internal

可以继续保存：

```text
reconstruction.json
reference_ranking.csv
reference_selection.yaml
effective config
logs
runtime provenance
internal mapping/debug outputs
```

这些不需要全部成为 Analyst 的 input contract。

---

# 26. 推荐的实际开发顺序

为了避免所有任务互相卡住，我建议 `revise-2.0` 按下面顺序推进。

## 第一阶段：打通统一 sc-SVC 输出

先：

```text
遍历所有 cell types
↓
random assembly
↓
merge
↓
whole-sample SVC.h5ad
```

这是最重要的第一目标。

---

## 第二阶段：生成 analysis-ready Raw

把：

```text
Level1
Level2
必要 confidence/posterior
```

映射回完整 Raw。

然后正式输出：

```text
raw.h5ad
```

---

## 第三阶段：完成新的 cluster assembly

在 whole-sample random pipeline 已经稳定的情况下，单独开发：

```text
cluster + spatial continuity
```

避免同时调试：

```text
多 cell-type
+
merge
+
cluster assignment
+
publication
```

四个问题。

---

## 第四阶段：比较 mean/random/cluster

统一输入和 output contract 后做正式 benchmark。

再决定最终 default。

---

## 第五阶段：Reference Selection

先建立：

```text
candidate input
→ ST confidence evaluation
→ ranking
→ selected reference artifact
```

然后再分别接：

```text
By Patient ID
By CELLxGENE
```

candidate discovery。

这部分与 assembly 相对独立，可以平行开发。

---

# 27. 后续真正需要你补材料/代码的地方

目前至少有下面几项，等实施时最好你提供现有材料，而不是重新设计。

### ① 遍历所有 cell types 的旧代码

你提到之前已经有类似实现。

开发 A3 时应优先复用这个思路。

### ② sc-SVC cluster assignment 的已有想法或实验

如果之前已经做过：

```text
cluster
→ spatial continuity
→ spatial assignment
```

相关代码、Notebook 或实验结果，开发 cluster mode 时需要一起看。

### ③ Reference pool / Patient ID 当前实际组织方式

决定：

```text
By Patient ID
```

怎么 discovery。

### ④ ST Unit Confidence 当前真实实现和输出

Reference ranking 不能只按抽象定义设计，要基于当前 GA posterior / confidence 的真实数值行为。

### ⑤ CELLxGENE 后续希望调用的数据入口

例如：

```text
CELLxGENE Census
下载 h5ad
预先缓存 reference pool
```

到底是哪种路线，需要之后单独定。

---

# 28. 当前可以直接交给 Codex 实现的部分

不需要等待更多科学设计的：

```text
1. 将 revise-2.0 reconstruction task ownership 从 iST per-cell-type 向 sample-level 重构。

2. 为 sc-SVC 建立“遍历 cell type 后 merge”的 orchestration 骨架；
   具体 eligible cell-type strategy 先保持简单或留接口。

3. 保留 mean/random assembly；
   第一阶段让 random 成为默认 unified output。

4. 所有 application route 最终 publication 收敛到 SVC.h5ad。

5. 新增 analysis-ready raw.h5ad publication。

6. 将已经存在的 Level1/Level2 annotation 写回完整 Raw observation axis。

7. 移除 Analyst 对 paired spatial/expr output 的依赖。

8. 保留 sST spot_name / cell_id 等已有真实 correspondence。

9. 不增加普通 gene-overlap metadata。

10. 为未来 Reference Selection 预留一个 reconstruction 输入点：
    reconstruction 接收“已经选好的 reference artifact/config”，
    而不是把 candidate discovery 逻辑直接耦合进 GA/LR。
```

---

# 29. 目前只应该建接口/TODO，不应该让 Codex自己猜算法的部分

```text
1. 所有 eligible Level1 cell types 的最终选择策略。
2. cluster assembly 的具体空间连续性算法。
3. mean/random/cluster 最终 benchmark 和默认选择。
4. out-of-spot reconstruction 与 provenance。
5. 最终保存哪些 posterior / uncertainty。
6. ST Unit Confidence → reference suitability score。
7. CELLxGENE candidate discovery / filtering / sampling。
8. patient reference pool 的真实组织与选择策略。
9. top-1 / top-k / ensemble reference。
10. 复杂 Raw↔SVC mapping 是否需要独立 sidecar object。
```

这些都应该：

> **等科学问题和现有代码更充分后再开发，不让 Codex自行补设计。**

---

# 30. 最终的 `revise-2.0` 目标状态

可以浓缩成：

```text
              ST
               │
               ▼
      Reference Selection
      ├─ Patient ID
      └─ CELLxGENE
               │
       selected reference
               │
               ▼
          Global Anchoring
               │
        Level1 / confidence
               │
               ▼
   per-broad-cell-type sc-SVC LR
               │
         Level2 / cluster
               │
               ▼
        Spatial Assembly
      mean / random / cluster
               │
               ▼
       whole-sample SVC.h5ad

同时：

original ST
    +
GA/LR昂贵 annotation
    ↓
analysis-ready raw.h5ad

最后：

raw.h5ad + SVC.h5ad
            │
            ▼
REVISE_Analysis_Agent
```

---

## 核心原则

整个这一轮 `revise-2.0` reconstruction 开发可以用一句话概括：

> **把 reconstruction 内部复杂性留在 REVISE 内部，只把 downstream 真正无法低成本重新获得的 inference 和一个统一的 whole-sample SVC 数据产品交给 Analyst；任何 Analysis 可以便宜推导的信息都不要为了“协议完整”额外持久化。**

这样既满足后续 `REVISE_Analysis_Agent`，又不会把 reconstruction output 做成一个越来越复杂的元数据仓库。
