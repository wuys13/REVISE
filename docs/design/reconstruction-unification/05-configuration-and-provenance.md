# Configuration and Provenance Contract

Parent index: [Reconstruction Unification Design Package](README.md)

## Single engine authority

Decided target：包内一个 typed authority 定义：

- `ENGINE_DEFAULTS`：所有正式、会影响结果的 engine 默认；
- `ROUTES`：Application/Benchmark route → profile/task/SVC kind/strategy；
- `LOCKED_KEYS`：用户层不得通过通用 override 改写的路线/solver/内部字段。

runner 构造前必须得到完整 effective config。pipeline 不再加载默认 YAML，adapter/runner dataclass defaults 不得成为隐藏的结果参数。Application preprocessing 在三层保持**有效行为一致**：函数默认是 spatial 60/100、reference None/100；官方 YAML 显式写出同样结果；public `run_application` override 默认 `None` 只表示“不覆盖 YAML”，而不是另一份结果默认。

## Compile and precedence

正式优先级：

```text
typed ENGINE_DEFAULTS
→ selected route/profile defaults
→ validated route YAML
→ explicitly supplied public Application overrides
→ locked-key validation
→ effective config + hash
```

Benchmark YAML 本身是唯一 route selector；删除 `--confounding`，launcher 为每类任务传对应 YAML，不做旧/新格式 auto-detection。

Public preprocessing overrides 不使用 sentinel。特别是 reference `min_transcript_counts: null` 由 YAML/config 承载并表示关闭该过滤；public override 的 `None` 不能同时表示“显式覆盖为 null”。因此自定义 YAML 若写了非空 threshold，调用者应修改 YAML，而不是用 public `None` 关闭它。

## Official YAML inventory

Application，恰好 5 个：

- `Xenium_T.yaml`
- `Xenium_Fib.yaml`
- `Xenium_Mono.yaml`
- `VisiumHD.yaml`
- `Visium.yaml`

Benchmark，恰好 6 个：

- `segmentation.yaml`
- `bin2cell.yaml`
- `batch_effect.yaml`
- `spot_size.yaml`
- `gene_panel.yaml`
- `gene_dropout.yaml`

不创建 `noise.yaml` 或 `imputation.yaml`。gene panel/dropout 分别包含运行自身所需的完整 imputation 配置。repo-facing templates 与 packaged templates 必须 byte-exact。

## Xenium application target

三个 Xenium YAML 使用统一 reference：

```yaml
inputs:
  reference:
    path: raw_data/Real_application/adata_sc_all_reanno.h5ad
    format: h5ad
    filter_column: Patient
    filter_value: P2CRC
```

Xenium_T 还固定：

- spatial path：`raw_data/Real_application/P2CRC_Xenium.h5ad`
- preprocessing：spatial 60/100，reference null/100
- GA broad column：`Level1`
- LR subtype/select：`Level2` / `T`
- graph alpha：`0.2`
- resolutions：`[0.6, 0.7, 0.8]`
- output dir：`output/P2CRC_Xenium/T`
- omit `output.name`

Xenium_Fib、Xenium_Mono 使用相同 filter/schema，output dir 分别为 `output/P2CRC_Xenium/Fib` 与 `output/P2CRC_Xenium/Mono`。Visium/VisiumHD 的既有行为除显式迁移的默认外不改变。

## `revise.yaml` migration

先逐字段建立旧参数 → typed authority/Application YAML/Benchmark YAML 的迁移表，再删除 `revise/revise.yaml`。禁止先删除再依靠 adapter/dataclass fallback 猜回参数。

需要迁移的结果参数至少包括：

- GA/LR solver 与 regularization；
- graph method/neighbors/alpha/resolutions/random state；
- local refinement strength；
- SR allocation/graph aggregation；
- imputation、prune 和 subcluster；
- preprocessing thresholds；
- runtime seed 与 route identity；
- I/O layout 中仍属于正式 contract 的字段。

## Provenance

每次 run 分开记录：

| Evidence | Required fields |
| --- | --- |
| Request config | source label/path、raw SHA-256、explicit overrides |
| Engine authority | schema/package version、authority SHA-256 |
| Effective config | canonical serialization、effective SHA-256、selected route/profile |
| Runtime | effective seed、resolved inputs、output roles、solver requested/completed |

配置 hash 必须包含所有影响结果的 Application preprocessing、alpha、resolutions 和 Benchmark leaf seed。不得把 output filename 当作算法 identity，也不得用 provenance metadata 反向修改结果。

Seed 范围冻结如下：runtime seed 只控制全局 RNG、显式抽样、SR allocation
和 Benchmark per-case seed。sc-SVC 的 PCA/neighbors/Leiden 使用深层
`graph.random_state: 0`，因为旧 `test_sc/test.py` 调用的三个 Scanpy API 默认值
均为 0；该值进入 typed authority 和 effective config hash，但不暴露为新的
Application YAML 用户参数。其他已有内部固定 seed 不在本轮迁移范围内。

## Distribution contract

wheel/sdist：

- 不包含 `revise.yaml`；
- 包含 5 个 Application 和 6 个 Benchmark templates；
- installed CLI 与 source checkout 使用相同 schema/authority；
- distribution test 对 repo/package template 做 byte comparison，并检查缺失/多余文件。

继续阅读：[Verification and P2CRC parity](06-verification-and-p2crc-parity.md)。
