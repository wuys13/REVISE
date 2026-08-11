# Application Entry, Preprocessing, Return, and Publication

Parent index: [Reconstruction Unification Design Package](README.md)

## Baseline

Current fact：Application preprocessing 分散在 `revise/backend/adapters.py`、`ApplicationSVC._adata_processing` 和 sc-SVC-sr constructor；adapter 还会在缺少 `transcript_counts` 时从 `X` 推导。`reconstruct.py` 返回 `ApplicationExecution`，并在入口内部实现 parser、override mapping 与 atomic H5AD publication。

## Three independent preprocessing functions

Decided target：包内只增加三个可独立调用的函数，不增加第四个 orchestration wrapper。

### `filter_reference`

- 输入：reference AnnData、`filter_column`、`filter_value`。
- 当两者非空时执行精确相等过滤并返回 copy；二者都为空时返回未过滤 copy。
- schema 保证二者成对出现；缺列或匹配零行直接失败。
- Xenium 正式值是 `Patient == P2CRC`，但接口不携带 patient-specific 命名。

### `preprocess_spatial`

- 默认 `min_transcript_counts=60`、`min_cell_counts=100`。
- transcript threshold 非空时要求 `obs["transcript_counts"]` 已存在，并执行 observation subset；不得从 `X` 推导。
- gene filtering 使用 `scanpy.pp.filter_genes(min_cells=min_cell_counts)`。
- 返回处理后的 AnnData，不依赖 runner config。

### `preprocess_reference`

- 默认 `min_transcript_counts=None`、`min_cell_counts=100`。
- 默认不做 sc transcript/cell count filtering；若调用者显式传非空 threshold，使用与 spatial 相同的已有列契约。
- gene filtering 使用 `scanpy.pp.filter_genes(min_cells=min_cell_counts)`。
- 返回处理后的 AnnData，不投影 obs 列、不做 slash normalization。

## Visible `reconstruct.py` flow

Decided target：完整 Application run 的 lifecycle callback 中必须直接看到并接收三个返回值：

```text
reference = filter_reference(reference, ...)
spatial = preprocess_spatial(spatial, ...)
reference = preprocess_reference(reference, ...)
```

之后由 route preparation 依次完成：

- standard sc-SVC：reference obs 投影到 configured Level1/Level2 → slash normalization → spatial-only overlap subset；
- sp-SVC：保留该路线既有 overlap 与 label normalization 语义；
- sc-SVC-sr：`ensure_all_cells_in_spot` 必须在 callback 前发生，然后再进入三个函数。

完整 run 未提供 Application preprocessing callback 时直接失败。Benchmark 不调用这三个函数；dry-run 不加载和处理完整矩阵。

## CLI and package ownership

Decided target：

- `revise/application/cli.py`：argument parser、legacy flag rejection、CLI payload。
- `revise/application/config.py`：strict schema、explicit overrides、compile effective Application config。
- `revise/application/publication.py`：output paths、H5AD metadata、staging 与 atomic replacement。
- `reconstruct.py`：可读高层线路与公共 `run_application`；不再暴露 `ApplicationExecution`。

公共 `run_application` 的 preprocessing override 参数默认都是 `None`，其唯一含义是“不覆盖 YAML”。三个预处理函数自身的默认是 spatial 60/100、reference None/100；官方 YAML 显式写出同样的有效行为，其中 reference `min_transcript_counts: null` 表示不做 transcript filtering。这里追求的是三层**有效行为一致**，不是三层参数字面默认相同。

由于 public `None` 已用于“无 override”，它不提供把一个自定义 YAML 的非空 reference threshold 显式改回关闭的第二种语义；本轮不为此引入 sentinel。最终实际采用的 threshold 必须进入 effective config hash/provenance。

## Return contract

| Route/action | Return |
| --- | --- |
| `sp-SVC` | primary reconstructed AnnData |
| `sc-SVC-sr` | final primary AnnData；graph aggregation enabled 时返回 graph-aggregated object |
| `sc-SVC` | `(spatial_adata, expression_adata)`，顺序固定 |
| dry-run | `None` |

Publication 与 return 必须使用同一批 AnnData 对象，不得写一个 copy 后返回另一个语义不同的对象。可以在写盘前就地附加正式 metadata，但不能改变轴或表达矩阵。

## Output naming

- standard sc-SVC 无 `output.name`：`spatial.h5ad`、`expr.h5ad`
- standard sc-SVC 有 `output.name=foo`：`foo_spatial.h5ad`、`foo_expr.h5ad`
- 单输出路线无 name：`svc.h5ad`
- 单输出路线有 name：`<name>.h5ad`

内部 run identity 不从 public filename 反推。未提供 name 时使用稳定 route identity 作为内部 sample/logger identity。

## Dry-run preflight

dry-run 不加载完整 AnnData 或执行 preprocessing/GA/LR，但必须：

1. 编译完整 effective config；
2. 检查输入路径/格式与 route-specific required fields；
3. 读取轻量 H5AD metadata，确认 reference filter column 存在且至少一行匹配 filter value；
4. 记录预期 output paths、配置来源与 hash；
5. 不创建 public H5AD。

## Code removal boundary

删除 adapter filters、`_ensure_transcript_counts`、Application `._adata_processing`、sc-SVC-sr constructor 中的重复 preprocessing，以及 public `ApplicationExecution` exports/docs/tests。Benchmark `._adata_processing` 不删除；slash label normalization 不删除。

继续阅读：[Configuration and provenance](05-configuration-and-provenance.md)。
