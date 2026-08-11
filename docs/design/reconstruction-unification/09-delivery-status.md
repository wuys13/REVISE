# Delivery Status

Parent index: [Reconstruction Unification Design Package](README.md)

**Status as of 2026-08-11: the agreed software boundaries and the P2CRC
Xenium_T/TACCO technical parity check are complete. Biological validity remains
unknown.**

本文是唯一把设计目标升级为当前交付事实的位置。

## Current fact

- 直接 TACCO/POT 求解均在 `OTKernel`；GA 经
  `GlobalAnchoringKernel` facade，LR 不引用它。
- `reconstruct.py` 显式顺序执行 reference filter、spatial preprocessing、
  reference preprocessing，并返回 AnnData 或 sc-SVC tuple。
- Application spatial preprocessing 现在分别支持 Xenium 的
  `min_transcript_counts` 和 sp-SVC 的 `min_counts`；reference 分别支持
  `min_transcript_counts` 和 sp-SVC 的 `min_genes`。VisiumHD 配置恢复为
  spatial `min_counts=20/min_cell_counts=30`、reference
  `min_genes=20/min_cell_counts=50`。
- Fib、Mono、T 的 Xenium spatial input 都是 P2CRC；repo/package 的 5+6
  YAML 模板镜像逐字节一致。
- Visium sc-SVC-sr 的正式 YAML 已声明预处理、PCA graph、
  `match_spot_sum` 与 `PM_on_cell.csv`。重命名后的 case notebook 只通过
  `reconstruct.py` 运行重建，再加载正式发布的 H5AD 做后续分析。
- Application 与 Benchmark 的显式 `seed: null` 统一解析为 42；兼容
  subsample 和 graph/Leiden 的固定 seed 0 不受此设置影响。
- 每个 Benchmark case 的 provenance 现在记录 route、已解析 io、算法
  override、runtime seed 的 effective request/hash，以及显式 CLI overrides。

## Verification ledger

| Gate | Status | Evidence |
| --- | --- | --- |
| Direct preprocessing/OT/config/provenance routes | Passed | 66 focused tests passed; one pre-existing AnnData view warning |
| OT static boundary | Passed | solver calls are confined to `OTKernel`; LocalAnchoring has no `OTKernel.couple` or `GlobalAnchoringKernel` dependency |
| P2CRC Xenium_T/TACCO checkpoints 1–11 | Passed (technical parity) | [compact parity evidence](evidence/p2crc-tacco-parity-2026-08-11.json); post-change new run compared with frozen legacy checkpoint artifacts |
| Full non-integration suite | Not run by scope decision | Only direct routes were run |
| Real POT dual run | Not run by scope decision | structural route test only |
| Biological validity | Scientific unknown | technical equality is not biological validation |

## P2CRC evidence boundary

The P2CRC check used Python 3.10, TACCO 0.5.0, seed 42, fixed one-thread
settings, the `Xenium_T.yaml` TACCO route, and fresh output under
`/private/tmp`. The old and new outputs have identical axes, three annotation
posteriors/labels/Confidence, graph inputs and outputs, Leiden results and
best resolution, reverse annotation, and final sparse `X` matrices. The
versioned evidence records only summary/hash/count facts, not H5AD or matrices.
This was rerun after removing implicit reference-category reordering and
zero-overlap reference-cell filtering; the same comparison remained exact.

The historical `max_score` column is intentionally not reconciled in this
round. `Confidence` is the parity comparison field. The spatial legacy column
can remain stale, as explicitly accepted for this delivery.

## Explicitly deferred

- Pairwise rollback if the second sc-SVC output-file replacement fails.
- dry-run full-file SHA I/O cost.
- inherited `max_score` cleanup and the generic `_build_svc` legacy fallback.
- real POT parity run, full non-integration suite, and biological validation.
- Treating a configured but missing `inputs.pm_on_cell.path` as random
  allocation. The current behavior is an input failure; only omission selects
  seeded random allocation.

## Delivery inclusion

The design package, this status file, its compact evidence file, and the main
implementation plan are intentionally force-tracked despite the broad
`docs/design/` and `docs/plans/` ignore rules. `raw_data/`, `test_sc/`, and
temporary P2CRC outputs are not included.

返回 [设计包入口](README.md)。
