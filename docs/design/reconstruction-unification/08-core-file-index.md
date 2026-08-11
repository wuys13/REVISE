# Core File and Test Index

Parent index: [Reconstruction Unification Design Package](README.md)

本索引按职责列出 baseline 文件、目标位置与验证入口。实施时以实际 checkout 为准；目标文件尚未存在时不得把本表误读为 Current fact。

## OT, GA, and LR

| Responsibility | Baseline/current entry | Decided target / tests |
| --- | --- | --- |
| Annotation GA | `revise/backend/kernels/global_anchoring.py` | 保留 facade；solver 下沉到 `revise/backend/kernels/ot.py` |
| sc annotation LR | `revise/backend/kernels/local_anchoring.py` | 只依赖 `OTKernel.annotate`；focused tests in `tests/backend/` |
| Matrix coupling | `revise/backend/ops/local_ot.py` 与多个 runner call sites | `OTKernel.couple`；迁移现有 local OT tests，增加 static boundary test |
| Runner-owned GA | `revise/backend/runners/base_svc_anchor.py` | 删除；Application/Benchmark bases 直接继承 `BaseSVC` |
| GA strategy | `revise/backend/adapters.py::RunnerBackedStrategy` | 唯一 stage-level GA facade call |
| sc-SVC LR | `revise/backend/runners/sc_svc_application.py` | 两次 LocalAnchoring + graph cluster trace |
| sp LR | `sp_svc_application.py`, `sp_svc_benchmark.py` | route-owned problem → `OTKernel.couple` |
| SR LR | `sc_svc_sr_application.py`, `sc_svc_sr_benchmark.py` | mandatory allocation + optional couple/graph aggregation |
| Imputation LR | `sc_svc_impute_benchmark.py` | in/all-panel route-owned problem → `OTKernel.couple` |

Primary existing tests to preserve/extend：

- `tests/backend/test_solver_routing.py`
- `tests/backend/test_tacco_global_freshness.py`
- `tests/backend/test_sc_local_ot.py`
- `tests/backend/test_scientific_contracts.py`
- `tests/integration/solvers/test_tacco_solver_smoke.py`
- `tests/integration/solvers/test_local_refinement_solver_smoke.py`
- new repository/static boundary coverage under `tests/repository/`

## Application entry and preprocessing

| Responsibility | Baseline/current entry | Decided target / tests |
| --- | --- | --- |
| Public entry | `reconstruct.py` | visible high-level flow; direct AnnData return |
| Parser | `reconstruct.py` | `revise/application/cli.py` |
| Config compiler | `revise/application/config.py` | strict schema + explicit overrides + preprocessing/filter fields |
| Three functions | `revise/backend/adapters.py`, `application_svc.py` | `revise/application/preprocess.py` |
| Publication | `reconstruct.py::_write_outputs` | `revise/application/publication.py` |
| Input/preflight | `revise/io/input_service.py` and Application tests | filter-column/value lightweight preflight |

Primary tests：

- `tests/application/test_request.py`
- `tests/application/test_preflight_selection.py`
- `tests/application/test_cli_dry_run.py`
- `tests/application/test_publication.py`
- `tests/application/test_service.py`
- `tests/backend/test_application_column_contract.py`
- new preprocessing and return-contract tests under `tests/application/`

## Config, Benchmark, and distribution

| Responsibility | Baseline/current entry | Decided target / tests |
| --- | --- | --- |
| Engine defaults/routes | `revise/revise.yaml`, `revise/config/loader.py` | one typed package authority |
| Application YAML | `configs/application/`, `revise/application/templates/` | 5 byte-exact templates |
| Benchmark selection | `revise/benchmark/cli.py --confounding` | 6 route YAMLs under repo/package template roots |
| Benchmark launcher | `revise/benchmark/launcher.py` | task YAML path is selector |
| Noise pollution | `revise/backend/adapters.py` SR Benchmark branch | remove noise-only symbols, protect graph aggregation |
| Packaging | `MANIFEST.in`, `pyproject.toml` | no `revise.yaml`; package 5+6 templates |

Primary tests：

- `tests/config/test_ot_config.py`
- `tests/benchmark/test_cli_contract.py`
- `tests/benchmark/test_launcher.py`
- `tests/application/test_templates.py`
- `tests/integration/distribution/test_artifacts.py`
- `tests/recon/test_provenance_identity.py`
- `tests/recon/test_run_manifest.py`

## Route-specific review order

```text
standard sc-SVC:
reconstruct.py → application preprocess/config → ScSvcApplicationStrategy
→ GlobalAnchoringKernel/OTKernel → ScSVC/LocalAnchoring/GraphCluster
→ publication → route/parity tests

sp-SVC:
route YAML → strategy → GA → sp runner problem construction
→ OTKernel.couple → GraphAggregate → output tests

sc-SVC-sr / SR Benchmark:
route YAML → ensure_all_cells_in_spot → GA → SpotSr/mandatory allocation
→ optional coupling/GraphAggregate → primary/evaluation tests

imputation:
benchmark YAML → prep/GA → uncertainty/subcluster
→ per-cell-type coupling → GeneImpute → paired outputs/evaluation
```

继续阅读：[Delivery status](09-delivery-status.md)。
