# Test architecture

Tests mirror production ownership. Run the fast synthetic suite during normal
development and run `tests/integration/` when validating a built distribution
or an optional solver. No test in this repository proves biological validity.

```bash
# Fast synthetic and repository contracts
pytest -q --ignore=tests/integration

# Distribution, installed CLI, and optional solver boundaries
pytest -q tests/integration
```

`tests/conftest.py` exposes the repository root as the `repo_root` fixture so
new tests do not depend on directory depth. Integration tests may require the
environment variables documented in their modules and CI jobs.

| Test file | Production owner | Detects / corresponding implementation | When used | Proof boundary |
| --- | --- | --- | --- | --- |
| `analysis/test_cli.py` | `revise/analysis/` | Biological-metrics console entry and delegation | Package/analysis edits | Static entry contract only |
| `application/test_cli_contract.py` | `revise/application/service.py` | Canonical publication, rollback, and result manifest | Application edits | Synthetic AnnData only |
| `application/test_cli_dry_run.py` | `revise/application/cli.py` | Source CLI preflight without reconstruction | Application/config edits | No scientific stages |
| `application/test_public_contract.py` | `revise/application/` | 1.x selector and internal route mapping | Every application change | Public vocabulary and dispatch only |
| `application/test_service.py` | `revise/application/service.py` | sc-SVC pair publication and rollback | Result-assembly edits | No kernels or real data |
| `backend/test_application_column_contract.py` | Application sp-SVC/sc-SVC-sr runners | Configured annotation columns reach route-specific local refinement | Application column changes | Stops at the configured downstream call |
| `backend/test_local_refinement_conditioning.py` | `revise/backend/ops/assignment.py`, `posterior_conditioning.py` | Strict global posterior axes and fixed local OT cost conditioning | Assignment/local-refinement edits | Synthetic matrices only |
| `backend/test_distance_contract.py` | `revise/backend/ops/distance.py` | Distance formulas and invalid inputs | Numerical-op edits | Small arrays |
| `backend/test_graph_cluster_spatial_score.py` | `revise/backend/` graph code | Graph aggregation, clustering, and scale smoke | Graph/backend edits | Component scale, not full pipeline |
| `backend/test_sc_local_ot.py` | `revise/backend/kernels/local_anchoring.py` | sc-SVC LR uses the selected POT/TACCO path | Local-OT edits | Synthetic local units |
| `backend/test_sc_graph_guidance.py` | sc-SVC graph route | Argmax cohort routing and ordinary GraphCluster invariance to soft-Q changes | Graph-route edits | Synthetic Graph and public-route proof only |
| `backend/test_sc_imputation_guidance.py` | sc-SVC imputation route | Guidance-free imputation, solver routing, and panel/dropout behavior | Imputation-route edits | Synthetic imputation proof |
| `backend/test_sc_sr_guidance.py` | sc-SVC-sr route | Mandatory allocation invariance, projected virtual Q, and conditioned local OT | SR-route edits | Synthetic SR proof |
| `backend/test_scientific_contracts.py` | Backend kernels/ops and analysis formulas | Cross-cutting scientific invariants | Scientific implementation edits | Formula/array contracts only |
| `backend/test_solver_routing.py` | `revise/backend/ops/` and route callers | Solver dispatch, event recording, and caller coverage | OT edits | Mocked/synthetic solver boundaries |
| `backend/test_spot_sr_assignment.py` | `revise/backend/kernels/spot_sr.py` | PM-on-cell and random cell-type allocation | Spot-SR edits | Assignment, not localization truth |
| `backend/test_spot_sr_min_cells.py` | `revise/backend/adapters.py` | Spot-SR gene filtering and overlap checks | Spot-SR preprocessing edits | Synthetic matrices |
| `backend/test_spot_sr_quota.py` | `revise/backend/kernels/spot_sr.py` | Rounded quota repair and large component smoke | Spot-SR allocation edits | Quota correctness only |
| `backend/test_sp_assignment_guidance.py` | sp-SVC route | Neighbor/replacement OT posterior conditioning and public application/benchmark routes | sp-route edits | Synthetic local-problem proof |
| `backend/test_tacco_global_freshness.py` | TACCO global anchoring | Freshness and failure behavior | TACCO edits | Synthetic solver interaction |
| `benchmark/test_cli_contract.py` | `revise/benchmark/cli.py` | Public benchmark options produce typed algorithm configuration without a generic override flag; process-scoped seeds become explicit per-leaf seeds | Benchmark CLI/config edits | Argument, seed, and wrapper contract only |
| `benchmark/test_local_refinement_cli.py` | `revise/benchmark/cli.py` | Unique strength flag, removed-flag migration error, and minimal leaf aggregation | Local-refinement CLI edits | Argument/report contract only |
| `benchmark/test_launcher.py` | `reproduce/benchmark_main.sh`, `revise/benchmark/launcher.py` | Bounded parallel launch and failure propagation | Benchmark launcher edits | Fake worker processes |
| `config/test_local_refinement_contract.py` | `revise/config/loader.py` | Route defaults, sc/imputation rejection, bounds, and removed YAML grammar | Local-refinement config edits | Resolved config only |
| `config/test_ot_config.py` | `revise/config/` plus OT application wiring | Schema, profiles, locked algorithm parameters, and solver mapping | Config/OT edits | Resolved configuration, not real runs |
| `io/test_h5ad_preflight.py` | `revise/io/` and input resolution | H5AD roles, axes, labels, overlap, and dependency preflight | Input-contract edits | Metadata and bounded synthetic values |
| `preprocess/test_cli.py` | `revise/preprocess/` | Histology-prior console entry | Preprocess packaging edits | Static entry contract only |
| `preprocess/test_histology_priors.py` | `revise/preprocess/histology_priors.py`, `revise/backend/ops/meta.py` | Segmentation centers persist in the standard cell-location table and missing centers use spot coordinates | Histology/SR input edits | Synthetic coordinates only |
| `recon/test_cross_process_determinism.py` | `revise/recon/` and deterministic utilities | Same-seed identities across processes | Determinism edits | Synthetic components |
| `recon/test_failure_provenance.py` | `revise/recon/context.py`, `revise/framework.py` | Failure/interruption manifests and publication compensation | Lifecycle edits | Catchable failures, not power loss |
| `recon/test_lifecycle_trace.py` | `revise/recon/pipeline.py` | Fixed stage order and trace transitions | Pipeline-stage edits | Synthetic strategy |
| `recon/test_local_refinement_record.py` | `revise/recon/context.py` | Minimal route/applied/strength record and monotonic applied state | Local-refinement provenance edits | In-memory/durable record only |
| `recon/test_provenance_identity.py` | `revise/utils/provenance.py` | Hashes, canonical identities, and deterministic paths | Provenance edits | Identity mechanics only |
| `recon/test_run_manifest.py` | Framework/provenance run envelope | Unique run directories, locks, and terminal manifests | Run-lifecycle edits | Local filesystem only |
| `recon/test_runtime_validation_errors.py` | Framework validation policies | Error typing and command failure semantics | Validation edits | Synthetic invalid requests |
| `repository/test_documentation_contract.py` | README/docs/package metadata | Public claims match code and tests | Documentation/interface edits | Claim consistency only |
| `repository/test_entrypoint_boundaries.py` | `reconstruct.py`, `reproduce/benchmark_main.*`, and `revise.benchmark` | The root exposes only reconstruction while benchmark launchers locate package-owned implementations | Entrypoint or package-boundary changes | Wrapper/path behavior only |
| `repository/test_import_side_effects.py` | Package imports and optional analysis imports | Host environment stays unchanged; lazy optional imports | Import/dependency edits | Subprocess import behavior |
| `repository/test_optional_dependencies.py` | `pyproject.toml` extras and optional imports | Missing extras fail with install guidance | Dependency edits | Dependency boundary only |
| `repository/test_repository_layout.py` | Repository root/package ownership | Removed maintenance surfaces and package-owned utilities | Repository cleanup | File/layout assertions |
| `repository/test_test_layout.py` | `tests/` | Ownership directories and this exhaustive index | Any test move/addition | Organization only |
| `integration/application/test_installed_cli.py` | Built wheel + application CLI | Wheel install, external-CWD preflight, small POT source/wheel parity | Release candidate | 52x52 synthetic run; no biology |
| `integration/distribution/test_artifacts.py` | Build metadata and artifacts | Exact wheel/sdist contents and clean installation | Release candidate | Packaging only |
| `integration/solvers/test_tacco_solver_smoke.py` | Optional TACCO installation | Real TACCO 0.5.0 GA/LR smoke and import order | TACCO/CI release gate | Tiny solver smoke; no parity claim |
| `integration/solvers/test_local_refinement_solver_smoke.py` | Installed POT/TACCO | Real posterior-conditioned local OT candidate | Local-refinement solver/release gate | Tiny matrices; skipped solvers are not evidence |

When adding a test module, place it under its production owner, add it to this
table, and state both what it detects and what it cannot prove.
