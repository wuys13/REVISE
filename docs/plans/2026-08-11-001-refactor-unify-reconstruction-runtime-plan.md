---
title: "refactor: Unify reconstruction runtime, OT boundaries, and route configuration"
type: refactor
date: 2026-08-11
---

# refactor: Unify reconstruction runtime, OT boundaries, and route configuration

## Summary

Unify all real TACCO/POT calls behind one bottom-level `OTKernel`, separate GA facade ownership from every LR path, expose the three Application preprocessing calls and direct AnnData returns in `reconstruct.py`, and replace `revise.yaml`/Benchmark `--confounding` with 5 Application plus 6 Benchmark route YAMLs. Verification is route-first and uses the legacy P2CRC sc-SVC implementation as a stepwise numerical baseline.

Behavioral authority and detailed contracts live in the [Reconstruction Unification Design Package](../design/reconstruction-unification/README.md). Delivery evidence is recorded only in [09-delivery-status.md](../design/reconstruction-unification/09-delivery-status.md).

## Frozen boundaries

- GA may use `GlobalAnchoringKernel`; LR/runners may not mention it.
- `OTKernel.annotate` and `OTKernel.couple` are sibling contracts and contain all direct TACCO/POT calls.
- Application uses exactly three independent preprocessing functions, directly sequenced by `reconstruct.py`; slash normalization remains route preparation.
- `run_application` returns AnnData, `(spatial, expression)` for standard sc-SVC, and `None` for dry-run.
- Benchmark has exactly six YAML routes and no noise/imputation meta-route.
- No solver fallback, category reorder, posterior renormalization, zero-overlap repair, or new slash collision policy.
- `raw_data/`, `test_sc/`, and generated outputs remain untracked.

## Progressive references

| Work | Read before implementation |
| --- | --- |
| Route/call order | [02-route-contracts.md](../design/reconstruction-unification/02-route-contracts.md) |
| OT/GA/LR | [03-ot-boundaries-and-contracts.md](../design/reconstruction-unification/03-ot-boundaries-and-contracts.md) |
| Application/preprocessing/publication | [04-application-entry-and-preprocessing.md](../design/reconstruction-unification/04-application-entry-and-preprocessing.md) |
| YAML/defaults/provenance | [05-configuration-and-provenance.md](../design/reconstruction-unification/05-configuration-and-provenance.md) |
| Tests/P2CRC | [06-verification-and-p2crc-parity.md](../design/reconstruction-unification/06-verification-and-p2crc-parity.md) |
| Deletion/risk/review decisions | [07-implementation-and-risk.md](../design/reconstruction-unification/07-implementation-and-risk.md) |
| Production/test navigation | [08-core-file-index.md](../design/reconstruction-unification/08-core-file-index.md) |

## Implementation units

### U1. Freeze documents and characterization baseline

Create the design package, enumerate 3 Application and 6 Benchmark route contracts, record current-vs-target language, and map production/test entrypoints. Verify internal links, `4/1/16/4/1/1` cardinality, and an initially honest “not implemented/not verified” delivery status.

### U2. Move matrix coupling into `OTKernel.couple`

Preserve route-owned masses/cost/support construction and POT/TACCO numerical behavior. Move only validation/solver execution and the existing hard-support/marginal continuation machinery. Verify dtype, support, reference-measure and marginal characterization before migrating every caller.

### U3. Move annotation into `OTKernel.annotate`

Make GA facade and sc-SVC LocalAnchoring call the same annotation contract. Remove LR→GA and runner-owned GA state/base layer. Verify fresh TACCO keys, reference copies, returned category order, and static import/call boundaries.

### U4. Centralize Application preprocessing and expose the entry flow

Add reference filtering, spatial preprocessing and reference preprocessing functions with function defaults 60/100 and None/100. Official YAMLs explicitly produce the same effective behavior; public override defaults remain `None` meaning only “do not override YAML”, without adding a sentinel. Call the functions directly in the `reconstruct.py` callback; preserve sc-SVC obs projection/slash/axis order and sc-SVC-sr pre-callback `ensure_all_cells_in_spot`. Delete duplicate adapter/class processing only after once-only route tests pass.

### U5. Return and publish AnnData directly

Move parser/config/publication mechanics into Application package modules, remove `ApplicationExecution`, and ensure publication uses the same final objects returned to Python callers. Verify single/pair atomic publication, graph-aggregated primary selection, dry-run `None`, and filenames: unnamed single output `svc.h5ad`, named single output `<name>.h5ad`, standard sc-SVC pair with optional `<name>_` prefix.

### U6. Establish typed config authority and Application YAMLs

Make defaults/routes/locked keys package-owned; compile complete effective config before runner construction. Migrate all result parameters, update five repo/package Application templates, and lock Xenium_T to Patient=P2CRC, 60/100, null/100, alpha 0.2 and resolutions 0.6/0.7/0.8.

### U7. Establish six Benchmark YAMLs and remove noise pollution

Make each Benchmark YAML the route selector, update launcher/CLI, preserve case expansion/evaluation, and delete only the specified SR noise code. After repo/package template and effective-config tests pass, remove `revise.yaml` and its distribution entries.

### U8. Route-driven verification and delivery evidence

Assert all nine route traces and static boundaries, build wheel/sdist, and run relevant suites. Then execute the legacy/package P2CRC eleven-checkpoint comparison in a frozen environment if the data/runtime is available. Update delivery status without conflating structural, synthetic, real numerical and scientific evidence.

## Acceptance gates

1. All direct `tacco.tl.annotate`, `tacco.utils.solve_OT` and `ot.unbalanced.sinkhorn_unbalanced` calls are inside `OTKernel`.
2. LR/runners contain no `GlobalAnchoringKernel`; standard sc-SVC records one GA and two LR annotations.
3. Three preprocessing functions run exactly once in the documented order and preserve spatial/reference output axes.
4. `reconstruct.py` visibly returns the reconstructed objects; output naming and atomic publication match the contract.
5. Effective config/provenance contains every result-affecting value; distribution contains 5+6 templates and no `revise.yaml`.
6. All Application/Benchmark route traces and Benchmark cardinalities pass.
7. P2CRC checkpoints report either the first precise divergence or all eleven technical parity checks; old `max_score` is compared by exact observation index and dtype-aware numeric tolerance to the corresponding new `Confidence` at the same annotation time, while the new public schema remains `Confidence`. Biological validity remains a separate scientific conclusion.

## Assumptions

- TACCO returned category order is authoritative.
- `test_sc/test.py` is the standard sc-SVC legacy semantics baseline, not a package/API shape baseline.
- Benchmark route selection is YAML-only after migration.
- Current checkout is the implementation input and source of code facts; delivery-status claims are owned only by `docs/design/reconstruction-unification/09-delivery-status.md`. No defensive checkout/cherry-pick or unrelated cleanup is part of this plan.
