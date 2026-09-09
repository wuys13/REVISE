---
title: Standard ST directories and hierarchical batch configuration
date: 2026-09-09
type: refactor
---

# Standard ST directories and hierarchical batch configuration

Approved implementation plan, 2026-09-09. This is a one-time work record;
maintained guidance lives in the batch reconstruction, protocol and framework
pages under [docs/design/batch/](../design/batch/README.md).

## Outcome

Replace experimental per-sample YAML preparation with standard `spatial.h5ad`
sample directories, project/group/sample `batch.yaml` inheritance and references
shared by path. Inputs remain unchanged. Keep existing output layout, scientific
single-run behavior, recovery and analysis adapters. Do not implement data
conversion or import reconstruction-impact algorithms.

## Protocol decisions

- Project schema version 2 declares separate input/output roots. Root keys are
  forbidden in overrides. Resolve paths against the declaring YAML before merging.
- Inherit project → input root → intermediate groups → sample. Recursively merge
  mappings, replace scalars/lists, and replace `inputs.reference` as a whole so
  old filters cannot leak across references. Validate nulls by field contract.
- Discover only directories with `spatial.h5ad`; IDs are relative paths. Reject
  nested samples; do not follow directory symlinks. `enabled` defaults to true
  and can be overridden; disabled results become inactive without deletion.
- Require standard H5AD X, spatial coordinates, configured labels and units.
  Check read-only; retain existing QC, reference filtering and gene alignment.
- Reuse checks inputs, effective settings, the complete configuration chain
  (including additions/deletions), code identity and artifact hashes. Analysis
  resolves the same chain before accepting reconstruction.
- Retain iST type defaults and paired/mean/random behavior. Remove the old prepare
  command and sample templates; explain v1 migration without moving/deleting data
  or trusting old cache entries as v2 successes.

## Implementation sequence

1. Add same-level usage, protocol and framework docs plus project/group/sample
   override templates, with stable anchors and a single navigation entry.
2. Add central configuration resolution and sample discovery; replace preparation
   with a read-only descriptor and checks.
3. Connect reconstruction and analysis to this shared resolution; preserve failure
   isolation, sequential tasks, handoff semantics and publication rollback.
4. Remove obsolete experimental interfaces and update packaging/documentation.
5. Verify regression behavior and bounded real-data parity; link maintained docs
   to the implementation and its tests.

## Acceptance

Test recursive inheritance, list/reference replacement, paths from other working
directories, disabled/invalid configurations, automatic discovery and duplicate
leaf names across groups. Assert shared input content and directory inventories
remain unchanged. Test successful reuse, configuration file addition/removal/edit
invalidation, failure continuation and stale-analysis rejection. Preserve iST
publication tests and compare bounded hST/iST/sST batch and single-run results.
Validate installation, commands, example templates and documentation navigation.

A new standard sample should usually need only one added `spatial.h5ad`; it
inherits configuration and uses an existing shared reference automatically.
