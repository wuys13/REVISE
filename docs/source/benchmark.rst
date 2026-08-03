Benchmark Reproduction
======================

Benchmark reproduction is a checkout-only paper path. It is separate from the
installed ``revise-reconstruct`` command and requires the corresponding paper
data.

Run one confounding family
--------------------------

From the repository root:

.. code-block:: bash

   python reproduce/benchmark_main.py \
     --confounding segmentation \
     --data-root raw_data/Sim2Real-ST \
     --dataset-task segmentation \
     --sample-name P2CRC/cut_part1 \
     --output-root output/benchmark

``--confounding`` accepts ``segmentation``, ``bin2cell``, ``batch_effect``,
``spot_size``, ``gene_panel``, or ``gene_dropout``. Route-specific parameters
determine the leaf runs inside that family:

- segmentation runs four segmentation leaves;
- bin2cell runs one leaf;
- batch effect runs four batch levels for every discovered spot size, falling
  back to the four configured spot sizes when none are discoverable;
- spot size runs four fixed spot-size leaves;
- gene panel and gene dropout each run one leaf.

Dedicated configuration controls
--------------------------------

Benchmark variants use named options rather than a generic key/value override.
These options are applied after the selected profile or custom ``--config``:

- ``--local-refinement-strength`` sets the non-negative finite OT conditioning
  strength for hST-SVC and sST-SVC routes;
- ``--sr-refinement-preset confidence_anchor`` selects the controlled graph
  refinement used by ``batch_effect`` and ``spot_size`` runs, while ``none``
  disables that refinement.

Omitting the strength creates no CLI override; route defaults remain
authoritative. Each result reports minimal ``local_refinement`` evidence with
``route``, ``applied``, and ``strength``. Removed policy and posterior flags are
rejected with a migration message rather than translated. sST-SVC always
performs composition, row assignment, and closed-form expression allocation;
the strength controls only posterior-conditioned local OT.

``--seed-scope process`` (the default) derives a reproducible effective seed
for each leaf run from ``--seed``. Every effective seed is stored as
``runtime.seed`` in that leaf's ``merged_config.json``, contributes to its
``config_hash``, and is included in the benchmark report. Use
``--seed-scope run`` to reuse the same ``--seed`` for every leaf instead.

For the example above, ``seg_1`` is one of four outputs. Its primary
compatibility metric table is shaped like:

.. code-block:: text

   output/benchmark/segmentation/P2CRC/cut_part1/seg_1/metrics_normalized.csv

The same run directory contains ``provenance.json`` and keyed tables below
``metrics/``. Treat a metric CSV as completed evidence only when the run
manifest records a successful evaluation and completed artifact hash.

Run the bounded launcher
------------------------

.. code-block:: bash

   bash reproduce/benchmark_main.sh

The shell entry stays in the foreground while ``revise.benchmark.launcher`` owns a
bounded set of child process groups. It writes ``launcher_status.json``,
forwards interruption, reaps children, and returns nonzero if any case fails.
The launcher does not convert a partial suite into success.

Metric preprocessing
--------------------

``compute_metric(ground_truth, prediction, ..., normalize=True)`` implements
the paper-compatibility table as follows:

1. total-normalize every observation in each input independently to ``1e4``;
2. select the requested genes or the shared gene names;
3. when requested, log-transform both inputs;
4. independently min-max normalize every gene in ground truth and prediction;
5. compute per-gene PCC, one-dimensional SSIM, MSE, and NRMSE.

The pipeline aligns the selected prediction and ground truth by shared
observation identifiers before this function is called. The metric function
then consumes those aligned arrays; it does not perform a second cell-ID join.

NRMSE orientation
-----------------

For each normalized gene, the implementation is:

.. code-block:: text

   NRMSE = sqrt(MSE) / mean(ground truth)

The denominator is the min-max-normalized ground-truth vector, not a symmetric
scale computed from both arrays. NRMSE is therefore directional: swapping
ground truth and prediction need not preserve the value.

SSIM row-order boundary
-----------------------

SSIM receives the aligned one-dimensional expression vectors in their current
row order with ``data_range=1`` for normalized metrics. It does not use spatial
coordinates, rasterize a tissue, or construct a spatial image. Even a common
permutation of both vectors can change local SSIM windows. This value is a
row-order compatibility metric, not coordinate-aware spatial-image evidence.

Undefined values
----------------

The implementation intentionally exposes constant and zero-mean mathematical
undefined cases:

- a constant normalized gene has undefined PCC and produces ``NaN``;
- when normalized ground-truth mean and error are both zero, NRMSE produces
  ``NaN``;
- when normalized ground-truth mean is zero but error is nonzero, NRMSE
  produces positive infinity.

REVISE does not replace, drop, or coerce these values to zero. Consumers
must retain the per-gene table and handle non-finite values explicitly when
aggregating.

Proof limit
-----------

Metric-contract tests prove formula orientation, identifier alignment, and row
order sensitivity on synthetic arrays. They do not prove biological
validation, real-tissue reconstruction quality, or comparability across
datasets with different preprocessing. The current real-data end-to-end suite
is deferred and is not part of this candidate's evidence.

Route-specific local-refinement tests and candidate-wheel solver smoke prove
configuration, axis alignment, Adapter invocation, solver dispatch, and the
minimal ``route/applied/strength`` record. The CI solver gate imports the
packaged candidate outside the source checkout, runs copied tests in pytest
importlib mode, asserts the test process is importing from the isolated wheel
environment, and exercises both POT and TACCO conditioned-cost paths. These
checks do not show that posterior conditioning improves any biological or
benchmark metric.

Benchmark notebooks
-------------------

The curated benchmark notebooks are tracked under ``reproduce/benchmark/``.
Their embedded historical outputs do not expand the proof boundary above and
do not establish that the current source checkout reran the paper datasets.
