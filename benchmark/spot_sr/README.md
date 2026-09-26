# Spot super-resolution benchmark handoff

This directory reproduces the **TESLA** and **iStar** baselines on the released
REVISE Sim2Real-ST spot-size cases and provides a separate REVISE spot-route
runner. The current scope is P1CRC `part1`–`part3`, spot sizes 50, 100, 150,
and 200 µm, with 324 genes in each case. The 24 baseline runs finished and
passed artifact validation on qz on 2026-09-26. Two REVISE single-case pilots
and the full four-size `part2` route passed with reconstructed inputs. See
[HANDOFF.zh-CN.md](HANDOFF.zh-CN.md) for the
Chinese step-by-step handoff.

## What is included

| File | Purpose |
| --- | --- |
| `prepare_inputs.py` | Register and crop H&E; convert pseudo-spots into TESLA/iStar inputs |
| `run_all.sh`, `run_tesla.py` | Prepare cases and run TESLA/iStar |
| `evaluate.py`, `summarize.py`, `validate_results.py` | Align predictions to held-out cells; score and audit all 24 runs |
| `istar_thread_limit.patch` | Cap BLAS/OpenMP threads in iStar on shared qz nodes |
| `prepare_revise_reference.py`, `prepare_revise_pm.py` | Reconstruct missing REVISE inputs from the released subsets for a **non-exact** pilot |
| `run_revise.sh` | Run the upstream REVISE spot-size route on selected regions |
| `results/` | Small CSV/JSON summaries; no raw data, checkpoints, images, or cell-level predictions |

### Input contract

Set `BENCHMARK_ROOT` to the **data/work directory**, separate from this code
checkout. The tested qz directory is
`/inspire/hdd/global_user/zhangyinghao-240107100021/REVISE_benchmark`.

```text
$BENCHMARK_ROOT/
  spot/part{1,2,3}/
    selected_xenium.h5ad
    real_sc_ref_part.h5ad
    spot_{50,100,150,200}/xenium_spot.h5ad
  Xenium_V1_Human_Colon_Cancer_P1_CRC_Add_on_FFPE_he_image.ome.tif
  Xenium_V1_Human_Colon_Cancer_P1_CRC_Add_on_FFPE_he_imagealignment.csv
  methods/{REVISE,istar,TESLA}/
  envs/{baselines,revise}/
```

The full-resolution H&E image and its 10x alignment CSV are essential for
TESLA/iStar. `prepare_inputs.py` fits the alignment keypoints, crops each
region with a 250 µm margin, and rescales it to 0.5 µm/pixel. The fitted
landmark residual on the tested file was about 2.87 original H&E pixels.
`image_transform.json` records the transform; verify it visually for any new
sample. Each case includes the 324-gene spot counts and spot centers.
`selected_xenium.h5ad` supplies the held-out cell coordinates **only during
scoring** for TESLA/iStar; its expression is never passed to these models.

### Tested source revisions and environment

| Repository | Commit |
| --- | --- |
| [REVISE](https://github.com/wuys13/REVISE) | `c83dc97d25b6d513b59cc301255e5bdc7e9c7cd9` |
| [iStar](https://github.com/daviddaiweizhang/istar) | `3cb0e5352a86df337f41c5841d0808da0003457b` |
| [TESLA](https://github.com/jianhuupenn/TESLA) | `8c4dcf896497dd4a34a6e720f03ec68c1502d571` |

The tested qz environments used Python 3.10, PyTorch 2.5.1, NumPy 1.26.4,
AnnData 0.11.4, Scanpy 1.11.5, pandas 2.3.3, SciPy 1.15.3, scikit-image
0.25.2, OpenCV, and pyvips 3.1.1. The baseline venv additionally had
`einops==0.6.1`, `pytorch-lightning==2.0.8`, `torchvision==0.20.1`,
`lightning-utilities==0.15.2`, and `torchmetrics==0.11.4`. The REVISE venv
additionally had `POT==0.9.5`, `squidpy==1.6.5`, and `scikit-misc==0.5.2`.
The qz venvs inherit `/opt/conda/envs/spacec`; creating a fresh environment
elsewhere requires the upstream method dependencies and a working `libvips`.
Pin NumPy to 1.26.4 on the tested qz image: installing Squidpy without a pin
upgraded it to NumPy 2.x and broke preinstalled binary packages.

Clone and pin the baseline methods into `$BENCHMARK_ROOT/methods`. Apply
`istar_thread_limit.patch` before running cases concurrently. The iStar
checkout also needs the official HIPT `vit256_small_dino.pth` and
`vit4k_xs_dino.pth` checkpoints under `methods/istar/checkpoints/`;
their tested SHA-256 values are
`6960cd5a8657dc8bb214671aa0c6dbd3f5b698e84386884955836487ddc89e24`
and `2b0bd9e9a602a35f2bb3f76da39d2b53a91f23fc3f115dc59a63267d95ad2b7b`.
The original iStar Box links were unavailable during this run; the checkpoints
were downloaded from the [HIPT repository](https://github.com/mahmoodlab/HIPT).

## Run TESLA and iStar

The scripts resolve their own code directory. Run from `$BENCHMARK_ROOT`, or
export `BENCHMARK_ROOT` and invoke them from elsewhere. Set `BASELINE_PYTHON`,
`ISTAR_REPO`, `TESLA_REPO`, `HE_IMAGE`, or `HE_ALIGNMENT` if your layout differs.

```bash
export BENCHMARK_ROOT=/inspire/hdd/global_user/zhangyinghao-240107100021/REVISE_benchmark
cd "$BENCHMARK_ROOT"
SCRIPT="$BENCHMARK_ROOT/methods/REVISE/benchmark/spot_sr"

bash "$SCRIPT/run_all.sh" prepare
bash "$SCRIPT/run_all.sh" tesla
GPU=0 bash "$SCRIPT/run_all.sh" istar

# Limit a rerun to one region/size; completed cases are skipped.
PARTS=part2 SIZES=200 bash "$SCRIPT/run_all.sh" tesla
PARTS=part2 SIZES=200 GPU=0 bash "$SCRIPT/run_all.sh" istar

envs/baselines/bin/python "$SCRIPT/summarize.py" \
  --prepared prepared --output results/baselines_summary.csv
envs/baselines/bin/python "$SCRIPT/validate_results.py" \
  --prepared prepared --output results/validation.json
```

TESLA uses `log1p` spot counts, a 10 µm superpixel grid, and 10 neighbors,
then maps each held-out cell to its nearest superpixel. iStar uses the
official shifted HIPT features once per region, then `impute.py` with 400
epochs and 5 states for each spot size. A held-out cell samples the
corresponding 16-pixel image patch. The released pseudo-spots are square,
whereas iStar expects circular capture disks; the adapter uses half the
square width as a radius proxy. That geometry is an interpretation limit.

Predictions and per-gene scores are written to `$BENCHMARK_ROOT/prepared/`.
The scoring code follows REVISE's normalized metric protocol: normalize each
cell's library to 10,000, min-max scale each gene separately, then compute
PCC, SSIM, MSE, and NRMSE per gene. `validate_results.py` checks the full
3 × 4 × 2 matrix, cell/gene counts, finite scores, and score summaries.

## Run REVISE's spot route

The downloaded spot archive lacked the required `real_sc_ref_all.h5ad` and
root-level `PM_on_cell.csv`, and used `partN` rather than the upstream
`P1CRC/cut_partN` path. If you obtain the **original** full reference and PM,
use those instead. To run an explicitly labeled **adapted pilot**, create a
reference from the union of the three released subsets (15,280 unique cells)
and a synthetic, reference-frequency PM prior as in REVISE's Visium example:

```bash
export BENCHMARK_ROOT=/inspire/hdd/global_user/zhangyinghao-240107100021/REVISE_benchmark
cd "$BENCHMARK_ROOT"
SCRIPT="$BENCHMARK_ROOT/methods/REVISE/benchmark/spot_sr"
envs/revise/bin/python "$SCRIPT/prepare_revise_reference.py" --spot-root spot
envs/revise/bin/python "$SCRIPT/prepare_revise_pm.py" --spot-root spot --seed 0
PARTS='part1 part2 part3' bash "$SCRIPT/run_revise.sh"
```

The PM creation script refuses to overwrite an existing file unless given
`--force` and writes `PM_on_cell.provenance.json`. The reference helper
creates the required directory links. The input reconstruction is **not**
proven identical to the authors' original full reference and PM, so these
REVISE outputs cannot claim exact paper reproduction. In the upstream
benchmark route, `selected_xenium.h5ad` also supplies cell IDs, spot
membership, and cell coordinates before reconstruction; true cell-type labels
are attached for a match-rate diagnostic. TESLA/iStar receive H&E and spot
counts instead, so runtimes and accuracy have different input contracts.

`run_revise.sh` runs the upstream four-size configuration for each selected
region and writes logs to `$BENCHMARK_ROOT/logs/`. The two completed pilot
cases used the separate single-size YAMLs in this directory.

## Current results and runtime accounting

`results/baselines_summary.csv` contains 12 TESLA and 12 iStar case rows;
`results/validation.json` reports 24/24 valid results. On qz, TESLA's
recorded inference time had a median of 85.5 s and a range of 26.4–1686 s.
iStar's recorded per-case imputation-plus-evaluation time had a median of
278.5 s and a range of 215–1286 s; this excludes shared HIPT extraction,
roughly 0.5–2 min per region. The first five iStar cases ran without the
thread cap, and one hit 1286 s amid CPU contention. `annotate_execution.py`
labels those historical results; do not use it for new runs.

The adapted REVISE `part2/spot_200` pilot took 34 s for 6,292 cells;
`part1/spot_50` took 25 s for 37,837 cells, including startup, writing, and
evaluation. The full upstream four-size `part2` route took 58 s and passed.
The `part1/spot_50` pilot selected no graph-aggregation updates, so these
times cannot be extrapolated to all regions or larger datasets. The
TESLA/iStar runtimes have different inclusion boundaries from REVISE's wall
times. Record GPU, CPU thread cap, cache state, input size, and timing start/
stop points in future large runs.

Raw data, H&E images, HIPT weights, virtual-cell predictions, and complete
per-gene outputs are intentionally outside this repository. Keep them in the
data/work directory and archive their checksums and exact code revisions for
the large-scale benchmark.
