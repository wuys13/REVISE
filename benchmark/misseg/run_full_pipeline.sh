#!/usr/bin/env bash
# Run the three P1CRC methods in order after the H&E nuclei mask is ready.
set -euo pipefail

ROOT=${ROOT:-/inspire/hdd/global_user/zhangyinghao-240107100021/REVISE_benchmark/misseg}
RESULTS="$ROOT/results"
MASK="$RESULTS/proseg/p1crc_stardist_masks.npy.gz"
PY="$ROOT/../envs/revise/bin/python"
RESOLVI_PY="$ROOT/envs/resolvi/bin/python"
R="$ROOT/envs/split/bin/Rscript"
mkdir -p "$RESULTS"
exec >> "$RESULTS/full_pipeline.log" 2>&1
trap 'status=$?; echo "$(date -Is) FAILED (exit $status) at: $BASH_COMMAND"; exit "$status"' ERR
echo "$(date -Is) Waiting for full H&E nuclei mask"
while [[ ! -s "$MASK" || ! -s "${MASK%.gz}.json" ]]; do
  sleep 60
done

run_stage() {
  local name=$1
  shift
  if [[ -e "$RESULTS/.$name.done" ]]; then
    echo "$(date -Is) SKIP $name (done marker exists)"
    return
  fi
  echo "$(date -Is) START $name"
  "$@"
  touch "$RESULTS/.$name.done"
  echo "$(date -Is) DONE $name"
}

run_stage proseg "$PY" "$ROOT/run_proseg.py" \
  --binary "$ROOT/methods/proseg/target/release/proseg" \
  --binned-outputs "$ROOT/raw_p1crc/binned_outputs" \
  --mask "$MASK" --output "$RESULTS/proseg/full" --nthreads 48

run_stage aggregate "$PY" "$ROOT/aggregate_stardist_cells.py" \
  --binned-outputs "$ROOT/raw_p1crc/binned_outputs" \
  --mask "$MASK" --output "$RESULTS/common/star_dist_cells_5um.h5ad" \
  --max-expansion-um 5

run_stage resolvi "$RESOLVI_PY" "$ROOT/run_resolvi.py" \
  --input "$RESULTS/common/star_dist_cells_5um.h5ad" \
  --output "$RESULTS/resolvi/full" --epochs 100 --batch-size 256

run_stage split_inputs "$PY" "$ROOT/prepare_split_inputs.py" \
  --spatial "$RESULTS/common/star_dist_cells_5um.h5ad" \
  --reference "$ROOT/adata_sc_all_reanno.h5ad" \
  --output "$RESULTS/split/full_input" \
  --max-spatial-cells 10000 --max-reference-per-type 500 --seed 42

run_stage split "$R" "$ROOT/run_split.R" \
  "$RESULTS/split/full_input" "$RESULTS/split/full" 16

run_stage split_h5ad "$PY" "$ROOT/convert_split_output.py" \
  --split-output "$RESULTS/split/full" \
  --split-input "$RESULTS/split/full_input" \
  --output "$RESULTS/split/full/split_purified.h5ad"

echo "$(date -Is) Full P1CRC pipeline finished"
