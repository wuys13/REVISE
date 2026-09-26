#!/usr/bin/env bash
# Start the ResolVI corrected-count export as soon as model training finishes.
set -euo pipefail

ROOT=${ROOT:-/inspire/hdd/global_user/zhangyinghao-240107100021/REVISE_benchmark/misseg}
RESULTS="$ROOT/results"
PY="$ROOT/envs/resolvi/bin/python"
OUTPUT="$RESULTS/resolvi/full/resolvi_corrected_counts.h5ad"
mkdir -p "$RESULTS/resolvi/full"
exec >> "$RESULTS/resolvi/full/corrected_export.log" 2>&1
trap 'status=$?; echo "$(date -Is) FAILED corrected export (exit $status)"; exit "$status"' ERR

echo "$(date -Is) Waiting for ResolVI training"
while [[ ! -e "$RESULTS/.resolvi.done" ]]; do
  if [[ -n "${PIPELINE_PID:-}" ]] && ! kill -0 "$PIPELINE_PID" 2>/dev/null; then
    echo "$(date -Is) Pipeline exited before ResolVI finished"
    exit 1
  fi
  sleep 30
done
if [[ -e "$RESULTS/.resolvi_corrected.done" ]]; then
  if [[ ! -s "$OUTPUT" ]]; then
    echo "Completion marker exists but corrected output is missing: $OUTPUT"
    exit 1
  fi
  echo "$(date -Is) Corrected output already complete"
  exit 0
fi

echo "$(date -Is) START corrected export"
CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-2} "$PY" "$ROOT/export_resolvi_corrected.py" \
  --input "$RESULTS/resolvi/full/resolvi_latent.h5ad" \
  --model "$RESULTS/resolvi/full/model" \
  --output "$OUTPUT" \
  --block-size 2048 --batch-size 256 --num-samples 3 --seed 42
touch "$RESULTS/.resolvi_corrected.done"
echo "$(date -Is) DONE corrected export"
