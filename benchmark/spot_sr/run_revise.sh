#!/usr/bin/env bash
# Run REVISE's upstream four-size spot route for each selected P1CRC region.
set -euo pipefail

ROOT=$(cd "${BENCHMARK_ROOT:-$PWD}" && pwd -P)
REVISE=${REVISE_REPO:-"$ROOT/methods/REVISE"}
PY=${REVISE_PYTHON:-"$ROOT/envs/revise/bin/python"}
PARTS=${PARTS:-'part1 part2 part3'}
OUTPUT=${REVISE_OUTPUT_ROOT:-"$ROOT/results/revise"}
CONFIG=${REVISE_CONFIG:-"$REVISE/configs/benchmark/spot_size.yaml"}

for required in "$ROOT/spot/PM_on_cell.csv" "$CONFIG" "$PY"; do
  [[ -e "$required" ]] || { echo "Missing required input: $required" >&2; exit 1; }
done
mkdir -p "$ROOT/logs" "$OUTPUT"
cd "$REVISE"

for part in $PARTS; do
  sample="P1CRC/cut_$part"
  for required in "$ROOT/spot/$sample/real_sc_ref_all.h5ad" \
                  "$ROOT/spot/$sample/selected_xenium.h5ad"; do
    [[ -e "$required" ]] || { echo "Missing required input: $required" >&2; exit 1; }
  done
  log="$ROOT/logs/revise_${part}.log"
  start=$(date +%s)
  echo "Running REVISE $sample (four spot sizes); log: $log"
  if env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=. \
    "$PY" -u reproduce/benchmark_main.py \
      --config "$CONFIG" --data-root "$ROOT/spot" \
      --sample-name "$sample" --output-root "$OUTPUT" \
      > "$log" 2>&1; then
    status=passed
  else
    status=failed
  fi
  end=$(date +%s)
  printf 'wall_seconds=%s status=%s\n' "$((end-start))" "$status" | tee -a "$log"
  [[ "$status" == passed ]] || { tail -40 "$log" >&2; exit 1; }
done
