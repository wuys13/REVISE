#!/usr/bin/env bash
# Verify the three full P1CRC outputs after both long-running branches finish.
set -euo pipefail

ROOT=${ROOT:-/inspire/hdd/global_user/zhangyinghao-240107100021/REVISE_benchmark/misseg}
RESULTS="$ROOT/results"
PY="$ROOT/../envs/revise/bin/python"
REPORT="$RESULTS/full_verification.json"
mkdir -p "$RESULTS"
exec >> "$RESULTS/full_verification.log" 2>&1
trap 'status=$?; echo "$(date -Is) FAILED full verification (exit $status)"; exit "$status"' ERR

echo "$(date -Is) Waiting for SPLIT and ResolVI corrected counts"
while [[ ! -e "$RESULTS/.split_h5ad.done" || ! -e "$RESULTS/.resolvi_corrected.done" ]]; do
  if [[ ! -e "$RESULTS/.split_h5ad.done" && -n "${PIPELINE_PID:-}" ]] \
      && ! kill -0 "$PIPELINE_PID" 2>/dev/null; then
    echo "$(date -Is) Pipeline exited before SPLIT finished"
    exit 1
  fi
  if [[ ! -e "$RESULTS/.resolvi_corrected.done" && -n "${CORRECTED_PID:-}" ]] \
      && ! kill -0 "$CORRECTED_PID" 2>/dev/null; then
    echo "$(date -Is) Corrected-count export exited before finishing"
    exit 1
  fi
  sleep 30
done

if [[ -e "$RESULTS/.verified.done" ]]; then
  test -s "$REPORT"
  echo "$(date -Is) Verification already complete"
  exit 0
fi

echo "$(date -Is) START full verification"
"$PY" "$ROOT/verify_full.py" --root "$ROOT"
test -s "$REPORT"
touch "$RESULTS/.verified.done"
echo "$(date -Is) DONE full verification"
