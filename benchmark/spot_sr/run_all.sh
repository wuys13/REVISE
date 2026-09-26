#!/usr/bin/env bash
# Run from the benchmark data root, or set BENCHMARK_ROOT explicitly.
# Usage: bash methods/REVISE/benchmark/spot_sr/run_all.sh prepare|tesla|istar
set -euo pipefail

ROOT=$(cd "${BENCHMARK_ROOT:-$PWD}" && pwd -P)
PY=${BASELINE_PYTHON:-"$ROOT/envs/baselines/bin/python"}
SCRIPT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
PREPARED="$ROOT/prepared"
ISTAR=${ISTAR_REPO:-"$ROOT/methods/istar"}
TESLA=${TESLA_REPO:-"$ROOT/methods/TESLA"}
HE=${HE_IMAGE:-"$ROOT/Xenium_V1_Human_Colon_Cancer_P1_CRC_Add_on_FFPE_he_image.ome.tif"}
ALIGN=${HE_ALIGNMENT:-"$ROOT/Xenium_V1_Human_Colon_Cancer_P1_CRC_Add_on_FFPE_he_imagealignment.csv"}
MODE=${1:?Expected prepare, tesla, or istar}
PARTS=${PARTS:-'part1 part2 part3'}
SIZES=${SIZES:-'50 100 150 200'}
GPU=${GPU:-0}

mkdir -p "$ROOT/logs" "$PREPARED"
cd "$ROOT"

prepare() {
  for part in $PARTS; do
    if [[ ! -s "$PREPARED/$part/image_transform.json" ]]; then
      "$PY" "$SCRIPT/prepare_inputs.py" part --root "$ROOT" --output "$PREPARED" \
        --part "$part" --image "$HE" --alignment "$ALIGN"
    fi
    for size in $SIZES; do
      if [[ ! -s "$PREPARED/$part/spot_$size/case_manifest.json" ]]; then
        "$PY" "$SCRIPT/prepare_inputs.py" case --root "$ROOT" --output "$PREPARED" \
          --part "$part" --size "$size"
      fi
    done
  done
}

run_tesla() {
  for part in $PARTS; do
    for size in $SIZES; do
      case_dir="$PREPARED/$part/spot_$size"
      label="${part}_spot_${size}"
      if [[ -s "$case_dir/tesla_metrics_summary.json" ]]; then
        echo "TESLA $label already evaluated"
        continue
      fi
      echo "Starting TESLA $label"
      "$PY" -u "$SCRIPT/run_tesla.py" --case "$case_dir" --tesla-repo "$TESLA" \
        --grid-microns 10 --neighbors 10 > "$ROOT/logs/tesla_${label}.log" 2>&1
      "$PY" "$SCRIPT/evaluate.py" --case "$case_dir" --method tesla \
        >> "$ROOT/logs/tesla_${label}.log" 2>&1
      echo "Completed TESLA $label"
    done
  done
}

make_image_cache() {
  local part=$1 cache="$PREPARED/$1/image_cache/"
  mkdir -p "$cache"
  if [[ ! -e "${cache}he-raw.jpg" ]]; then
    ln -s ../he-raw.jpg "${cache}he-raw.jpg"
  fi
  printf '0.5\n' > "${cache}pixel-size-raw.txt"
  printf '0.5\n' > "${cache}pixel-size.txt"
  cd "$ISTAR"
  if [[ ! -s "${cache}he.jpg" ]]; then
    "$PY" rescale.py "$cache" --image
    "$PY" preprocess.py "$cache" --image
  fi
  if [[ ! -s "${cache}embeddings-hist.pickle" ]]; then
    CUDA_VISIBLE_DEVICES="$GPU" "$PY" -u extract_features.py "$cache" --device=cuda \
      > "$ROOT/logs/istar_features_${part}.log" 2>&1
  fi
  if [[ ! -s "${cache}mask-small.png" ]]; then
    "$PY" get_mask.py "${cache}embeddings-hist.pickle" "${cache}mask-small.png"
  fi
  cd "$ROOT"
}

run_istar() {
  for part in $PARTS; do
    make_image_cache "$part"
    cache="$PREPARED/$part/image_cache"
    for size in $SIZES; do
      case_dir="$PREPARED/$part/spot_$size"
      label="${part}_spot_${size}"
      if [[ -s "$case_dir/istar_metrics_summary.json" && -s "$case_dir/istar_run.json" ]]; then
        echo "iStar $label already evaluated"
        continue
      fi
      ln -sfn ../image_cache/he.jpg "$case_dir/he.jpg"
      ln -sfn ../image_cache/embeddings-hist.pickle "$case_dir/embeddings-hist.pickle"
      ln -sfn ../image_cache/mask-small.png "$case_dir/mask-small.png"
      cd "$ISTAR"
      "$PY" rescale.py "$case_dir/" --locs --radius \
        > "$ROOT/logs/istar_${label}.log" 2>&1
      echo "Starting iStar $label"
      start=$(date +%s)
      CUDA_VISIBLE_DEVICES="$GPU" "$PY" -u impute.py "$case_dir/" \
        --epochs=400 --n-states=5 --device=cuda \
        >> "$ROOT/logs/istar_${label}.log" 2>&1
      cd "$ROOT"
      "$PY" "$SCRIPT/evaluate.py" --case "$case_dir" --method istar \
        >> "$ROOT/logs/istar_${label}.log" 2>&1
      end=$(date +%s)
      "$PY" -c 'import json,sys; from pathlib import Path; p=Path(sys.argv[1]); p.write_text(json.dumps({"method":"iStar","epochs":400,"n_states":5,"shifted_histology_features":True,"runtime_seconds":int(sys.argv[2])-int(sys.argv[3])},indent=2)+"\n")' \
        "$case_dir/istar_run.json" "$end" "$start"
      echo "Completed iStar $label"
    done
  done
}

case "$MODE" in
  prepare) prepare ;;
  tesla) run_tesla ;;
  istar) run_istar ;;
  *) echo "Unknown mode: $MODE" >&2; exit 2 ;;
esac
