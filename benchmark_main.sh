#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="${script_dir}${PYTHONPATH:+:${PYTHONPATH}}"
exec "${BENCHMARK_LAUNCHER_PYTHON:-python}" -u "${script_dir}/revise/benchmark/launcher.py"
