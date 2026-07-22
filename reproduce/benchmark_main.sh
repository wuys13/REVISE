#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository_root="$(cd "${script_dir}/.." && pwd)"
export PYTHONPATH="${repository_root}${PYTHONPATH:+:${PYTHONPATH}}"
exec "${BENCHMARK_LAUNCHER_PYTHON:-python}" -u "${repository_root}/revise/benchmark/launcher.py"
