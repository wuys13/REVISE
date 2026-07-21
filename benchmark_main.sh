#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "${BENCHMARK_LAUNCHER_PYTHON:-python}" -u "${script_dir}/benchmark_launcher.py"
