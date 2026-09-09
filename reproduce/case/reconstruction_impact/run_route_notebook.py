"""Execute one route source and materialize its sole final notebook artifact.

The source lives in ``reproduce/case/reconstruction_impact``.  The final
executed notebook always lives at
``output/reconstruction_impact/<sample_id>/notebook/<route>.ipynb``.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

import yaml


ROOT = Path(__file__).resolve().parents[3]
ROUTES = {
    "visiumhd": (
        "configs/analysis/reconstruction_impact_visiumhd_p1crc.yaml",
        "VisiumHD_sp_SVC_Reconstruction_Impact.ipynb",
    ),
    "xenium": (
        "configs/analysis/reconstruction_impact_xenium_p2crc_fibroblast.yaml",
        "Xenium_sc_SVC_Fibroblast_Reconstruction_Impact.ipynb",
    ),
}


def route_paths(route: str, output_root: str | None = None) -> tuple[Path, Path, Path]:
    """Return source, output root, and the configured sole final notebook path."""
    config_relative, source_name = ROUTES[route]
    with (ROOT / config_relative).open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    configured_output = Path(output_root or config["output"]["dir"])
    if not configured_output.is_absolute():
        configured_output = (ROOT / configured_output).resolve()
    final_relative = Path(config["output"]["final_notebook"])
    if final_relative.parent != Path("notebook") or final_relative.suffix != ".ipynb":
        raise ValueError("output.final_notebook must be notebook/<route>.ipynb")
    return ROOT / "reproduce/case/reconstruction_impact" / source_name, configured_output, configured_output / final_relative


def verify_executed_copy(source: Path, executed: Path) -> None:
    """Reject an output notebook whose executable source differs from its input."""
    source_notebook = json.loads(source.read_text(encoding="utf-8"))
    executed_notebook = json.loads(executed.read_text(encoding="utf-8"))
    source_code = [cell["source"] for cell in source_notebook["cells"] if cell["cell_type"] == "code"]
    executed_code = [cell["source"] for cell in executed_notebook["cells"] if cell["cell_type"] == "code"]
    if source_code != executed_code:
        raise RuntimeError("Executed notebook code cells do not match the route source")
    errors = [
        output
        for cell in executed_notebook["cells"]
        if cell["cell_type"] == "code"
        for output in cell.get("outputs", [])
        if output.get("output_type") == "error"
    ]
    if errors:
        raise RuntimeError("Executed notebook contains error output")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("route", choices=sorted(ROUTES))
    parser.add_argument("--output-root", help="Override the configured output root")
    parser.add_argument("--timeout", type=int, default=0)
    parser.add_argument("--checkpoint", help="Restore a bounded workflow checkpoint before rendering")
    args = parser.parse_args()

    source, output_root, final = route_paths(args.route, args.output_root)
    final.parent.mkdir(parents=True, exist_ok=True)
    environment = os.environ | {
        "REVISE_REPOSITORY_ROOT": str(ROOT),
        "REVISE_ANALYSIS_OUTPUT_ROOT": str(output_root),
    }
    if args.checkpoint:
        environment["RECONSTRUCTION_IMPACT_CACHE_PATH"] = str(Path(args.checkpoint).resolve())
    command = [
        sys.executable,
        "-m",
        "jupyter",
        "nbconvert",
        "--to",
        "notebook",
        "--execute",
        str(source),
        "--output",
        final.stem,
        "--output-dir",
        str(final.parent),
        f"--ExecutePreprocessor.timeout={args.timeout}",
    ]
    subprocess.run(command, cwd=ROOT, env=environment, check=True)
    verify_executed_copy(source, final)
    print(final)


if __name__ == "__main__":
    main()
