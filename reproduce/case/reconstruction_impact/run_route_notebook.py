"""Execute one route source and materialize its sole final notebook artifact.

The source lives in ``reproduce/case/reconstruction_impact``.  The final
executed notebook always lives at
``output/reconstruction_impact/<sample_id>/notebook/<route>.ipynb``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import pandas as pd
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
    configured_output = Path(output_root) if output_root else Path(config["output"]["dir"])
    if not configured_output.is_absolute():
        configured_output = (ROOT / configured_output).resolve()
    final_relative = Path(config["output"]["final_notebook"])
    if final_relative.parent != Path("notebook") or final_relative.suffix != ".ipynb":
        raise ValueError("output.final_notebook must be notebook/<route>.ipynb")
    return ROOT / "reproduce/case/reconstruction_impact" / source_name, configured_output, configured_output / final_relative


def verify_executed_copy(source: Path, executed: Path) -> dict:
    """Reject an output notebook whose executable source differs from its input."""
    source_notebook = json.loads(source.read_text(encoding="utf-8"))
    executed_notebook = json.loads(executed.read_text(encoding="utf-8"))
    if [c["cell_type"] for c in source_notebook["cells"]] != [c["cell_type"] for c in executed_notebook["cells"]]:
        raise RuntimeError("Executed notebook cell types or positions differ from source")
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
    source_ids = [cell.get("id") for cell in source_notebook["cells"]]
    executed_ids = [cell.get("id") for cell in executed_notebook["cells"]]
    if not all(source_ids) or len(set(source_ids)) != len(source_ids) or source_ids != executed_ids:
        raise RuntimeError("Executed notebook cell identities differ from source")
    code_cells = [cell for cell in executed_notebook["cells"] if cell["cell_type"] == "code"]
    if any(cell.get("execution_count") is None for cell in code_cells if cell["source"]):
        raise RuntimeError("Executed notebook contains unexecuted code cells")
    stderr = []
    for cell in code_cells:
        for output in cell.get("outputs", []):
            if output.get("output_type") == "stream" and output.get("name") == "stderr":
                stderr.append({"cell_id": cell["id"], "text": "".join(output.get("text", []))})
    return {"executed_code_cells": len(code_cells), "stderr": stderr}


def append_execution_artifacts(output_root: Path, final: Path, audit_path: Path) -> Path:
    """Append the final notebook and execution audit to the artifact index."""
    index_path = output_root / "analysis" / "artifacts.csv"
    index_path.parent.mkdir(parents=True, exist_ok=True)
    required_columns = ["path", "size_bytes", "sha256"]
    rows = []
    fieldnames = list(required_columns)
    if index_path.exists():
        with index_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            fieldnames = list(reader.fieldnames or required_columns)
            for column in required_columns:
                if column not in fieldnames:
                    fieldnames.append(column)
            rows = [
                {column: row.get(column, "") for column in fieldnames}
                for row in reader
            ]

    additions = [final, audit_path]
    report = output_root / "report.html"
    if report.is_file():
        additions.append(report)
    additions.extend(sorted((output_root / "analysis" / "conclusions").glob("*")))
    replacements = {p.relative_to(output_root).as_posix() for p in additions}
    rows = [row for row in rows if row.get("path") not in replacements]
    for artifact in additions:
        if not artifact.is_file():
            raise FileNotFoundError(f"Cannot index missing execution artifact: {artifact}")
        try:
            relative_path = artifact.relative_to(output_root).as_posix()
        except ValueError as exc:
            raise ValueError(f"Execution artifact is outside output root: {artifact}") from exc
        rows.append({
            **{column: "" for column in fieldnames},
            "path": relative_path,
            "size_bytes": str(artifact.stat().st_size),
            "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest(),
        })

    with index_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return index_path


def attach_conclusions(final: Path, package: dict, *, formatter=None) -> None:
    """Populate presentation placeholders from this run without editing code."""
    if formatter is None:
        from reproduce.case.reconstruction_impact.report_records import format_notebook_records
        formatter = format_notebook_records
    notebook = json.loads(final.read_text(encoding="utf-8"))
    for cell in notebook["cells"]:
        metadata = cell.get("metadata", {})
        node_id = metadata.get("impact_node")
        section = metadata.get("impact_records")
        if node_id:
            try:
                cell["source"] = formatter(package, node_id=node_id)
            except TypeError:
                # Keep compatibility with the pre-node formatter while the
                # report reader is upgraded.
                cell["source"] = formatter(package, layer=None)
        elif section:
            layer = None if section == "overview" else int(section)
            cell["source"] = formatter(package, layer=layer)
        if cell["cell_type"] == "markdown":
            source = cell["source"]
            text = source if isinstance(source, str) else "".join(source)
            method_root = os.path.relpath(ROOT / "docs/design/reconstruction-impact", final.parent)
            cell["source"] = text.replace(
                "](../../../docs/design/reconstruction-impact/",
                f"]({method_root}/",
            )
    final.write_text(json.dumps(notebook, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")


def publish_report(output_root: Path, final: Path) -> dict:
    """Build reading artifacts exclusively from the completed scientific output."""
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))
    from reproduce.case.reconstruction_impact.report_records import build_records, save_records, render_report
    scientific_audit = json.loads((output_root / "analysis" / "audit.json").read_text())
    package = build_records(output_root, scientific_audit["sample_id"], scientific_audit["route"])
    records_dir = output_root / "analysis" / "conclusions"
    review_path = records_dir / "review.json"
    review = json.loads(review_path.read_text()) if review_path.is_file() else None
    save_records(package, records_dir, review_manifest=review)
    package = json.loads((records_dir / "report_records.json").read_text(encoding="utf-8"))
    render_report(package, output_root / "report.html")
    attach_conclusions(final, package)
    return package


def local_source_digests() -> dict:
    return {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(Path(__file__).parent.glob("*.py"))
            if path.name not in {"report_records.py", "report_html.py", "report_interpretations.py"}}


def report_source_digests() -> dict:
    """Include the actual report template alongside its readers in provenance."""
    paths = [
        "reproduce/case/reconstruction_impact/report_records.py",
        "reproduce/case/reconstruction_impact/report_html.py",
        "reproduce/case/reconstruction_impact/report_interpretations.py",
        "reproduce/case/reconstruction_impact/content_contract.py",
        "docs/design/reconstruction-impact/templates/report.html",
        "docs/design/reconstruction-impact/templates/notebook-block.md",
    ]
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            for name in paths if (ROOT / name).is_file()}



def route_source_digests(route: str) -> dict:
    """Identify the configuration and scientific implementations used by this route."""
    paths = [ROUTES[route][0], *[
        "revise/" + name + ".py" for name in (
            "analysis/reconstruction_impact", "analysis/basic/spatial_region",
            "analysis/basic/partition_change", "analysis/paired_moran",
            "analysis/advanced/aucell", "analysis/basic/gene_set_scoring",
            "application/ist_assembly",
        )
    ]]
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            if (ROOT / name).is_file() else "missing" for name in paths}


def write_run_status(output_root: Path, status: str, reason: str = "") -> None:
    """Do not present a previous report as the result of a running or failed attempt."""
    (output_root / "report.html").write_text(
        '<!doctype html><meta charset="utf-8"><title>Impact run status</title>'
        f'<h1>Report status: {html.escape(status)}</h1><p>{html.escape(reason)}</p>'
        '<p>Previous conclusions are not current for this attempt.</p>', encoding="utf-8")
    directory = output_root / "analysis/conclusions"
    path = directory / "report_records.json"
    if path.is_file():
        package = json.loads(path.read_text())
        review = {"status": "pending", "not_current": True, "reason": status}
        package["review"] = review
        for record in package.get("records", []):
            record["review"] = review
        path.write_text(json.dumps(package, ensure_ascii=False, indent=2) + "\n")
        (directory / "review.json").write_text(json.dumps(review) + "\n")
        csv_path = directory / "report_records.csv"
        if csv_path.is_file():
            with csv_path.open(newline="") as handle:
                reader = csv.DictReader(handle); fields = reader.fieldnames; rows = list(reader)
            for row in rows:
                row["review_status"] = "pending"
            with csv_path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(rows)


def verify_fresh_evidence(output_root: Path, started_ns: int) -> list[str]:
    """Reject retained scientific files that were not regenerated by this full run."""
    paths = []
    for directory in (output_root / "analysis", output_root / "figures"):
        for path in directory.rglob("*"):
            if not path.is_file():
                continue
            relative = path.relative_to(output_root)
            if relative.parts[:2] in (("analysis", "conclusions"), ("analysis", "integration")) or relative.as_posix() == "analysis/artifacts.csv":
                continue
            paths.append(path)
    stale = [str(p.relative_to(output_root)) for p in paths if p.stat().st_mtime_ns < started_ns]
    if not paths or stale:
        raise RuntimeError(f"Scientific evidence was not regenerated: {stale or ['no evidence']}")
    return [str(p.relative_to(output_root)) for p in paths]


def _dataframe_digest(frame: pd.DataFrame) -> str:
    """Hash a CSV after the same pandas read used by the baseline contract."""

    values = pd.util.hash_pandas_object(frame, index=True).values.tobytes()
    return hashlib.sha256(values).hexdigest()


def validate_baseline(output_root: Path, baseline_path: str | os.PathLike[str] | None) -> dict:
    """Validate old scientific tables before any report publication."""

    if baseline_path is None:
        return {"status": "skipped"}
    path = Path(baseline_path)
    try:
        baseline = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"baseline could not be read: {path}") from exc
    sample_id = output_root.name
    if sample_id not in baseline:
        raise RuntimeError(f"baseline has no sample entry for {sample_id}")
    expected_tables = baseline[sample_id].get("tables", {})
    if not isinstance(expected_tables, dict):
        raise RuntimeError(f"baseline tables for {sample_id} are not a mapping")
    checked = []
    for relative_name, expected in expected_tables.items():
        relative = Path(relative_name)
        if relative.is_absolute() or ".." in relative.parts:
            raise RuntimeError(f"baseline table path is unsafe: {relative_name}")
        table_path = output_root / relative
        if not table_path.is_file():
            raise RuntimeError(f"baseline table is missing: {relative_name}")
        try:
            frame = pd.read_csv(table_path, low_memory=False)
        except (OSError, EOFError, pd.errors.EmptyDataError, pd.errors.ParserError, UnicodeError) as exc:
            raise RuntimeError(f"baseline table could not be read: {relative_name}") from exc
        actual_columns = list(frame.columns)
        actual_rows = len(frame)
        actual_digest = _dataframe_digest(frame)
        if (
            actual_columns != expected.get("columns")
            or actual_rows != expected.get("rows")
            or actual_digest != expected.get("digest")
        ):
            raise RuntimeError(
                f"baseline mismatch for {sample_id}/{relative_name}: "
                f"columns={actual_columns!r}, rows={actual_rows}, digest={actual_digest}"
            )
        checked.append(relative_name)
    return {
        "status": "passed",
        "sample_id": sample_id,
        "baseline": str(path),
        "checked_tables": checked,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("route", choices=sorted(ROUTES))
    parser.add_argument("--output-root", help="Override the configured output root")
    parser.add_argument("--timeout", type=int, default=0)
    parser.add_argument("--kernel", help="Installed Jupyter kernel; defaults to notebook metadata")
    parser.add_argument("--baseline", help="Optional scientific table baseline JSON")
    args = parser.parse_args()

    source, output_root, final = route_paths(args.route, args.output_root)
    source_digest = hashlib.sha256(source.read_bytes()).hexdigest()
    helper_digests = local_source_digests()
    dependency_digests = route_source_digests(args.route)
    started_ns = time.time_ns()
    final.parent.mkdir(parents=True, exist_ok=True)
    environment = os.environ | {
        "REVISE_REPOSITORY_ROOT": str(ROOT),
        "REVISE_ANALYSIS_OUTPUT_ROOT": str(output_root),
        "TQDM_DISABLE": "1",
        "MPLBACKEND": "module://matplotlib_inline.backend_inline",
    }
    environment.pop("RECONSTRUCTION_IMPACT_CACHE_PATH", None)
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
    if args.kernel:
        command.append(f"--ExecutePreprocessor.kernel_name={args.kernel}")
    audit_path = final.parent / "execution.json"
    audit_path.write_text(json.dumps({"status": "running", "source_sha256": source_digest}) + "\n")
    try:
        write_run_status(output_root, "running")
        subprocess.run(command, cwd=ROOT, env=environment, check=True)
        if hashlib.sha256(source.read_bytes()).hexdigest() != source_digest:
            raise RuntimeError("Source notebook changed during execution; rerun the final source")
        if local_source_digests() != helper_digests:
            raise RuntimeError("Local notebook helpers changed during execution; rerun frozen sources")
        if route_source_digests(args.route) != dependency_digests:
            raise RuntimeError("Route configuration or scientific implementation changed during execution")
        written_paths = verify_fresh_evidence(output_root, started_ns)
        audit = verify_executed_copy(source, final)
        baseline_audit = validate_baseline(output_root, args.baseline)
        report_digest = report_source_digests()
        publish_report(output_root, final)
        if report_source_digests() != report_digest:
            raise RuntimeError("Report implementation changed during publication")
        audit["report_source_sha256"] = report_digest
        verify_executed_copy(source, final)
        audit.update({
            "status": "completed",
            "source": str(source.relative_to(ROOT)),
            "source_sha256": source_digest,
            "executed_sha256": hashlib.sha256(final.read_bytes()).hexdigest(),
            "checkpoint_used": False,
            "local_source_sha256": helper_digests,
            "scientific_source_sha256": dependency_digests,
            "run_started_ns": started_ns,
            "fresh_evidence": written_paths,
            "baseline": baseline_audit,
        })
        audit_path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
        append_execution_artifacts(output_root, final, audit_path)
    except (Exception, KeyboardInterrupt) as exc:
        reason = str(exc) or "interrupted"
        audit_path.write_text(json.dumps({"status": "failed", "source_sha256": source_digest, "run_started_ns": started_ns, "reason": reason}) + "\n")
        write_run_status(output_root, "failed", reason)
        raise
    print(final)


if __name__ == "__main__":
    main()
