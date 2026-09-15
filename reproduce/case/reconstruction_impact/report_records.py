"""Standardized records and a small static report for reconstruction-impact outputs.

This module is deliberately a reader.  It consumes the CSV/JSON tables and
figures already published by the route notebooks; it does not load H5AD files,
rerun analysis, or introduce a second scientific calculation path.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any, Iterable, Mapping

import pandas as pd

from .content_contract import NODE_SPECS, QUESTION_LABELS, QUESTION_LAYERS, SCOPES
from .report_interpretations import build_node_narratives


SCHEMA_VERSION = "reconstruction-impact-records/v1"
QUESTIONS = tuple(QUESTION_LAYERS)
ROUTE_ALIASES = {
    "visiumhd": "sp_svc",
    "hd": "sp_svc",
    "sp_svc": "sp_svc",
    "xenium": "sc_svc",
    "sc_svc": "sc_svc",
}


def _jsonable(value: Any) -> Any:
    """Convert pandas/numpy values to deterministic JSON values."""

    if value is None or value is pd.NA or value is pd.NaT:
        return None
    if isinstance(value, Mapping):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(item) for item in value]
    if hasattr(value, "item") and not isinstance(value, (str, bytes)):
        try:
            return _jsonable(value.item())
        except (TypeError, ValueError):
            pass
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    return value


def _canonical(value: Any) -> bytes:
    return json.dumps(_jsonable(value), sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")


def _digest(value: Any) -> str:
    return hashlib.sha256(_canonical(value)).hexdigest()


def _without_review(package: Mapping[str, Any]) -> dict[str, Any]:
    """Return the part covered by the review digest.

    Review metadata is intentionally excluded so an approval can be attached
    without changing the digest that the reviewer signed.
    """

    def strip_review(value: Any) -> Any:
        if isinstance(value, Mapping):
            return {key: strip_review(item) for key, item in value.items() if key != "review"}
        if isinstance(value, list):
            return [strip_review(item) for item in value]
        if isinstance(value, tuple):
            return [strip_review(item) for item in value]
        return value

    copied = {key: value for key, value in package.items() if key not in {"review", "record_digest"}}
    return strip_review(copied)


def _as_number(value: Any) -> int | float | None:
    if value is None or value is pd.NA:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return int(number) if number.is_integer() else number


def _fact(value: Any) -> str:
    """Format a saved fact compactly for conclusions while keeping numbers intact in results."""

    value = _jsonable(value)
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def _as_bool(value: Any) -> bool | None:
    if value is None or value is pd.NA:
        return None
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"true", "1", "yes", "y"}:
        return True
    if text in {"false", "0", "no", "n"}:
        return False
    return None


def _first_non_none(*values: Any) -> Any:
    """Return the first value that is present, preserving valid zeroes."""

    for value in values:
        if value is not None and value is not pd.NA:
            try:
                if pd.isna(value):
                    continue
            except (TypeError, ValueError):
                pass
            return value
    return None


def _value(row: Mapping[str, Any], key: str, *, number: bool = False, boolean: bool = False) -> Any:
    value = row.get(key)
    if number:
        return _as_number(value)
    if boolean:
        return _as_bool(value)
    return _jsonable(value)


def _read_json(path: Path) -> tuple[Any, str | None]:
    if not path.is_file():
        return {}, "missing"
    try:
        return json.loads(path.read_text(encoding="utf-8")), None
    except (OSError, json.JSONDecodeError) as exc:
        return {}, f"invalid: {exc}"


def _read_csv(path: Path) -> tuple[pd.DataFrame, str | None]:
    if not path.is_file():
        return pd.DataFrame(), "missing"
    try:
        return pd.read_csv(path), None
    except (OSError, EOFError, pd.errors.EmptyDataError, pd.errors.ParserError, UnicodeError) as exc:
        return pd.DataFrame(), f"invalid: {exc}"


def _row(frame: pd.DataFrame, scope: str, *, scope_column: str = "scope") -> tuple[dict[str, Any], str | None]:
    if frame.empty or scope_column not in frame.columns:
        return {}, "missing"
    matches = frame.loc[frame[scope_column].astype(str).eq(scope)]
    if matches.empty:
        return {}, "missing"
    issue = "duplicate" if len(matches) > 1 else None
    return matches.iloc[0].to_dict(), issue


def _relative(path: Path, output_root: Path) -> str:
    try:
        return path.resolve().relative_to(output_root.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def _artifact(path: Path, output_root: Path) -> dict[str, Any]:
    artifact = {"path": _relative(path, output_root), "exists": path.is_file(), "sha256": None}
    if path.is_file():
        digest = hashlib.sha256()
        try:
            with path.open("rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
            artifact["sha256"] = digest.hexdigest()
        except OSError:
            artifact["sha256"] = None
    return artifact


def _condition(*issues: str | None, **values: Any) -> dict[str, Any]:
    status = "valid"
    clean_issues = [item for item in issues if item]
    if any(str(item).startswith("invalid") for item in clean_issues):
        status = "invalid"
    elif any(item == "missing" for item in clean_issues):
        status = "missing"
    elif clean_issues:
        status = "review"
    result = {"status": status}
    if clean_issues:
        result["issues"] = clean_issues
    result.update({key: _jsonable(value) for key, value in values.items()})
    return result


def _figure_candidates(output_root: Path, scope: str | None = None, *names: str) -> list[str]:
    found: list[str] = []
    for name in names:
        path = output_root / "figures" / name
        if path.is_file():
            found.append(_relative(path, output_root))
    if scope:
        slug = scope.lower().replace("_", "_")
        for name in names:
            path = output_root / "figures" / name.format(scope=slug)
            if path.is_file():
                found.append(_relative(path, output_root))
    return list(dict.fromkeys(found))


def _evidence(output_root: Path, tables: Iterable[str] = (), figures: Iterable[str] = (), condition: Mapping[str, Any] | None = None) -> dict[str, Any]:
    table_items = [_artifact(output_root / path, output_root) for path in tables]
    figure_items = [_relative(output_root / path, output_root) for path in figures]
    return {"tables": table_items, "figures": list(dict.fromkeys(figure_items)), "condition": dict(condition or {})}


def _make_record(
    *,
    sample_id: str,
    route_kind: str,
    run_identity: Mapping[str, Any],
    scope: str,
    question: str,
    results: Mapping[str, Any],
    judgment: str,
    condition: Mapping[str, Any],
    conclusion_zh: str,
    conclusion_en: str,
    limitations_zh: Iterable[str] = (),
    limitations_en: Iterable[str] = (),
    evidence: Mapping[str, Any] | None = None,
    headline_eligible: bool | None = None,
) -> dict[str, Any]:
    record_id = f"{sample_id}::{route_kind}::{scope}::{question}"
    rule_id = f"{SCHEMA_VERSION}/{question}"
    performance = {"judgment": judgment, "headline_eligible": headline_eligible}
    baseline = {"local_vs_raw_leiden": "raw_leiden", "local_vs_raw_level2": "raw_level2"}.get(question)
    return {
        "identity": {
            "record_id": record_id,
            "sample_id": sample_id,
            "route_kind": route_kind,
            "scope": scope,
            "task_cell_type": scope,
            "question": question,
            "baseline": baseline,
            "layer": QUESTION_LAYERS[question],
            "rule_id": rule_id,
            "spec_version": SCHEMA_VERSION,
        },
        "spec_version": SCHEMA_VERSION,
        "rule_id": rule_id,
        "results": _jsonable(dict(results)),
        "performance": performance,
        "performance_judgment": judgment,
        "evidence_condition": _jsonable(dict(condition)),
        "conclusion": {"zh": conclusion_zh, "en": conclusion_en},
        "limitations": {"zh": list(limitations_zh), "en": list(limitations_en)},
        "evidence": _jsonable(dict(evidence or {})),
        "review": {"status": "pending"},
        "run_identity": _jsonable(dict(run_identity)),
    }


def _resolve_roots(output_dir: str | os.PathLike[str], sample_id: str) -> tuple[Path, Path]:
    given = Path(output_dir).expanduser().resolve()
    candidates = [given]
    if given.name != "analysis":
        candidates.extend([given / "analysis", given / sample_id / "analysis"])
    for candidate in candidates:
        if candidate.is_dir() and candidate.name == "analysis":
            return candidate.parent, candidate
    raise ValueError(f"output_dir does not contain an analysis directory: {given}")


def _run_identity(output_root: Path, analysis: Path, sample_id: str, route_kind: str, audit: Mapping[str, Any], parameters: Mapping[str, Any]) -> dict[str, Any]:
    input_records = []
    for item in audit.get("inputs", []) if isinstance(audit.get("inputs", []), list) else []:
        if not isinstance(item, Mapping):
            continue
        input_records.append({"role": item.get("role"), "path": Path(str(item.get("path", ""))).name, "sha256": item.get("sha256")})
    return {
        "sample_id": sample_id,
        "route_kind": route_kind,
        "route": audit.get("route") or parameters.get("route"),
        "seed": _as_number(audit.get("seed", parameters.get("seed"))),
        "output_root": str(output_root),
        "cohort": audit.get("cohort"),
        "expression": audit.get("expression"),
        "region": audit.get("region"),
        "inputs": input_records,
        "record_builder": _builder_identity(),
        "analysis_files": {
            "audit": _artifact(analysis / "audit.json", output_root),
            "parameters": _artifact(analysis / "parameters.json", output_root),
        },
    }


def _summary_row(summary: pd.DataFrame, scope: str) -> tuple[dict[str, Any], str | None]:
    return _row(summary, scope)


def _moran_result(distribution: pd.DataFrame, summary_row: Mapping[str, Any], scope: str, gene_set: str) -> tuple[dict[str, Any], str | None]:
    subset = distribution.loc[
        distribution.get("scope", pd.Series(dtype=str)).astype(str).eq(scope)
        & distribution.get("gene_set", pd.Series(dtype=str)).astype(str).eq(gene_set)
    ] if not distribution.empty else pd.DataFrame()
    sides: dict[str, dict[str, Any]] = {}
    for side in ("Raw", "Reconstruction"):
        rows = subset.loc[subset["side"].astype(str).eq(side)] if not subset.empty and "side" in subset.columns else pd.DataFrame()
        row = rows.iloc[0].to_dict() if not rows.empty else {}
        sides[side.lower()] = {
            "n_valid": _value(row, "n_valid", number=True),
            "median": _value(row, "median", number=True),
            "q1": _value(row, "q1", number=True),
            "q75": _value(row, "q75", number=True),
        }
    raw_median = sides["raw"].get("median")
    reconstruction_median = sides["reconstruction"].get("median")
    median_difference = None if raw_median is None or reconstruction_median is None else _as_number(reconstruction_median - raw_median)
    delta = _value(summary_row, "paired_delta_median_moran", number=True)
    if delta is None and gene_set == "shared_valid":
        delta = _value(summary_row, "paired_delta_median", number=True)
    paired_n = _value(summary_row, "paired_delta_n", number=True) if gene_set == "shared_valid" else None
    paired_median = delta if gene_set == "shared_valid" else None
    results = {
        "gene_set": gene_set,
        "comparison_kind": "paired_shared_valid" if gene_set == "shared_valid" else "marginal_side_specific",
        "raw": sides["raw"],
        "reconstruction": sides["reconstruction"],
        "paired_delta": {"n": paired_n, "median": paired_median},
        "median_difference": median_difference if gene_set == "all_valid" else None,
        "missing_status_counts": {
            "raw_unmeasured": _value(summary_row, "raw_unmeasured_n", number=True),
            "raw_not_computable": _value(summary_row, "raw_not_computable_n", number=True),
            "raw_insufficient_support": _value(summary_row, "raw_insufficient_support_n", number=True),
            "reconstruction_unmeasured": _value(summary_row, "reconstruction_unmeasured_n", number=True),
            "reconstruction_not_computable": _value(summary_row, "reconstruction_not_computable_n", number=True),
            "reconstruction_insufficient_support": _value(summary_row, "reconstruction_insufficient_support_n", number=True),
        },
    }
    status = "valid" if sides["raw"]["n_valid"] is not None and sides["reconstruction"]["n_valid"] is not None else "missing"
    return results, status


def _threshold(path: Path, window_metrics_path: Path, n_bootstrap: Any) -> dict[str, Any]:
    value, issue = _read_json(path)
    if isinstance(value, list):
        value = value[0] if value else {}
    if not isinstance(value, Mapping):
        value = {}
        issue = issue or "invalid threshold payload"
    status = value.get("status")
    threshold = _as_number(value.get("threshold"))
    ci_lower = _as_number(value.get("ci_lower"))
    ci_upper = _as_number(value.get("ci_upper"))
    ci_width = None if ci_lower is None or ci_upper is None else _as_number(ci_upper - ci_lower)
    relative = None
    neff_range = None
    neff_min = None
    neff_max = None
    window_issue = None
    if window_metrics_path.is_file():
        try:
            metrics = pd.read_csv(window_metrics_path, usecols=["valid_window", "neff_recon"])
            valid_mask = metrics["valid_window"].map(lambda item: _as_bool(item) is True)
            valid = pd.to_numeric(metrics.loc[valid_mask, "neff_recon"], errors="coerce").dropna()
            if not valid.empty:
                neff_min = _as_number(float(valid.min()))
                neff_max = _as_number(float(valid.max()))
                neff_range = _as_number(float(neff_max - neff_min))
        except (OSError, ValueError, KeyError, pd.errors.EmptyDataError, pd.errors.ParserError) as exc:
            window_issue = f"invalid: {exc}"
    else:
        window_issue = "missing"
    valid_bootstrap = _as_number(value.get("n_valid_bootstrap"))
    total_bootstrap = _as_number(n_bootstrap)
    valid_fraction = None
    if valid_bootstrap is not None and total_bootstrap not in {None, 0}:
        valid_fraction = _as_number(float(valid_bootstrap) / float(total_bootstrap))
    if ci_width is not None and neff_range not in {None, 0}:
        relative = _as_number(ci_width / abs(float(neff_range)))
    reasons: list[str] = []
    if status in {"ok", "no_stable_threshold"}:
        if relative is None:
            reasons.append("threshold_not_identifiable")
        elif relative > 0.25:
            reasons.append("ci_too_wide")
        if valid_fraction is None or valid_fraction < 0.8:
            reasons.append("bootstrap_insufficient")
        analysis_status = "stable" if status == "ok" and not reasons and threshold is not None else "valid_unstable"
    else:
        analysis_status = "not_computable"
    reason = "; ".join(reasons) if reasons else None
    return {
        "status": status,
        "analysis_status": analysis_status,
        "reason": reason,
        "reasons": reasons,
        "threshold": threshold,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper,
        "ci_width": ci_width,
        "relative_ci_width": relative,
        "neff_recon_min": neff_min,
        "neff_recon_max": neff_max,
        "neff_recon_range": neff_range,
        "n_windows": _as_number(value.get("n_windows")),
        "n_valid_bootstrap": valid_bootstrap,
        "n_total_bootstrap": total_bootstrap,
        "valid_bootstrap_fraction": valid_fraction,
        "issue": window_issue or issue,
    }


def _first_existing(paths: Iterable[Path]) -> Path | None:
    for path in paths:
        if path.is_file():
            return path
    return None


def _observation_profile(path: Path, output_root: Path | None = None, scope: str | None = None) -> dict[str, Any]:
    """Read only the saved observation carrier needed for foundation metadata."""

    display_path = _relative(path, output_root) if output_root is not None else path.as_posix()
    if not path.is_file():
        return {"status": "missing", "path": display_path}
    try:
        frame = pd.read_csv(path)
    except (OSError, EOFError, pd.errors.EmptyDataError, pd.errors.ParserError, UnicodeError) as exc:
        return {"status": "invalid", "path": display_path, "error": str(exc)}
    id_column = _first_existing_column(frame.columns, ("unit_id", "observation_id", "cell_id"))
    coordinate_columns = [column for column in ("x", "y") if column in frame.columns]
    included = frame["included"] if "included" in frame.columns else pd.Series(dtype=bool)
    excluded_n = int((~included.fillna(False).astype(bool)).sum()) if not included.empty else None
    return {
        "status": "valid",
        "path": display_path,
        "n_rows": int(len(frame)),
        "included_n": int(included.map(_as_bool).eq(True).sum()) if not included.empty else None,
        "id_column": id_column,
        "coordinate_columns": coordinate_columns,
        "excluded_n": excluded_n,
        "reason_values": sorted(frame["reason"].dropna().astype(str).unique().tolist()) if "reason" in frame.columns else [],
        "selection": _selection_profile(frame, scope or "All"),
    }


def _selection_profile(frame: pd.DataFrame, scope: str) -> dict[str, Any]:
    """Summarize parent selection from the saved full-context membership table."""

    result: dict[str, Any] = {
        "status": "missing",
        "scope": scope,
        "scope_basis": "raw_level1_exact",
        "raw_parent_label": "Mono/Macro" if scope == "Mono_Macro" else scope,
        "raw_parent_n": None,
        "reconstruction_covered_parent_n": None,
        "eligible_n": None,
        "sampled_n": None,
        "stage_exclusions": {"outside_raw_level1": None, "not_reconstructed": None, "not_sampled": None},
    }
    if frame.empty or "included" not in frame.columns:
        return result
    if scope == "All" or "raw_level1" not in frame.columns:
        parent = frame
        result["scope_basis"] = "all_scopes"
        result["raw_parent_label"] = None
        result["outside_raw_level1_n"] = 0
    else:
        label = result["raw_parent_label"]
        parent_mask = frame["raw_level1"].astype(str).eq(str(label))
        parent = frame.loc[parent_mask]
        result["outside_raw_level1_n"] = int((~parent_mask).sum())
    result["raw_parent_n"] = int(len(parent))
    if "reconstruction_covered" in parent.columns:
        covered = parent["reconstruction_covered"].map(_as_bool)
    elif "reason" in parent.columns:
        covered = ~parent["reason"].astype(str).eq("not_reconstructed")
    else:
        return result
    result["reconstruction_covered_parent_n"] = int(covered.eq(True).sum())
    result["eligible_n"] = result["reconstruction_covered_parent_n"]
    result["sampled_n"] = int(parent["included"].map(_as_bool).eq(True).sum())
    result["stage_exclusions"] = {
        "outside_raw_level1": result.get("outside_raw_level1_n"),
        "not_reconstructed": int(covered.eq(False).sum()),
        "not_sampled": int((covered.eq(True) & ~parent["included"].map(_as_bool).eq(True)).sum()),
    }
    result["status"] = "valid"
    return result


def _table_status(frame: pd.DataFrame, issue: str | None) -> str:
    if issue == "missing":
        return "missing"
    if issue:
        return "invalid"
    return "empty" if frame.empty else "valid"


def _table_rows(frame: pd.DataFrame) -> list[dict[str, Any]]:
    return [{str(key): _jsonable(value) for key, value in row.items()} for row in frame.to_dict("records")]


def _table_facts(frame: pd.DataFrame, issue: str | None, path: Path, output_root: Path) -> dict[str, Any]:
    return {
        "status": _table_status(frame, issue),
        "path": _relative(path, output_root),
        "n_rows": None if issue else int(len(frame)),
        "rows": _table_rows(frame),
    }


def _scope_frame(frame: pd.DataFrame, scope: str) -> pd.DataFrame:
    if frame.empty or "scope" not in frame.columns:
        return frame
    return frame.loc[frame["scope"].astype(str).eq(scope)]


def _anatomy_summary(frame: pd.DataFrame, issue: str | None, scope: str, baseline: str) -> list[dict[str, Any]]:
    """Keep baseline-specific anatomy medians from the saved summary table."""

    if issue or frame.empty:
        return []
    scoped = _scope_frame(frame, scope)
    if scoped.empty:
        return []
    suffix = "raw" if baseline == "raw_leiden" else "level2"
    metrics = ("k_obs", "neff", "evenness")
    rows: list[dict[str, Any]] = []
    for source in scoped.to_dict("records"):
        row: dict[str, Any] = {
            "level1_region": _value(source, "level1_region"),
            "n_valid_windows": _value(source, "n_valid_windows", number=True),
            "scale_um": _value(source, "scale_um", number=True),
        }
        for metric in metrics:
            row[f"median_{metric}_{suffix}"] = _value(source, f"median_{metric}_{suffix}", number=True)
            row[f"median_{metric}_recon"] = _value(source, f"median_{metric}_recon", number=True)
            row[f"median_delta_{metric}_vs_{baseline}"] = _value(source, f"median_delta_{metric}_vs_{baseline}", number=True)
        rows.append(_jsonable(row))
    return rows


def _gene_side_facts(frame: pd.DataFrame, side: str) -> dict[str, Any]:
    side_frame = frame.loc[frame["side"].astype(str).str.lower().eq(side)] if not frame.empty and "side" in frame.columns else pd.DataFrame()
    status_counts: dict[str, int] = {}
    if "status" in side_frame.columns:
        for value in side_frame["status"].dropna().astype(str):
            status_counts[value] = status_counts.get(value, 0) + 1
    elif "provided" in side_frame.columns:
        for value in side_frame["provided"].map(_as_bool):
            label = "provided" if value is True else "not_provided" if value is False else "unknown"
            status_counts[label] = status_counts.get(label, 0) + 1
    expression_views: dict[str, int] = {}
    if "expression_view" in side_frame.columns:
        for value in side_frame["expression_view"].dropna().astype(str):
            expression_views[value] = expression_views.get(value, 0) + 1
    provided_n = None
    if "provided" in side_frame.columns:
        provided_n = int(side_frame["provided"].map(_as_bool).eq(True).sum())
    return {
        "n_rows": int(len(side_frame)),
        "n_genes": int(side_frame["gene_id"].dropna().astype(str).nunique()) if "gene_id" in side_frame.columns else None,
        "provided_n": provided_n,
        "status_counts": dict(sorted(status_counts.items())),
        "expression_view_counts": dict(sorted(expression_views.items())),
    }


def _gene_facts(
    frame: pd.DataFrame,
    issue: str | None,
    scope: str,
    path: Path,
    output_root: Path,
    *,
    gene_space: Any,
    union_gene_n: Any,
    partition_audit_values: Mapping[str, Any],
) -> dict[str, Any]:
    scoped = _scope_frame(frame, scope)
    status = _table_status(frame, issue)
    if issue is None and "scope" in frame.columns and scoped.empty:
        status = "missing"
    return {
        "status": status,
        "path": _relative(path, output_root),
        "scope": scope,
        "n_rows": None if issue or status == "missing" else int(len(scoped)),
        "gene_space": gene_space,
        "union_gene_n": union_gene_n,
        "n_shared_genes": _as_number(partition_audit_values.get("n_shared_genes")),
        "input_shared_genes": _as_number(partition_audit_values.get("input_shared_genes")),
        "excluded_genes": _as_number(partition_audit_values.get("excluded_raw_qc_genes")),
        "n_features": _as_number(partition_audit_values.get("n_features")),
        "feature_selection": _jsonable(partition_audit_values.get("feature_selection")),
        "by_side": {side: _gene_side_facts(scoped, side) for side in ("raw", "reconstruction")},
    }


_PAIRING_FIELDS = (
    "coordinate_source",
    "coordinate_unit",
    "rtol",
    "atol",
    "n_checked",
    "id_unique",
    "id_subset",
    "coordinates_match",
    "coordinate_finite",
)


def _pairing_facts(frame: pd.DataFrame, issue: str | None, scope: str, path: Path, output_root: Path) -> dict[str, Any]:
    result: dict[str, Any] = {
        "status": _table_status(frame, issue),
        "path": _relative(path, output_root),
        "scope": scope,
        "missing_columns": [],
    }
    result.update({field: None for field in _PAIRING_FIELDS})
    if issue:
        return result
    if frame.empty or "scope" not in frame.columns:
        result["status"] = "missing" if frame.empty else "invalid"
        if "scope" not in frame.columns:
            result["missing_columns"] = ["scope", *_PAIRING_FIELDS]
        return result
    missing_columns = [field for field in _PAIRING_FIELDS if field not in frame.columns]
    result["missing_columns"] = missing_columns
    matches = frame.loc[frame["scope"].astype(str).eq(scope)]
    if matches.empty:
        result["status"] = "missing"
        return result
    row = matches.iloc[0].to_dict()
    result["status"] = "review" if len(matches) > 1 else "valid"
    result["n_rows"] = int(len(matches))
    if missing_columns:
        result["status"] = "invalid"
        return result
    result["coordinate_source"] = _value(row, "coordinate_source")
    result["coordinate_unit"] = _value(row, "coordinate_unit")
    result["rtol"] = _value(row, "rtol", number=True)
    result["atol"] = _value(row, "atol", number=True)
    result["n_checked"] = _value(row, "n_checked", number=True)
    result["id_unique"] = _value(row, "id_unique", boolean=True)
    result["id_subset"] = _value(row, "id_subset", boolean=True)
    result["coordinates_match"] = _value(row, "coordinates_match", boolean=True)
    result["coordinate_finite"] = _value(row, "coordinate_finite", boolean=True)
    if "checked_basis" in row:
        result["checked_basis"] = _value(row, "checked_basis")
    return result


def _saved_preprocessing_fact(
    row: Mapping[str, Any],
    issue: str | None,
    path: Path,
    output_root: Path,
    fields: Iterable[str],
    *,
    status_source: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    source = row if status_source is None else status_source
    status = "valid" if source else "missing" if issue == "missing" else "invalid" if issue else "empty"
    return {
        "status": status,
        "path": _relative(path, output_root),
        "fields": {field: _jsonable(row.get(field)) for field in fields},
    }


def _foundation_facts(
    *,
    scope: str,
    canonical_route: str,
    audit: Mapping[str, Any],
    observation_profile: Mapping[str, Any],
    partition_audit_values: Mapping[str, Any],
    carriers: pd.DataFrame,
    carriers_issue: str | None,
    source_files: pd.DataFrame,
    source_files_issue: str | None,
    gene_availability: pd.DataFrame,
    gene_availability_issue: str | None,
    pairing_audit: pd.DataFrame,
    pairing_audit_issue: str | None,
    output_root: Path,
    analysis: Path,
    input_units: Any,
    paired_units: Any,
    excluded_units: Mapping[str, Any],
    gene_space: Any,
    union_gene_n: Any,
    preprocessing: Mapping[str, Any],
    analysis_denominators: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "carrier": _table_facts(carriers, carriers_issue, analysis / "inputs/carriers.csv", output_root),
        "source_files": _table_facts(source_files, source_files_issue, analysis / "inputs/source_files.csv", output_root),
        "cohort_audit": {
            "status": "valid" if audit.get("cohort") is not None else "missing",
            "label": _jsonable(audit.get("cohort")),
            "scope": scope,
            "route": canonical_route,
            "observation": _jsonable(dict(observation_profile)),
        },
        "denominator": {
            "input_units": input_units,
            "paired_units": paired_units,
            "excluded_units": _jsonable(dict(excluded_units)),
            "observation_rows": _jsonable(observation_profile.get("n_rows")),
            "included_units": _jsonable(observation_profile.get("included_n")),
            "moran_units": _jsonable(analysis_denominators.get("moran_units")),
            "emt_raw_units": _jsonable(analysis_denominators.get("emt_raw_units")),
            "emt_reconstruction_units": _jsonable(analysis_denominators.get("emt_reconstruction_units")),
            "emt_paired_units": _jsonable(analysis_denominators.get("emt_paired_units")),
        },
        "selection": _jsonable(observation_profile.get("selection", {})),
        "preprocessing": _jsonable(dict(preprocessing)),
        "gene": _gene_facts(
            gene_availability,
            gene_availability_issue,
            scope,
            analysis / "inputs/gene_availability.csv",
            output_root,
            gene_space=gene_space,
            union_gene_n=union_gene_n,
            partition_audit_values=partition_audit_values,
        ),
        "pairing_audit": _pairing_facts(
            pairing_audit,
            pairing_audit_issue,
            scope,
            analysis / "inputs/pairing_audit.csv",
            output_root,
        ),
    }


def _section_summaries(records: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    record_list = list(records)
    summaries: list[dict[str, Any]] = []
    for node in NODE_SPECS:
        node_id = str(node["node_id"])
        if node_id not in {"1.6", "2.6", "3.9", "4.5", "summary"}:
            continue
        questions = set(node["questions"])
        summaries.append({
            "node_id": node_id,
            "conclusion": {"zh": "", "en": ""},
            "record_ids": [
                str(record.get("identity", {}).get("record_id"))
                for record in record_list
                if str(record.get("identity", {}).get("question", "")) in questions
            ],
            "review": {"status": "pending"},
        })
    return summaries


def _builder_identity() -> dict[str, str]:
    """Identify the reader script so review invalidates after implementation edits."""

    path = Path(__file__).resolve()
    try:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError:
        digest = "unavailable"
    return {"path": str(path), "sha256": digest, "schema_version": SCHEMA_VERSION}


def _first_existing_column(columns: Iterable[Any], candidates: Iterable[str]) -> str | None:
    known = {str(column) for column in columns}
    return next((candidate for candidate in candidates if candidate in known), None)


def build_records(output_dir: str | os.PathLike[str], sample_id: str, route_kind: str) -> dict[str, Any]:
    """Build fixed-order records from an already-published route output."""

    if route_kind not in ROUTE_ALIASES:
        raise ValueError(f"unknown route_kind: {route_kind}")
    output_root, analysis = _resolve_roots(output_dir, sample_id)
    canonical_route = ROUTE_ALIASES[route_kind]
    audit, audit_issue = _read_json(analysis / "audit.json")
    parameters, parameters_issue = _read_json(analysis / "parameters.json")
    if not isinstance(audit, Mapping):
        audit = {}
    if not isinstance(parameters, Mapping):
        parameters = {}
    audit_sample = audit.get("sample_id")
    if audit_sample and str(audit_sample) != sample_id:
        raise ValueError(f"sample_id mismatch: requested {sample_id}, audit has {audit_sample}")
    saved_route = audit.get("route") or parameters.get("route")
    if saved_route and str(saved_route) != canonical_route:
        raise ValueError(f"route_kind mismatch: requested {route_kind}, output has {saved_route}")
    run_identity = _run_identity(output_root, analysis, sample_id, route_kind, audit, parameters)

    summary, summary_issue = _read_csv(analysis / "cross_parent_summary.csv")
    moran_summary, moran_summary_issue = _read_csv(analysis / "moran/moran_summary.csv")
    moran_distribution, moran_distribution_issue = _read_csv(analysis / "moran/moran_distribution_summary.csv")
    pathway_summary, pathway_summary_issue = _read_csv(analysis / "pathway_activity/summary.csv")
    pathway_availability, pathway_availability_issue = _read_csv(analysis / "pathway_activity/availability.csv")
    cutoff_audit, cutoff_audit_issue = _read_csv(analysis / "pathway_activity/cutoff_audit.csv")
    pathway_resource, pathway_resource_issue = _read_json(analysis / "pathway_activity/resource.json")
    if not isinstance(pathway_resource, Mapping):
        pathway_resource = {}
    moran_graph_audit, moran_graph_audit_issue = _read_csv(analysis / "moran/moran_graph_audit.csv")
    carriers, carriers_issue = _read_csv(analysis / "inputs/carriers.csv")
    source_files, source_files_issue = _read_csv(analysis / "inputs/source_files.csv")
    gene_availability, gene_availability_issue = _read_csv(analysis / "inputs/gene_availability.csv")
    pairing_audit, pairing_audit_issue = _read_csv(analysis / "inputs/pairing_audit.csv")
    anatomy, anatomy_issue = _read_csv(analysis / "anatomy/anatomy_context_summary.csv")
    anatomy_support, anatomy_support_issue = _read_csv(analysis / "anatomy/support_sensitivity.csv")
    local_diversity, local_diversity_issue = _read_csv(analysis / "cross_parent_local_diversity.csv")
    extent, extent_issue = _read_csv(analysis / "cross_parent_state_region_extent.csv")
    complexity, complexity_issue = _read_csv(analysis / "reconstruction_impact/partition/complexity_summary.csv")
    matched_k, matched_k_issue = _read_csv(analysis / "reconstruction_impact/partition/matched_k_summary.csv")
    spatial_parameters = parameters.get("spatial", {}) if isinstance(parameters.get("spatial", {}), Mapping) else {}
    threshold_bootstraps = spatial_parameters.get("threshold_bootstraps")
    if canonical_route == "sc_svc":
        carrier_condition_zh = "Xenium 重建表达由 cluster-mean 投影到空间单元；这不是独立的空间生物学验证。"
        carrier_condition_en = "The Xenium reconstructed expression is a cluster-mean projection onto spatial units; this is not independent spatial biological validation."
        emt_field_condition = "cluster_mean_projection"
        carrier_metadata = {"raw_native": True, "reconstruction_cluster_mean_projection": True, "reconstruction_expression_spatial": False}
    else:
        carrier_condition_zh = "HD 使用保存的重建表达空间载体，结果描述该载体上的变化。"
        carrier_condition_en = "HD uses the saved reconstructed-expression spatial carrier; results describe changes on that carrier."
        emt_field_condition = "reconstructed_expression_spatial"
        carrier_metadata = {"raw_native": True, "reconstruction_cluster_mean_projection": False, "reconstruction_expression_spatial": True}

    scopes = ("All", *SCOPES) if canonical_route == "sp_svc" else SCOPES
    records: list[dict[str, Any]] = []
    for scope in scopes:
        row, row_issue = _summary_row(summary, scope)
        mrow, mrow_issue = _summary_row(moran_summary, scope)
        complexity_row, complexity_row_issue = _row(complexity, scope)
        matched_row, matched_row_issue = _row(matched_k, scope)
        partition_audit, partition_audit_issue = _read_json(analysis / "reconstruction_impact/partition" / f"{scope}_audit.json")
        if not isinstance(partition_audit, Mapping):
            partition_audit = {}
        partition_audit_values = partition_audit.get("audit", {}) if isinstance(partition_audit.get("audit", {}), Mapping) else {}
        observation_profile = _observation_profile(analysis / scope / "inputs" / "observations.csv.gz", output_root, scope)
        observation_status = observation_profile.get("status")
        observation_issue = observation_status if observation_status not in {None, "valid"} else None
        base_condition = _condition(audit_issue, parameters_issue, summary_issue, row_issue, partition_audit_issue, observation_issue, scope=scope, route=canonical_route, observation_status=observation_status)
        figures_common: list[str] = []
        foundation_paired_units = _as_number(partition_audit_values.get("n_units"))
        foundation_input_units = _as_number(_first_non_none(partition_audit_values.get("input_units"), observation_profile.get("included_n")))
        foundation_excluded_units = {
            "raw_qc": _as_number(partition_audit_values.get("excluded_raw_qc_units")),
            "post_gene_filter": _as_number(partition_audit_values.get("excluded_post_gene_filter_units")),
            "observation_carrier": observation_profile.get("excluded_n"),
        }
        foundation_gene_space = "full" if audit.get("expression") == "complete gene space per side" else "unknown"
        foundation_gene_union_n = _as_number(row.get("union_gene_n"))
        pathway_row, pathway_row_issue = _row(pathway_summary, scope)
        graph_row, graph_row_issue = _row(moran_graph_audit, scope)
        partition_parameters = parameters.get("partition", {}) if isinstance(parameters.get("partition", {}), Mapping) else {}
        partition_fields = {
            key: partition_parameters.get(key)
            for key in ("mode", "within_level1_resolution", "random_state", "n_top_genes", "raw_qc_min_genes", "raw_qc_min_cells", "sample_n_units", "parent_sample_n_units", "feature_selection")
        }
        partition_fields.update({key: partition_audit_values.get(key) for key in ("n_units", "n_shared_genes", "n_features")})
        partition_source = {key: value for key, value in partition_fields.items() if value is not None}
        moran_fields = {
            key: graph_row.get(key)
            for key in ("graph_id", "graph_status", "coordinate_source", "coordinate_unit", "coordinate_to_microns", "coord_type_actual", "n_neighbors_actual", "n_units", "n_edges_row_normalized", "n_isolates", "n_components", "support_status", "weight_transformation", "actual_graph_parameters", "min_units")
        }
        moran_fields.update({key: mrow.get(key) for key in ("raw_nvalid", "reconstruction_nvalid", "matched_gene_n", "paired_delta_n")})
        cutoff_scope = cutoff_audit.loc[cutoff_audit["scope"].astype(str).eq(scope)] if not cutoff_audit.empty and "scope" in cutoff_audit.columns else pd.DataFrame()
        cutoff_fields = {}
        for side in ("raw", "reconstruction"):
            side_matches = cutoff_scope.loc[cutoff_scope["side"].astype(str).str.lower().eq(side)] if not cutoff_scope.empty and "side" in cutoff_scope.columns else pd.DataFrame()
            side_row = side_matches.iloc[0].to_dict() if not side_matches.empty else {}
            cutoff_fields[side] = {
                key: _jsonable(side_row.get(key))
                for key in ("n_observations", "n_genes", "detected_count_quantile", "provider_auc_threshold", "effective_rank_length", "rank_cutoff_zero_based", "cutoff_formula", "rank_length_formula", "status", "reason")
            }
        resource_provider = pathway_resource.get("provider", {}) if isinstance(pathway_resource.get("provider", {}), Mapping) else {}
        emt_fields = {
            "pathway": pathway_resource.get("pathway"),
            "scorer": pathway_resource.get("scorer"),
            "cutoff_policy": pathway_resource.get("cutoff_policy"),
            "cutoff_formula": pathway_resource.get("cutoff_formula"),
            "provider": {"name": resource_provider.get("name"), "scorer": resource_provider.get("scorer"), "version": resource_provider.get("version")},
            "cutoff_by_side": cutoff_fields,
        }
        analysis_denominators = {
            "moran_units": _as_number(graph_row.get("n_units")),
            "emt_raw_units": _as_number(cutoff_fields["raw"].get("n_observations")),
            "emt_reconstruction_units": _as_number(cutoff_fields["reconstruction"].get("n_observations")),
            "emt_paired_units": _value(pathway_row, "paired_n", number=True),
        }
        foundation_preprocessing = {
            "partition": _saved_preprocessing_fact(
                partition_fields,
                partition_audit_issue or ("missing" if not partition_fields else None),
                analysis / "reconstruction_impact/partition" / f"{scope}_audit.json",
                output_root,
                tuple(partition_fields),
                status_source=partition_source,
            ),
            "moran": _saved_preprocessing_fact(
                moran_fields,
                moran_graph_audit_issue or graph_row_issue or moran_summary_issue or mrow_issue,
                analysis / "moran/moran_graph_audit.csv",
                output_root,
                tuple(moran_fields),
                status_source=graph_row,
            ),
            "emt": _saved_preprocessing_fact(
                emt_fields if any(item.get("status") is not None for item in cutoff_fields.values()) else {},
                cutoff_audit_issue or pathway_resource_issue or pathway_summary_issue or pathway_row_issue,
                analysis / "pathway_activity/cutoff_audit.csv",
                output_root,
                tuple(emt_fields),
                status_source=cutoff_fields if any(item.get("status") is not None for item in cutoff_fields.values()) else {},
            ),
        }
        foundation_facts = _foundation_facts(
            scope=scope,
            canonical_route=canonical_route,
            audit=audit,
            observation_profile=observation_profile,
            partition_audit_values=partition_audit_values,
            carriers=carriers,
            carriers_issue=carriers_issue,
            source_files=source_files,
            source_files_issue=source_files_issue,
            gene_availability=gene_availability,
            gene_availability_issue=gene_availability_issue,
            pairing_audit=pairing_audit,
            pairing_audit_issue=pairing_audit_issue,
            output_root=output_root,
            analysis=analysis,
            input_units=foundation_input_units,
            paired_units=foundation_paired_units,
            excluded_units=foundation_excluded_units,
            gene_space=foundation_gene_space,
            union_gene_n=foundation_gene_union_n,
            preprocessing=foundation_preprocessing,
            analysis_denominators=analysis_denominators,
        )
        foundation_tables = [
            "analysis/audit.json",
            "analysis/parameters.json",
            f"analysis/reconstruction_impact/partition/{scope}_audit.json",
            f"analysis/{scope}/inputs/observations.csv.gz",
        ]
        for optional_path, optional_issue in (
            ("analysis/inputs/carriers.csv", carriers_issue),
            ("analysis/inputs/source_files.csv", source_files_issue),
            ("analysis/inputs/gene_availability.csv", gene_availability_issue),
            ("analysis/inputs/pairing_audit.csv", pairing_audit_issue),
        ):
            if optional_issue is None:
                foundation_tables.insert(2, optional_path)

        records.append(_make_record(
            sample_id=sample_id, route_kind=route_kind, run_identity=run_identity, scope=scope, question="foundation",
            results={"scope": scope, "paired_units": foundation_paired_units, "input_units": foundation_input_units, "excluded_units": foundation_excluded_units, "observation_axis": {"id_column": observation_profile.get("id_column"), "coordinate_columns": observation_profile.get("coordinate_columns"), "carrier": observation_profile.get("path"), "reason_values": observation_profile.get("reason_values", [])}, "gene_carrier": {"n_shared_genes": _as_number(_first_non_none(partition_audit_values.get("n_shared_genes"), row.get("union_gene_n"))), "input_shared_genes": _as_number(_first_non_none(partition_audit_values.get("input_shared_genes"), row.get("union_gene_n"))), "excluded_genes": _as_number(partition_audit_values.get("excluded_raw_qc_genes")), "n_features": _as_number(partition_audit_values.get("n_features")), "feature_selection": partition_audit_values.get("feature_selection")}, "gene_space": foundation_gene_space, "union_gene_n": foundation_gene_union_n, "expression": audit.get("expression"), "cohort": audit.get("cohort"), "route": canonical_route, "seed": _as_number(_first_non_none(audit.get("seed"), parameters.get("seed"))), **foundation_facts},
            judgment="foundation available" if partition_audit_values or observation_profile.get("status") == "valid" else "foundation incomplete", condition=base_condition,
            conclusion_zh=f"{scope} 的 partition QC 后配对 n={foundation_paired_units if foundation_paired_units is not None else 'NA'}；表达/Moran/EMT cohort 保留 input_units={foundation_input_units if foundation_input_units is not None else 'NA'}。ID、坐标和排除来自已保存 observation carrier，基因空间为 {audit.get('expression') or 'NA'}。",
            conclusion_en=f"For {scope}, partition QC leaves paired n={foundation_paired_units if foundation_paired_units is not None else 'NA'}; the expression/Moran/EMT cohort retains input_units={foundation_input_units if foundation_input_units is not None else 'NA'}. IDs, coordinates, and exclusions come from the saved observation carrier, with expression contract {audit.get('expression') or 'NA'}.",
            limitations_zh=["paired units 是分析范围，不等同于原始输入单元总数。"], limitations_en=["Paired units are the analyzed scope and are not the raw input total."],
            evidence=_evidence(output_root, foundation_tables, figures_common, base_condition),
        ))

        complexity_ari = _value(complexity_row, "ARI", number=True)
        complexity_status = "available" if complexity_row else "unavailable"
        complexity_condition = _condition(complexity_issue, complexity_row_issue, scope=scope, status=complexity_status, comparison="complexity_summary")
        if canonical_route == "sc_svc":
            complexity_raw_k = _value(complexity_row, "Raw K at reference resolution", number=True)
            complexity_fixed_k = _value(complexity_row, "Fixed final-cluster K", number=True)
            complexity_results = {
                "comparison_kind": "raw_reference_vs_fixed_final_clusters",
                "raw_k": complexity_raw_k,
                "raw_k_reference_resolution": complexity_raw_k,
                "fixed_final_cluster_k": complexity_fixed_k,
                "reconstruction_k_at_raw_resolution": None,
                "ari": complexity_ari,
                "status": complexity_status,
                "source": "complexity_summary",
            }
            complexity_conclusion_zh = f"复杂度诊断比较 reference-resolution Raw K={_fact(complexity_raw_k)} 与固定最终 cluster K={_fact(complexity_fixed_k)}，ARI={_fact(complexity_ari)}；这不是同分辨率重建 K 比较。"
            complexity_conclusion_en = f"The complexity diagnostic compares reference-resolution Raw K={_fact(complexity_raw_k)} with fixed final-cluster K={_fact(complexity_fixed_k)}, ARI={_fact(complexity_ari)}; it is not a same-resolution reconstruction-K comparison."
        else:
            complexity_raw_k = _value(complexity_row, "Raw K", number=True)
            complexity_recon_k = _value(complexity_row, "Recon K at Raw resolution", number=True)
            complexity_results = {
                "comparison_kind": "raw_vs_reconstruction_at_raw_resolution",
                "raw_k": complexity_raw_k,
                "reconstruction_k_at_raw_resolution": complexity_recon_k,
                "fixed_final_cluster_k": None,
                "ari": complexity_ari,
                "status": complexity_status,
                "source": "complexity_summary",
            }
            complexity_conclusion_zh = f"复杂度诊断为 Raw K={_fact(complexity_raw_k)}、重建 Raw-resolution K={_fact(complexity_recon_k)}、ARI={_fact(complexity_ari)}；它与 matched-K 结果分开。"
            complexity_conclusion_en = f"The complexity diagnostic is Raw K={_fact(complexity_raw_k)}, reconstruction-at-raw-resolution K={_fact(complexity_recon_k)}, ARI={_fact(complexity_ari)}; it is separate from matched-K."
        records.append(_make_record(
            sample_id=sample_id, route_kind=route_kind, run_identity=run_identity, scope=scope, question="complexity", results=complexity_results,
            judgment="complexity diagnostic available" if complexity_row else "complexity diagnostic unavailable", condition=complexity_condition,
            conclusion_zh=complexity_conclusion_zh,
            conclusion_en=complexity_conclusion_en,
            limitations_zh=["聚类复杂度是路线范围内的描述，不直接证明生物学改善。"], limitations_en=["Cluster complexity is route-scoped and does not directly prove biological improvement."],
            evidence=_evidence(output_root, ["analysis/reconstruction_impact/partition/complexity_summary.csv", "analysis/reconstruction_impact/partition/summary.csv"], _figure_candidates(output_root, None, f"{scope.lower()}_reconstructed_clusters.png"), complexity_condition),
        ))
        matched_status = str(matched_row.get("status")) if matched_row.get("status") is not None else None
        matched_eligible = _value(matched_row, "headline_eligible", boolean=True)
        matched = matched_eligible is True and matched_status == "ok"
        matched_condition = _condition(matched_k_issue, matched_row_issue, scope=scope, status=matched_status, comparison="matched_k_summary")
        records.append(_make_record(
            sample_id=sample_id, route_kind=route_kind, run_identity=run_identity, scope=scope, question="matched_k",
            results={"raw_k": _value(matched_row, "Raw K", number=True), "reconstruction_k": _value(matched_row, "Recon K", number=True), "matched": matched, "status": matched_status, "unit_change_fraction": _value(matched_row, "ST-unit change", number=True), "balanced_change": _value(matched_row, "balanced change", number=True), "ari": _value(matched_row, "ARI", number=True), "headline_eligible": matched_eligible},
            judgment="matched-K eligible" if matched else "matched-K unavailable",
            condition=matched_condition,
            conclusion_zh="Raw 与重建 K 相同，可在该范围内使用匹配诊断。" if matched else "Raw 与重建 K 不一致，匹配诊断不作为 headline。",
            conclusion_en="Raw and reconstructed K match, so the matching diagnostic is eligible for this scope." if matched else "Raw and reconstructed K do not match; the matching diagnostic is not a headline result.",
            limitations_zh=["匹配失败时，下游 assignment 结论受限。"] if not matched else ["匹配诊断不是 Level1 生物身份变化的直接测量。"],
            limitations_en=["Downstream assignment conclusions are limited when matched-K is unavailable."] if not matched else ["The matching diagnostic is not a direct measurement of Level1 biological identity change."],
            evidence=_evidence(output_root, ["analysis/reconstruction_impact/partition/matched_k_summary.csv", "analysis/cross_parent_summary.csv"], _figure_candidates(output_root, None, "matched_k_contingency.png", f"{scope.lower()}_matched_clusters.png"), matched_condition), headline_eligible=matched,
        ))

        for gene_set, question in (("all_valid", "moran_all_valid"), ("shared_valid", "moran_shared_valid")):
            moran_results, moran_status = _moran_result(moran_distribution, mrow if mrow else row, scope, gene_set)
            delta = moran_results["paired_delta"].get("median") if gene_set == "shared_valid" else moran_results.get("median_difference")
            direction = "higher" if delta is not None and delta > 0 else "lower" if delta is not None and delta < 0 else "unchanged or unavailable"
            raw_moran = moran_results.get("raw", {})
            recon_moran = moran_results.get("reconstruction", {})
            delta_label = "paired median ΔI" if gene_set == "shared_valid" else "marginal median difference"
            distribution_condition_zh = "两侧基因集合不同，此差值仅描述各自分布，不是逐基因变化。" if gene_set == "all_valid" else "差值来自共同有效同名基因。"
            distribution_condition_en = "The gene sets differ; this contrasts side-specific distributions, not paired gene changes." if gene_set == "all_valid" else "The delta uses shared valid genes."
            moran_condition = _condition(moran_distribution_issue, moran_summary_issue, mrow_issue, scope=scope, gene_set=gene_set)
            records.append(_make_record(
                sample_id=sample_id, route_kind=route_kind, run_identity=run_identity, scope=scope, question=question, results=moran_results,
                judgment=("side-specific distributions" if gene_set == "all_valid" else "paired directional shift") if moran_status == "valid" else "unavailable", condition=moran_condition,
                conclusion_zh=f"{gene_set}：Raw 中位数={_fact(raw_moran.get('median'))}、Q75={_fact(raw_moran.get('q75'))}；重建中位数={_fact(recon_moran.get('median'))}、Q75={_fact(recon_moran.get('q75'))}；{delta_label}={_fact(delta)}，方向为{('升高' if direction == 'higher' else '降低' if direction == 'lower' else '未确定')}。{distribution_condition_zh} {carrier_condition_zh}",
                conclusion_en=f"{gene_set}: Raw median={_fact(raw_moran.get('median'))}, Q75={_fact(raw_moran.get('q75'))}; reconstructed median={_fact(recon_moran.get('median'))}, Q75={_fact(recon_moran.get('q75'))}; {delta_label}={_fact(delta)}, direction {direction}. {distribution_condition_en} No amplitude binning or biological claim is assigned. {carrier_condition_en}",
                limitations_zh=["all_valid 与 shared_valid 分母分开；缺失状态保持 NA。", carrier_condition_zh], limitations_en=["All-valid and shared-valid denominators are separate; missing statuses remain NA.", carrier_condition_en],
                evidence=_evidence(output_root, ["analysis/moran/moran_summary.csv", "analysis/moran/moran_distribution_summary.csv", "analysis/moran/moran_graph_audit.csv"], _figure_candidates(output_root, None, *("moran_all_gene_distribution.png", "moran_q75_heatmap.png") if gene_set == "all_valid" else ("moran_shared_gene_delta_distribution.png", "moran_shared_gene_scatter.png")), moran_condition),
            ))

        pathway_row, pathway_row_issue = _row(pathway_summary, scope)
        emt_figure_slug = "overall" if scope == "All" else scope.lower()
        emt_figure_names = [f"aucell_{emt_figure_slug}_distribution_and_delta.png"]
        if scope != "All":
            emt_figure_names.append("aucell_overall_distribution_and_delta.png")
        availability = pathway_availability.loc[pathway_availability.get("scope", pd.Series(dtype=str)).astype(str).eq(scope)] if not pathway_availability.empty and "scope" in pathway_availability.columns else pd.DataFrame()
        availability_values: dict[str, dict[str, Any]] = {}
        for side in ("raw", "reconstruction"):
            matches = availability.loc[availability["side"].astype(str).str.lower().eq(side)] if not availability.empty and "side" in availability.columns else pd.DataFrame()
            av = matches.iloc[0].to_dict() if not matches.empty else {}
            availability_values[side] = {"resource_gene_count": _value(av, "resource_gene_count", number=True), "available_gene_count": _value(av, "available_gene_count", number=True), "coverage": _value(av, "coverage", number=True), "detected_signature_gene_count": _value(av, "detected_signature_gene_count", number=True), "status": av.get("status")}
        raw_coverage = availability_values["raw"].get("coverage")
        recon_coverage = availability_values["reconstruction"].get("coverage")
        coverage_delta = None if raw_coverage is None or recon_coverage is None else _as_number(recon_coverage - raw_coverage)
        coverage_direction = "increased" if coverage_delta is not None and coverage_delta > 0 else "decreased" if coverage_delta is not None and coverage_delta < 0 else "unchanged or unavailable"
        coverage_condition = _condition(pathway_availability_issue, pathway_row_issue, scope=scope, pathway=availability_values["raw"].get("status"))
        records.append(_make_record(
            sample_id=sample_id, route_kind=route_kind, run_identity=run_identity, scope=scope, question="emt_coverage", results={"pathway": _first_non_none(pathway_row.get("pathway"), "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"), "raw": availability_values["raw"], "reconstruction": availability_values["reconstruction"], "carrier_condition": emt_field_condition, "carrier_metadata": carrier_metadata},
            judgment="coverage differs; comparability limited" if availability_values["raw"]["coverage"] is not None and availability_values["reconstruction"]["coverage"] is not None and availability_values["raw"]["coverage"] != availability_values["reconstruction"]["coverage"] else "coverage recorded", condition=coverage_condition,
            conclusion_zh=f"EMT 基因覆盖率由 Raw {_fact(raw_coverage)} 变为重建 {_fact(recon_coverage)}，差值={_fact(coverage_delta)}（{coverage_direction}）；覆盖率与分数分开报告。{carrier_condition_zh}", conclusion_en=f"EMT coverage changed from Raw {_fact(raw_coverage)} to reconstructed {_fact(recon_coverage)} (delta={_fact(coverage_delta)}, {coverage_direction}); coverage and score are reported separately. {carrier_condition_en}",
            limitations_zh=["覆盖率差异反映可用基因集合，不等同于生物学变化。", carrier_condition_zh], limitations_en=["Coverage reflects the available gene set and is not itself a biological change.", carrier_condition_en],
            evidence=_evidence(output_root, ["analysis/pathway_activity/availability.csv", "analysis/pathway_activity/cutoff_audit.csv", "analysis/pathway_activity/resource.json"], _figure_candidates(output_root, None, *emt_figure_names), coverage_condition),
        ))
        score_raw = {"n_valid": _value(pathway_row, "raw_n_valid", number=True), "mean": _first_non_none(_value(pathway_row, "raw_mean", number=True), _value(row, "raw_mean_aucell", number=True)), "median": _first_non_none(_value(pathway_row, "raw_median", number=True), _value(row, "raw_median_aucell", number=True)), "q1": _value(pathway_row, "raw_q1", number=True), "q3": _value(pathway_row, "raw_q3", number=True)}
        score_recon = {"n_valid": _value(pathway_row, "reconstruction_n_valid", number=True), "mean": _first_non_none(_value(pathway_row, "reconstruction_mean", number=True), _value(row, "reconstruction_mean_aucell", number=True)), "median": _first_non_none(_value(pathway_row, "reconstruction_median", number=True), _value(row, "reconstruction_median_aucell", number=True)), "q1": _value(pathway_row, "reconstruction_q1", number=True), "q3": _value(pathway_row, "reconstruction_q3", number=True)}
        cutoff_scope = cutoff_audit.loc[cutoff_audit["scope"].astype(str).eq(scope)] if not cutoff_audit.empty and "scope" in cutoff_audit.columns else pd.DataFrame()
        cutoff_values: dict[str, dict[str, Any]] = {}
        for side in ("raw", "reconstruction"):
            cutoff_matches = cutoff_scope.loc[cutoff_scope["side"].astype(str).str.lower().eq(side)] if not cutoff_scope.empty and "side" in cutoff_scope.columns else pd.DataFrame()
            cutoff_row = cutoff_matches.iloc[0].to_dict() if not cutoff_matches.empty else {}
            cutoff_values[side] = {"n_observations": _value(cutoff_row, "n_observations", number=True), "n_genes": _value(cutoff_row, "n_genes", number=True), "detected_count_quantile": _value(cutoff_row, "detected_count_quantile", number=True), "provider_auc_threshold": _value(cutoff_row, "provider_auc_threshold", number=True), "effective_rank_length": _value(cutoff_row, "effective_rank_length", number=True), "rank_cutoff_zero_based": _value(cutoff_row, "rank_cutoff_zero_based", number=True), "cutoff_formula": cutoff_row.get("cutoff_formula"), "status": cutoff_row.get("status"), "reason": cutoff_row.get("reason")}
        carrier_scope = carriers
        if not carriers.empty and "role" in carriers.columns and scope != "All":
            scope_text = " ".join(scope.replace("_", " ").split()).casefold()
            roles = carriers["role"].astype(str).str.replace("_", " ", regex=False).str.replace("/", " ", regex=False).str.split().str.join(" ").str.casefold()
            expected_roles = {"raw full context", f"{scope_text} expression carrier", f"{scope_text} spatial carrier"}
            if canonical_route == "sp_svc":
                expected_roles.add("reconstruction full carrier")
            carrier_scope = carriers.loc[roles.isin(expected_roles)]
        carrier_values = [_jsonable(item) for item in carrier_scope.to_dict("records")] if not carrier_scope.empty else []
        score_results = {"pathway": _first_non_none(pathway_row.get("pathway"), "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"), "raw": score_raw, "reconstruction": score_recon, "paired": {"n": _first_non_none(_value(pathway_row, "paired_n", number=True), _value(row, "paired_n", number=True)), "mean_delta": _first_non_none(_value(pathway_row, "paired_delta_mean", number=True), _value(row, "paired_delta_mean_aucell", number=True)), "median_delta": _first_non_none(_value(pathway_row, "paired_delta_median", number=True), _value(row, "paired_delta_median_aucell", number=True))}, "cutoff_audit": cutoff_values, "expression_carriers": carrier_values, "carrier_condition": emt_field_condition, "carrier_metadata": carrier_metadata}
        score_delta = score_results["paired"].get("median_delta")
        score_direction = "increased" if score_delta is not None and score_delta > 0 else "decreased" if score_delta is not None and score_delta < 0 else "unchanged or unavailable"
        score_condition = _condition(pathway_summary_issue, pathway_row_issue, cutoff_audit_issue, carriers_issue, scope=scope, score_status=pathway_row.get("comparison_status"))
        records.append(_make_record(
            sample_id=sample_id, route_kind=route_kind, run_identity=run_identity, scope=scope, question="emt_score", results=score_results,
            judgment="descriptive score shift" if score_results["paired"].get("median_delta") is not None else "unavailable", condition=score_condition,
            conclusion_zh=f"EMT 分数 Raw 中位数={_fact(score_raw.get('median'))}，重建中位数={_fact(score_recon.get('median'))}，paired median Δ={_fact(score_delta)}（{score_direction}）；cutoff 与表达载体条件见结果，分数不单独证明生物学改善。{carrier_condition_zh}", conclusion_en=f"EMT score median changed from Raw {_fact(score_raw.get('median'))} to reconstructed {_fact(score_recon.get('median'))}; paired median delta={_fact(score_delta)} ({score_direction}). Cutoff and expression-carrier conditions are retained in the result; the score alone does not prove biological improvement. {carrier_condition_en}",
            limitations_zh=["排名 cutoff、覆盖率和投影层限制分数解释。", carrier_condition_zh], limitations_en=["Rank cutoffs, coverage, and the projection layer limit score interpretation.", carrier_condition_en],
            evidence=_evidence(output_root, ["analysis/pathway_activity/summary.csv", "analysis/pathway_activity/cutoff_audit.csv"], _figure_candidates(output_root, None, *emt_figure_names), score_condition),
        ))

        anatomy_rows = anatomy.to_dict("records") if not anatomy.empty else []
        anatomy_results = {"regions": [_jsonable(item) for item in anatomy_rows], "n_regions": len(anatomy_rows), "support_scales": [_jsonable(item) for item in anatomy_support.to_dict("records")]} 
        anatomy_condition = _condition(anatomy_issue, anatomy_support_issue, scope=scope)
        records.append(_make_record(
            sample_id=sample_id, route_kind=route_kind, run_identity=run_identity, scope=scope, question="anatomy", results=anatomy_results,
            judgment="anatomy context available" if anatomy_rows else "anatomy context unavailable", condition=anatomy_condition,
            conclusion_zh="解剖区域与窗口支持来自已保存的几何上下文。", conclusion_en="Anatomy labels and window support come from the saved geometric context.",
            limitations_zh=["区域上下文不等于状态区域的生物学解释。"], limitations_en=["Anatomy context does not by itself establish a biological state-region interpretation."],
            evidence=_evidence(output_root, ["analysis/anatomy/anatomy_context_summary.csv", "analysis/anatomy/support_sensitivity.csv"], _figure_candidates(output_root, None, "anatomy_regions.png", "anatomy_window_decision.png"), anatomy_condition),
        ))

        if scope == "All" and canonical_route == "sp_svc":
            change_path = analysis / "reconstruction_impact/partition/level1_summary.csv"
            change, change_issue = _read_csv(change_path)
            change_rows = change.loc[change["scope"].astype(str).eq("All")].to_dict("records") if not change.empty and "scope" in change.columns else []
            overall_change = next((item for item in change_rows if str(item.get("level1")) == "Overall"), None)
            change_evidence_tables = ["analysis/reconstruction_impact/partition/level1_summary.csv"]
            change_evidence_figures = _figure_candidates(output_root, None, "level1_change_fraction.png", "changed_units.png")
            change_conclusion_en = "The HD global Level1 change fraction is reported by saved Level1 stratum, with Overall as its separate denominator."
            change_conclusion_zh = "HD 全局 Level1 改变比例按已保存的 Level1 分层报告，并单独保留 Overall 分母。"
        else:
            change_path = analysis / "local_state" / scope / "change_by_anatomy.csv"
            change, change_issue = _read_csv(change_path)
            change_rows = change.to_dict("records") if not change.empty else []
            overall_change = next((item for item in change_rows if str(item.get("level1_region")) == "Overall"), None)
            change_evidence_tables = [f"analysis/local_state/{scope}/change_by_anatomy.csv", "analysis/cross_parent_summary.csv"]
            change_evidence_figures = _figure_candidates(output_root, None, "changed_units.png", "level1_change_fraction.png")
            change_conclusion_en = "Changed-unit fractions describe assignment differences and are not direct measurements of cell-identity change."
            change_conclusion_zh = "改变单元比例用于描述重建前后的分配差异，不直接等同于细胞身份变化。"
        change_condition = _condition(change_issue, summary_issue, scope=scope, source=_relative(change_path, output_root), matched_k_status=matched_status)
        records.append(_make_record(
            sample_id=sample_id, route_kind=route_kind, run_identity=run_identity, scope=scope, question="changed_units", results={"overall": _jsonable(overall_change), "by_region": [_jsonable(item) for item in change_rows], "stratification": "level1" if scope == "All" and canonical_route == "sp_svc" else "anatomy"},
            judgment="descriptive changed-unit fraction" if change_rows else "changed-unit result unavailable", condition=change_condition,
            conclusion_zh=change_conclusion_zh, conclusion_en=change_conclusion_en,
            limitations_zh=["匹配 K 不可用时，改变单元的下游解释受限。"], limitations_en=["Downstream interpretation is limited when matched-K is unavailable."],
            evidence=_evidence(output_root, change_evidence_tables, change_evidence_figures, change_condition),
        ))

        fields_path = analysis / "spatial_fields" / "values.csv.gz"
        fields, fields_issue = _read_csv(fields_path)
        field_rows = fields.loc[fields["scope"].astype(str).eq(scope)] if not fields.empty and "scope" in fields.columns else pd.DataFrame()
        field_results = {"artifact": _relative(fields_path, output_root), "available": fields_path.is_file(), "n_units": int(len(field_rows)), "features": sorted(field_rows["feature"].dropna().astype(str).unique().tolist()) if not field_rows.empty and "feature" in field_rows.columns else [], "expression_layers": sorted(field_rows["expression_layer"].dropna().astype(str).unique().tolist()) if not field_rows.empty and "expression_layer" in field_rows.columns else [], "coordinate_units": sorted(field_rows["coordinate_unit"].dropna().astype(str).unique().tolist()) if not field_rows.empty and "coordinate_unit" in field_rows.columns else [], "carrier_condition": emt_field_condition, "carrier_metadata": carrier_metadata}
        fields_condition = _condition(fields_issue, scope=scope, field_scope="EMT")
        records.append(_make_record(
            sample_id=sample_id, route_kind=route_kind, run_identity=run_identity, scope=scope, question="emt_spatial_fields", results=field_results,
            judgment="EMT spatial field available" if field_results["available"] else "EMT spatial field unavailable", condition=fields_condition,
            conclusion_zh=("Xenium EMT 空间场是 cluster-mean 表达投影到空间单元的字段，用于定位而不是独立空间验证。" if canonical_route == "sc_svc" else "HD EMT 空间场来自保存的重建表达空间载体，用于定位而不是追加分析。"), conclusion_en=("The Xenium EMT spatial field is a cluster-mean expression projection onto spatial units; it supports localization rather than independent spatial validation." if canonical_route == "sc_svc" else "The HD EMT spatial field uses the saved reconstructed-expression spatial carrier; it supports localization without adding an analysis."),
            limitations_zh=["空间字段是分数投影；受覆盖率和 cutoff 限制。"], limitations_en=["The spatial field is a score projection and inherits coverage and cutoff limits."],
            evidence=_evidence(output_root, ["analysis/spatial_fields/values.csv.gz", "analysis/pathway_activity/summary.csv"], _figure_candidates(output_root, None, f"{scope.lower()}_hallmark_epithelial_mesenchymal_transition_spatial.png"), fields_condition),
        ))

        window_path = analysis / "local_state" / scope / "window_decision.json"
        window_decision, window_decision_issue = _read_json(window_path)
        if not isinstance(window_decision, Mapping):
            window_decision = {}
        selection = window_decision.get("selection", {}) if isinstance(window_decision.get("selection", {}), Mapping) else {}
        scale = window_decision.get("scale", {}) if isinstance(window_decision.get("scale", {}), Mapping) else {}
        support_path = analysis / "local_state" / scope / "support_sensitivity.csv"
        support, support_issue = _read_csv(support_path)
        window_results = {"main_window_side_um": _as_number(_first_non_none(selection.get("window_side_length"), scale.get("main_window_side_um"), row.get("scale_um"))), "min_parent_units": _as_number(_first_non_none(selection.get("min_parent_units"), scale.get("min_parent_units"), row.get("min_parent_units"))), "rarefaction_draws": _as_number(_first_non_none(scale.get("rarefaction_draws"), row.get("rarefaction_draws"))), "total_windows": _value(row, "total_windows", number=True), "n_valid_windows": _value(row, "n_valid_windows", number=True), "valid_units": _value(row, "valid_units", number=True), "support_sensitivity": support.to_dict("records")}
        window_condition = _condition(window_decision_issue, support_issue, scope=scope)
        records.append(_make_record(
            sample_id=sample_id, route_kind=route_kind, run_identity=run_identity, scope=scope, question="window_support", results=window_results,
            judgment="window support available" if window_results["n_valid_windows"] is not None else "window support unavailable", condition=window_condition,
            conclusion_zh="窗口主尺度、最小父群支持和有效窗口分母来自已保存的窗口选择。", conclusion_en="Main window scale, minimum parent support, and valid-window denominators come from the saved window selection.",
            limitations_zh=["局部指标只在有效窗口和指定支持条件下解释。"], limitations_en=["Local metrics are interpreted only on valid windows under the specified support condition."],
            evidence=_evidence(output_root, [f"analysis/local_state/{scope}/window_decision.json", f"analysis/local_state/{scope}/support_sensitivity.csv", "analysis/cross_parent_summary.csv"], _figure_candidates(output_root, None, f"{scope.lower()}_window_decision.png"), window_condition),
        ))

        if not local_diversity.empty and "scope" in local_diversity.columns and "metric" in local_diversity.columns:
            local_scope = local_diversity.loc[local_diversity["scope"].astype(str).eq(scope)]
        else:
            local_scope = pd.DataFrame()
        local_metrics: dict[str, dict[str, Any]] = {}
        for metric, label in (("k_obs", "Kobs"), ("neff", "Neff"), ("evenness", "evenness")):
            matches = local_scope.loc[local_scope["metric"].astype(str).eq(metric)] if not local_scope.empty else pd.DataFrame()
            values: dict[str, Any] = {}
            for baseline in ("raw_leiden", "raw_level2"):
                baseline_matches = matches.loc[matches["baseline"].astype(str).eq(baseline)] if not matches.empty and "baseline" in matches.columns else pd.DataFrame()
                local_row = baseline_matches.iloc[0].to_dict() if not baseline_matches.empty else {}
                values[baseline] = {
                    "baseline_n_observations": _value(local_row, "baseline_n_observations", number=True),
                    "baseline_median": _value(local_row, "baseline_median", number=True),
                    "baseline_q1": _value(local_row, "baseline_q1", number=True),
                    "baseline_q3": _value(local_row, "baseline_q3", number=True),
                    "reconstruction_n_observations": _value(local_row, "reconstruction_n_observations", number=True),
                    "reconstruction_median": _value(local_row, "reconstruction_median", number=True),
                    "delta_n_observations": _value(local_row, "delta_n_observations", number=True),
                    "delta_median": _value(local_row, "delta_median", number=True),
                    "delta_q1": _value(local_row, "delta_q1", number=True),
                    "delta_q3": _value(local_row, "delta_q3", number=True),
                }
            local_metrics[label] = values
        anatomy_diversity_path = analysis / "local_state" / scope / "diversity_by_anatomy.csv"
        anatomy_diversity, anatomy_diversity_issue = _read_csv(anatomy_diversity_path)
        local_dependency = [] if matched else ["Matched-K is unavailable for this scope; downstream local assignment interpretation is limited."]
        for baseline, question in (("raw_leiden", "local_vs_raw_leiden"), ("raw_level2", "local_vs_raw_level2")):
            metric_values = {label: values.get(baseline, {}) for label, values in local_metrics.items()}
            anatomy_summary = _anatomy_summary(anatomy_diversity, anatomy_diversity_issue, scope, baseline)
            neff_delta = metric_values["Neff"].get("delta_median")
            if neff_delta is None:
                judgment = "direction unknown"
                direction_en = "the saved Neff delta is unavailable"
                direction_zh = "Neff 方向未知"
            elif neff_delta == 0:
                judgment = "no directional support"
                direction_en = "the saved Neff delta is zero and gives no directional support"
                direction_zh = "保存的 Neff delta 为 0，不支持方向性判断"
            elif baseline == "raw_leiden" and neff_delta < 0:
                judgment = "supports concentrated state direction"
                direction_en = "the saved Neff delta supports a concentrated-state direction"
                direction_zh = "保存的 Neff delta 支持集中状态方向"
            elif baseline == "raw_level2" and neff_delta > 0:
                judgment = "supports finer state direction"
                direction_en = "the saved Neff delta supports a finer-state direction"
                direction_zh = "保存的 Neff delta 支持更细状态方向"
            else:
                judgment = "does not support requested direction"
                direction_en = "the saved Neff delta does not support the requested direction"
                direction_zh = "Neff delta 未支持局部状态更集中方向" if baseline == "raw_leiden" else "Neff delta 未支持更细内部状态方向"
            scale_value = _value(row, "scale_um", number=True)
            window_count = _value(row, "n_valid_windows", number=True)
            baseline_count = metric_values["Neff"].get("baseline_n_observations")
            support_valid = all(value is not None and value > 0 for value in (scale_value, window_count, baseline_count))
            local_condition = _condition(local_diversity_issue, None if support_valid else "missing", scope=scope, baseline=baseline,
                                         baseline_status="available" if baseline_count is not None and baseline_count > 0 else "unavailable",
                                         scale_um=scale_value, n_valid_windows=window_count,
                                         interpretation_condition="supported" if support_valid else "insufficient_support", matched_k_status=matched_status)
            if not support_valid:
                judgment = "direction unknown"
                direction_zh = "baseline 或窗口支持不足，方向暂不可判读"
                direction_en = "baseline or window support is insufficient for directional interpretation"
            limitation_en = ["Evenness is descriptive; no higher-is-better interpretation is assigned.", "These are conditioned directions, not improvement proof.", *local_dependency]
            limitation_zh = ["均匀度仅作描述，不作越高越好的解释。", "这些是有条件的方向，不是改善证明。", *( ["matched-K 不可用，下游局部 assignment 解释受限。"] if not matched else [])]
            records.append(_make_record(
                sample_id=sample_id, route_kind=route_kind, run_identity=run_identity, scope=scope, question=question,
                results={"baseline": baseline, "scale_um": _value(row, "scale_um", number=True), "n_valid_windows": _value(row, "n_valid_windows", number=True), "metrics": metric_values, "anatomy_summary": anatomy_summary},
                judgment=judgment, condition=local_condition,
                conclusion_zh=f"相对 {('Raw Leiden' if baseline == 'raw_leiden' else 'Raw Level2')}，{direction_zh}；Neff delta={_fact(neff_delta)}。",
                conclusion_en=f"Against {('Raw Leiden' if baseline == 'raw_leiden' else 'Raw Level2')}, {direction_en} (Neff delta={_fact(neff_delta)}).",
                limitations_zh=limitation_zh, limitations_en=limitation_en,
                evidence=_evidence(output_root, ["analysis/cross_parent_local_diversity.csv", f"analysis/local_state/{scope}/scale_sensitivity.csv", f"analysis/local_state/{scope}/diversity_by_anatomy.csv"], _figure_candidates(output_root, None, f"{scope.lower()}_kobs_matrix.png", f"{scope.lower()}_neff_matrix.png", f"{scope.lower()}_evenness_matrix.png"), local_condition),
            ))

        threshold_state_path = analysis / "local_state" / scope / "state_threshold.json"
        window_metrics_path = analysis / "local_state" / scope / "window_metrics.csv"
        threshold_state = _threshold(threshold_state_path, window_metrics_path, threshold_bootstraps)
        threshold_condition = _condition(threshold_state.get("issue"), scope=scope, state_status=threshold_state.get("status"), window_metrics= _relative(window_metrics_path, output_root))
        relative_width = threshold_state.get("relative_ci_width")
        valid_fraction = threshold_state.get("valid_bootstrap_fraction")
        reliable = threshold_state.get("status") == "ok" and threshold_state.get("threshold") is not None and relative_width is not None and relative_width <= 0.25 and valid_fraction is not None and valid_fraction >= 0.8
        threshold_judgment = "reliable" if reliable else "not computable / NA" if threshold_state.get("analysis_status") == "not_computable" else "unstable / extent NA"
        threshold_reason = threshold_state.get("reason") or ""
        threshold_status_text = _fact(threshold_state.get("status"))
        threshold_reason_text = threshold_reason or "not computable"
        threshold_conclusion_zh = "先报告阈值稳定性，再解释区域范围；阈值可靠。" if reliable else f"先报告阈值稳定性；阈值状态为 {threshold_status_text}，原因：{threshold_reason_text}；区域范围保持 NA。"
        threshold_conclusion_en = "Threshold reliability is assessed before extent; the threshold is reliable." if reliable else f"Threshold reliability is assessed before extent; the threshold status is {threshold_status_text} ({threshold_reason_text}), so extent remains NA."
        records.append(_make_record(
            sample_id=sample_id, route_kind=route_kind, run_identity=run_identity, scope=scope, question="threshold_reliability", results={"state": threshold_state},
            judgment=threshold_judgment, condition=threshold_condition,
            conclusion_zh=threshold_conclusion_zh, conclusion_en=threshold_conclusion_en,
            limitations_zh=["CI 宽度和 bootstrap 有效数来自已保存阈值输出。"], limitations_en=["CI width and valid-bootstrap counts are read from the saved threshold output."],
            evidence=_evidence(output_root, [f"analysis/local_state/{scope}/state_threshold.json", f"analysis/local_state/{scope}/window_metrics.csv", f"analysis/local_state/{scope}/state_threshold_bootstrap.csv"], _figure_candidates(output_root, None, "threshold_reliability.png", f"{scope.lower()}_state_region.png"), threshold_condition),
        ))

        extent_scope_issue = None
        extent_scope = pd.DataFrame()
        if extent.empty or "scope" not in extent.columns:
            extent_scope_issue = "missing"
        else:
            extent_scope = extent.loc[extent["scope"].astype(str).eq(scope)]
            if extent_scope.empty:
                extent_scope_issue = "missing"
        overall_matches = extent_scope.loc[extent_scope["level1_region"].astype(str).eq("Overall")] if not extent_scope.empty and "level1_region" in extent_scope.columns else pd.DataFrame()
        overall_extent = overall_matches.iloc[0].to_dict() if not overall_matches.empty else {}
        overall_extent = {key: _jsonable(value) for key, value in overall_extent.items()}
        if "region_area_fraction" not in overall_extent:
            overall_extent["region_area_fraction"] = overall_extent.get("area_fraction")
        if "region_unit_fraction" not in overall_extent:
            overall_extent["region_unit_fraction"] = overall_extent.get("unit_fraction")
        extent_status = str(_first_non_none(overall_extent.get("threshold_status"), threshold_state.get("status")) or "")
        extent_usable = reliable and extent_status == "ok"
        region_available = _as_bool(overall_extent.get("region_available"))
        if not extent_usable or region_available is not True:
            for key in ("region_windows", "region_area_um2", "region_area_mm2", "area_fraction", "region_area_fraction", "region_units", "unit_fraction", "region_unit_fraction"):
                if key in overall_extent:
                    overall_extent[key] = None
        extent_condition = _condition(extent_issue, extent_scope_issue, scope=scope, threshold_status=extent_status, overall_rows=int(len(overall_matches)), region_available=region_available, matched_k_status=matched_status)
        by_region = extent_scope.to_dict("records") if not extent_scope.empty else []
        for item in by_region:
            for key, value in list(item.items()):
                item[key] = _jsonable(value)
        if not matched:
            downstream_extent_limitation = "Matched-K is unavailable for this scope; downstream region assignment interpretation is limited."
        else:
            downstream_extent_limitation = ""
        records.append(_make_record(
            sample_id=sample_id, route_kind=route_kind, run_identity=run_identity, scope=scope, question="region_extent", results={"overall": overall_extent, "by_region": by_region, "threshold_status": extent_status},
            judgment="extent available" if extent_usable and overall_extent and region_available is True else "extent NA", condition=extent_condition,
            conclusion_zh="区域范围只在稳定阈值和可用区域条件下报告。" if extent_usable and region_available is True else "阈值不满足可靠性条件或区域不可用，区域范围为 NA。", conclusion_en="Region extent is reported only when the threshold is reliable and the region is available." if extent_usable and region_available is True else "The threshold is not reliable or the region is unavailable; region extent is NA.",
            limitations_zh=["区域面积与单元比例不是生物学改善的直接证明。", *( ["matched-K 不可用，下游区域 assignment 解释受限。"] if downstream_extent_limitation else [])], limitations_en=["Area and unit fractions are not direct proof of biological improvement.", *( [downstream_extent_limitation] if downstream_extent_limitation else [])],
            evidence=_evidence(output_root, ["analysis/cross_parent_state_region_extent.csv", f"analysis/local_state/{scope}/state_region_extent_by_anatomy.csv"], _figure_candidates(output_root, None, f"{scope.lower()}_state_region.png"), extent_condition),
        ))

        scale_path = analysis / "local_state" / scope / "scale_sensitivity.csv"
        scale_frame, scale_issue = _read_csv(scale_path)
        scale_condition = _condition(scale_issue, scope=scope, interpretation="conditional")
        records.append(_make_record(
            sample_id=sample_id, route_kind=route_kind, run_identity=run_identity, scope=scope, question="scale_sensitivity", results={"rows": scale_frame.to_dict("records"), "main_scale_um": _value(row, "scale_um", number=True)},
            judgment="conditional sensitivity" if not scale_frame.empty else "sensitivity unavailable", condition=scale_condition,
            conclusion_zh="尺度敏感性用于检查条件稳定性；它不是改善证明。", conclusion_en="Scale sensitivity checks conditional stability; it is not proof of improvement.",
            limitations_zh=["解释受窗口支持、最小父群和随机种子条件约束。"], limitations_en=["Interpretation is conditioned on window support, minimum parent units, and the saved random seed."],
            evidence=_evidence(output_root, [f"analysis/local_state/{scope}/scale_sensitivity.csv", "analysis/anatomy/support_sensitivity.csv"], _figure_candidates(output_root, None, f"{scope.lower()}_window_decision.png"), scale_condition),
        ))

    global_questions = {"foundation", "complexity", "matched_k", "moran_all_valid", "moran_shared_valid", "changed_units"}
    if canonical_route == "sp_svc":
        records = [
            record for record in records
            if record.get("identity", {}).get("scope") != "All" or record.get("identity", {}).get("question") in global_questions
        ]
    _enrich_saved_conclusions(records)
    package: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "sample_id": sample_id,
        "route_kind": route_kind,
        "run_identity": run_identity,
        "records": records,
        "node_narratives": build_node_narratives(records),
        "section_summaries": _section_summaries(records),
    }
    digest = _digest(_without_review(package))
    for record in records:
        record["review"] = {"status": "pending", "record_digest": digest}
    for summary_item in package["section_summaries"]:
        summary_item["review"] = {"status": "pending", "record_digest": digest}
    for narrative in package["node_narratives"]:
        narrative["review"] = {"status": "pending", "record_digest": digest}
    package["record_digest"] = digest
    package["review"] = {"status": "pending", "record_digest": digest}
    return package


def _enrich_saved_conclusions(records):
    """Add short factual conclusions from already-saved record results.

    The function mutates the supplied package or record list and returns the
    same object.  It deliberately touches only the four saved-result
    questions that this enrichment owns.
    """

    def number(value):
        if value is None or isinstance(value, bool):
            return None
        try:
            result = float(value)
        except (TypeError, ValueError):
            return None
        return result if math.isfinite(result) else None

    def same_number(left, right):
        left_value = number(left)
        right_value = number(right)
        return left_value is not None and right_value is not None and abs(left_value - right_value) <= 1e-9

    def fmt_number(value):
        value = number(value)
        if value is None:
            return "NA"
        if value == 0:
            return "0"
        if value.is_integer():
            return f"{int(value):,}"
        return f"{value:.4g}"

    def fmt_fraction(value):
        value = number(value)
        return "NA" if value is None else f"{value * 100:.4g}%"

    def fmt_bool(value):
        if value is True:
            return "true"
        if value is False:
            return "false"
        return "NA"

    def set_conclusion(record, zh, en):
        conclusion = record.get("conclusion")
        if not isinstance(conclusion, dict):
            conclusion = {}
            record["conclusion"] = conclusion
        conclusion["zh"] = zh
        conclusion["en"] = en

    for record in records:
        if not isinstance(record, dict):
            continue
        identity = record.get("identity")
        results = record.get("results")
        if not isinstance(identity, dict) or not isinstance(results, dict):
            continue
        question = identity.get("question")
        if question not in {"window_support", "scale_sensitivity", "region_extent", "changed_units"}:
            continue
        if results.get("not_applicable") is True:
            continue

        if question == "window_support":
            main_scale = results.get("main_window_side_um")
            support_rows = results.get("support_sensitivity")
            selected = None
            if isinstance(support_rows, list):
                selected = next(
                    (row for row in support_rows if isinstance(row, dict) and same_number(row.get("window_side_length"), main_scale)),
                    None,
                )
            if selected is None:
                continue
            total_windows = selected.get("n_tissue_windows")
            valid_windows = selected.get("n_valid_windows")
            valid_fraction = selected.get("valid_window_fraction")
            retained_units = selected.get("retained_parent_units")
            if total_windows is None or valid_windows is None or retained_units is None:
                continue
            min_units = results.get("min_parent_units")
            draws = results.get("rarefaction_draws")
            set_conclusion(
                record,
                f"选定 {fmt_number(main_scale)} μm 窗口下，{fmt_number(valid_windows)}/{fmt_number(total_windows)} 个窗口有效（{fmt_fraction(valid_fraction)}），保留 {fmt_number(retained_units)} 个 parent units；最小支持为 {fmt_number(min_units)}，rarefaction draws={fmt_number(draws)}。",
                f"At the selected {fmt_number(main_scale)} μm window, {fmt_number(valid_windows)}/{fmt_number(total_windows)} windows were valid ({fmt_fraction(valid_fraction)}) and {fmt_number(retained_units)} parent units were retained; minimum support={fmt_number(min_units)}, rarefaction draws={fmt_number(draws)}.",
            )
            continue

        if question == "scale_sensitivity":
            rows = results.get("rows")
            if not isinstance(rows, list):
                continue
            usable = [row for row in rows if isinstance(row, dict) and number(row.get("window_side_length")) is not None]
            leiden = [number(row.get("median_delta_neff_vs_raw_leiden")) for row in usable]
            leiden = [value for value in leiden if value is not None]
            level2 = [number(row.get("median_delta_neff_vs_raw_level2")) for row in usable]
            level2 = [value for value in level2 if value is not None]
            scales = [number(row.get("window_side_length")) for row in usable]
            scales = [value for value in scales if value is not None]
            if not scales or not leiden or not level2:
                continue
            set_conclusion(
                record,
                f"在已保存的 {fmt_number(min(scales))}–{fmt_number(max(scales))} μm 候选窗口中，Neff 相对 Raw Leiden 的 median delta 范围为 {fmt_number(min(leiden))} 至 {fmt_number(max(leiden))}，相对 Raw Level2 为 {fmt_number(min(level2))} 至 {fmt_number(max(level2))}；选定尺度为 {fmt_number(results.get('main_scale_um'))} μm。",
                f"Across the saved {fmt_number(min(scales))}–{fmt_number(max(scales))} μm candidate windows, the median Neff delta ranged from {fmt_number(min(leiden))} to {fmt_number(max(leiden))} versus Raw Leiden and from {fmt_number(min(level2))} to {fmt_number(max(level2))} versus Raw Level2; the selected scale was {fmt_number(results.get('main_scale_um'))} μm.",
            )
            continue

        if question == "region_extent":
            overall = results.get("overall")
            if not isinstance(overall, dict):
                continue
            threshold_status = results.get("threshold_status") or overall.get("threshold_status") or "NA"
            region_available = overall.get("region_available")
            valid_windows = overall.get("valid_windows")
            if threshold_status == "ok" and region_available is True:
                region_windows = overall.get("region_windows")
                valid_units = overall.get("valid_units")
                region_units = overall.get("region_units")
                area_mm2 = overall.get("region_area_mm2")
                area_fraction = overall.get("region_area_fraction", overall.get("area_fraction"))
                unit_fraction = overall.get("region_unit_fraction", overall.get("unit_fraction"))
                if all(value is not None for value in (region_windows, valid_windows, region_units, valid_units, area_mm2, area_fraction, unit_fraction)):
                    set_conclusion(
                        record,
                        f"State Region 包含 {fmt_number(region_windows)}/{fmt_number(valid_windows)} 个有效窗口，面积为 {fmt_number(area_mm2)} mm²（{fmt_fraction(area_fraction)}），单元为 {fmt_number(region_units)}/{fmt_number(valid_units)}（{fmt_fraction(unit_fraction)}）；threshold_status={threshold_status}。",
                        f"The State Region contains {fmt_number(region_windows)}/{fmt_number(valid_windows)} valid windows, covers {fmt_number(area_mm2)} mm² ({fmt_fraction(area_fraction)}), and contains {fmt_number(region_units)}/{fmt_number(valid_units)} units ({fmt_fraction(unit_fraction)}); threshold_status={threshold_status}.",
                    )
            elif valid_windows is not None:
                set_conclusion(
                    record,
                    f"threshold_status={threshold_status}；已评估 {fmt_number(valid_windows)} 个有效窗口，但 Region extent 保持 NA。",
                    f"threshold_status={threshold_status}; {fmt_number(valid_windows)} valid windows were evaluated, but Region extent remains NA.",
                )
            continue

        overall = results.get("overall")
        if not isinstance(overall, dict):
            continue
        total_units = overall.get("paired_units", overall.get("total_units"))
        changed_units = overall.get("changed_units")
        change_fraction = overall.get("change_fraction")
        if total_units is None or changed_units is None or change_fraction is None:
            continue
        conditions = []
        stratification = results.get("stratification")
        matched_status = overall.get("matched_cluster_status")
        headline_eligible = overall.get("headline_eligible")
        if stratification:
            conditions.append(f"stratification={stratification}")
        if matched_status is not None:
            conditions.append(f"matched_cluster_status={matched_status}")
        if headline_eligible is not None:
            conditions.append(f"headline_eligible={fmt_bool(headline_eligible)}")
        suffix = f"; {', '.join(conditions)}" if conditions else ""
        set_conclusion(
            record,
            f"Overall 有 {fmt_number(changed_units)}/{fmt_number(total_units)} 个单元发生 assignment 差异（{fmt_fraction(change_fraction)}）{suffix}。该结果描述 assignment 差异，不等同于细胞身份变化。",
            f"Overall, {fmt_number(changed_units)}/{fmt_number(total_units)} units showed assignment differences ({fmt_fraction(change_fraction)}){suffix}. This describes assignment differences and is not a direct measurement of cell identity.",
        )

    return records


def _normalise_package(records: Mapping[str, Any] | list[Mapping[str, Any]]) -> dict[str, Any]:
    if isinstance(records, Mapping) and "records" in records:
        package = dict(records)
        package["records"] = list(package["records"])
        if "node_narratives" in package:
            narratives = package["node_narratives"]
            package["node_narratives"] = [dict(item) for item in narratives if isinstance(item, Mapping)] if isinstance(narratives, list) else []
        else:
            package["node_narratives"] = build_node_narratives(package["records"])
        if "section_summaries" in package:
            summaries = package["section_summaries"]
            package["section_summaries"] = [dict(item) for item in summaries if isinstance(item, Mapping)] if isinstance(summaries, list) else []
        return package
    if isinstance(records, list):
        sample_id = str(records[0].get("identity", {}).get("sample_id", "unknown")) if records else "unknown"
        route_kind = str(records[0].get("identity", {}).get("route_kind", "unknown")) if records else "unknown"
        package = {"schema_version": SCHEMA_VERSION, "sample_id": sample_id, "route_kind": route_kind, "run_identity": {}, "records": list(records), "node_narratives": build_node_narratives(records)}
        digest = _digest(_without_review(package))
        package["record_digest"] = digest
        package["review"] = {"status": "pending", "record_digest": digest}
        return package
    raise TypeError("records must be a package mapping or list")


def save_records(records: Mapping[str, Any] | list[Mapping[str, Any]], destination: str | os.PathLike[str], *, review_manifest: Mapping[str, Any] | None = None) -> Path:
    """Save stable records and a digest-bound review manifest.

    A record set remains ``pending`` unless an explicit manifest says
    ``status=reviewed`` and carries the exact record digest.
    """

    package = _normalise_package(records)
    destination_path = Path(destination)
    if destination_path.suffix.lower() == ".json":
        destination_path.parent.mkdir(parents=True, exist_ok=True)
        records_path = destination_path
        review_path = destination_path.with_name("review.json") if destination_path.name == "report_records.json" else destination_path.with_suffix(".review.json")
    else:
        destination_path.mkdir(parents=True, exist_ok=True)
        records_path = destination_path / "report_records.json"
        review_path = destination_path / "review.json"
    digest = _digest(_without_review(package))
    package["record_digest"] = digest
    manifest = dict(review_manifest or {})
    if manifest.get("status") == "reviewed" and manifest.get("record_digest") == digest:
        review = {**manifest, "status": "reviewed", "record_digest": digest}
    else:
        review = {"status": "pending", "record_digest": digest}
        if review_manifest:
            review["rejected_manifest"] = {"status": review_manifest.get("status"), "record_digest": review_manifest.get("record_digest")}
    package["review"] = review
    if "section_summaries" in package:
        normalised_summaries = []
        for item in package["section_summaries"]:
            if not isinstance(item, Mapping):
                continue
            summary = dict(item)
            conclusion = summary.get("conclusion") if isinstance(summary.get("conclusion"), Mapping) else {}
            summary["conclusion"] = {
                "zh": _jsonable(conclusion.get("zh", "")),
                "en": _jsonable(conclusion.get("en", "")),
            }
            review_value = summary.get("review") if isinstance(summary.get("review"), Mapping) else {}
            reviewed_text = bool(str(summary["conclusion"]["zh"] or "").strip() and str(summary["conclusion"]["en"] or "").strip())
            if review["status"] == "reviewed" and reviewed_text:
                summary_review = {**dict(review_value), **{key: value for key, value in review.items() if key != "rejected_manifest"}}
                summary["review"] = {**summary_review, "status": "reviewed", "record_digest": digest}
            else:
                summary["review"] = {"status": "pending", "record_digest": digest}
            normalised_summaries.append(summary)
        package["section_summaries"] = normalised_summaries
    if "node_narratives" in package:
        normalised_narratives = []
        for item in package["node_narratives"]:
            if not isinstance(item, Mapping):
                continue
            narrative = dict(item)
            lead = narrative.get("lead") if isinstance(narrative.get("lead"), Mapping) else {}
            interpretations = narrative.get("interpretations") if isinstance(narrative.get("interpretations"), Mapping) else {}
            reviewed_text = bool(
                str(lead.get("zh", "") or "").strip()
                and str(lead.get("en", "") or "").strip()
                and isinstance(interpretations.get("zh"), list)
                and isinstance(interpretations.get("en"), list)
                and interpretations.get("zh")
                and interpretations.get("en")
                and all(str(value or "").strip() for value in interpretations.get("zh", []))
                and all(str(value or "").strip() for value in interpretations.get("en", []))
            )
            if review["status"] == "reviewed" and reviewed_text:
                narrative_review = {key: value for key, value in review.items() if key != "rejected_manifest"}
                narrative["review"] = {**narrative_review, "status": "reviewed", "record_digest": digest}
            else:
                narrative["review"] = {"status": "pending", "record_digest": digest}
            normalised_narratives.append(narrative)
        package["node_narratives"] = normalised_narratives
    package["records"] = [dict(record, review={"status": review["status"], "record_digest": digest}) for record in package.get("records", [])]
    if isinstance(records, dict):
        records.clear()
        records.update(package)
    records_path.write_text(json.dumps(_jsonable(package), indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    csv_rows = []
    for record in package.get("records", []):
        identity = record.get("identity", {})
        conclusion = record.get("conclusion", {})
        csv_rows.append({
            "record_id": identity.get("record_id"),
            "sample_id": identity.get("sample_id"),
            "route_kind": identity.get("route_kind"),
            "scope": identity.get("scope"),
            "task_cell_type": identity.get("task_cell_type"),
            "question": identity.get("question"),
            "layer": identity.get("layer"),
            "baseline": identity.get("baseline"),
            "performance_judgment": record.get("performance_judgment"),
            "evidence_status": record.get("evidence_condition", {}).get("status"),
            "conclusion_en": conclusion.get("en"),
            "conclusion_zh": conclusion.get("zh"),
            "results_json": json.dumps(_jsonable(record.get("results", {})), ensure_ascii=False, sort_keys=True),
            "evidence_json": json.dumps(_jsonable(record.get("evidence", {})), ensure_ascii=False, sort_keys=True),
            "review_status": record.get("review", {}).get("status"),
            "record_digest": digest,
        })
    pd.DataFrame(csv_rows).to_csv(records_path.with_suffix(".csv"), index=False)
    review_path.write_text(json.dumps(_jsonable(review), indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    return records_path


def _load_package(records_or_path: Mapping[str, Any] | list[Mapping[str, Any]] | str | os.PathLike[str]) -> dict[str, Any]:
    if isinstance(records_or_path, (str, os.PathLike)):
        path = Path(records_or_path)
        if path.is_dir():
            path = path / "report_records.json"
        return _normalise_package(json.loads(path.read_text(encoding="utf-8")))
    return _normalise_package(records_or_path)


def render_report(records_or_path: Mapping[str, Any] | list[Mapping[str, Any]] | str | os.PathLike[str], destination: str | os.PathLike[str]) -> Path:
    """Render through the dedicated static HTML module."""

    from .report_html import render_report as delegated

    return delegated(records_or_path, destination)


def format_notebook_records(records_or_path: Mapping[str, Any] | list[Mapping[str, Any]] | str | os.PathLike[str], layer: int | str | None = None, *, node_id: str | None = None) -> str:
    """Render the saved English result matrix through the HTML module."""

    from .report_html import format_notebook_records as delegated

    return delegated(_load_package(records_or_path), layer=layer, node_id=node_id)
