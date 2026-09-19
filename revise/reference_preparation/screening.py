"""Rank whole-file reference candidates using host-provided Global Anchoring."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

from .config import Candidate
from .evidence import stable_file_identity
from .ga_contract import GAResponse, GARunner, validate_global_anchoring_result
from .scoring import SCORE_METHODS, row_metrics, summarize_scores


_RUN_IDENTITY_FIELDS = (
    "host_adapter", "backend", "route", "reconstruction_config_sha256",
    "st_input_sha256", "st_axis_sha256", "effective_parameters_sha256",
    "effective_seed", "effective_solver",
)


def _failure_row(candidate: Candidate, error: Exception) -> dict[str, Any]:
    return {
        "candidate_id": candidate.candidate_id,
        "path": str(candidate.path),
        "status": "failed",
        "scores": {method: None for method in SCORE_METHODS},
        "ranks": {method: None for method in SCORE_METHODS},
        "failure_reason": f"{type(error).__name__}: {error}",
    }


def _comparison(value: object, baseline: object) -> str:
    if value is None or baseline is None:
        return "unavailable"
    return "equal" if value == baseline else "different"


def screen_candidates(
    candidates: Sequence[Candidate],
    reconstruction_config: Path,
    ga_runner: GARunner,
) -> dict[str, Any]:
    """Try every candidate and return compact, deterministically ranked evidence."""
    rows: list[dict[str, Any]] = []
    canonical_st_axis: list[str] | None = None
    baseline_labels: list[str] | None = None
    baseline_gene_digest: str | None = None
    baseline_gene_count: int | None = None
    baseline_run_identity: dict[str, Any] | None = None
    candidate_states: dict[str, tuple[int, int, int, int]] = {}

    for candidate in candidates:
        try:
            candidate_path = Path(candidate.path)
            if not candidate_path.is_file() or candidate_path.suffix.lower() != ".h5ad":
                raise FileNotFoundError(
                    f"candidate is not an existing H5AD file: {candidate_path}"
                )
            before = candidate_path.stat()
            candidate_states[candidate.candidate_id] = (
                before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns,
            )
            response = ga_runner(candidate_path, reconstruction_config)
            after = candidate_path.stat()
            after_state = (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
            if after_state != candidate_states[candidate.candidate_id]:
                raise RuntimeError("candidate file changed during Global Anchoring")
            if not isinstance(response, GAResponse):
                raise TypeError("GA runner must return a GAResponse")
            result = validate_global_anchoring_result(
                response.result, response.expected_st_unit_ids,
            )
            scores = summarize_scores(row_metrics(result.distribution))
            run_identity = {
                field: result.metadata.get(field) for field in _RUN_IDENTITY_FIELDS
            }
            if canonical_st_axis is None:
                canonical_st_axis = list(result.st_unit_ids)
                baseline_labels = list(result.cell_type_labels)
                baseline_gene_digest = result.metadata.get("scoring_genes_sha256")
                baseline_gene_count = result.metadata.get("scoring_n_vars")
                baseline_run_identity = run_identity
                comparability = {
                    "st_unit_ids": "baseline",
                    "cell_type_labels": "baseline",
                    "scoring_genes": "baseline" if baseline_gene_digest is not None else "unavailable",
                    "execution_context": "baseline",
                }
            else:
                if result.st_unit_ids != canonical_st_axis:
                    raise ValueError(
                        "Global Anchoring ST-unit IDs differ from the first valid candidate"
                    )
                if run_identity != baseline_run_identity:
                    changed = sorted(
                        field for field in _RUN_IDENTITY_FIELDS
                        if run_identity[field] != baseline_run_identity[field]
                    )
                    raise ValueError(
                        "Global Anchoring execution context differs from the first valid "
                        f"candidate: {changed}"
                    )
                gene_value = (
                    result.metadata.get("scoring_genes_sha256"),
                    result.metadata.get("scoring_n_vars"),
                )
                gene_baseline = (baseline_gene_digest, baseline_gene_count)
                comparability = {
                    "st_unit_ids": "equal",
                    "cell_type_labels": _comparison(result.cell_type_labels, baseline_labels),
                    "scoring_genes": _comparison(gene_value, gene_baseline)
                    if None not in gene_value and None not in gene_baseline else "unavailable",
                    "execution_context": "equal",
                }
            rows.append({
                "candidate_id": candidate.candidate_id,
                "path": str(candidate.path),
                "status": "success",
                "scores": scores,
                "ranks": {method: None for method in SCORE_METHODS},
                "failure_reason": None,
                "matrix_shape": [int(value) for value in result.distribution.shape],
                "cell_type_labels": list(result.cell_type_labels),
                "comparability": comparability,
                "ga_metadata": dict(result.metadata),
            })
        except Exception as error:
            rows.append(_failure_row(candidate, error))
        finally:
            response = None
            result = None

    for method in SCORE_METHODS:
        ranked_rows = sorted(
            (row for row in rows if row["status"] == "success"),
            key=lambda row: (-float(row["scores"][method]), row["candidate_id"]),
        )
        for rank, row in enumerate(ranked_rows, start=1):
            row["ranks"][method] = rank

    selected = next((row for row in rows if row["ranks"]["max_median"] == 1), None)
    selected_hash = None
    if selected is not None:
        selected_path = Path(selected["path"])
        current = selected_path.stat()
        current_state = (current.st_dev, current.st_ino, current.st_size, current.st_mtime_ns)
        if current_state != candidate_states[selected["candidate_id"]]:
            raise RuntimeError("selected reference changed after Global Anchoring")
        selected_hash = stable_file_identity(selected_path)["sha256"]
        callback_hash = selected["ga_metadata"].get("reference_sha256")
        if callback_hash is not None and callback_hash != selected_hash:
            raise RuntimeError(
                "selected reference SHA-256 does not match Global Anchoring metadata"
            )

    return {
        "status": "selected" if selected is not None else "no_reference",
        "default_score_method": "max_median",
        "selected_candidate_id": selected["candidate_id"] if selected is not None else None,
        "selected_reference_path": selected["path"] if selected is not None else None,
        "selected_reference_sha256": selected_hash,
        "candidates": rows,
    }
