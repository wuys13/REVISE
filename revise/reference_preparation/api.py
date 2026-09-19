"""Public preparation entry: produce a reference, never run reconstruction."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .config import (
    PairedConfig,
    atomic_write_json,
    load_preparation_config,
    load_reference_config,
    write_reference_config,
)
from .evidence import REPORT_KIND, stable_file_identity
from .ga_contract import GARunner
from .paired import extract_pair
from .screening import screen_candidates


@dataclass(frozen=True)
class PreparationResult:
    reference_path: Path
    reference_config_path: Path
    report_path: Path


class NoUsableReferenceError(RuntimeError):
    """No candidate succeeded; report_path contains per-candidate outcomes."""

    def __init__(self, report_path: Path):
        self.report_path = report_path
        super().__init__(f"No usable reference; inspect {report_path}")


def prepare_reference(
    config_path: str | Path,
    *,
    ga_runner: GARunner | None = None,
) -> PreparationResult:
    """Prepare one reference in a new output directory."""
    config_path = Path(config_path).expanduser().resolve()
    config_identity = stable_file_identity(config_path)
    config = load_preparation_config(config_path)
    if stable_file_identity(config_path) != config_identity:
        raise RuntimeError(f"Preparation configuration changed while resolving: {config_path}")
    paired = isinstance(config, PairedConfig)
    if not paired and not callable(ga_runner):
        raise ValueError("screen mode requires a callable ga_runner; no GA backend is bundled")

    config.output_dir.mkdir(parents=True, exist_ok=False)
    report_path = config.output_dir / "report.json"
    reference_config_path = config.output_dir / "reference.yaml"
    report = {
        "schema_version": 1,
        "report_kind": REPORT_KIND,
        "mode": "paired" if paired else "screen",
        "preparation_config": str(config_path),
        "preparation_config_sha256": config_identity["sha256"],
    }
    try:
        if paired:
            reference_path = config.output_dir / "reference.h5ad"
            source_identity = stable_file_identity(config.source)
            details = extract_pair(
                config.source, config.pair_column, config.pair_key, reference_path,
            )
            if stable_file_identity(config.source) != source_identity:
                raise RuntimeError("Paired source changed during extraction")
            output_hash = stable_file_identity(reference_path)["sha256"]
            report.update({
                "status": "selected",
                "selected_reference_path": str(reference_path),
                "selected_reference_sha256": output_hash,
                "paired_output_sha256": output_hash,
                "source_sha256": source_identity["sha256"],
                "source_size_bytes": source_identity["size_bytes"],
                **details,
            })
        else:
            report["reconstruction_config"] = str(config.reconstruction_config)
            report.update(screen_candidates(
                config.candidates, config.reconstruction_config, ga_runner,
            ))
            selected_path = report["selected_reference_path"]
            reference_path = Path(selected_path) if selected_path is not None else None
        if stable_file_identity(config_path) != config_identity:
            raise RuntimeError("Preparation configuration changed during execution")
    except Exception as error:
        report.update({"status": "failed", "failure_reason": f"{type(error).__name__}: {error}"})
        atomic_write_json(report_path, report)
        raise

    atomic_write_json(report_path, report)
    if reference_path is None:
        raise NoUsableReferenceError(report_path)
    try:
        if stable_file_identity(config_path) != config_identity:
            raise RuntimeError("Preparation configuration changed before publication")
        if (
            stable_file_identity(reference_path)["sha256"]
            != report["selected_reference_sha256"]
        ):
            raise RuntimeError("Selected reference changed before publication")
        write_reference_config(reference_config_path, reference_path)
    except Exception as error:
        report.update({"status": "failed", "failure_reason": f"{type(error).__name__}: {error}"})
        try:
            atomic_write_json(report_path, report, replace=True)
        except Exception as reporting_error:
            raise error from reporting_error
        raise
    return PreparationResult(reference_path, reference_config_path, report_path)


__all__ = [
    "prepare_reference", "load_reference_config", "PreparationResult", "NoUsableReferenceError",
]
