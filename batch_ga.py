#!/usr/bin/env python3
"""Run TACCO global anchoring for one or more references from a YAML file."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile

import pandas as pd
from anndata import AnnData, read_h5ad

from revise.application.batch_ga_config import BatchGAConfig, load_batch_config
from revise.application.preprocess import (
    prepare_sc_svc_pair,
    preprocess_reference,
    preprocess_spatial,
)
from revise.backend.kernels.ot import OTKernel
from revise.backend.ops.assignment import GlobalAssignment
from revise.utils import set_global_seed, write_json


@dataclass
class BatchGAResult:
    assignments: dict[str, GlobalAssignment]
    failures: dict[str, str]
    status_path: Path


def global_anchor(
    spatial_adata: AnnData,
    reference_adata: AnnData,
    broad_column: str,
    seed: int,
) -> GlobalAssignment:
    set_global_seed(seed=seed, deterministic=True)
    annotated = OTKernel.annotate(
        spatial_adata,
        reference_adata,
        method="tacco",
        annotation_key=broad_column,
        confidence_key="Confidence",
        multi_center=1,
        lamb=0.001,
    )
    return GlobalAssignment(
        labels=annotated.obs[broad_column].copy(),
        posterior=annotated.obsm[broad_column].copy(),
    )


def _reference_ids(reference: AnnData, config: BatchGAConfig) -> list[str]:
    if config.filter_column not in reference.obs:
        raise KeyError(f"Missing reference filter column: {config.filter_column}")
    values = reference.obs[config.filter_column]
    if config.filter_value == "all":
        return list(dict.fromkeys(values.tolist()))
    return [config.filter_value]


def _write_posterior(path: Path, posterior: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = None
    try:
        with NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            posterior.to_csv(handle)
        temporary_path.replace(path)
    except BaseException:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
        raise


def run_batch_ga(config_path: str | Path) -> BatchGAResult:
    config = load_batch_config(config_path)
    spatial = preprocess_spatial(read_h5ad(config.spatial_path), **config.spatial_preprocessing)
    combined_reference = read_h5ad(config.reference_path)
    reference_ids = _reference_ids(combined_reference, config)
    config.output_dir.mkdir(parents=True, exist_ok=True)
    status_path = config.output_dir / "ga_status.json"
    tasks: list[dict] = []
    assignments: dict[str, GlobalAssignment] = {}
    failures: dict[str, str] = {}
    write_json(status_path, {"schema_version": 1, "status": "running", "tasks": tasks, "summary": None})

    for reference_id in reference_ids:
        task = {
            "reference_id": reference_id,
            "status": "failed",
            "posterior_path": None,
            "shape": None,
            "error": None,
        }
        try:
            current_reference = combined_reference[
                combined_reference.obs[config.filter_column] == reference_id, :
            ].copy()
            current_reference = preprocess_reference(
                current_reference,
                **config.reference_preprocessing,
            )
            current_spatial, current_reference = prepare_sc_svc_pair(
                spatial,
                current_reference,
                broad_column=config.broad_column,
                subtype_column=None,
            )
            assignment = global_anchor(
                current_spatial,
                current_reference,
                config.broad_column,
                config.seed,
            )
            output_path = config.output_dir / f"{reference_id}.csv"
            _write_posterior(output_path, assignment.posterior)
            assignments[reference_id] = assignment
            task.update(
                status="succeeded",
                posterior_path=output_path.name,
                shape=list(assignment.posterior.shape),
            )
        except Exception as exc:
            failures[reference_id] = f"{type(exc).__name__}: {exc}"
            task["error"] = failures[reference_id]
        tasks.append(task)
        write_json(
            status_path,
            {
                "schema_version": 1,
                "status": "running",
                "tasks": tasks,
                "summary": None,
            },
        )

    summary = {
        "total": len(reference_ids),
        "succeeded": len(assignments),
        "failed": len(failures),
    }
    write_json(
        status_path,
        {"schema_version": 1, "status": "completed", "tasks": tasks, "summary": summary},
    )
    return BatchGAResult(assignments=assignments, failures=failures, status_path=status_path)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Run batch TACCO global anchoring")
    parser.add_argument("--config", required=True, help="batch GA YAML")
    args = parser.parse_args(argv)
    result = run_batch_ga(args.config)
    print(f"Succeeded: {len(result.assignments)}")
    print(f"Failed: {len(result.failures)}")
    print(f"Status: {result.status_path}")


__all__ = ["BatchGAResult", "global_anchor", "load_batch_config", "run_batch_ga"]


if __name__ == "__main__":
    main()
