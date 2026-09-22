#!/usr/bin/env python3
"""Measure full P2 Xenium load/preprocess memory without running reconstruction."""

from __future__ import annotations

import argparse
import copy
from datetime import datetime, timezone
import json
import hashlib
import os
from pathlib import Path
import resource
import subprocess
import sys
import traceback

import numpy as np
from scipy import sparse
import yaml

from revise.application.config import compile_application_config, load_application_yaml
from revise.application.expression import bind_expression_sources
from revise.application.inputs import load_data, preprocess_data


def now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat()


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def hash_jsonable(value) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def scientific_config(document: dict) -> dict:
    normalized = copy.deepcopy(document)
    normalized["inputs"]["st"].pop("path")
    normalized["inputs"]["reference"].pop("path")
    normalized["output"].pop("dir")
    return normalized


def read_integer(path: Path) -> int | None:
    try:
        value = path.read_text().strip()
    except OSError:
        return None
    if value == "max":
        return None
    try:
        return int(value)
    except ValueError:
        return None


def read_key_values(path: Path) -> dict[str, int] | None:
    try:
        lines = path.read_text().splitlines()
        return {key: int(value) for key, value in (line.split() for line in lines)}
    except (OSError, ValueError):
        return None


def current_rss_bytes() -> int | None:
    try:
        for line in Path("/proc/self/status").read_text().splitlines():
            if line.startswith("VmRSS:"):
                return int(line.split()[1]) * 1024
    except (OSError, ValueError):
        pass
    return None


def peak_rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def resources() -> dict:
    memory_current = read_integer(Path("/sys/fs/cgroup/memory.current"))
    memory_stat = read_key_values(Path("/sys/fs/cgroup/memory.stat"))
    inactive_file = memory_stat.get("inactive_file") if memory_stat is not None else None
    return {
        "process_rss_bytes": current_rss_bytes(),
        "process_peak_rss_bytes": peak_rss_bytes(),
        "cgroup_memory_current_bytes": memory_current,
        "cgroup_memory_working_set_bytes": (
            max(0, memory_current - inactive_file)
            if memory_current is not None and inactive_file is not None
            else None
        ),
        "cgroup_memory_peak_bytes": read_integer(Path("/sys/fs/cgroup/memory.peak")),
        "cgroup_memory_limit_bytes": read_integer(Path("/sys/fs/cgroup/memory.max")),
        "cgroup_memory_stat": (
            {key: memory_stat.get(key) for key in ("anon", "file", "inactive_file")}
            if memory_stat is not None
            else None
        ),
        "cgroup_memory_events": read_key_values(Path("/sys/fs/cgroup/memory.events")),
    }


def matrix_record(matrix) -> dict:
    result = {
        "shape": [int(value) for value in matrix.shape],
        "dtype": str(matrix.dtype),
        "sparse": bool(sparse.issparse(matrix)),
    }
    if sparse.issparse(matrix):
        value = matrix.tocsr(copy=False)
        result.update(
            {
                "format": value.format,
                "nnz": int(value.nnz),
                "storage_bytes": int(value.data.nbytes + value.indices.nbytes + value.indptr.nbytes),
            }
        )
    else:
        result["storage_bytes"] = int(np.asarray(matrix).nbytes)
    return result


def adata_record(adata) -> dict:
    return {
        "shape": [int(adata.n_obs), int(adata.n_vars)],
        "x": matrix_record(adata.X),
        "obs_columns": [str(value) for value in adata.obs.columns],
        "var_names_unique": bool(adata.var_names.is_unique),
        "obs_names_unique": bool(adata.obs_names.is_unique),
    }


def atomic_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")
    os.replace(temporary, path)


def run_text(command: list[str], *, cwd: Path) -> str:
    return subprocess.run(
        command,
        cwd=cwd,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    ).stdout.strip()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-st-sha256", required=True)
    parser.add_argument("--expected-reference-sha256", required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    config_path = (args.config if args.config.is_absolute() else repo / args.config).resolve()
    output = (args.output if args.output.is_absolute() else repo / args.output).resolve()
    if output.exists():
        raise SystemExit(f"refusing to overwrite capacity evidence: {output}")
    git_dirty = run_text(["git", "status", "--short"], cwd=repo).splitlines()
    if git_dirty:
        raise SystemExit(f"refusing to probe from a dirty worktree: {git_dirty}")
    source, document = load_application_yaml(config_path)
    base_config_path = repo / "configs/acceptance/full/P2CRC_Xenium-random.yaml"
    base_document = yaml.safe_load(base_config_path.read_text())
    normalized_science = scientific_config(document)
    if normalized_science != scientific_config(base_document):
        raise SystemExit("probe config changes frozen scientific parameters")
    if document["paths"]["root_dir"] != ".":
        raise SystemExit("capacity probe requires paths.root_dir: .")
    st_path = (repo / document["inputs"]["st"]["path"]).resolve()
    reference_path = (repo / document["inputs"]["reference"]["path"]).resolve()
    record = {
        "schema_version": 1,
        "kind": "P2CRC Xenium full load and preprocessing capacity probe",
        "status": "running",
        "reconstruction_started": False,
        "config": {
            "path": str(config_path),
            "sha256": digest(config_path),
            "base_path": str(base_config_path),
            "scientific_config_sha256": hash_jsonable(normalized_science),
        },
        "code_commit": run_text(["git", "rev-parse", "HEAD"], cwd=repo),
        "branch": run_text(["git", "branch", "--show-current"], cwd=repo),
        "git_dirty": git_dirty,
        "started_at": now(),
        "stages": {"baseline": {"resources": resources()}},
    }
    atomic_json(output, record)
    try:
        input_identity = {
            "spatial": {"path": str(st_path), "size": st_path.stat().st_size, "sha256": digest(st_path)},
            "reference": {
                "path": str(reference_path),
                "size": reference_path.stat().st_size,
                "sha256": digest(reference_path),
            },
        }
        if input_identity["spatial"]["sha256"] != args.expected_st_sha256:
            raise AssertionError("spatial input hash does not match the capacity gate")
        if input_identity["reference"]["sha256"] != args.expected_reference_sha256:
            raise AssertionError("reference input hash does not match the capacity gate")
        record["input_identity"] = input_identity
        record["observed_at"] = now()
        atomic_json(output, record)
        config = compile_application_config(document, source=source, cwd=repo)
        if (
            config.svc_type != "sc-SVC"
            or config.mode != "cluster"
            or config.ist_mapping != "random"
            or config.seed != 42
            or config.ot_method != "tacco"
            or (config.reference_filter_column, config.reference_filter_value) != ("Patient", "P2CRC")
        ):
            raise AssertionError("probe config is not the frozen full P2 Xenium random request")

        spatial, reference = load_data(config)
        record["stages"]["loaded"] = {
            "resources": resources(),
            "spatial": adata_record(spatial),
            "reference": adata_record(reference),
        }
        record["observed_at"] = now()
        atomic_json(output, record)
        if "Patient" not in reference.obs:
            raise KeyError("loaded reference is missing the configured Patient filter column")
        original_reference_patient = reference.obs["Patient"].copy()

        config = bind_expression_sources(config, spatial, reference)
        spatial, reference = preprocess_data(spatial, reference, config)
        record["stages"]["preprocessed"] = {
            "resources": resources(),
            "spatial": adata_record(spatial),
            "reference": adata_record(reference),
        }
        record["observed_at"] = now()
        atomic_json(output, record)

        surviving_reference_patient = original_reference_patient.reindex(reference.obs_names)
        if surviving_reference_patient.isna().any():
            missing = surviving_reference_patient.index[surviving_reference_patient.isna()].tolist()
            raise AssertionError(
                "preprocessed reference observations are not all traceable to loaded Patient metadata: "
                f"{missing[:10]}"
            )
        unexpected_patient = surviving_reference_patient.astype(str) != "P2CRC"
        if unexpected_patient.any():
            raise AssertionError(
                "preprocessed reference contains observations outside Patient=P2CRC: "
                f"{surviving_reference_patient[unexpected_patient].value_counts(dropna=False).to_dict()}"
            )
        record["stages"]["preprocessed"].update({
            "reference_patient_counts": {
                str(key): int(value)
                for key, value in surviving_reference_patient.value_counts(dropna=False).items()
            },
            "spatial_broad_counts": {
                str(key): int(value)
                for key, value in spatial.obs[config.broad_column].value_counts(dropna=False).items()
            },
            "reference_broad_counts": {
                str(key): int(value)
                for key, value in reference.obs[config.broad_column].value_counts(dropna=False).items()
            },
        })
        record["status"] = "capacity_observed"
        record["observed_at"] = now()
        atomic_json(output, record)
    except BaseException as exc:
        record["status"] = "failed"
        record["error"] = f"{type(exc).__name__}: {exc}"
        record["traceback"] = traceback.format_exc()
        record["observed_at"] = now()
        record["resources_at_failure"] = resources()
        atomic_json(output, record)
        raise


if __name__ == "__main__":
    main()
