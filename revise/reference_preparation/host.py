"""REVISE host adapter for portable reference screening."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from revise.application.config import (
    _compile_engine_config,
    compile_application_config,
    load_application_yaml,
)
from revise.application.inputs import load_data, preprocess_data
from revise.application.publication import application_metadata
from revise.framework import REVISEPipeline
from revise.reference_preparation.ga_contract import (
    GAResponse,
    GlobalAnchoringResult,
    validate_global_anchoring_result,
)
from revise.utils import hash_jsonable, input_identities


def _path_digest(path: Path, role: str) -> str:
    records = input_identities([SimpleNamespace(role=role, path=str(path))])
    return records[0]["sha256"]


def run_global_anchoring(
    reference_path: Path,
    reconstruction_config: Path,
) -> GAResponse:
    """Run the configured Application route through Global Anchoring only.

    Every invocation recompiles the request and reloads both inputs. The
    candidate replaces the configured reference before its former filter is
    parsed, so whole-file candidates are compared under one unchanged route,
    preprocessing contract, seed, and solver.
    """
    candidate = Path(reference_path).resolve()
    config_path = Path(reconstruction_config).resolve()
    source, document = load_application_yaml(config_path)
    config = compile_application_config(
        document,
        source=source,
        reference_override=candidate,
    )
    input_digests = {
        "config": _path_digest(config_path, "reconstruction_config"),
        "st": _path_digest(config.st_path, "st"),
        "reference": _path_digest(candidate, "reference"),
    }
    if input_digests["config"] != config.config_sha256:
        raise RuntimeError("Reconstruction config changed while it was being compiled")
    spatial_adata, reference_adata = load_data(config)
    spatial_adata, reference_adata = preprocess_data(
        spatial_adata,
        reference_adata,
        config,
    )
    runtime, io, algorithm = _compile_engine_config(config)
    execution = REVISEPipeline().run_global_anchoring(
        svc_type=config.svc_type,
        application_mode=config.mode,
        runtime_overrides=runtime,
        io_overrides=io,
        algorithm_overrides=algorithm,
        st_adata=spatial_adata,
        sc_ref_adata=reference_adata,
        application_config_metadata=application_metadata(config, paths={}),
    )
    final_digests = {
        "config": _path_digest(config_path, "reconstruction_config"),
        "st": _path_digest(config.st_path, "st"),
        "reference": _path_digest(candidate, "reference"),
    }
    changed = sorted(
        role for role, digest in input_digests.items()
        if final_digests[role] != digest
    )
    if changed:
        raise RuntimeError(
            "Global Anchoring input changed during execution: " + ", ".join(changed)
        )

    posterior = execution.annotated_spatial.obsm.get(config.broad_column)
    if not isinstance(posterior, pd.DataFrame):
        raise TypeError(
            "Global Anchoring must publish its broad posterior as a pandas "
            f"DataFrame in obsm[{config.broad_column!r}]"
        )
    result_ids = [str(value) for value in posterior.index.tolist()]
    labels = [str(value) for value in posterior.columns.tolist()]
    expected_ids = list(execution.expected_st_unit_ids)
    scoring_genes = list(execution.scoring_genes)
    solver = str(execution.merged_config["ot"]["ga"]["solver"])
    metadata = {
        "host_adapter": "revise.reference_preparation.host.run_global_anchoring",
        "backend": "revise.backend.kernels.GlobalAnchoringKernel",
        "route": execution.route_key,
        "normalization": "adapter_declared_probability",
        "reconstruction_config_sha256": input_digests["config"],
        "st_input_sha256": input_digests["st"],
        "reference_sha256": input_digests["reference"],
        "effective_seed": int(execution.merged_config["runtime"]["seed"]),
        "effective_solver": solver,
        "effective_parameters_sha256": execution.algorithm_config_hash,
        "st_axis_sha256": hash_jsonable(expected_ids),
        "cell_type_labels": labels,
        "reference_n_obs": int(reference_adata.n_obs),
        "reference_n_vars": int(reference_adata.n_vars),
        "scoring_n_vars": len(scoring_genes),
        "scoring_genes_sha256": hash_jsonable(scoring_genes),
    }
    result = GlobalAnchoringResult(
        distribution=posterior.to_numpy(dtype=np.float64, copy=True),
        st_unit_ids=result_ids,
        cell_type_labels=labels,
        metadata=metadata,
    )
    validated = validate_global_anchoring_result(result, expected_ids)
    return GAResponse(result=validated, expected_st_unit_ids=expected_ids)


__all__ = ["run_global_anchoring"]
