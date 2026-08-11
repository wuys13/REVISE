"""Application output naming, metadata, and H5AD publication."""

from __future__ import annotations

import os
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any, Mapping

from revise.utils.provenance import hash_jsonable

from .config import ApplicationConfig


def output_paths(config: ApplicationConfig) -> dict[str, Path]:
    if config.svc_type == "sc-SVC":
        prefix = f"{config.output_name}_" if config.output_name else ""
        return {
            "spatial": config.output_dir / f"{prefix}spatial.h5ad",
            "expression": config.output_dir / f"{prefix}expr.h5ad",
        }
    filename = f"{config.output_name}.h5ad" if config.output_name else "svc.h5ad"
    return {"svc": config.output_dir / filename}


def application_metadata(
    config: ApplicationConfig,
    *,
    cli_overrides: Mapping[str, Any],
    paths: Mapping[str, Path],
    dry_run: bool,
) -> dict[str, Any]:
    if config.svc_type == "sc-SVC":
        local_refinement = {
            "subtype_column": config.subtype_column,
            "select_cell_type": config.select_cell_type,
            "alpha": config.local_refinement_alpha,
            "resolutions": list(config.local_refinement_resolutions or ()),
        }
    elif config.svc_type == "sc-SVC-sr":
        local_refinement = {
            "strength": config.local_refinement_strength,
            "graph": (
                {
                    "method": config.local_refinement_graph_method,
                    "alpha": config.local_refinement_graph_alpha,
                    "n_neighbors": config.local_refinement_graph_n_neighbors,
                    "exp_neighbors": config.local_refinement_graph_exp_neighbors,
                    "spatial_neighbors": config.local_refinement_graph_spatial_neighbors,
                }
                if config.local_refinement_graph_method is not None
                else None
            ),
            "match_spot_sum": config.local_refinement_match_spot_sum,
        }
    else:
        local_refinement = {"strength": config.local_refinement_strength}
    effective_request = {
        "svc_type": config.svc_type,
        "algorithm": {"ot_method": config.ot_method},
        "inputs": {
            "st_format": config.st_format,
            "reference_filter": {
                "column": config.reference_filter_column,
                "value": config.reference_filter_value,
            },
        },
        "preprocessing": {
            "spatial": {
                "min_transcript_counts": config.spatial_min_transcript_counts,
                "min_counts": config.spatial_min_counts,
                "min_cell_counts": config.spatial_min_cell_counts,
            },
            "reference": {
                "min_transcript_counts": config.reference_min_transcript_counts,
                "min_genes": config.reference_min_genes,
                "min_cell_counts": config.reference_min_cell_counts,
            },
        },
        "global_anchoring": {"broad_column": config.broad_column},
        "local_refinement": local_refinement,
        "output": {
            "dir": str(config.output_dir),
            "name": config.output_name,
        },
        "execution": {"seed": config.seed},
    }
    return {
        "source_path": config.source_path,
        "source_sha256": config.config_sha256,
        "cli_overrides": dict(cli_overrides),
        "declared_root": config.declared_root,
        "resolved_root": str(config.resolved_root),
        "cwd": str(config.cwd),
        "resolved_inputs": config.resolved_inputs,
        "output_name": config.output_name,
        "output_paths": {key: str(path) for key, path in paths.items()},
        "effective_action": "preflight" if dry_run else "run",
        "effective_request": effective_request,
        "effective_request_hash": hash_jsonable(effective_request),
    }


def _published_artifacts(config: ApplicationConfig, svc) -> list[tuple[str, Any]]:
    outputs = dict(svc.artifacts.get("outputs", {}))
    if config.svc_type == "sp-SVC":
        required = (("svc", "sp_svc"),)
    elif config.svc_type == "sc-SVC":
        required = (
            ("spatial", "sc_svc_spatial"),
            ("expression", "sc_svc_expr"),
        )
    else:
        required = (("svc", svc.provenance.get("primary_output_key")),)

    missing = [str(key) for _, key in required if key not in outputs]
    if missing:
        raise RuntimeError(
            f"{config.svc_type} did not produce required output(s): {', '.join(missing)}"
        )
    return [(role, outputs[key]) for role, key in required]


def publish_outputs(config: ApplicationConfig, paths: Mapping[str, Path], ctx):
    """Publish the exact artifact objects and return those same references."""
    artifacts = _published_artifacts(config, ctx.svc)
    config.output_dir.mkdir(parents=True, exist_ok=True)
    metadata = {
        "svc_type": config.svc_type,
        "output_name": config.output_name,
        "profile": ctx.profile,
        "run_manifest": str(Path(ctx.run_dir) / "provenance.json"),
        "selected_cell_type": config.select_cell_type,
        "ot": ctx.merged_config.get("ot"),
    }
    temporary: list[tuple[Path, Path]] = []
    try:
        for role, adata in artifacts:
            adata.uns["revise_reconstruction"] = dict(metadata, output_role=role)
            target = paths[role]
            with NamedTemporaryFile(
                dir=config.output_dir,
                prefix=f".{target.name}.",
                suffix=".tmp.h5ad",
                delete=False,
            ) as handle:
                temporary_path = Path(handle.name)
            temporary.append((temporary_path, target))
            adata.write_h5ad(temporary_path)
        for temporary_path, target in temporary:
            os.replace(temporary_path, target)
    finally:
        for temporary_path, _ in temporary:
            temporary_path.unlink(missing_ok=True)

    values = tuple(adata for _, adata in artifacts)
    return values if config.svc_type == "sc-SVC" else values[0]


__all__ = ["application_metadata", "output_paths", "publish_outputs"]
