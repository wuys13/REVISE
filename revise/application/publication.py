"""Application output naming, metadata, and H5AD publication."""

from __future__ import annotations

import os
from dataclasses import replace
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any, Mapping

from revise.utils.provenance import completed_artifact, hash_jsonable

from .config import ApplicationConfig


def output_paths(config: ApplicationConfig) -> dict[str, Path]:
    if config.mode == "cluster" and getattr(config, "ist_mapping", "paired") != "paired":
        filename = f"{config.output_name}.h5ad" if config.output_name else "SVC.h5ad"
        return {"svc": config.output_dir / filename}
    if config.mode == "cluster":
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
    paths: Mapping[str, Path],
) -> dict[str, Any]:
    if config.mode == "cluster":
        local_refinement = {
            "subtype_column": config.subtype_column,
            "select_cell_type": config.select_cell_type,
            "alpha": config.local_refinement_alpha,
            "resolutions": list(config.local_refinement_resolutions or ()),
        }
    elif config.mode == "sr":
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
        "mode": config.mode,
        "application_route": config.svc_type,
        "application_mode": config.mode,
        "selected_cell_type": config.select_cell_type,
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
            "root": str(config.output_root),
            "dir": str(config.output_dir),
            "name": config.output_name,
        },
        "execution": {"seed": config.seed},
    }
    if config.mode == "cluster" and getattr(config, "ist_mapping", "paired") != "paired":
        effective_request["output"]["ist_mapping"] = config.ist_mapping
    return {
        "source_path": config.source_path,
        "source_sha256": config.config_sha256,
        "declared_root": config.declared_root,
        "resolved_root": str(config.resolved_root),
        "cwd": str(config.cwd),
        "resolved_inputs": config.resolved_inputs,
        "output_root": str(config.output_root),
        "output_dir": str(config.output_dir),
        "output_name": config.output_name,
        "output_paths": {key: str(path) for key, path in paths.items()},
        "effective_request": effective_request,
        "effective_request_hash": hash_jsonable(effective_request),
    }


def _published_artifacts(config: ApplicationConfig, svc) -> list[tuple[str, Any]]:
    outputs = dict(svc.artifacts.get("outputs", {}))
    if config.svc_type == "sp-SVC":
        required = (("svc", "sp_svc"),)
    elif config.mode == "cluster":
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
    if config.mode == "cluster" and getattr(config, "ist_mapping", "paired") != "paired":
        from .ist_assembly import assemble_ist

        assembled = assemble_ist(
            outputs["sc_svc_spatial"], outputs["sc_svc_expr"],
            mapping=config.ist_mapping, seed=config.seed,
        )
        return [("svc", assembled)]
    return [(role, outputs[key]) for role, key in required]


def _owned_alternatives(config, paths):
    """Find only same-name application outputs superseded by a mode switch."""
    if config.mode != "cluster" or not isinstance(config, ApplicationConfig):
        return []
    import h5py

    alternatives = {}
    for mapping in ("paired", "mean"):
        alternatives.update(output_paths(replace(config, ist_mapping=mapping)))
    owned = []
    for path in set(alternatives.values()) - set(paths.values()):
        if not path.is_file() or path.is_symlink():
            continue
        try:
            with h5py.File(path, "r") as handle:
                metadata = handle["uns/revise_reconstruction"]
                mode = metadata["application_mode"].asstr()[()]
                selected = metadata["selected_cell_type"].asstr()[()]
                if mode != "cluster" or selected != config.select_cell_type:
                    continue
                if Path(metadata["output_dir"].asstr()[()]).resolve() != config.output_dir.resolve():
                    continue
                if Path(metadata["output_root"].asstr()[()]).resolve() != config.output_root.resolve():
                    continue
                role = metadata["output_role"].asstr()[()]
                recorded_path = metadata["output_paths"][role].asstr()[()]
                if Path(recorded_path).resolve() != path.resolve():
                    continue
                recorded_inputs = metadata.get("resolved_inputs")
                if recorded_inputs is not None and any(
                    name in recorded_inputs
                    and Path(recorded_inputs[name].asstr()[()]).resolve() != Path(value).resolve()
                    for name, value in config.resolved_inputs.items()
                ):
                    continue
                owned.append(path)
        except (OSError, KeyError, TypeError, AttributeError):
            continue
    return owned


def publish_outputs(config: ApplicationConfig, paths: Mapping[str, Path], ctx):
    """Publish pipeline carriers or assembled iST output transactionally."""
    artifacts = _published_artifacts(config, ctx.svc)
    config.output_dir.mkdir(parents=True, exist_ok=True)
    metadata = dict(ctx.application_config_metadata)
    metadata.update({
        "svc_type": config.svc_type,
        "mode": config.mode,
        "application_route": config.svc_type,
        "application_mode": config.mode,
        "output_name": config.output_name,
        "profile": ctx.profile,
        "run_manifest": str(Path(ctx.run_dir) / "provenance.json"),
        "selected_cell_type": config.select_cell_type,
        "ot": ctx.merged_config.get("ot"),
    })
    temporary: list[tuple[Path, Path]] = []
    try:
        for role, adata in artifacts:
            assembly = (
                adata.uns.get("revise_reconstruction", {})
                if config.mode == "cluster" and getattr(config, "ist_mapping", "paired") != "paired"
                else {}
            )
            adata.uns["revise_reconstruction"] = dict(metadata, **assembly, output_role=role)
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
        backups: dict[Path, Path] = {}
        published: set[Path] = set()
        publication_records: list[dict[str, Any]] = []

        def commit() -> None:
            for backup in backups.values():
                backup.unlink(missing_ok=True)
            backups.clear()

        def rollback() -> None:
            for target in published:
                target.unlink(missing_ok=True)
            published.clear()
            for target, backup in tuple(backups.items()):
                if backup.exists():
                    os.replace(backup, target)
            backups.clear()
            if publication_records:
                ctx.artifact_records[:] = [
                    record
                    for record in ctx.artifact_records
                    if record not in publication_records
                ]

        manages_publication = hasattr(ctx, "set_pending_publication")
        if manages_publication:
            ctx.set_pending_publication(commit=commit, rollback=rollback)
        try:
            targets = [target for _, target in temporary] + _owned_alternatives(config, paths)
            for target in targets:
                if target.exists():
                    with NamedTemporaryFile(
                        dir=config.output_dir,
                        prefix=f".{target.name}.",
                        suffix=".backup",
                        delete=False,
                    ) as handle:
                        backup = Path(handle.name)
                    os.replace(target, backup)
                    backups[target] = backup
            for temporary_path, target in temporary:
                os.replace(temporary_path, target)
                published.add(target)
            publication_records.extend(
                completed_artifact(f"publication:{role}", target)
                for role, target in paths.items()
            )
            for artifact in publication_records:
                ctx.record_artifact(artifact)
        except BaseException:
            rollback()
            raise
        else:
            if not manages_publication:
                commit()
    finally:
        for temporary_path, _ in temporary:
            temporary_path.unlink(missing_ok=True)

    values = tuple(adata for _, adata in artifacts)
    return values if config.mode == "cluster" and getattr(config, "ist_mapping", "paired") == "paired" else values[0]


__all__ = ["application_metadata", "output_paths", "publish_outputs"]
