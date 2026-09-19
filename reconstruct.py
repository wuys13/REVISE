#!/usr/bin/env python3
"""Readable Application reconstruction entrypoint."""

from __future__ import annotations

import argparse
from contextlib import redirect_stdout
from dataclasses import replace
from pathlib import Path
import sys
from anndata import AnnData

from revise.application.config import (
    ApplicationConfig,
    ApplicationConfigError,
    compile_application_config,
    _compile_engine_config,
    load_application_yaml,
    override_select_cell_type,
)
from revise.application.inputs import load_data, preprocess_data
from revise.application.publication import application_metadata, output_paths, publish_outputs
from revise.application.delivery import is_sample_delivery, source_identity, validate_raw
from revise.application.expression import bind_expression_sources
from revise.framework import REVISEPipeline
from revise.io import REVISEInputService
from revise.utils.spot_sr_input import SR_CELL_COUNT_PROVENANCE_KEY

def reconstruct(
    spatial_adata: AnnData,
    reference_adata: AnnData,
    config: ApplicationConfig,
    *,
    raw_adata: AnnData | None = None,
    raw_source_identity: dict | None = None,
) -> AnnData | tuple[AnnData, AnnData]:
    """Run GA/LR and publish the exact returned Application artifact(s)."""
    config = bind_expression_sources(
        config, raw_adata if raw_adata is not None else spatial_adata, reference_adata
    )
    if config.reference_preparation is not None:
        from revise.reference_preparation.evidence import assert_reference_unchanged

        assert_reference_unchanged(config.reference_preparation)
    raw_snapshot = None
    if is_sample_delivery(config):
        if raw_adata is not None:
            validate_raw(raw_adata)
            # The caller may pass the same object as the working input. Freeze
            # the explicit original before any runner can normalize in place.
            raw_snapshot = raw_adata.copy()
        elif raw_source_identity is None:
            # A direct in-memory call may use its configured original source.
            raw_source_identity = source_identity(config)
    runtime, io, algorithm = _compile_engine_config(config)
    paths = output_paths(config)
    metadata = application_metadata(config, paths=paths)
    cell_count = spatial_adata.uns.get(SR_CELL_COUNT_PROVENANCE_KEY)
    if cell_count is not None:
        metadata["cell_count"] = dict(cell_count)
    result: AnnData | tuple[AnnData, AnnData] | None = None

    def finalize(ctx) -> None:
        nonlocal result
        if config.reference_preparation is not None:
            assert_reference_unchanged(config.reference_preparation)
        original = raw_snapshot
        if is_sample_delivery(config) and original is None:
            if source_identity(config) != raw_source_identity:
                raise ValueError("Original spatial input changed during reconstruction")
            original = REVISEInputService(io_config=io).read_st_adata(config.st_path)
            if source_identity(config) != raw_source_identity:
                raise ValueError("Original spatial input changed while loading Raw")
        if is_sample_delivery(config):
            result = publish_outputs(config, paths, ctx, raw=original, raw_source=raw_source_identity,
                                     raw_is_owned=True)
        else:
            result = publish_outputs(config, paths, ctx)

    REVISEPipeline().run(
        svc_type=config.svc_type,
        application_mode=config.mode,
        cf=None,
        runtime_overrides=runtime,
        io_overrides=io,
        algorithm_overrides=algorithm,
        st_adata=spatial_adata,
        sc_ref_adata=reference_adata,
        finalize_callback=finalize,
        application_config_metadata={
            **metadata,
            "reference_filter": {
                "column": config.reference_filter_column,
                "value": config.reference_filter_value,
            },
        },
    )
    if result is None:
        raise RuntimeError("Application reconstruction completed without publication")
    return result


def run_application(
    config_path: str | Path,
    *,
    select_ct: str | None = None,
    reference_config: str | Path | None = None,
) -> AnnData | tuple[AnnData, AnnData]:
    """Compile YAML, load inputs, preprocess, reconstruct, and return AnnData."""
    source, document = load_application_yaml(config_path)
    reference_path = reference_evidence = None
    if reference_config is not None:
        from revise.reference_preparation.evidence import resolve_reference_input

        reference_path, reference_evidence = resolve_reference_input(reference_config)
    if reference_path is None:
        config = compile_application_config(document, source=source)
    else:
        config = compile_application_config(
            document, source=source, reference_override=reference_path
        )
        config = replace(config, reference_preparation=reference_evidence)
    config = override_select_cell_type(config, select_ct)
    original_identity = source_identity(config) if is_sample_delivery(config) else None
    spatial_adata, reference_adata = load_data(config)
    if reference_evidence is not None:
        from revise.reference_preparation.evidence import assert_reference_unchanged

        assert_reference_unchanged(reference_evidence)
    if original_identity is not None and source_identity(config) != original_identity:
        raise ValueError("Original spatial input changed while loading")
    config = bind_expression_sources(config, spatial_adata, reference_adata)
    spatial_adata, reference_adata = preprocess_data(
        spatial_adata,
        reference_adata,
        config,
    )
    if is_sample_delivery(config):
        return reconstruct(spatial_adata, reference_adata, config, raw_source_identity=original_identity)
    return reconstruct(spatial_adata, reference_adata, config)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Reconstruct one SVC from an Application YAML")
    parser.add_argument("--config", required=True, help="Application YAML")
    parser.add_argument(
        "--reference-config", help="Optional prepared reference.yaml; replaces reference and filters"
    )
    parser.add_argument(
        "--select-ct",
        help="Override local_refinement.select_cell_type for sc-SVC cluster mode",
    )
    args = parser.parse_args(argv)
    try:
        with redirect_stdout(sys.stderr):
            kwargs = {"select_ct": args.select_ct}
            if args.reference_config is not None:
                kwargs["reference_config"] = args.reference_config
            result = run_application(args.config, **kwargs)
    except ApplicationConfigError as exc:
        parser.error(str(exc))
    print("Finished")
    if isinstance(result, tuple):
        print(f"spatial: {result[0]}")
        print(f"expression: {result[1]}")
    else:
        print(result)


__all__ = [
    "ApplicationConfigError",
    "load_data",
    "preprocess_data",
    "reconstruct",
    "run_application",
]


if __name__ == "__main__":
    main()
