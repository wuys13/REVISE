#!/usr/bin/env python3
"""Readable Application reconstruction entrypoint."""

from __future__ import annotations

import argparse
from contextlib import redirect_stdout
from pathlib import Path
import sys
from anndata import AnnData

from revise.application.config import (
    ApplicationConfig,
    ApplicationConfigError,
    compile_application_config,
    _compile_engine_config,
    load_application_yaml,
)
from revise.application.preprocess import (
    filter_reference,
    normalize_reference_labels,
    prepare_sc_svc_pair,
    preprocess_reference,
    preprocess_spatial,
)
from revise.application.publication import application_metadata, output_paths, publish_outputs
from revise.framework import REVISEPipeline
from revise.io import REVISEInputService
from revise.utils.spot_sr_input import ensure_all_cells_in_spot

def load_data(config: ApplicationConfig) -> tuple[AnnData, AnnData]:
    """Load Application inputs through the package-owned I/O adapters."""
    _, io, _ = _compile_engine_config(config)
    input_service = REVISEInputService(io_config=io)
    return (
        input_service.read_st_adata(config.st_path),
        input_service.read_sc_ref_adata(config.reference_path),
    )


def preprocess_data(
    spatial_adata: AnnData,
    reference_adata: AnnData,
    config: ApplicationConfig,
) -> tuple[AnnData, AnnData]:
    """Apply the visible Application preprocessing flow."""
    if config.svc_type == "sc-SVC-sr":
        ensure_all_cells_in_spot(spatial_adata)
    reference_adata = filter_reference(
        reference_adata,
        config.reference_filter_column,
        config.reference_filter_value,
    )
    spatial_adata = preprocess_spatial(
        spatial_adata,
        config.spatial_min_transcript_counts,
        config.spatial_min_cell_counts,
        min_counts=config.spatial_min_counts,
    )
    reference_adata = preprocess_reference(
        reference_adata,
        config.reference_min_transcript_counts,
        config.reference_min_cell_counts,
        min_genes=config.reference_min_genes,
    )
    if config.svc_type == "sc-SVC":
        return prepare_sc_svc_pair(
            spatial_adata,
            reference_adata,
            broad_column=config.broad_column,
            subtype_column=config.subtype_column,
        )
    return spatial_adata, normalize_reference_labels(
        reference_adata,
        (config.broad_column, config.subtype_column),
    )


def reconstruct(
    spatial_adata: AnnData,
    reference_adata: AnnData,
    config: ApplicationConfig,
) -> AnnData | tuple[AnnData, AnnData]:
    """Run GA/LR and publish the exact returned Application artifact(s)."""
    runtime, io, algorithm = _compile_engine_config(config)
    paths = output_paths(config)
    metadata = application_metadata(config, paths=paths)
    result: AnnData | tuple[AnnData, AnnData] | None = None

    def finalize(ctx) -> None:
        nonlocal result
        result = publish_outputs(config, paths, ctx)

    REVISEPipeline().run(
        svc_type=config.svc_type,
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


def run_application(config_path: str | Path) -> AnnData | tuple[AnnData, AnnData]:
    """Compile YAML, load inputs, preprocess, reconstruct, and return AnnData."""
    source, document = load_application_yaml(config_path)
    config = compile_application_config(document, source=source)
    spatial_adata, reference_adata = load_data(config)
    spatial_adata, reference_adata = preprocess_data(
        spatial_adata,
        reference_adata,
        config,
    )
    return reconstruct(spatial_adata, reference_adata, config)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Reconstruct one SVC from an Application YAML")
    parser.add_argument("--config", required=True, help="Application YAML")
    args = parser.parse_args(argv)
    try:
        with redirect_stdout(sys.stderr):
            result = run_application(args.config)
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
