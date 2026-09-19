"""Package-owned loading and preprocessing for Application inputs."""

from __future__ import annotations

from anndata import AnnData

from revise.application.config import ApplicationConfig, _compile_engine_config
from revise.application.preprocess import (
    filter_reference,
    normalize_reference_labels,
    prepare_sc_svc_pair,
    preprocess_reference,
    preprocess_spatial,
)
from revise.io import REVISEInputService
from revise.utils.spot_sr_input import (
    ALL_CELLS_IN_SPOT_KEY,
    SR_CELL_COUNT_PROVENANCE_KEY,
    ensure_all_cells_in_spot,
)


def load_data(config: ApplicationConfig) -> tuple[AnnData, AnnData]:
    """Load fresh Application inputs through the package-owned I/O adapters."""
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
    """Apply the public Application preprocessing contract."""
    if config.mode == "sr":
        provided_mapping = (
            ALL_CELLS_IN_SPOT_KEY in spatial_adata.uns
            and spatial_adata.uns[ALL_CELLS_IN_SPOT_KEY] is not None
        )
        ensure_all_cells_in_spot(
            spatial_adata,
            cell_count_method=config.sr_cell_count_method,
        )
        if provided_mapping:
            cell_count = {
                "method": "provided_mapping",
                "source": "input.uns['all_cells_in_spot']",
            }
        else:
            cell_count = {
                "method": config.sr_cell_count_method,
                "source": (
                    "raw_full_gene_X"
                    if config.sr_cell_count_method == "cyto_linear_v1"
                    else "spot_transcript_counts"
                ),
            }
        spatial_adata.uns[SR_CELL_COUNT_PROVENANCE_KEY] = cell_count
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
    if config.mode == "cluster":
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


__all__ = ["load_data", "preprocess_data"]
