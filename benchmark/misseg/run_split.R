#!/usr/bin/env Rscript
# Run RCTD doublet deconvolution followed by SPLIT on matched P1CRC inputs.

suppressPackageStartupMessages({
  library(Matrix)
  library(spacexr)
  library(SPLIT)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript run_split.R <prepared-input-dir> <output-dir> <cores>")
}
input <- normalizePath(args[[1]], mustWork = TRUE)
output <- args[[2]]
cores <- as.integer(args[[3]])
dir.create(output, recursive = TRUE, showWarnings = FALSE)

read_lines <- function(name) readLines(file.path(input, name), warn = FALSE)
genes <- read_lines("genes.tsv")
spatial_cells <- read_lines("spatial_cells.tsv")
reference_cells <- read_lines("reference_cells.tsv")
coords <- read.csv(file.path(input, "spatial_coordinates.csv"), row.names = 1)
labels <- read.csv(file.path(input, "reference_labels.csv"), row.names = 1)

read_counts <- function(name, cells) {
  matrix <- as(readMM(gzfile(file.path(input, name))), "dgCMatrix")
  stopifnot(nrow(matrix) == length(genes), ncol(matrix) == length(cells))
  rownames(matrix) <- genes
  colnames(matrix) <- cells
  matrix
}
spatial_counts <- read_counts("spatial_genes_by_cells.mtx.gz", spatial_cells)
reference_counts <- read_counts("reference_genes_by_cells.mtx.gz", reference_cells)
stopifnot(identical(rownames(coords), spatial_cells))
stopifnot(identical(rownames(labels), reference_cells))
stopifnot(all(is.finite(as.matrix(coords[, c("x_um", "y_um")]))))

reference <- Reference(
  counts = reference_counts,
  cell_types = setNames(factor(labels$cell_type), reference_cells),
  require_int = TRUE
)
spatial <- SpatialRNA(
  coords = coords[, c("x_um", "y_um")],
  counts = spatial_counts,
  require_int = TRUE
)
message("RCTD: ", ncol(spatial_counts), " spatial cells, ", ncol(reference_counts),
        " reference cells, ", nrow(spatial_counts), " genes")

rctd_path <- file.path(output, "rctd.rds")
if (file.exists(rctd_path)) {
  rctd <- readRDS(rctd_path)
} else {
  rctd <- create.RCTD(spatial, reference, max_cores = cores)
  rctd <- run.RCTD(rctd, doublet_mode = "doublet")
  saveRDS(rctd, rctd_path, compress = "xz")
}
rctd <- SPLIT::run_post_process_RCTD(rctd)
message("RCTD status: ", paste(capture.output(print(table(rctd@results$results_df$spot_class))), collapse = " "))

converted <- SPLIT::convert_rctd_result_to_purify_input(rctd)
result <- SPLIT::purify(
  counts = spatial_counts,
  reference = Matrix::t(converted$reference),
  primary_cell_type = converted$primary_cell_type,
  deconvolution_weights = converted$deconvolution_weights,
  DO_output_sce = FALSE,
  DO_run_in_chunks = TRUE,
  chunk_size = 10000
)
stopifnot(!is.null(result$purified_counts), !is.null(result$cell_meta))
saveRDS(result, file.path(output, "split_result.rds"), compress = "xz")

purified <- as(result$purified_counts, "dgCMatrix")
writeMM(purified, file.path(output, "purified_counts.mtx"))
writeLines(rownames(purified), file.path(output, "genes.tsv"))
writeLines(colnames(purified), file.path(output, "cells.tsv"))
write.csv(result$cell_meta, file.path(output, "cell_metadata.csv"))
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
message("SPLIT output: ", ncol(purified), " cells, ", nrow(purified), " genes")
