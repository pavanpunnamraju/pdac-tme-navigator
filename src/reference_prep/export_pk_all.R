# Convert Zenodo 6024273's pk_all.rds (Seurat object, all cells across 5
# integrated PDAC cohorts, $Cell_type = label-transferred cell-type calls per
# Code2 "Data integration.R" in the record's Code.zip) into a CellRanger-style
# mtx triplet + metadata CSV that Python/anndata can load without needing
# rpy2/SeuratDisk.
#
# RNA assay's "data" slot (log-normalized) is exported, not "counts": the
# 5 source cohorts were profiled on different protocols/depths, so raw counts
# aren't on a comparable scale across cohorts anyway, and the classifier only
# needs expression values, not depth-aware modeling.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
in_path <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pk_all <- readRDS(in_path)
DefaultAssay(pk_all) <- "RNA"

mat <- GetAssayData(pk_all, layer = "data")
cat("Expression matrix:", nrow(mat), "genes x", ncol(mat), "cells\n")

writeMM(mat, file.path(out_dir, "matrix.mtx"))
writeLines(rownames(mat), file.path(out_dir, "genes.tsv"))
writeLines(colnames(mat), file.path(out_dir, "barcodes.tsv"))

meta <- pk_all@meta.data
write.csv(meta, file.path(out_dir, "meta.csv"), row.names = TRUE)

cat("meta.data columns:", paste(colnames(meta), collapse = ", "), "\n")
cat("Done.\n")
