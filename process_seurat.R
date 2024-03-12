#!/usr/bin/env Rscript

library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
filtered_matrix <- args[1]
sample_name <- args[2]

sc_10x <- Read10X(data.dir = paste0(filtered_matrix, "/matrix"))
get_seurat_obj <- function(sc_count, sample_name) {
    sc_obj <- CreateSeuratObject(counts = sc_count, project = sample_name, min.cells = 3, min.features = 200)
    sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")
    return(sc_obj)
}

sc_obj <- get_seurat_obj(sc_10x, sample_name)
