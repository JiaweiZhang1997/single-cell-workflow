#!/usr/bin/env Rscript

library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
seurat_obj_file <- args[1]
sample_name <- args[2]
mt_thresh <- as.numeric(args[3])
max_features <- as.numeric(args[4])

load(seurat_obj_file)
sc_4261 <- readRDS(seurat_obj_file)

sc_qc_analyze <- function(seurat_obj, sample_name = "test", mt_thresh = 10, max_features = NULL, min_nCount_RNA = NULL, max_nCount_RNA = NULL) {
  #QC befor
  plot <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste0(sample_name, "_QC_befor.pdf"), width = 6, height = 6, plot)
  
  if (!is.null(mt_thresh)) {
    seurat_obj <- subset(seurat_obj, subset = percent.mt <= mt_thresh)
  }
  if (!is.null(max_features)) {
    seurat_obj <- subset(seurat_obj, subset = nFeature_RNA <= max_features)
  }
  if (!is.null(min_nCount_RNA)) {
    seurat_obj <- subset(seurat_obj, subset = nCount_RNA >= min_nCount_RNA)
  }
  if (!is.null(max_nCount_RNA)) {
    seurat_obj <- subset(seurat_obj, subset = nCount_RNA <= max_nCount_RNA)
  }
  #QC after
  plot_ <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste0(sample_name, "_QC_after.pdf"), width = 6, height = 6, plot_)
  return(seurat_obj)
}

sc_qc_analyze(sc_4261, sample_name, mt_thresh, max_features)