#!/usr/bin/env Rscript

library(dplyr)
library(Seurat)
library(patchwork)
library(tibble)
library(limma)
library(scRepertoire)
library(tidyr)
library(gsubfn)
library(ggplot2)
library(celldex)
library(SingleR)

args <- commandArgs(trailingOnly = TRUE)
workdir <- args[1]
sample_name <- args[2]
mt_thresh <- as.numeric(args[3])
max_features <- as.numeric(args[4])
tcr_react <- args[5]
tcr_dir <- args[6]

sc_10x <- Read10X(data.dir = paste0(workdir, "/filtered_feature_bc_matrix"))
get_seurat_obj <- function(sc_count, sample_name) {
    sc_obj <- CreateSeuratObject(counts = sc_count, project = sample_name, min.cells = 3, min.features = 200)
    sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")
    return(sc_obj)
}

sc_obj <- get_seurat_obj(sc_10x, sample_name)

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

sc_qc <- sc_qc_analyze(sc_obj, sample_name, mt_thresh, max_features)

tcr <- read.csv(tcr_dir) %>% select(cdr3_aa1, cdr3_aa2, CTaa, barcode_raw)

#tumor reactive TCR
tcr_react_list <- read.table(tcr_react, sep = "\t", header = TRUE) %>% select(CD4.CD8, CDR3A_B, Archival.Prospective, Tumor.ID)

metadata_add_TCR <- function (seurat_obj, tcr_df, tcr_react_df){
  patient_id <- seurat_obj@meta.data$orig.ident[1]
  temp_mtx <- tcr_df %>% remove_rownames %>% column_to_rownames(var = "barcode_raw")
  temp_combined <- merge(seurat_obj@meta.data, temp_mtx, by = 0, sort = FALSE, all.x = TRUE)
  temp_combined <- merge(temp_combined, subset(tcr_react_df, Tumor.ID == patient_id), by.x = "CTaa", by.y = "CDR3A_B", sort = FALSE, all.x = TRUE) %>% remove_rownames %>% column_to_rownames(var = "Row.names")
  raw_barcode <- row.names(seurat_obj@meta.data)
  temp_combined <- temp_combined[raw_barcode, ]
  seurat_obj@meta.data <- temp_combined
  return(seurat_obj)
}

sc_qc_TCR <- metadata_add_TCR(sc_qc, tcr, tcr_react_list)
summary(sc_qc_TCR)


