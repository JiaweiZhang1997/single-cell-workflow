library(dplyr)
library(Seurat)
library(patchwork)
library(tcintr)
library(tibble)
library(limma)
library(scRepertoire)
library(tidyr)
library(gsubfn)
library(ggplot2)
library(celldex)
library(SingleR)

args <- commandArgs(trailingOnly = TRUE)
tcr_react <- args[1]
tcr_dir <- args[2]

tcr_4421 <- read.csv(tcr_dir) %>% select(cdr3_aa1, cdr3_aa2, CTaa, barcode_raw)

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

sc_4261_qc_TCR <- metadata_add_TCR(sc_4261_qc, tcr_4261, tcr_react_list)

