#!/usr/bin/env Rscript

library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
filtered_matrix <- args[1]
sample_name <- args[2]

sc_10x <- Read10X(data.dir = paste0(filtered_matrix, "/matrix"))

# sc_genes_filter <- function (sc_count){
#   raw_genes_list <- rownames(sc_count)
#   TRAV_genes_list <- raw_genes_list[grep("^TRAV", raw_genes_list)]
#   TRAD_genes_list <- raw_genes_list[grep("^TRAD", raw_genes_list)]
#   TRAJ_genes_list <- raw_genes_list[grep("^TRAJ", raw_genes_list)]
#   TRAC_genes_list <- raw_genes_list[grep("^TRAC", raw_genes_list)]
  
#   TRBV_genes_list <- raw_genes_list[grep("^TRBV", raw_genes_list)]
#   TRBD_genes_list <- raw_genes_list[grep("^TRBD", raw_genes_list)]
#   TRBJ_genes_list <- raw_genes_list[grep("^TRBJ", raw_genes_list)]
#   TRBC_genes_list <- raw_genes_list[grep("^TRBC", raw_genes_list)]
  
#   TRDV_genes_list <- raw_genes_list[grep("^TRDV", raw_genes_list)]
#   TRDD_genes_list <- raw_genes_list[grep("^TRDD", raw_genes_list)]
#   TRDJ_genes_list <- raw_genes_list[grep("^TRDJ", raw_genes_list)]
#   TRDC_genes_list <- raw_genes_list[grep("^TRDC", raw_genes_list)]
  
#   TRGV_genes_list <- raw_genes_list[grep("^TRGV", raw_genes_list)]
#   TRGD_genes_list <- raw_genes_list[grep("^TRGD", raw_genes_list)]
#   TRGJ_genes_list <- raw_genes_list[grep("^TRGJ", raw_genes_list)]
#   TRGC_genes_list <- raw_genes_list[grep("^TRGC", raw_genes_list)]
  
#   # MT_genes_list <- raw_genes_list[grep("^MT-", raw_genes_list)]
  
#   # LINC_genes_list <- raw_genes_list[grep("^LINC", raw_genes_list)]
  
#   TR_genes_list <- c(TRAV_genes_list, TRAD_genes_list, TRAJ_genes_list, TRAC_genes_list,
#                      TRBV_genes_list, TRBD_genes_list, TRBJ_genes_list, TRBC_genes_list,
#                      TRDV_genes_list, TRDD_genes_list, TRDJ_genes_list, TRDC_genes_list,
#                      TRGV_genes_list, TRGD_genes_list, TRGJ_genes_list, TRGC_genes_list)
#   raw_genes_list_filter <- setdiff(raw_genes_list, TR_genes_list)
#   # sc_count_filter <- sc_count[raw_genes_list_filter, ]
#   return(raw_genes_list_filter)
# }

get_seurat_obj <- function(sc_count, sample_name) {
    sc_obj <- CreateSeuratObject(counts = sc_count, project = sample_name, min.cells = 3, min.features = 200)
    sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")
    return(sc_obj)
}

sc_obj <- get_seurat_obj(sc_10x, sample_name)
