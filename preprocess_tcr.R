args <- commandArgs(trailingOnly = TRUE)
file_in <- args[1]

library(dplyr)
library(Seurat)
library(patchwork)
library(tibble)
library(scRepertoire)
library(gridExtra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# print(file_in)
#conda activate R_4.0.5
#setwd("/share/home/xxwang/project/2023.05.25_CVAE_tumor_react")
#sc_data_10x_1 <- Read10X(data.dir = "/share/home/xxwang/project/2022.03.04_single_cell/PRJEB37120/SampleID_1_UM10")
#sc_data_1 <- CreateSeuratObject(counts = sc_data_10x_1, project = "s1", min.cells = 10, min.features = 200)
tcr_path <- read.csv(file_in)
combined <- combineTCR(tcr_path, samples = 1, ID = 1, cells = "T-AB", removeNA = F, removeMulti = T, filterMulti = F)
barcode_raw <- str_split_fixed(combined[[1]]$barcode, "_", 3)[,3]
combined[[1]]$barcode_raw <- barcode_raw
write.csv(combined[[1]], file = paste0(file_in, ".format.csv"), row.names = FALSE)

