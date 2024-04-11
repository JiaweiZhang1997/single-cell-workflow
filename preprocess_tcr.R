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

tcr_path <- read.csv(file_in)
combined <- combineTCR(tcr_path, samples = 1, ID = 1, cells = "T-AB", removeNA = F, removeMulti = T, filterMulti = F)
barcode_raw <- str_split_fixed(combined[[1]]$barcode, "_", 3)[,3]
combined[[1]]$barcode_raw <- barcode_raw
write.csv(combined[[1]], file = paste0(file_in, ".format.csv"), row.names = FALSE)

