# 按照要求提取对应的rds
# 因为要做monocle2，需要提取对应的count数据
# 最终的rds数据来自上次已经整理好的 "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Main_mapping_count.rds"
# 需要提取的MP的信息来自于"/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou11/Save_data/MP_result_obs.csv"
rm(list = ls())
library(reticulate)
library(sceasy)
library(Seurat)
library(Matrix)
library(tidyverse)
setwd("/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou13/")

sce <- readRDS("/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Main_mapping_count.rds")
meta <- read.csv("/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou11/Save_data/MP_result_obs.csv",check.names = F,sep = '\t')
colnames(meta)[1] <- 'barcode'

subrds <- subset(sce,cells = meta$barcode)
mp_vector <- setNames(meta$MP_group, meta$barcode)
subrds$MP_group <- mp_vector[rownames(subrds@meta.data)]

saveRDS(object = subrds,file = "MP.rds")


# 
meta1 <- meta[meta$MP_group != "unresolved",]
subrds1 <- subset(sce,cells = meta1$barcode)
mp_vector <- setNames(meta1$MP_group, meta1$barcode)
subrds1$MP_group <- mp_vector[rownames(subrds1@meta.data)]

saveRDS(object = subrds1,file = "MP_remove_unsolved.rds")


library(jpeg)
library(pdftools,lib.loc = "/mnt/NexData_file/share/scRNA-tools/envs/pip_cellchat/lib/R/library/")

pdf_files <- list.files(path = "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou13/Monocle_result/MP/Pseudotime/", pattern = "*.pdf", full.names = TRUE)
for (i in seq_along(pdf_files)) {
  file_name <- pdf_files[i]
  pdf_image <- pdf_render_page(file_name, dpi = 300)
  writeJPEG(pdf_image, str_replace(string = file_name,pattern = "pdf",replacement = "png"))
}
cat("Processed", length(pdf_files), "PDF files.\n")
