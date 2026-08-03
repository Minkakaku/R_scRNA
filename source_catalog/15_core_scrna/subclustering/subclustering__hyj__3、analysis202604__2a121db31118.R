all_data <- readRDS("R_OUTPUT/all_data.rds")
all_data@meta.data %>% head()

#一、加入圈选区域------
all_data_ROI <- all_data
all_data_ROI@meta.data$ROI_old <- all_data_ROI@meta.data$ROI
all_data_ROI@meta.data$ROI <- NULL
all_data_ROI$cell_id <- gsub("TMA_","",rownames(all_data_ROI@meta.data))


##读取文件夹内所有CSV文件到列表-----
library(tidyverse)
folder_path <- "./20260417ROI/"

csv_files <- list.files(
  path = folder_path,
  pattern = "\\.csv$",
  full.names = TRUE
)

# 读取所有CSV文件到命名列表
dataframes <- map(csv_files, ~ read_csv(.x, skip = 2)) %>%
  set_names(tools::file_path_sans_ext(basename(csv_files)))

print(head(dataframes[["WT_2_Tumor"]]))

##更新ROI信息-----
all_data_ROI@meta.data$ROI <- all_data_ROI@meta.data$sample

obs <- all_data_ROI@meta.data %>% 
  mutate(
    ROI = as.character(ROI),
    cell_id = as.character(cell_id)
  )

walk2(names(dataframes), dataframes, ~ {
  cell_ids <- .y[["Cell ID"]] %>% as.character()
  obs[obs$cell_id %in% cell_ids, "ROI"] <<- .x
})

obs <- obs %>%
  mutate(
    ROI = factor(ROI),
    cell_id = factor(cell_id)
  )
obs$ROI %>% table()

all_data_ROI@meta.data <- obs

##排除未在ROI内的样本,并更新对象------
exclude_samples <- unique(all_data_ROI@meta.data$sample)
all_data_ROI <- all_data_ROI[, !all_data_ROI@meta.data$ROI %in% exclude_samples]

all_data_ROI$ROI <- as.character(all_data_ROI$ROI)
all_data_ROI$ROI <- as.factor(all_data_ROI$ROI)
all_data_ROI$ROI %>% table()

##构建区域分组-----
all_data_ROI$area <- case_when(
  grepl("Tumor", all_data_ROI$ROI)      ~ "Tumor",
  grepl("Boundary", all_data_ROI$ROI)   ~ "Boundary",
  grepl("Adjacent", all_data_ROI$ROI)   ~ "Adjacent",
  TRUE ~ "Other"
)
table(all_data_ROI$area,all_data_ROI$ROI)

all_data_ROI$ROI_group <- paste0(all_data_ROI$group, "_",all_data_ROI$area)
table(all_data_ROI$ROI_group,all_data_ROI$ROI)

table(all_data_ROI$sample,all_data_ROI$ROI)
saveRDS(all_data_ROI,"result/all_data_ROI.rds")
#all_data_ROI <- readRDS("result/all_data_ROI.rds")

result_dir <- 'result' 
dir.create(result_dir)

#二、all_data中add imagesFOV----
source("/mnt/NexData_file/home/hyj/h5ad2rds(1).R")
library(Seurat)
library(SeuratDisk)
library(sceasy)
library(reticulate)
library(dplyr)
library(ggplot2)
library(argparser)
library(patchwork)
library(RColorBrewer)
library(arrow)

all_data2 <-  JoinLayers(all_data,layers = 'count')
rownames(all_data2@assays$RNA@layers$count) <- rownames(all_data2)
colnames(all_data2@assays$RNA@layers$count) <- colnames(all_data2)
xenium <- CreateSeuratObject(
  all_data2@assays$RNA@layers$count,
  assay = "Xenium")

#saveRDS(xenium,"xenium_pro.rds")
#xenium <- readRDS("xenium_pro.rds")

xenium@meta.data <- all_data@meta.data
###改行名
rownames(xenium@meta.data) <- gsub("TMA_","",rownames(xenium@meta.data))
xenium@meta.data$cell_id <- rownames(xenium@meta.data)
rownames(xenium@meta.data) <- paste0(xenium@meta.data$sample,"_",xenium@meta.data$cell_id)

umap <- all_data@reductions$umap@cell.embeddings
rownames(umap) <- rownames(xenium@meta.data)
metadata <- xenium@meta.data
xenium <- Add_meta_umapData_Xenium(object = xenium,umapfile = umap,metedata = metadata)
colnames(xenium@assays$Xenium) <-  colnames(xenium)
row.names(xenium@reductions$umap@cell.embeddings) <- colnames(xenium)
xenium <- RenameCells(xenium,new.names = colnames(xenium))

#saveRDS(xenium,"xenium_renames.rds")
#xenium_renames <- readRDS("xenium_renames.rds")

##所有fov在一个文件夹内-----
has_dt <- requireNamespace("data.table", quietly = TRUE) && requireNamespace("R.utils", quietly = TRUE)
data.dir = "/data/xenium/FWSP20255039/output-XETG00305__0077158__TMA__20251017__072720"
xenium@meta.data$label = xenium@meta.data$sample
centroids = {
  if (has_dt) {
    centroids <- as.data.frame(data.table::fread(file.path(data.dir,
                                                           "cells.csv.gz")))
  } else {
    centroids <- read.csv(file.path(data.dir, "cells.csv.gz"))
  }
  centroids <- data.frame(x = centroids$x_centroid,
                          y = centroids$y_centroid, cell = centroids$cell_id,
                          stringsAsFactors = FALSE)}

segmentations = {
  if (has_dt) {
    segmentations <- as.data.frame(data.table::fread(file.path(data.dir,
                                                               "cell_boundaries.csv.gz")))
  } else {
    segmentations <- read.csv(file.path(data.dir,
                                        "cell_boundaries.csv.gz"), stringsAsFactors = FALSE)
  }
  segmentations <- data.frame(x = segmentations$vertex_x,y = segmentations$vertex_y,
                              cell = segmentations$cell_id)}
mols.qv.threshold = 20
# microns = {
#   if (has_dt) {
#     microns <- data.table::fread(file.path(data.dir, "transcripts.csv"))
#     microns <- subset(microns, qv >= mols.qv.threshold)
#   } else {
#     microns <- read.csv(file.path(data.dir, "transcripts.csv"))
#     microns <- subset(microns, qv >= mols.qv.threshold)
#   }
#   microns <- data.frame(x = microns$x_location, y = microns$y_location,
#                         cell = microns$cell_id,
#                         gene = microns$feature_name, stringsAsFactors = FALSE)
# }

microns = {
  microns <- read_parquet(file.path(data.dir, "transcripts.parquet"))
  microns <- subset(microns, qv >= mols.qv.threshold)
  microns <- data.frame(x = microns$x_location, y = microns$y_location,
                        cell = microns$cell_id,
                        gene = microns$feature_name, stringsAsFactors = FALSE)
}

object <- xenium
assay = "Xenium"
xenium_fov <- function(object,label,assay){
  cells_data <- data.frame(cell_id = object@meta.data$cell_id,sample = object@meta.data$sample)
  cells <- subset(cells_data,sample == label)$cell_id
  features <- rownames(object)
  
  segmentations_sub <- subset(segmentations,cell %in% cells)
  centroids_sub <- subset(centroids,cell %in% cells)
  microns_sub <- subset(microns,gene %in% features & cell %in% cells)
  
  segmentations_sub$cell <- paste0(label,"_",segmentations_sub$cell)
  centroids_sub$cell <- paste0(label,"_",centroids_sub$cell)
  microns_sub$cell <- paste0(label,"_",microns_sub$cell)
  
  segmentations.data <- list(centroids = CreateCentroids(centroids_sub),
                             segmentation = CreateSegmentation(segmentations_sub))
  
  coords <- CreateFOV(coords = segmentations.data, type = c("segmentation",
                                                            "centroids"), molecules = microns_sub, assay = assay)
  return(coords)
  
}

sample_label <- unique(xenium@meta.data$sample)
for (i in sample_label) {
  fov <- xenium_fov(object = xenium,label = i,assay = "Xenium")
  xenium@images[[i]] <- fov
}

future::plan("sequential")
xenium <- NormalizeData(xenium)

names(xenium@images)
###测试一个样本----
DefaultBoundary(xenium@images[["WT_2_2"]]) <- "segmentation"
p <- ImageFeaturePlot(
  object = xenium,
  fov = "WT_2_2",
  nmols = 20000,       # 每个图像中显示的最大分子数（随机子抽样）
  border.color = "white",  # 细胞边界的颜色（深蓝色）
  border.size = 0.1,    # 细胞边界的线宽
  features = 'INHBA',
  max.cutoff = "q90",
  axes = FALSE,
  cols = c('#003366', "#E7475E"),
  size = 1.5,
  dark.background = T
)
ggsave(p,filename = paste0("feature_spatial.png"),width = 12 ,height = 12,dpi = 300)


###Imagedimplot_annotation----
all_colors <- c(
  "#A6CEE3", "#5B9EC9", "#2D82AF", "#7EBA98", "#98D277",
  "#52AF43", "#6F9E4C", "#DD9A88", "#FDBB69", "#FE982C",
  "#F78620", "#D9A295", "#feee7d", "#02c9c9", "#600473",
  "#B294C7", "#7D54A5", "#9E8099", "#F0EB99", "#DBB466",
  "#B15928", # 原始21色
  "#FF6B6B", "#4ECDC4", "#FF9F1C", "#7209B7", "#4361EE", 
  "#F15BB5", "#00BBF9", "#00F5D4", "#9B5DE5", "#00F5A7",
  "#FEE440", "#9C89B8", "#06D6A0"  # 新增13色
)

all_colors8 <- c(
  "#FF6B6B","#A6CEE3","#00F5A7", "#FDBB69","#feee7d", "#9E8099",
  "#00BBF9","#7209B7"
)

sample_label <- unique(xenium@meta.data$sample)
for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "black", 
                    border.size = 0.1,
                    cols = all_colors8,
                    #shuffle.cols = F, #随机打乱颜色顺序
                    #size = 1.2,
                    group.by = "annotation"
                    #molecules = c("INHBA"),
                    #mols.size = 1,
                    #mols.cols = "red",
                    #mols.alpha = 0.5, #分子的透明度值
                    
  )
  ggsave(p,filename = paste0("result/Imagedimplot/",i,"_anno_spatial.png"),width = 12 ,height = 12,dpi = 300)
  ggsave(p,filename = paste0("result/Imagedimplot/",i,"_anno_spatial.pdf"),width = 12 ,height = 12)
}

saveRDS(xenium, "result/xenium.rds")
#xenium <- readRDS("result/xenium.rds")


xenium@meta.data$ROI_old <- xenium@meta.data$ROI
xenium@meta.data$ROI <- NULL
##添加ROI信息到xenium@meta.data----
xenium@meta.data %>% head()
all_data_ROI@meta.data %>% head()

roi_meta <- all_data_ROI@meta.data[, c("cell_id", "ROI", "area", "ROI_group")]

roi_meta_sorted <- roi_meta[match(xenium@meta.data$cell_id, roi_meta$cell_id), ]

xenium$ROI <- roi_meta_sorted$ROI
xenium$area <- roi_meta_sorted$area
xenium$ROI_group <- roi_meta_sorted$ROI_group

head(cbind(colnames(xenium), roi_meta_sorted$cell_id, xenium$ROI, xenium$area, xenium$ROI_group))

table(xenium$ROI_group,xenium$sample)

#xenium_metata <- xenium@meta.data[, c("cell_id", "ROI", "area", "ROI_group")]

xenium$ROI <- as.character(xenium$ROI)
xenium$ROI[is.na(xenium$ROI)] <- "Other"
xenium$area[is.na(xenium$area)] <- "Other"
xenium$ROI_group[is.na(xenium$ROI_group)] <- "Other"
table(is.na(xenium$ROI))
table(is.na(xenium$area))
table(is.na(xenium$ROI_group))



###Imagedimplot_area----
xenium$area <- factor(xenium$area,levels = c("Tumor",   
                                             "Boundary",
                                             "Adjacent",
                                             "Other"))
all_colors3 <- c(
  "Tumor"    = "#7209B7",
  "Boundary" = "#9E8099",
  "Adjacent" = "#A6CEE3",
  "Other"    = "gray90"
)
sample_label <- unique(xenium@meta.data$sample)
for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "black", 
                    border.size = 0.1,
                    cols = all_colors3,
                    #shuffle.cols = F, #随机打乱颜色顺序
                    #size = 1.2,
                    group.by = "area"
                    #molecules = c("INHBA"),
                    #mols.size = 1,
                    #mols.cols = "red",
                    #mols.alpha = 0.5, #分子的透明度值
                    
  )
  ggsave(p,filename = paste0("result/Imagedimplot/",i,"_area_spatial.png"),width = 12 ,height = 12,dpi = 300)
  ggsave(p,filename = paste0("result/Imagedimplot/",i,"_area_spatial.pdf"),width = 12 ,height = 12)
}

saveRDS(xenium, "result/xenium.rds")

xenium <- readRDS("result/xenium.rds")

#三、xenium对应all_data_ROI中细胞--------

xenium_ROI <- subset(xenium, area != "Other")
##ROI在样本中-----
###Imagedimplot_annotation----
# all_colors8 <- c(
#   "#FF6B6B","#A6CEE3","#00F5A7", "#FDBB69","#feee7d", "#9E8099",
#   "#00BBF9","#7209B7"
# )
# 
# sample_label <- unique(xenium_ROI@meta.data$sample)
# for( i in sample_label){
#   #j = gsub("_","-",i)
#   DefaultBoundary(xenium_ROI@images[[i]]) <- "segmentation"
#   p <- ImageDimPlot(xenium_ROI,
#                     fov = i,
#                     nmols = 20000,
#                     dark.background = F,
#                     border.color = "black", 
#                     border.size = 0.1,
#                     cols = all_colors8,
#                     #shuffle.cols = F, #随机打乱颜色顺序
#                     #size = 1.2,
#                     group.by = "annotation"
#                     #molecules = c("INHBA"),
#                     #mols.size = 1,
#                     #mols.cols = "red",
#                     #mols.alpha = 0.5, #分子的透明度值
#                     
#   )
#   ggsave(p,filename = paste0("result/Imagedimplot_ROI/",i,"_anno_spatial.png"),width = 12 ,height = 12,dpi = 300)
#   ggsave(p,filename = paste0("result/Imagedimplot_ROI/",i,"_anno_spatial.pdf"),width = 12 ,height = 12)
# }


###Imagedimplot_area----
# all_colors3 <- c(
#   "Tumor"    = "#7209B7",
#   "Boundary" = "#9E8099",
#   "Adjacent" = "#A6CEE3"#,
#   #"Other"    = "gray90"
# )
# sample_label <- unique(xenium_ROI@meta.data$sample)
# for( i in sample_label){
#   #j = gsub("_","-",i)
#   DefaultBoundary(xenium_ROI@images[[i]]) <- "segmentation"
#   p <- ImageDimPlot(xenium_ROI,
#                     fov = i,
#                     nmols = 20000,
#                     dark.background = F,
#                     border.color = "black", 
#                     border.size = 0.1,
#                     cols = all_colors3,
#                     #shuffle.cols = F, #随机打乱颜色顺序
#                     #size = 1.2,
#                     group.by = "area"
#                     #molecules = c("INHBA"),
#                     #mols.size = 1,
#                     #mols.cols = "red",
#                     #mols.alpha = 0.5, #分子的透明度值
#                     
#   )
#   ggsave(p,filename = paste0("result/Imagedimplot_ROI/",i,"_area_spatial.png"),width = 12 ,height = 12,dpi = 300)
#   ggsave(p,filename = paste0("result/Imagedimplot_ROI/",i,"_area_spatial.pdf"),width = 12 ,height = 12)
# }

saveRDS(xenium_ROI, "result/xenium_ROI.rds")

#xenium_ROI <- readRDS("result/xenium_ROI.rds")


##单独ROI-------
###Imagedimplot_annotation----
all_colors8 <- c(
  "Endothelial"     = "#FF6B6B",    
  "Fibroblast"      = "#A6CEE3",    
  "Neuron"          = "#00F5A7",    
  "Oligodendrocyte" = "#FDBB69",    
  "Pericyte"        = "#feee7d",    
  "T cell"          = "#9E8099",    
  "TAM"             = "#00BBF9",    
  "Tumor"           = "#7209B7"     
)

table(xenium_ROI$ROI,xenium_ROI$sample)

Idents(xenium_ROI) <- "annotation"
Idents(xenium_ROI) %>% unique()

sample_label <- unique(xenium_ROI@meta.data$sample)

for( i in sample_label){
  
  dat1 <- subset(xenium_ROI,sample == i)
  
  for (r in unique(dat1$ROI)) {
    
    dat2 <- subset(dat1,ROI == r)
    
    DefaultBoundary(dat2@images[[i]]) <- "segmentation"
    p <- ImageDimPlot(dat2,
                      fov = i,
                      nmols = 20000,
                      dark.background = F,
                      border.color = "black", 
                      border.size = 0.1,
                      cols = all_colors8
    )+coord_equal()
    
    ggsave(p,filename = paste0("result/F1C Imagedimplot_ROI/",r,"_anno_spatial.png"),width = 15 ,height = 15,dpi = 300)
    ggsave(p,filename = paste0("result/F1C Imagedimplot_ROI/",r,"_anno_spatial.pdf"),width = 15 ,height = 15)
  }
  
}


###Imagedimplot_area----
table(xenium_ROI$ROI,xenium_ROI$sample)

Idents(xenium_ROI) <- "area"
all_colors3 <- c(
  "Tumor"    = "#7209B7",
  "Boundary" = "#9E8099",
  "Adjacent" = "#A6CEE3"#,
  #"Other"    = "gray90"
)

for( i in sample_label){
  
  dat1 <- subset(xenium_ROI,sample == i)
  
  for (r in unique(dat1$ROI)) {
    
    dat2 <- subset(dat1,ROI == r)
    
    DefaultBoundary(dat2@images[[i]]) <- "segmentation"
    p <- ImageDimPlot(dat2,
                      fov = i,
                      nmols = 20000,
                      dark.background = F,
                      border.color = "black", 
                      border.size = 0.1,
                      cols = all_colors3
    )+coord_equal()
    
    ggsave(p,filename = paste0("result/F1C Imagedimplot_ROI/",r,"_area_spatial.png"),width = 15 ,height = 15,dpi = 300)
    ggsave(p,filename = paste0("result/F1C Imagedimplot_ROI/",r,"_area_spatial.pdf"),width = 15 ,height = 15)
  }
  
}





