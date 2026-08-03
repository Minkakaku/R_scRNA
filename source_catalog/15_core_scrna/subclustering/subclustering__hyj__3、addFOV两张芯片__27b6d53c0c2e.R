result_dir <- 'result' 
dir.create(result_dir)

#add FOV画ImageDimPlot----
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

all_data <- readRDS("R_OUTPUT/all_data.rds") 

all_data2 <-  JoinLayers(all_data2)


###改行名
T1<- data_frame(rownames(all_data2@meta.data),all_data2@meta.data[c("orig.ident", "cell_id")])
# all_data2@meta.data$cell_id <- gsub("^.+?_(?=[a-z]{8}-\\d+$)", "",
#                                  rownames(all_data2@meta.data),
#                                  perl = TRUE)

# rownames(xenium@meta.data) <- gsub("TMA_","",rownames(xenium@meta.data))
all_data2@meta.data$cell_id_row <- rownames(all_data2@meta.data)
rownames(all_data2@meta.data) <- paste0(all_data2@meta.data$sample,"_",all_data2@meta.data$cell_id)

rownames(all_data2@assays$RNA@layers$count) <- rownames(all_data2)
colnames(all_data2@assays$RNA@layers$count) <- colnames(all_data2)
xenium <- CreateSeuratObject(
  all_data2@assays$RNA@layers$count,
  assay = "Xenium")

#saveRDS(xenium,"xenium_pro.rds")
#xenium <- readRDS("xenium_pro.rds")

xenium@meta.data <- all_data2@meta.data
umap <- all_data2@reductions$umap@cell.embeddings
rownames(umap) <- rownames(xenium@meta.data)
metadata <- xenium@meta.data
xenium <- Add_meta_umapData_Xenium(object = xenium,umapfile = umap,metedata = metadata)
colnames(xenium@assays$Xenium) <-  colnames(xenium)
row.names(xenium@reductions$umap@cell.embeddings) <- colnames(xenium)
xenium <- RenameCells(xenium,new.names = colnames(xenium))


#所有fov在一个文件夹内-----
#第一张芯片------
has_dt <- requireNamespace("data.table", quietly = TRUE) && requireNamespace("R.utils", quietly = TRUE)
data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/XE2505092086-01_yangrui/XE2505092086-01/output-XETG00305__0097541__tma1__20260318__080246"
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

tab <- table(xenium$sample, xenium$project)
samples_tma1 <- rownames(tab)[tab[, "tma1"] > 0]
samples_tma2 <- rownames(tab)[tab[, "tma2"] > 0]

print(samples_tma1)
print(samples_tma2)

#sample_label <- unique(xenium@meta.data$sample)
for (i in samples_tma1) {
  fov <- xenium_fov(object = xenium,label = i,assay = "Xenium")
  xenium@images[[i]] <- fov
}

#第二张芯片--------
has_dt <- requireNamespace("data.table", quietly = TRUE) && requireNamespace("R.utils", quietly = TRUE)
data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/XE2505092086-01_yangrui/XE2505092086-01/output-XETG00305__0097549__tma2__20260318__080247"
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

tab <- table(xenium$sample, xenium$project)
samples_tma1 <- rownames(tab)[tab[, "tma1"] > 0]
samples_tma2 <- rownames(tab)[tab[, "tma2"] > 0]

print(samples_tma1)
print(samples_tma2)

#sample_label <- unique(xenium@meta.data$sample)
for (i in samples_tma2) {
  fov <- xenium_fov(object = xenium,label = i,assay = "Xenium")
  xenium@images[[i]] <- fov
}



xenium <- NormalizeData(xenium)

saveRDS(xenium,"result/xenium.rds")

xenium <- readRDS("result/xenium.rds")

###测试一个样本----
DefaultBoundary(xenium@images[["IRDL3"]]) <- "segmentation"
p <- ImageFeaturePlot(
  object = xenium,
  fov = "IRDL3",
  nmols = 20000,       # 每个图像中显示的最大分子数（随机子抽样）
  border.color = "white",  # 细胞边界的颜色（深蓝色）
  border.size = 0.1,    # 细胞边界的线宽
  features = 'Cd3e',
  max.cutoff = "q90",
  axes = FALSE,
  cols = c('#003366', "#E7475E"),
  size = 1.5,
  dark.background = T
)
ggsave(p,filename = paste0("IRDL3","_feature_spatial.png"),width = 12 ,height = 12,dpi = 300)


###Imagedimplot----

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
                    dark.background = T,
                    border.color = "white", 
                    border.size = 0.1,
                    cols = all_colors,
                    shuffle.cols = F, #随机打乱颜色顺序
                    #size = 1.2,
                    group.by = "annotation"
                    #molecules = c("INHBA"),
                    #mols.size = 1,
                    #mols.cols = "red",
                    #mols.alpha = 0.5, #分子的透明度值
                    
  )
  ggsave(p,filename = paste0("result/Imagedimplot/",i,"_anno_spatial.png"),width = 12 ,height = 12,dpi = 300)
}


###SegImagefeature----
SegImagefeature_dir <- 'result/SegImagefeature'
dir.create(SegImagefeature_dir)

gene_list <- c('INHBA')
valid_gene <- gene_list[gene_list %in% rownames(xenium)]
cat("有效基因:", length(valid_gene), "/", length(gene_list))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_gene <- setdiff(gene_list, valid_gene)

gene_cus <- valid_gene

for (sample_id in unique(xenium@meta.data$sample)) {
  # 智能匹配样本ID
  j <- names(xenium@images)[grepl(sample_id, names(xenium@images))]
  
  if (length(j) == 0) {
    warning("跳过未匹配样本: ", sample_id)
    next
  }
  
  for (gene in gene_cus) {
    tryCatch({
      # 空值检查
      if (!is.null(xenium@images[[j]])) {
        DefaultBoundary(xenium@images[[j]]) <- "segmentation"
      }
      
      p <- ImageFeaturePlot(
        object = xenium,
        fov = j,
        nmols = 20000,       # 每个图像中显示的最大分子数（随机子抽样）
        border.color = "white",  # 细胞边界的颜色（深蓝色）
        border.size = 0.1,    # 细胞边界的线宽
        features = gene,
        max.cutoff = "q90",
        axes = FALSE,
        cols = c('#003366', "#E7475E"),
        size = 0.8,
        dark.background = TRUE
      )
      
      ggsave(
        filename = file.path(SegImagefeature_dir, paste0(sample_id, "_", gene, "_spatial.png")),
        plot = p,
        width = 12,
        height = 12,
        dpi = 300
      )
      
      message("成功保存: ", sample_id, "_", gene)
      
    }, error = function(e) {
      warning("失败[", sample_id, ":", gene, "]: ", e$message)
    })
  }
}


xenium@meta.data %>% head()

selected_data <- data.frame(
  cell_id = xenium@meta.data$cell_id,
  group = xenium@meta.data$annotation
)

# 检查数据是否包含缺失值
cat("cell_id列缺失值数量:", sum(is.na(selected_data$cell_id)), "\n")
cat("group列缺失值数量:", sum(is.na(selected_data$group)), "\n")

# 导出为CSV文件
write.csv(
  selected_data,
  file = "cell_annotation.csv",  # 文件名也相应修改
  row.names = FALSE
)



