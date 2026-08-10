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

all_data2@meta.data$sample %>% table()
all_data3 <-  JoinLayers(all_data2,layers = 'counts')
all_data3@meta.data %>% head()

T1<- data_frame(rownames(all_data3@meta.data),all_data3@meta.data[c("orig.ident", "cell_id")])
###改行名
all_data3@meta.data$cell_id <- gsub("^.+?_(?=[a-z]{8}-\\d+$)", "", 
                                    rownames(all_data3@meta.data), 
                                    perl = TRUE)

rownames(all_data3@meta.data) <- paste0(all_data3@meta.data$sample,"_",all_data3@meta.data$cell_id)

rownames(all_data3@assays$RNA@layers$count) <- rownames(all_data3)
colnames(all_data3@assays$RNA@layers$count) <- colnames(all_data3)
xenium <- CreateSeuratObject(
  all_data3@assays$RNA@layers$count,
  assay = "Xenium")

xenium@meta.data <- all_data3@meta.data

umap <- all_data3@reductions$umap@cell.embeddings
rownames(umap) <- rownames(xenium@meta.data)
metadata <- xenium@meta.data
xenium <- Add_meta_umapData_Xenium(object = xenium,umapfile = umap,metedata = metadata)
colnames(xenium@assays$Xenium) <-  colnames(xenium)
row.names(xenium@reductions$umap@cell.embeddings) <- colnames(xenium)

xenium <- RenameCells(xenium,new.names = colnames(xenium))
label <- unique(xenium@meta.data$sample)

options(future.globals.maxSize = 20 * 1024^3)
xenium <- NormalizeData(xenium)

#多个文件夹，一个文件夹内包含多个fov###----

##加入1个测试----
library(arrow)
coord1 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0097339__XENIUM3__20260306__090428",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "N16546-3-1",object = xenium,cov = NULL)
xenium@images[["N16546-3-1"]] = coord1

DefaultBoundary(xenium@images[["N16546-3-1"]]) <- "segmentation"
p <- ImageFeaturePlot(xenium, fov = "N16546-3-1",features = c("CD3E"),
                      max.cutoff = "q90",axes = FALSE,
                      border.color = "white", 
                      border.size = 0.1,
                      cols = c('#003366', '#FFD700'))
ggsave(p,filename = paste0("result/","_feature_spatial.png"),width = 12 ,height = 12,dpi = 300)


##加入12个样本FOV(正式)----
# coord1 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0097339__XENIUM3__20260306__090428",assay = 'Xenium',
#                            mols.qv.threshold = 20,label = "N16546-3-1",object = xenium,cov = NULL)
# xenium@images[["N16546-3-1"]] = coord1

coord2 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0097339__XENIUM3__20260306__090428",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "N16546-3-2",object = xenium,cov = NULL)
xenium@images[["N16546-3-2"]] = coord2

coord3 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0097339__XENIUM3__20260306__090428",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "N16546-3-3",object = xenium,cov = NULL)
xenium@images[["N16546-3-3"]] = coord3

coord4 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "N18717_G2_3",object = xenium,cov = NULL)
xenium@images[["N18717_G2_3"]] = coord4

coord5 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "N18717_G2_2",object = xenium,cov = NULL)
xenium@images[["N18717_G2_2"]] = coord5

coord6 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "N18717_G2_1",object = xenium,cov = NULL)
xenium@images[["N18717_G2_1"]] = coord6

coord7 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "N18530_G1_3",object = xenium,cov = NULL)
xenium@images[["N18530_G1_3"]] = coord7

coord8 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "N18530_G1_1",object = xenium,cov = NULL)
xenium@images[["N18530_G1_1"]] = coord8

coord9 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "N16326_G1_3",object = xenium,cov = NULL)
xenium@images[["N16326_G1_3"]] = coord9

coord10 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "N16326_G1_2",object = xenium,cov = NULL)
xenium@images[["N16326_G1_2"]] = coord10

coord11 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "N16326_G1_1",object = xenium,cov = NULL)
xenium@images[["N16326_G1_1"]] = coord11



###绘图一个样本测试----
DefaultBoundary(xenium@images[["N16326_G1_1"]]) <- "segmentation"
p <- ImageFeaturePlot(
  object = xenium,
  fov = "N16326_G1_1",
  nmols = 20000,       # 每个图像中显示的最大分子数（随机子抽样）
  border.color = "white",  # 细胞边界的颜色（深蓝色）
  border.size = 0.1,    # 细胞边界的线宽
  features = 'CD3E',
  max.cutoff = "q90",
  axes = FALSE,
  cols = c('#003366', "#E7475E"),
  size = 1.5,
  dark.background = T
)
ggsave(p,filename = paste0("result/sham1","_feature_spatial.png"),width = 12 ,height = 12,dpi = 300)


#saveRDS(xenium,"result/xenium_FOV.rds")
xenium <- readRDS("result/xenium_FOV.rds")

##Imagedimplot----
xenium@meta.data$annotation %>% table()
xenium@meta.data$annotation<- as.factor(xenium@meta.data$annotation)

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


all_colors10 <- c(
  "B_cell" = "#00F5A7", 
  "Endothelial" = "#98D277",
  "Epithelial" = "#9B5DE5",
  "Fibroblast" = "#FDBB69",
  "Myeloid" = "#feee7d", #"#9E8099",
  "Pericyte" = "#00BBF9",
  "Plasma" =  "#52AF43",
  "SMC" ="#A6CEE3",
  "T_cell" = "#D9A295",
  "Uncertain" = "gray90"
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
                    cols = all_colors11,
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



source("ImageFeature-plot.R")
##空间featureplt----
Module_Imagefeature_dir <- 'result/Module_Imagefeature' 
dir.create(Module_Imagefeature_dir)

module_names <- new_names

###只画cells细胞,其他不显示--------
for (sample_id in unique(xenium_score@meta.data$sample)) {
  for (module in module_names) {
    
    p <- ImageFeaturePlot(
      object = xenium_score,
      fov = sample_id,
      nmols = 20000,       
      border.color = "white",  
      border.size = 0.01,    
      features = module,
      max.cutoff = "q90",
      axes = FALSE,
      cols = c('#003366', '#FFD700'),
      # cols = c('blue', "red"),
      size = 0.8,
      dark.background = TRUE#,
      #cells = WhichCells(xenium_score, idents = "Macrophage" )
    )
    
    ggsave(
      filename = file.path(Module_Imagefeature_dir, paste0(sample_id, "_", module,"_spatial.png")),
      plot = p,
      width = 12,
      height = 12,
      dpi = 300
    )
    ggsave(
      filename = file.path(Module_Imagefeature_dir, paste0(sample_id, "_", module,"_spatial.pdf")),
      plot = p,
      width = 12,
      height = 12
    )
    
    message("成功保存: ", sample_id, "_", module)
  }
}

###只画cells.use 细胞,其他灰色--------
for (sample_id in unique(xenium_score@meta.data$sample)){
  for (module in module_names) {
    
    p <- ImageFeaturePlot_wyt(xenium_score, 
                              features = module, 
                              size = 0.5, 
                              cols = c('#003366', '#FFD700'),
                              #cols = c('blue', "red"),
                              fov =sample_id,
                              nmols = 20000,
                              border.color = all_colors8, 
                              border.size = 0.2,
                              border.by = "annotation",
                              #cells.use = c("Macrophage"),
                              axes = FALSE,
                              dark.background = F
    )
    
    ggsave(
      filename = file.path(Module_Imagefeature_dir, paste0(sample_id, "_", module,"_spatial.png")),
      plot = p,
      width = 12,
      height = 12,
      dpi = 300
    )
    ggsave(
      filename = file.path(Module_Imagefeature_dir, paste0(sample_id, "_", module,"_spatial.pdf")),
      plot = p,
      width = 12,
      height = 12
    )
    
    message("成功保存: ", sample_id, "_", module)
  }
}


##把五个基因高表达的区域空间上映射出来------
name_mapping <- c(
  "N16546-3-1" = "N16546-3-1",  
  "N16546-3-2" = "N16546-3-2",  
  "N16546-3-3" = "N16546-3-3",  
  "N16326_G1_1" = "N16326-G1-1", 
  "N16326_G1_2" = "N16326-G1-2", 
  "N16326_G1_3" = "N16326-G1-3", 
  "N18530_G1_1" = "N18530-G1-1", 
  "N18530_G1_3" = "N18530-G1-3", 
  "N18717_G2_1" = "N18717-G2-1", 
  "N18717_G2_2" = "N18717-G2-2", 
  "N18717_G2_3" = "N18717-G2-3"
)

names(xenium@images) <- name_mapping[names(xenium@images)]
names(xenium@images)
#sample_label <- unique(xenium@meta.data$sample)
sample_label <- names(xenium@images)

xenium@meta.data$annotation %>% unique()
xenium@meta.data$annotation<- as.factor(xenium@meta.data$annotation)


all_colors10 <- c(
  "B_cell" = "#00F5A7", 
  "Endothelial" = "#98D277",
  "Epithelial" = "#9B5DE5",
  "Fibroblast" = "#FDBB69",
  "Myeloid" = "#feee7d", #"#9E8099",
  "Pericyte" = "#00BBF9",
  "Plasma" =  "#52AF43",
  "SMC" ="#A6CEE3",
  "T_cell" = "#D9A295",
  "Uncertain" = "gray90"
) 


all_colors10 <- c(
  "B_cell" = "gray90", 
  "Endothelial" = "#98D277",
  "Epithelial" = "#9B5DE5",
  "Fibroblast" = "#FDBB69",
  "Myeloid" = "gray90", #"#9E8099",
  "Pericyte" = "gray90",
  "Plasma" = "gray90",
  "SMC" = "gray90", 
  "T_cell" = "gray90",
  "Uncertain" = "gray90"
) 

genelist <- c("EFNA1","EPCAM","EPHA3","SMA","COL5A1","PECAM1")
valid_gene_sets <- genelist[genelist %in% rownames(xenium)]
cat("有效基因:", length(valid_gene_sets), "/", length(genelist))
missing_genes <- setdiff(genelist, valid_gene_sets)

Idents(xenium) <- "annotation"
for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "black", 
                    border.size = 0.1,
                    cols = all_colors10,
                    #cells = WhichCells(xenium, idents = c("Endothelial","Epithelial","Fibroblast")),
                    molecules = genelist,
                    mols.size = 0.005,
                    mols.alpha = 0.7,
                    mols.cols    = c(
                      "EFNA1"   = "#FF6B6B",         
                      "EPCAM"   = "#FEE440",        
                      "EPHA3"   = "#2D82AF",       
                      "COL5A1"  = "#02c9c9",   
                      "PECAM1"  = "#7209B7"  
                    )
  )
  
  ggsave(p,filename = paste0("result/Imagedimplot/",i,"_alltype_5gene_spatial.png"),width = 15 ,height = 15,dpi = 300)
  ggsave(p,filename = paste0("result/Imagedimplot/",i,"_alltype_5gene_spatial.pdf"),width = 15 ,height = 15)
}




#20260408analysis------
##加入老师圈选的ROI-------
xenium@meta.data %>% head()

#1.读取文件夹内所有CSV文件到列表
library(tidyverse)

folder_path <- "./result/ROI/"

csv_files <- list.files(
  path = folder_path,
  pattern = "\\.csv$",
  full.names = TRUE
)

# 读取所有CSV文件到命名列表（跳过前两行）
dataframes <- map(csv_files, ~ read_csv(.x, skip = 2)) %>%
  set_names(tools::file_path_sans_ext(basename(csv_files)))

# 验证
print(head(dataframes[["XENIUM3_1"]]))

names(dataframes)

dataframes <- imap(dataframes, function(df, sample_name) {
  
  # 处理样本名：去掉最后一个 _ 及其后面的所有内容
  sample_prefix <- sub("_[^_]*$", "", sample_name)
  
  df %>%
    mutate(
      names = sample_prefix,
      cell_id_names = paste(names, `Cell ID`, sep = "_"),  
                                          
    )
})


#2. 更新观测数据（obs）的样本信息
xenium@meta.data$ROI <- xenium@meta.data$sample
xenium@meta.data$project %>% table()
xenium@meta.data$names <- gsub("XENIUM-6","XENIUM6",xenium@meta.data$project)
xenium@meta.data$names %>% table()
xenium@meta.data$cell_id <- paste0(xenium@meta.data$names,"_",xenium@meta.data$cell_id_raw)
xenium@meta.data$cell_id %>% head()

obs <- xenium@meta.data %>% 
  mutate(
    ROI = as.character(ROI),
    cell_id = as.character(cell_id)
  )

for(roi_name in names(dataframes)){
  
  df <- dataframes[[roi_name]]
  cell_ids_in_roi <- df$cell_id_names

  obs$ROI[obs$cell_id %in% cell_ids_in_roi] <- roi_name
}

obs %>% head()

xenium@meta.data <- obs
xenium@meta.data %>% head()
Idents(xenium) <- "annotation"
xenium@meta.data$ROI %>% table()


#3. 排除指定样本并更新对象
exclude_samples <- unique(xenium@meta.data$sample)
# 过滤Seurat对象
xenium_ROI <- xenium[,!xenium@meta.data$ROI %in% exclude_samples]
xenium_ROI$ROI <- as.character(xenium_ROI$ROI)

saveRDS(xenium_ROI,"result/xenium_ROI.rds")
xenium_ROI <- readRDS("result/xenium_ROI.rds")

###ROI中细胞映射-----
all_colors10 <- c(
  "B_cell" = "#00F5A7", 
  "Endothelial" = "#98D277",
  "Epithelial" = "#9B5DE5",
  "Fibroblast" = "#FDBB69",
  "Myeloid" = "#feee7d", #"#9E8099",
  "Pericyte" = "#00BBF9",
  "Plasma" =  "#52AF43",
  "SMC" ="#A6CEE3",
  "T_cell" = "#D9A295",
  "Uncertain" = "gray90"
) 


all_colors10 <- c(
  "B_cell" = "gray90", 
  "Endothelial" = "#98D277",
  "Epithelial" = "#9B5DE5",
  "Fibroblast" = "#FDBB69",
  "Myeloid" = "gray90", #"#9E8099",
  "Pericyte" = "gray90",
  "Plasma" = "gray90",
  "SMC" = "gray90", 
  "T_cell" = "gray90",
  "Uncertain" = "gray90"
) 

genelist <- c("EFNA1","EPCAM","EPHA3","SMA","COL5A1","PECAM1")
valid_gene_sets <- genelist[genelist %in% rownames(xenium_ROI)]
cat("有效基因:", length(valid_gene_sets), "/", length(genelist))
missing_genes <- setdiff(genelist, valid_gene_sets)


for(i in sample_label){
  xenium_ROI@images[[i]] <- NULL
}

name_mapping <- c(
  "N16546.3.1" = "N16546-3-1",  
  "N16546.3.2" = "N16546-3-2",  
  "N16546.3.3" = "N16546-3-3",  
  "N16326.G1.1" = "N16326-G1-1", 
  "N16326.G1.2" = "N16326-G1-2", 
  "N16326.G1.3" = "N16326-G1-3", 
  "N18530.G1.1" = "N18530-G1-1", 
  "N18530.G1.3" = "N18530-G1-3", 
  "N18717.G2.1" = "N18717-G2-1", 
  "N18717.G2.2" = "N18717-G2-2", 
  "N18717.G2.3" = "N18717-G2-3"
)

names(xenium_ROI@images) <- name_mapping[names(xenium_ROI@images)]
names(xenium_ROI@images)
ROI_label <- names(xenium_ROI@images)

Idents(xenium_ROI) <- "annotation"

# for( i in ROI_label){
#   #j = gsub("_","-",i)
#   DefaultBoundary(xenium_ROI@images[[i]]) <- "segmentation"
#   p <- ImageDimPlot(xenium_ROI,
#                     fov = i,
#                     nmols = 20000,
#                     dark.background = F,
#                     border.color = "black", 
#                     border.size = 0.1,
#                     cols = all_colors10#,
#                     #cells = WhichCells(xenium_ROI, idents = c("Endothelial","Epithelial","Fibroblast")),
#                     # molecules = genelist,
#                     # mols.size = 0.005,
#                     # mols.alpha = 0.7,
#                     # mols.cols    = c(
#                     #   "EFNA1"   = "#FF6B6B",         
#                     #   "EPCAM"   = "#FEE440",        
#                     #   "EPHA3"   = "#2D82AF",       
#                     #   "COL5A1"  = "#02c9c9",   
#                     #   "PECAM1"  = "#7209B7"  
#                     )
#   
#   
#   ggsave(p,filename = paste0("result/Imagedimplot_ROI/",i,"_alltype_spatial.png"),width = 15 ,height = 15,dpi = 300)
#   ggsave(p,filename = paste0("result/Imagedimplot_ROI/",i,"_alltype_spatial.pdf"),width = 15 ,height = 15)
# }

xenium_ROI$sample%>% unique()
xenium_ROI$sample <- gsub("_","-",xenium_ROI$sample)

table(xenium_ROI$ROI,xenium_ROI$sample)

for( i in ROI_label){
  
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
                    cols = all_colors10
  )

  ggsave(p,filename = paste0("result/Imagedimplot_ROI/",r,"_spatial.png"),width = 15 ,height = 15,dpi = 300)
  ggsave(p,filename = paste0("result/Imagedimplot_ROI/",r,"_spatial.pdf"),width = 15 ,height = 15)
}}


##sample中映射ROI内细胞-----
xenium@meta.data %>% head()
xenium@meta.data$ROI %>% table()
ROI_names <- names(dataframes)
xenium@meta.data <- xenium@meta.data %>% 
  mutate(
      type3_in_ROI = case_when( 
      ROI %in% ROI_names ~ annotation,  
      TRUE ~ "other"                   
    )
  )
table(xenium@meta.data$type3_in_ROI, useNA = "ifany")

all_colors10 <- c(
  "B_cell" = "gray90", 
  "Endothelial" = "#98D277",
  "Epithelial" = "#9B5DE5",
  "Fibroblast" = "#FDBB69",
  "Myeloid" = "gray90", #"#9E8099",
  "Pericyte" = "gray90",
  "Plasma" = "gray90",
  "SMC" = "gray90", 
  "T_cell" = "gray90",
  "Uncertain" = "gray90",
  "other" = "gray90"
) 

genelist <- c("EFNA1","EPCAM","EPHA3","SMA","COL5A1","PECAM1")
valid_gene_sets <- genelist[genelist %in% rownames(xenium)]
cat("有效基因:", length(valid_gene_sets), "/", length(genelist))
missing_genes <- setdiff(genelist, valid_gene_sets)

Idents(xenium) <- "type3_in_ROI"
for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "black", 
                    border.size = 0.1,
                    cols = all_colors10#,
                    #cells = WhichCells(xenium, idents = c("Endothelial","Epithelial","Fibroblast")),
                    # molecules = genelist,
                    # mols.size = 0.005,
                    # mols.alpha = 0.7,
                    # mols.cols    = c(
                    #   "EFNA1"   = "#FF6B6B",         
                    #   "EPCAM"   = "#FEE440",        
                    #   "EPHA3"   = "#2D82AF",       
                    #   "COL5A1"  = "#02c9c9",   
                    #   "PECAM1"  = "#7209B7"  
                    # )
  )
  
  ggsave(p,filename = paste0("result/Imagedimplot/",i,"_type3_in_ROI_spatial.png"),width = 15 ,height = 15,dpi = 300)
  ggsave(p,filename = paste0("result/Imagedimplot/",i,"_type3_in_ROI_spatial.pdf"),width = 15 ,height = 15)
}


#20260421analysis-----
###ROI中细胞映射-----
all_colors10 <- c(
  "B_cell" = "#00F5A7", 
  "Endothelial" = "#98D277",
  "Epithelial" = "#9B5DE5",
  "Fibroblast" = "#FDBB69",
  "Myeloid" = "#feee7d", #"#9E8099",
  "Pericyte" = "#00BBF9",
  "Plasma" =  "#52AF43",
  "SMC" ="#A6CEE3",
  "T_cell" = "#D9A295",
  "Uncertain" = "gray90"
) 

genelist <- c("EFNA1","EPCAM","EPHA3","SMA","COL5A1","PECAM1")
valid_gene_sets <- genelist[genelist %in% rownames(xenium_ROI)]
cat("有效基因:", length(valid_gene_sets), "/", length(genelist))
missing_genes <- setdiff(genelist, valid_gene_sets)


for(i in sample_label){
  xenium_ROI@images[[i]] <- NULL
}

name_mapping <- c(
  "N16546.3.1" = "N16546-3-1",  
  "N16546.3.2" = "N16546-3-2",  
  "N16546.3.3" = "N16546-3-3",  
  "N16326.G1.1" = "N16326-G1-1", 
  "N16326.G1.2" = "N16326-G1-2", 
  "N16326.G1.3" = "N16326-G1-3", 
  "N18530.G1.1" = "N18530-G1-1", 
  "N18530.G1.3" = "N18530-G1-3", 
  "N18717.G2.1" = "N18717-G2-1", 
  "N18717.G2.2" = "N18717-G2-2", 
  "N18717.G2.3" = "N18717-G2-3"
)

names(xenium_ROI@images) <- name_mapping[names(xenium_ROI@images)]
names(xenium_ROI@images)
ROI_label <- names(xenium_ROI@images)

Idents(xenium_ROI) <- "annotation"

xenium_ROI$sample%>% unique()
xenium_ROI$sample <- gsub("_","-",xenium_ROI$sample)

table(xenium_ROI$ROI,xenium_ROI$sample)

all_colors10 <- c(
  "B_cell" = "gray90", 
  "Endothelial" = "#fdcbc0",
  "Epithelial" = "#2694ab",
  "Fibroblast" = "#96ceb4",
  "Myeloid" = "gray90", #"#9E8099",
  "Pericyte" = "gray90",
  "Plasma" = "gray90",
  "SMC" = "gray90", 
  "T_cell" = "gray90",
  "Uncertain" = "gray90"
) 

for( i in ROI_label){
  
  dat1 <- subset(xenium_ROI,sample == i)
  
  for (r in unique(dat1$ROI)) {
    
    dat2 <- subset(dat1,ROI == r)
    
    #dat2@images$`N16326-G1-1`@molecules$molecules
    
    DefaultBoundary(dat2@images[[i]]) <- "segmentation"
    p <- ImageDimPlot(dat2,
                      fov = i,
                      nmols = 20000,
                      dark.background = F,
                      border.color = "black", 
                      border.size = 0.1,
                      cols = all_colors10,
                      molecules = genelist,
                      mols.size = 0.005,
                      mols.alpha = 0.8,
                      mols.cols    = c(
                        "EFNA1"   = "#FF9F1C",         
                        "EPCAM"   = "#FEE440",        
                        "EPHA3"   = "#FF6B6B",       
                        "COL5A1"  = "#5B9EC9",   
                        "PECAM1"  = "#9B5DE5")
    )
    
    ggsave(p,filename = paste0("result/Imagedimplot_ROI_0421/",r,"_spatial.png"),width = 15 ,height = 15,dpi = 300)
    ggsave(p,filename = paste0("result/Imagedimplot_ROI_0421/",r,"_spatial.pdf"),width = 15 ,height = 15)
  }}



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


##xenium_ROI中根据圈选区域加入fov------

####add FOV----
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


xenium_ROI@meta.data$sample %>% table()

#all_data2 <-  JoinLayers(xenium_ROI,layers = 'counts')

all_data2 <- xenium_ROI

all_data2@meta.data %>% head()

T1<- data_frame(rownames(all_data2@meta.data),all_data2@meta.data[c("orig.ident", "cell_id")])
# ###改行名
# all_data2@meta.data$cell_id <- gsub("^.+?_(?=[a-z]{8}-\\d+$)", "", 
#                                     rownames(all_data2@meta.data), 
#                                     perl = TRUE)

all_data2$ROI %>% table()
all_data2$ROI <- gsub("_","-",all_data2$ROI )
all_data2@meta.data$cell_id_sample <- rownames(all_data2@meta.data)
rownames(all_data2@meta.data) <- paste0(all_data2@meta.data$ROI,"_",all_data2@meta.data$cell_id_raw)
all_data2@meta.data$cell_id <- all_data2@meta.data$cell_id_raw

all_data2$sample_raw <- all_data2$sample
all_data2$sample <- all_data2$ROI

T1<- data_frame(rownames(all_data2@meta.data),all_data2@meta.data[c("sample_raw","sample", "cell_id")])

rownames(all_data2@assays$Xenium@layers$counts) <- rownames(all_data2)
colnames(all_data2@assays$Xenium@layers$counts) <- colnames(all_data2)
xenium_ROI_FOV <- CreateSeuratObject(
  all_data2@assays$Xenium@layers$counts,
  assay = "Xenium")

xenium_ROI_FOV@meta.data <- all_data2@meta.data

umap <- all_data2@reductions$umap@cell.embeddings
rownames(umap) <- rownames(xenium_ROI_FOV@meta.data)
metadata <- xenium_ROI_FOV@meta.data
xenium_ROI_FOV <- Add_meta_umapData_Xenium(object = xenium_ROI_FOV,umapfile = umap,metedata = metadata)
colnames(xenium_ROI_FOV@assays$Xenium) <-  colnames(xenium_ROI_FOV)
row.names(xenium_ROI_FOV@reductions$umap@cell.embeddings) <- colnames(xenium_ROI_FOV)

xenium_ROI_FOV <- RenameCells(xenium_ROI_FOV,new.names = colnames(xenium_ROI_FOV))
label <- unique(xenium_ROI_FOV@meta.data$sample)

options(future.globals.maxSize = 20 * 1024^3)
xenium_ROI_FOV <- NormalizeData(xenium_ROI_FOV)

####多个文件夹，一个文件夹内包含多个fov###----

#####加入1个测试----
library(arrow)
# coord1 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0097339__XENIUM3__20260306__090428",assay = 'Xenium',
#                            mols.qv.threshold = 20,label = "XENIUM3_1",object = xenium_ROI_FOV,cov = NULL)
# xenium_ROI_FOV@images[["XENIUM3_1"]] = coord1

DefaultBoundary(xenium_ROI_FOV@images[["XENIUM3-1"]]) <- "segmentation"
p <- ImageFeaturePlot(xenium_ROI_FOV, fov = "XENIUM3-1",features = c("CD3E"),
                      max.cutoff = "q90",axes = FALSE,
                      border.color = "white", 
                      border.size = 0.1,
                      cols = c('#003366', '#FFD700'))
ggsave(p,filename = paste0("result/","_feature_spatial.png"),width = 12 ,height = 12,dpi = 300)


#####加入圈选区域FOV(正式)----
library(arrow)
coord1 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0097339__XENIUM3__20260306__090428",
                           assay = 'Xenium',
                           mols.qv.threshold = 20,
                           label = "XENIUM3-1",
                           object = xenium_ROI_FOV,
                           cov = NULL)
xenium_ROI_FOV@images[["XENIUM3-1"]] = coord1

#####test===========
# data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0097339__XENIUM3__20260306__090428"
# assay = 'Xenium'
# mols.qv.threshold = 20
# label = "XENIUM3-1"
# object = xenium_ROI_FOV
# cov = NULL
# 
# cells_data <- data.frame(cell_id = paste0(paste(label,collapse = ":"),"_",object@meta.data$cell_id),sample = object@meta.data$sample)
# 
# cells_data_sub <- subset(cells_data,sample %in% label)
# 
# cells <- as.character(cells_data_sub$cell_id)
# features <- rownames(object)
# 
# has_dt <- requireNamespace("data.table", quietly = TRUE) && 
#   requireNamespace("R.utils", quietly = TRUE)
# 
# centroids <- as.data.frame(data.table::fread(file.path(data.dir, 
#                                                        "cells.csv.gz")))
# centroids <- data.frame(x = centroids$x_centroid, 
#                         y = centroids$y_centroid, cell = paste0(paste(label,collapse = ":"),"_",centroids$cell_id), 
#                         stringsAsFactors = FALSE)
# 
# 
# segmentations <- as.data.frame(data.table::fread(file.path(data.dir, 
#                                                            "cell_boundaries.csv.gz")))
# 
# segmentations <- data.frame(x = segmentations$vertex_x,y = segmentations$vertex_y,
#                             cell = paste0(paste(label,collapse = ":"),"_",segmentations$cell_id))
# 
# microns <- read_parquet(file.path(data.dir, "transcripts.parquet"))
# microns <- subset(microns, qv >= mols.qv.threshold)
# microns <- data.frame(x = microns$x_location, y = microns$y_location,
#                       cell = paste0(paste(label,collapse = ":"),"_",microns$cell_id),
#                       gene = microns$feature_name, stringsAsFactors = FALSE)
# microns %>% head()
# 
# segmentations <- subset(segmentations,cell %in% cells)
# centroids <- subset(centroids,cell %in% cells)
# microns <- subset(microns,gene %in% features & cell %in% cells)
# 
# segmentations.data <- list(centroids = CreateCentroids(centroids), 
#                            segmentation = CreateSegmentation(segmentations))
# 
# coords <- CreateFOV(coords = segmentations.data, type = c("segmentation", 
#                                                           "centroids"), molecules = microns, assay = assay)
# 
# coord1 <- coords
# 
# 



coord2 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0097339__XENIUM3__20260306__090428",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "XENIUM3-2",object = xenium_ROI_FOV,cov = NULL)
xenium_ROI_FOV@images[["XENIUM3-2"]] = coord2

coord3 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0097339__XENIUM3__20260306__090428",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "XENIUM3-3",object = xenium_ROI_FOV,cov = NULL)
xenium_ROI_FOV@images[["XENIUM3-3"]] = coord3



# 1
coord4 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "XENIUM6-2",object = xenium_ROI_FOV,cov = NULL)
xenium_ROI_FOV@images[["XENIUM6-2"]] = coord4

# 2
coord5 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "XENIUM6-1",object = xenium_ROI_FOV,cov = NULL)
xenium_ROI_FOV@images[["XENIUM6-1"]] = coord5

# 3
coord6 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "XENIUM6-11",object = xenium_ROI_FOV,cov = NULL)
xenium_ROI_FOV@images[["XENIUM6-11"]] = coord6

# 4
coord7 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "XENIUM6-7",object = xenium_ROI_FOV,cov = NULL)
xenium_ROI_FOV@images[["XENIUM6-7"]] = coord7

# 5
coord8 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "XENIUM6-3",object = xenium_ROI_FOV,cov = NULL)
xenium_ROI_FOV@images[["XENIUM6-3"]] = coord8

# 6
coord9 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "XENIUM6-6",object = xenium_ROI_FOV,cov = NULL)
xenium_ROI_FOV@images[["XENIUM6-6"]] = coord9

# 7
coord10 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "XENIUM6-5",object = xenium_ROI_FOV,cov = NULL)
xenium_ROI_FOV@images[["XENIUM6-5"]] = coord10

# 8
coord11 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "XENIUM6-4",object = xenium_ROI_FOV,cov = NULL)
xenium_ROI_FOV@images[["XENIUM6-4"]] = coord11

# 9
coord12 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "XENIUM6-10",object = xenium_ROI_FOV,cov = NULL)
xenium_ROI_FOV@images[["XENIUM6-10"]] = coord12

# 10
coord13 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "XENIUM6-8",object = xenium_ROI_FOV,cov = NULL)
xenium_ROI_FOV@images[["XENIUM6-8"]] = coord13

# 11
coord14 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "XENIUM6-13",object = xenium_ROI_FOV,cov = NULL)
xenium_ROI_FOV@images[["XENIUM6-13"]] = coord14

# 12
coord15 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "XENIUM6-17",object = xenium_ROI_FOV,cov = NULL)
xenium_ROI_FOV@images[["XENIUM6-17"]] = coord15

# 13
coord16 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "XENIUM6-9",object = xenium_ROI_FOV,cov = NULL)
xenium_ROI_FOV@images[["XENIUM6-9"]] = coord16

# 14
coord17 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "XENIUM6-14",object = xenium_ROI_FOV,cov = NULL)
xenium_ROI_FOV@images[["XENIUM6-14"]] = coord17

# 15
coord18 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "XENIUM6-15",object = xenium_ROI_FOV,cov = NULL)
xenium_ROI_FOV@images[["XENIUM6-15"]] = coord18

# 16
coord19 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "XENIUM6-12",object = xenium_ROI_FOV,cov = NULL)
xenium_ROI_FOV@images[["XENIUM6-12"]] = coord19

# 17
coord20 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2601290804-01/AE2601290804-01/XE2601290804-02/output-XETG00305__0065784__XENIUM-6__20260202__100927",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "XENIUM6-16",object = xenium_ROI_FOV,cov = NULL)
xenium_ROI_FOV@images[["XENIUM6-16"]] = coord20



#####绘图一个样本测试----
DefaultBoundary(xenium_ROI_FOV@images[["XENIUM6-16"]]) <- "segmentation"
p <- ImageFeaturePlot(
  object = xenium_ROI_FOV,
  fov = "XENIUM6-16",
  nmols = 20000,       # 每个图像中显示的最大分子数（随机子抽样）
  border.color = "white",  # 细胞边界的颜色（深蓝色）
  border.size = 0.1,    # 细胞边界的线宽
  features = 'CD3E',
  max.cutoff = "q90",
  axes = FALSE,
  cols = c('#003366', "#E7475E"),
  size = 1.5,
  dark.background = F
)

ggsave(p,filename = paste0("result/","_feature_spatial.png"),width = 12 ,height = 12,dpi = 300)


###ROI中细胞映射-----
genelist <- c("EFNA1",
              #"EPCAM",
              "EPHA3"#,"SMA","COL5A1","PECAM1"
              )
valid_gene_sets <- genelist[genelist %in% rownames(xenium_ROI_FOV)]
cat("有效基因:", length(valid_gene_sets), "/", length(genelist))
missing_genes <- setdiff(genelist, valid_gene_sets)

names(xenium_ROI_FOV@images)
ROI_label <- names(xenium_ROI_FOV@images)

Idents(xenium_ROI_FOV) <- "annotation"

xenium_ROI_FOV$sample%>% unique()

table(xenium_ROI_FOV$ROI,xenium_ROI_FOV$sample)

all_colors10 <- c(
  "B_cell" = "gray90", 
  "Endothelial" = "#fdcbc0",
  "Epithelial" = "#2694ab",
  "Fibroblast" = "#96ceb4",
  "Myeloid" = "gray90", #"#9E8099",
  "Pericyte" = "gray90",
  "Plasma" = "gray90",
  "SMC" = "gray90", 
  "T_cell" = "gray90",
  "Uncertain" = "gray90"
) 

for( i in ROI_label){
    DefaultBoundary(xenium_ROI_FOV@images[[i]]) <- "segmentation"
    p <- ImageDimPlot(xenium_ROI_FOV,
                      fov = i,
                      nmols = 20000,
                      dark.background = F,
                      border.color = "black", 
                      border.size = 0.1,
                      cols = all_colors10,
                      molecules = genelist,
                      mols.size = 1.5,
                      mols.alpha = 0.9,
                      mols.cols    = c(
                        "EFNA1"   = "red",         
                        #"EPCAM"   = "#FE982C",        
                        "EPHA3"   = "#7209B7"#,       
                        #"COL5A1"  = "#5B9EC9",   
                        #"PECAM1"  = "#9B5DE5"
                        )
    )
    p
    ggsave(p,filename = paste0("result/Imagedimplot_ROI_0423_size1.5/",i,"_spatial.png"),width = 15 ,height = 15,dpi = 300)
    ggsave(p,filename = paste0("result/Imagedimplot_ROI_0423_size1.5/",i,"_spatial.pdf"),width = 15 ,height = 15)
  }


for( i in ROI_label){
  DefaultBoundary(xenium_ROI_FOV@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium_ROI_FOV,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "black", 
                    border.size = 0.1,
                    cols = all_colors10,
                    molecules = genelist,
                    mols.size = 1.2,
                    mols.alpha = 0.9,
                    mols.cols    = c(
                      "EFNA1"   = "red",         
                      #"EPCAM"   = "#FE982C",        
                      "EPHA3"   = "#7209B7"#,       
                      #"COL5A1"  = "#5B9EC9",   
                      #"PECAM1"  = "#9B5DE5"
                    )
  )
  p
  ggsave(p,filename = paste0("result/Imagedimplot_ROI_0423_size1.2/",i,"_spatial.png"),width = 15 ,height = 15,dpi = 300)
  ggsave(p,filename = paste0("result/Imagedimplot_ROI_0423_size1.2/",i,"_spatial.pdf"),width = 15 ,height = 15)
}


#saveRDS(xenium_ROI_FOV,"result/xenium_ROI_FOV.rds")

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


#只输出13号ROI-------
all_colors10 <- c(
  "B_cell" = "gray90", 
  "Endothelial" = "#fdcbc0",
  "Epithelial" = "#2694ab",
  "Fibroblast" = "#96ceb4",
  "Myeloid" = "gray90", #"#9E8099",
  "Pericyte" = "gray90",
  "Plasma" = "gray90",
  "SMC" = "gray90", 
  "T_cell" = "gray90",
  "Uncertain" = "gray90"
) 


i <- "XENIUM6-13"

##ImageDimPlot------
DefaultBoundary(xenium_ROI_FOV@images[[i]]) <- "segmentation"
p <- ImageDimPlot(xenium_ROI_FOV,
                  fov = i,
                  nmols = 20000,
                  dark.background = F,
                  border.color = "black", 
                  border.size = 0.1,
                  cols = all_colors10#,
                  # molecules = genelist,
                  # mols.size = 1.5,
                  # mols.alpha = 0.9,
                  # mols.cols    = c(
                  #   "EFNA1"   = "red",         
                  #   #"EPCAM"   = "#FE982C",        
                  #   "EPHA3"   = "#7209B7"#,       
                  #   #"COL5A1"  = "#5B9EC9",   
                  #   #"PECAM1"  = "#9B5DE5"
                  # )
)
p
ggsave(p,filename = paste0("result/Imagedimplot_ROI_0423_xenium6_13/",i,"_spatial1.png"),width = 15 ,height = 15,dpi = 300)
ggsave(p,filename = paste0("result/Imagedimplot_ROI_0423_xenium6_13/",i,"_spatial1.pdf"),width = 15 ,height = 15)


DefaultBoundary(xenium_ROI_FOV@images[[i]]) <- "segmentation"
p <- ImageDimPlot(xenium_ROI_FOV,
                  fov = i,
                  nmols = 20000,
                  dark.background = F,
                  border.color = "black", 
                  border.size = 0.1,
                  cols = all_colors10,
                  molecules = genelist,
                  mols.size = 1.5,
                  mols.alpha = 0.9,
                  mols.cols    = c(
                    "EFNA1"   = "red",         
                    #"EPCAM"   = "#FE982C",        
                    "EPHA3"   = "#7209B7"#,       
                    #"COL5A1"  = "#5B9EC9",   
                    #"PECAM1"  = "#9B5DE5"
                  )
)
p
ggsave(p,filename = paste0("result/Imagedimplot_ROI_0423_xenium6_13/",i,"_spatial2.png"),width = 15 ,height = 15,dpi = 300)
ggsave(p,filename = paste0("result/Imagedimplot_ROI_0423_xenium6_13/",i,"_spatial2.pdf"),width = 15 ,height = 15)


##ImageFeaturePlot-----
DefaultBoundary(xenium_ROI_FOV@images[[i]]) <- "segmentation"
p <- ImageFeaturePlot(
  object = xenium_ROI_FOV,
  fov = i,
  nmols = 20000,       # 每个图像中显示的最大分子数（随机子抽样）
  border.color = "black",  # 细胞边界的颜色（深蓝色）
  border.size = 0.1,    # 细胞边界的线宽
  features = 'EFNA1',
  max.cutoff = "q90",
  axes = FALSE,
  cols = c("#5B9EC9", '#FFD700'),
  size = 1.5,
  dark.background = F
)

p
ggsave(p,filename = paste0("result/Imagedimplot_ROI_0423_xenium6_13/",i,"_spatial_EFNA1.png"),width = 15 ,height = 15,dpi = 300)
ggsave(p,filename = paste0("result/Imagedimplot_ROI_0423_xenium6_13/",i,"_spatial_EFNA1.pdf"),width = 15 ,height = 15)



DefaultBoundary(xenium_ROI_FOV@images[[i]]) <- "segmentation"
p <- ImageFeaturePlot(
  object = xenium_ROI_FOV,
  fov = i,
  nmols = 20000,       # 每个图像中显示的最大分子数（随机子抽样）
  border.color = "black",  # 细胞边界的颜色（深蓝色）
  border.size = 0.1,    # 细胞边界的线宽
  features = 'EPHA3',
  max.cutoff = "q90",
  axes = FALSE,
  cols = c("#5B9EC9", '#FFD700'),
  size = 1.5,
  dark.background = F
)

p
ggsave(p,filename = paste0("result/Imagedimplot_ROI_0423_xenium6_13/",i,"_spatial_EPHA3.png"),width = 15 ,height = 15,dpi = 300)
ggsave(p,filename = paste0("result/Imagedimplot_ROI_0423_xenium6_13/",i,"_spatial_EPHA3.pdf"),width = 15 ,height = 15)


## 双基因共表达分析（只显示基因点）------
# DefaultBoundary(xenium_ROI_FOV@images[[i]]) <- "segmentation"
# p <- ImageFeaturePlot(
#   object = xenium_ROI_FOV,
#   fov = i,
#   nmols = 20000,
#   border.color = "black",
#   border.size = 0.1,
#   blend = TRUE, #双基因共表达分析
#   features = c('EFNA1', 'EPHA3'),
#   max.cutoff = "q90",
#   axes = FALSE,
#   cols = #c("lightgrey", "firebrick1"),
#   c("#5B9EC9", "#00ff00", "#ff0000"),
#   size = 0.5,
#   dark.background = F,
#   boundaries = "segmentation",
#   overlap = FALSE  # 防止边界被分子覆盖或错位
# )
# 
# p
# ggsave(p,filename = paste0("result/Imagedimplot_ROI_0423_xenium6_13/",i,"_spatial3.png"),width = 15 ,height = 15,dpi = 300)
# ggsave(p,filename = paste0("result/Imagedimplot_ROI_0423_xenium6_13/",i,"_spatial3.pdf"),width = 15 ,height = 15)
# 

##共表达作图-------
#help_function
mergeColor <- function(..., color = NULL, base_color = "white", use_alpha = FALSE) {
  
  if ('list' %in% class(...elt(1))) {
    dat_ls <- ...elt(1)
  } else {
    dat_ls <- as.list(...)
  }
  
  # length of dat
  n_gene <- length(dat_ls)
  if (n_gene > 5) {
    stop("Up to 5 genes are supported.")
  }
  
  len_ls <- dat_ls %>% map_dbl(length)
  if (!all(len_ls[1] == len_ls)) {
    stop("The lengths of all the data must be the same")
  }
  
  # color
  if (is.null(color)) {
    color <- c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0")
    color <- color[1:n_gene]
  } else {
    stopifnot(length(color) == n_gene)
  }
  color <- map(color, ~ c(base_color, .x))
  
  # set color
  color_ls <-
    map2(dat_ls, color, function(dat, col) {
      col_numeric(col, domain = range(dat))(dat)
    })
  
  # weight
  weight_df <-
    dat_ls %>%
    map_dfc(rescale, to = c(1e-10, 1))
  
  weight_df_scaled <-
    weight_df %>%
    sweep(
      MARGIN = 1,
      STATS = rowSums(.),
      FUN = "/"
    )
  
  # merge color
  color_mat <- matrix(0, nrow = nrow(weight_df_scaled), ncol = 3)
  for (i in seq_along(color_ls)) {
    color_used <- color_ls[[i]]
    color_used_weighted <- t(col2rgb(color_used)) * weight_df_scaled[[i]]
    color_mat <- color_mat + color_used_weighted
  }
  
  color_merged <- rgb(color_mat, maxColorValue = 255)
  
  # alpha
  if (use_alpha) {
    color_merged <- alpha(color_merged, rowSums(weight_df))
  }
  
  # return
  list(
    mergeColor   = color_merged,
    colorPalette = color_ls,
    Palette      = color
  )
}

###按照ROI划分cells_boundary_ls-----
cells_boundary_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_boundary_ls.rds")

xenium_ROI_FOV$sample <- as.character(xenium_ROI_FOV$sample)

sample_vec <- intersect(unique(na.omit(xenium_ROI_FOV$sample)), names(xenium_ROI_FOV@images))

cells_boundary_list <- list()

for (i in sample_vec) {
  
  img <- xenium_ROI_FOV@images[[i]]
  
  if (!"segmentation" %in% names(img@boundaries)) {
    stop("Cannot find 'segmentation' boundary in image: ", i)
  }
  
  DefaultBoundary(img) <- "segmentation"
  
  cells_boundary <- GetTissueCoordinates(
    img@boundaries[["segmentation"]]
  ) %>% 
    as_tibble()
  
  if (!"ids" %in% colnames(cells_boundary)) {
    id_col <- intersect(c("cell_id", "cell", "Cell", "barcode", "name"), colnames(cells_boundary))
    
    if (length(id_col) == 0) {
      stop("Cannot find cell id column in segmentation coordinates: ", i)
    }
    
    cells_boundary <- cells_boundary %>%
      dplyr::rename(ids = all_of(id_col[1]))
  }
  
  if (!all(c("vertex_x", "vertex_y") %in% colnames(cells_boundary))) {
    x_col <- intersect(c("vertex_x", "x", "X", "imagecol"), colnames(cells_boundary))
    y_col <- intersect(c("vertex_y", "y", "Y", "imagerow"), colnames(cells_boundary))
    
    if (length(x_col) == 0 || length(y_col) == 0) {
      stop("Cannot find x/y columns in segmentation coordinates: ", i)
    }
    
    cells_boundary <- cells_boundary %>%
      dplyr::rename(
        vertex_x = all_of(x_col[1]),
        vertex_y = all_of(y_col[1])
      )
  }
  
  cells_boundary_list[[i]] <- cells_boundary
}

saveRDS(cells_boundary_list,"R_OUTPUT/data_for_analysis/cells_boundary_ls_ROI.rds")


###两基因共表达作图-------
library(parallel)
library(ggnewscale)
library(scales)

gene  <- c('EFNA1', 'EPHA3')

data_exp <- FetchData(xenium_ROI_FOV,vars = c(gene,'sample',"ROI"))

outdir <- "result/Imagedimplot_ROI_0423_xenium6_13/"
dir.create(outdir)

rownames(data_exp) %>% head()

table(xenium_ROI_FOV$sample_raw,xenium_ROI_FOV$sample)


samples <- unique(data_exp$sample)

i <- "XENIUM6-13"
plot_one_sample <- function(i) { browser
  tmp_boundary <- cells_boundary_list[[i]]
  tmp_exp <- data_exp %>% filter(sample == i) %>% rownames_to_column("cell_id")
  dat_merge_color1 <- mergeColor(tmp_exp[, gene], color = c("red", "#7209B7"), base_color = "gray90")
  tmp_exp$merge_col <- dat_merge_color1$mergeColor
  plot_data <- tmp_exp %>% left_join(tmp_boundary,by = c("cell_id"= "ids"))
  
  p <- ggplot(plot_data, aes(x = vertex_y,           
                             y = -vertex_x,          
                             )) +
    geom_polygon(aes(group = cell_id, fill = merge_col#,color = subtype
                     ),
                     linewidth = 0.1, show.legend = FALSE) +
    scale_fill_identity() +
    #geom_point(aes(color = subtype), size = 3, alpha = 0, show.legend = TRUE) +
    scale_color_manual(values = colors, name = "Cell Type",
                       guide = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    coord_equal() + scale_y_reverse() + theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9))
  
  for (genes in seq_along(gene)) {
    g <- gene[genes]
    pal <- dat_merge_color1$Palette[[genes]]
    p <- p + ggnewscale::new_scale_fill() +
      geom_polygon(data = plot_data,
                   aes(x = vertex_y,           
                       y = -vertex_x, group = cell_id, fill = .data[[g]]),
                   linewidth = 0, color = NA, alpha = 0) +
      scale_fill_gradientn(colours = pal, name = g)#+
      #coord_flip()
      
  }
  
  ggsave(filename = file.path(outdir, paste0(i, "_polygon.pdf")),
         plot = p, width = 15, height = 15)
  ggsave(filename = file.path(outdir, paste0(i, "_polygon.png")),
         plot = p, width = 15, height = 15, dpi = 300)
  
  return(i)
}
# 并行执行
mclapply(samples, plot_one_sample, mc.cores = 8)




###两基因共表达加细胞颜色边界作图-------
gene  <- c('EFNA1', 'EPHA3')

data_exp <- FetchData(xenium_ROI_FOV,vars = c(gene,'sample',"ROI","annotation"))

outdir <- "result/Imagedimplot_ROI_0423_xenium6_13/"
dir.create(outdir)

#colors <- setNames(object = color3,nm = unique(data_exp$subtype))
rownames(data_exp) %>% head()

table(xenium_ROI_FOV$sample_raw,xenium_ROI_FOV$sample)


samples <- unique(data_exp$sample)

i <- "XENIUM6-13"
plot_one_sample <- function(i) { browser
  tmp_boundary <- cells_boundary_list[[i]]
  tmp_exp <- data_exp %>% filter(sample == i) %>% rownames_to_column("cell_id")
  dat_merge_color1 <- mergeColor(tmp_exp[, gene], color = c("red", "#7209B7"), base_color = "gray90")
  tmp_exp$merge_col <- dat_merge_color1$mergeColor
  plot_data <- tmp_exp %>% left_join(tmp_boundary,by = c("cell_id"= "ids"))
  
  p <- ggplot(plot_data, aes(x = vertex_y,           
                             y = -vertex_x,          
  )) +
    geom_polygon(aes(group = cell_id, fill = merge_col,color = annotation
    ),
    linewidth = 0.5, show.legend = FALSE) +
    scale_fill_identity() +
    geom_point(aes(color = annotation), size = 3, alpha = 0, show.legend = TRUE) +
    scale_color_manual(values = all_colors10, name = "Celltype",
                       guide = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    coord_equal() + scale_y_reverse() + theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9))
  
  for (genes in seq_along(gene)) {
    g <- gene[genes]
    pal <- dat_merge_color1$Palette[[genes]]
    p <- p + ggnewscale::new_scale_fill() +
      geom_polygon(data = plot_data,
                   aes(x = vertex_y,           
                       y = -vertex_x, group = cell_id, fill = .data[[g]]),
                   linewidth = 0, color = NA, alpha = 0) +
      scale_fill_gradientn(colours = pal, name = g)#+
    #coord_flip()
    
  }
  
  ggsave(filename = file.path(outdir, paste0(i, "_anno_polygon.pdf")),
         plot = p, width = 15, height = 15)
  ggsave(filename = file.path(outdir, paste0(i, "_anno_polygon.png")),
         plot = p, width = 15, height = 15, dpi = 300)
  
  return(i)
}
# 并行执行
mclapply(samples, plot_one_sample, mc.cores = 8)


xenium_ROI_FOV <- readRDS("result/xenium_ROI_FOV.rds")
table(xenium_ROI_FOV$sample,xenium_ROI_FOV$sample_raw)

xenium_ROI_FOV$sample_cell_id <- paste0(xenium_ROI_FOV$sample,"_",xenium_ROI_FOV$cell_id)
###添加新分组------
xenium_ROI_FOV$sample %>% table()

xenium_ROI_FOV$group_old <- xenium_ROI_FOV$group
xenium_ROI_FOV@meta.data <- xenium_ROI_FOV@meta.data %>%
  select(-group)

sample_type_df <- data.frame(
  sample_raw = c(
    "N16326-G1-1","N16326-G1-2","N16326-G1-3",
    "N18530-G1-1","N18530-G1-3",
    "N18717-G2-1","N18717-G2-2","N18717-G2-3",
    "N16546-3-1","N16546-3-2","N16546-3-3"
  ),
  group = c(
    "ADC","ADC","ADC",
    "ADC","ADC",
    "SCC","SCC","SCC",
    "SCC","SCC","SCC"
  )
)

xenium_ROI_FOV@meta.data <- xenium_ROI_FOV@meta.data %>%
  left_join(sample_type_df, by = "sample_raw")
rownames(xenium_ROI_FOV@meta.data) <- xenium_ROI_FOV$sample_cell_id

xenium_ROI_FOV$group %>% table()
table(xenium_ROI_FOV$sample_raw, xenium_ROI_FOV$group)
table(xenium_ROI_FOV$sample, xenium_ROI_FOV$group)


