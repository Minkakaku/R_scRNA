
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


all_data2 <-  JoinLayers(all_data_rctd,layers = 'counts')

all_data2$sample %>% unique()
all_data2$sample_old <- all_data2$sample
all_data2$sample <- gsub("_","-",all_data2$sample) 
all_data2$sample %>% unique()

###改行名
# all_data2@meta.data$cell_id <- gsub("^.+?_(?=[a-z]{8}-\\d+$)", "", 
#                                     rownames(all_data2@meta.data), 
#                                     perl = TRUE)
T1<- data_frame(rownames(all_data2@meta.data),all_data2@meta.data[c("orig.ident", "cell_id")])
rownames(all_data2@meta.data) <- paste0(all_data2@meta.data$sample,"_",all_data2@meta.data$cell_id)

rownames(all_data2@assays$RNA@layers$count) <- rownames(all_data2)
colnames(all_data2@assays$RNA@layers$count) <- colnames(all_data2)
xenium <- CreateSeuratObject(
  all_data2@assays$RNA@layers$count,
  assay = "Xenium")

xenium@meta.data <- all_data2@meta.data

umap <- all_data2@reductions$umap@cell.embeddings
rownames(umap) <- rownames(xenium@meta.data)
metadata <- xenium@meta.data
xenium <- Add_meta_umapData_Xenium(object = xenium,umapfile = umap,metedata = metadata)
colnames(xenium@assays$Xenium) <-  colnames(xenium)
row.names(xenium@reductions$umap@cell.embeddings) <- colnames(xenium)

xenium <- RenameCells(xenium,new.names = colnames(xenium))
label <- unique(xenium@meta.data$sample)

options(future.globals.maxSize = 11 * 1024^3)
xenium <- NormalizeData(xenium)

Idents(xenium) <- "first_type"
#多个文件夹，一个文件夹内包含多个fov###----

##加入1个测试----
library(arrow)
coord_01 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0029706__D24-45111-1__20250402__040839",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "low-1-1",object = xenium,cov = NULL)
xenium@images[["low-1-1"]] = coord_01

DefaultBoundary(xenium@images[["low-1-1"]]) <- "segmentation"
p <- ImageFeaturePlot(xenium, fov = "low-1-1",features = c("CD3E"),
                      max.cutoff = "q90",axes = FALSE,
                      border.color = "white", 
                      border.size = 0.1,
                      cols = c('#003366', '#FFD700'))
ggsave(p,filename = paste0("result/","feature_spatial.png"),width = 12 ,height = 12,dpi = 300)


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

xenium$first_type <- as.character(xenium$first_type)
xenium$first_type[is.na(xenium$first_type)] <- "others"
xenium$first_type <- factor(xenium$first_type)

DefaultBoundary(xenium@images[["low-1-1"]]) <- "segmentation"
p <- ImageDimPlot(xenium,
                  fov = "low-1-1",
                  nmols = 20000,
                  dark.background = F,
                  border.color = "black", 
                  border.size = 0.1,
                  cols = all_colors,
                  #shuffle.cols = F, #随机打乱颜色顺序
                  #size = 1.2,
                  group.by = "first_type"
                  #molecules = c("INHBA"),
                  #mols.size = 1,
                  #mols.cols = "red",
                  #mols.alpha = 0.5, #分子的透明度值
                  
)
ggsave(p,filename = paste0("result/Imagedimplot/","_anno_spatial.png"),width = 12 ,height = 12,dpi = 300)

xenium@meta.data$first_type %>% unique()
table(xenium$sample,xenium$first_type)
xenium$sample %>% table()


##加入12个样本FOV(正式)----
# 2
coord_02 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0029706__D24-45111-1__20250402__040839",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "low-1-2",object = xenium,cov = NULL)
xenium@images[["low-1-2"]] = coord_02

# 3
coord_03 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0029706__D24-45111-1__20250402__040839",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "low-1-3",object = xenium,cov = NULL)
xenium@images[["low-1-3"]] = coord_03

# 4
coord_04 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0029708__D24-13981-2__20250422__063812",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "low-2-1",object = xenium,cov = NULL)
xenium@images[["low-2-1"]] = coord_04

# 5
coord_05 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0029708__D24-13981-2__20250422__063812",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "low-2-2",object = xenium,cov = NULL)
xenium@images[["low-2-2"]] = coord_05

# 6
coord_06 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0029708__D24-13981-2__20250422__063812",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "low-2-3",object = xenium,cov = NULL)
xenium@images[["low-2-3"]] = coord_06

# 7
coord_07 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0029708__D24-49895-2__20250422__063812",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "Her2-0-1-1",object = xenium,cov = NULL)
xenium@images[["Her2-0-1-1"]] = coord_07

# 8
coord_08 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0029708__D24-49895-2__20250422__063812",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "Her2-0-1-2",object = xenium,cov = NULL)
xenium@images[["Her2-0-1-2"]] = coord_08

# 9
coord_09 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0029708__D24-49895-2__20250422__063812",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "Her2-0-1-3",object = xenium,cov = NULL)
xenium@images[["Her2-0-1-3"]] = coord_09

# 10
coord_10 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0057846__D24-27540__20250514__092136",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "Her2-0-10",object = xenium,cov = NULL)
xenium@images[["Her2-0-10"]] = coord_10

# 11
coord_11 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0057846__D24-37766__20250514__092136",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "Her2-0-13-1",object = xenium,cov = NULL)
xenium@images[["Her2-0-13-1"]] = coord_11

# 12
coord_12 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0057846__D24-37766__20250514__092136",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "Her2-0-13-2",object = xenium,cov = NULL)
xenium@images[["Her2-0-13-2"]] = coord_12

# 13
coord_13 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0057846__D24-36283__20250514__092136",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "low-12",object = xenium,cov = NULL)
xenium@images[["low-12"]] = coord_13

# 14
coord_14 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0057846__D24-27520__20250514__092136",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "Her2-0-14",object = xenium,cov = NULL)
xenium@images[["Her2-0-14"]] = coord_14

# 15
coord_15 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0057846__D24-25697__20250514__092136",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "Her2-0-15-1",object = xenium,cov = NULL)
xenium@images[["Her2-0-15-1"]] = coord_15

# 16
coord_16 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0057846__D24-25697__20250514__092136",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "Her2-0-15-2",object = xenium,cov = NULL)
xenium@images[["Her2-0-15-2"]] = coord_16

# 17
coord_17 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0057847__D24-24957__20250514__092136",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "Her2-0-16",object = xenium,cov = NULL)
xenium@images[["Her2-0-16"]] = coord_17

# 18
coord_18 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0060608__D25-00191-6__20250620__060536",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "Her2-030",object = xenium,cov = NULL)
xenium@images[["Her2-030"]] = coord_18

# 19
coord_19 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0060626__D24-65708-1__20250620__060537",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "Her2-028",object = xenium,cov = NULL)
xenium@images[["Her2-028"]] = coord_19

# 20
coord_20 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0066062__D24-54475-1__20250703__074541",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "low-26",object = xenium,cov = NULL)
xenium@images[["low-26"]] = coord_20

# 21
coord_21 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0066051__D24-65747-1__20250717__081753",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "Her2-029",object = xenium,cov = NULL)
xenium@images[["Her2-029"]] = coord_21

# 22
coord_22 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0066051__D25-00003-7-2__20250717__081753",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "Her2-026-2",object = xenium,cov = NULL)
xenium@images[["Her2-026-2"]] = coord_22

# 23
coord_23 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0066051__D25-00003-7-1-D25-05658-1__20250717__081753",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "Her2-026-1",object = xenium,cov = NULL)
xenium@images[["Her2-026-1"]] = coord_23

# 24
coord_24 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0066051__D25-00003-7-1-D25-05658-1__20250717__081753",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "low-24",object = xenium,cov = NULL)
xenium@images[["low-24"]] = coord_24

# 25
coord_25 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0066051__D25-09086-1-D25-21921-2__20250717__081753",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "Her2-031",object = xenium,cov = NULL)
xenium@images[["Her2-031"]] = coord_25

# 26
coord_26 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0066051__D25-13349-3-1-D25-13349-3-2__20250717__081753",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "low-25-1",object = xenium,cov = NULL)
xenium@images[["low-25-1"]] = coord_26

# 27
coord_27 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0066051__D25-13349-3-1-D25-13349-3-2__20250717__081753",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "low-25-2",object = xenium,cov = NULL)
xenium@images[["low-25-2"]] = coord_27

# 28
coord_28 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0066051__D23-48937-4-D25-12733-1-2__20250717__081753",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "Her2-027-2",object = xenium,cov = NULL)
xenium@images[["Her2-027-2"]] = coord_28

# 29
coord_29 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0066051__D25-12733-1-1__20250717__081753",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "Her2-027-1",object = xenium,cov = NULL)
xenium@images[["Her2-027-1"]] = coord_29

# 30
coord_30 <- Read_fov_spatial(data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2508274589-02_mateng/rawdata/output-XETG00305__0065812__tma__20251115__033004",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "Her2-034",object = xenium,cov = NULL)
xenium@images[["Her2-034"]] = coord_30




##绘图一个样本测试----
DefaultBoundary(xenium@images[["Her2-034"]]) <- "segmentation"
p <- ImageFeaturePlot(
  object = xenium,
  fov = "Her2-034",
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
ggsave(p,filename = paste0("result/Her2-034","_feature_spatial.png"),width = 12 ,height = 12,dpi = 300)


saveRDS(xenium,"result/xenium.rds")
xenium <- readRDS("result/xenium.rds")

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

all_colors11 <- c(
  "#FF6B6B","#A6CEE3","#00F5A7", "#FDBB69","#feee7d", "#9E8099",
  "#00BBF9","#7209B7","#B294C7", "#52AF43","#2D82AF"
)

all_colors <- c(
  # B细胞&浆细胞 蓝色系
  "B cells"          = "#A6CEE3",
  "IgA Plasma"       = "#5B9EC9",
  "IgG Plasma"       = "#2D82AF",
  "IgM Plasma"       = "#02c9c9",
  
  # T/NK细胞 绿色系
  "Cytotoxic CD8T"   = "#7EBA98",
  "Effector CD8T"    = "#98D277",
  "Exhausted CD4T"   = "#52AF43",
  "Naive CD4T"       = "#6F9E4C",
  "Naive CD8T"       = "#00F5A7",
  "Treg"             = "#82B74B",
  "NK"               = "#78866b",
  
  # 髓系免疫 黄橙棕系
  "Myeloid cells"    = "#FDBB69",
  "Mast cells"       = "#FE982C",
  "Neutrophil"       = "#F78620",
  "ILC"              = "#DD9A88",
  "pDCs"             = "#D9A295",
  "others"           = "gray90",
  
  # 基质内皮肌上皮 紫粉系
  "Fibroblasts"      =   "#B15928", 
  "Endothelial cells"=   "#DBB466", 
  "Myoepithelial cells" ="#9E8099", 

  # 肿瘤亚群 红调渐变
  "Tumor-like C1"    = "#FF6B6B",
  "Tumor-like C2"    = "#F15BB5",
  "Tumor-like C3"    = "#9B5DE5",
  "Tumor-like C4"    = "#600473",
  "Tumor-like C5"    = "#B294C7",
  "Tumor-like C6"    = "#7D54A5"
)


unique(xenium$first_type)

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
                    cols = all_colors,
                    shuffle.cols = F, #随机打乱颜色顺序
                    #size = 1.2,
                    group.by = "first_type"
                    #molecules = c("INHBA"),
                    #mols.size = 1,
                    #mols.cols = "red",
                    #mols.alpha = 0.5, #分子的透明度值
                    
  )
  ggsave(p,filename = paste0("result/Imagedimplot/",i,"_anno_spatial.png"),width = 12 ,height = 12,dpi = 300)
  ggsave(p,filename = paste0("result/Imagedimplot/",i,"_anno_spatial.pdf"),width = 12 ,height = 12)
}

Idents(xenium) <- "first_type"
###Imagedimplot(关注细胞类型)----
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
                    cols = all_colors,
                    #shuffle.cols = F, #随机打乱颜色顺序
                    #size = 1.2,
                    group.by = "first_type",
                    #molecules = c("INHBA"),
                    #mols.size = 1,
                    #mols.cols = "red",
                    #mols.alpha = 0.5, #分子的透明度值
                    cells = WhichCells(xenium, idents = c("IgA Plasma", "IgG Plasma", "IgM Plasma"))
                    
  )
  ggsave(p,filename = paste0("result/Imagedimplot/",i,"_anno_Plasma_spatial.png"),width = 12 ,height = 12,dpi = 300)
  ggsave(p,filename = paste0("result/Imagedimplot/",i,"_anno_Plasma_spatial.pdf"),width = 12 ,height = 12)
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

#修改分组-------
xenium$sample %>% table()
xenium$group %>% table()
xenium$annotation %>% table()
xenium$first_type %>% table()

table(xenium$sample,xenium$group)

xenium$group_raw <- xenium$group
xenium$group_raw %>% table()
##group: 4cyc 8cyc+病人序号-----
sample_to_group <- c(
  "low-1-1" = "4cyc-pCP-1",
  "low-1-2" = "4cyc-pCP-1",
  "low-1-3" = "4cyc-pCP-1",
  "low-2-1" = "8cyc-pCP-1",
  "low-2-2" = "8cyc-pCP-1",
  "low-2-3" = "8cyc-pCP-1",
  "Her2-0-1-1" = "4cyc-pCP-2",
  "Her2-0-1-2" = "4cyc-pCP-2",
  "Her2-0-1-3" = "4cyc-pCP-2",
  "Her2-0-10" = "8cyc-pCP-2",
  "Her2-0-13-1" = "8cyc-pCP-3",
  "Her2-0-13-2" = "8cyc-pCP-3",
  "low-12" = "8cyc-pCP-4",
  "Her2-0-14" = "non-pCR-2",
  "Her2-0-15-1" = "4cyc-pCP-3",
  "Her2-0-15-2" = "4cyc-pCP-3",
  "Her2-0-16" = "non-pCR-2",
  "Her2-030" = "Surgery",
  "Her2-028" = "Surgery",
  "low-26" = "Surgery",
  "Her2-029" = "Surgery",
  "Her2-026-2" = "Surgery",
  "Her2-026-1" = "Surgery",
  "low-24" = "Surgery",
  "Her2-031" = "Surgery",
  "low-25-1" = "Surgery",
  "low-25-2" = "Surgery",
  "Her2-027-2" = "Surgery",
  "Her2-027-1" = "Surgery",
  "Her2-034" = "4cyc-pCP-5"
)

xenium$group <- plyr::mapvalues(
  x = xenium$sample,
  from = names(sample_to_group),
  to = as.character(sample_to_group)
)

table(xenium$group)
table(xenium$sample,xenium$group)

##group1：4cyc 8cyc-----
sample_to_group1 <- c(
  "low-1-1" = "4cyc-pCP",
  "low-1-2" = "4cyc-pCP",
  "low-1-3" = "4cyc-pCP",
  "low-2-1" = "8cyc-pCP",
  "low-2-2" = "8cyc-pCP",
  "low-2-3" = "8cyc-pCP",
  "Her2-0-1-1" = "4cyc-pCP",
  "Her2-0-1-2" = "4cyc-pCP",
  "Her2-0-1-3" = "4cyc-pCP",
  "Her2-0-10" = "8cyc-pCP",
  "Her2-0-13-1" = "8cyc-pCP",
  "Her2-0-13-2" = "8cyc-pCP",
  "low-12" = "8cyc-pCP",
  "Her2-0-14" = "non-pCR",
  "Her2-0-15-1" = "4cyc-pCP",
  "Her2-0-15-2" = "4cyc-pCP",
  "Her2-0-16" = "non-pCR",
  "Her2-030" = "Surgery",
  "Her2-028" = "Surgery",
  "low-26" = "Surgery",
  "Her2-029" = "Surgery",
  "Her2-026-2" = "Surgery",
  "Her2-026-1" = "Surgery",
  "low-24" = "Surgery",
  "Her2-031" = "Surgery",
  "low-25-1" = "Surgery",
  "low-25-2" = "Surgery",
  "Her2-027-2" = "Surgery",
  "Her2-027-1" = "Surgery",
  "Her2-034" = "4cyc-pCP"
)

xenium$group1 <- plyr::mapvalues(
  x = xenium$sample,
  from = names(sample_to_group1),
  to = as.character(sample_to_group1)
)

table(xenium$group1)
table(xenium$sample,xenium$group1)
saveRDS(xenium,"result/xenium_group.rds")

#xenium <- readRDS("result/xenium_group.rds")


data_frame(table(xenium$group,xenium$group1))

#组间相关性比较4cyc vs 8cyc------
plotdata <- xenium[[]] %>% 
  filter(group1%in%c("4cyc-pCP","8cyc-pCP") )%>% 
  filter(!first_type %in% c("others"))

table(plotdata$first_type)
table(plotdata$group,plotdata$sample)

# sample_order <- plotdata %>%
#   select(sample, group) %>%
#   distinct() %>%
#   arrange(group) %>%  # 按 group 排序
#   pull(sample)        # 提取排好的样本顺序

sample_order <- plotdata %>%
  select(sample, group) %>%
  distinct() %>%
  arrange(group)

p_col_proj <-
  plotdata %>% 
  group_by(sample,group, first_type) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = first_type)) + 
  geom_col() +
  #facet_wrap(~group,nrow = 1)+
  #ggsci::scale_fill_igv() +
  scale_fill_manual(values = all_colors) +
  scale_x_discrete(
    limits = sample_order$sample,
    labels = paste0(sample_order$sample,"\n",sample_order$group )
  ) +
  
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

# export
ggsave(
  p_col_proj,
  filename = file.path(result_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 12
)



p_col_group <-
  plotdata %>% 
  group_by(group, first_type) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = first_type)) + 
  geom_col() +
  #ggsci::scale_fill_igv() +
  scale_fill_manual(values = all_colors) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_group,
  filename = file.path(result_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 10
)




library(tidyverse)
library(Seurat)
library(ggalluvial)
library(ggthemes)
library(ggpubr)


dat_props_immune <- xenium[[]] %>% 
  filter(group1%in%c("4cyc-pCP","8cyc-pCP")) %>% 
  group_by(group, first_type) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

dat_props_immune %>% head
# # A tibble: 6 × 4
# # Groups:   group [1]
# group annotation  count   freq
# <fct> <fct>       <int>  <dbl>
#   1 Ctrl  B cell       1437 0.153 
# 2 Ctrl  Plasma cell  1515 0.162 
# 3 Ctrl  T cell       4222 0.451 
# 4 Ctrl  NK cell       152 0.0162
# 5 Ctrl  Macrophage   1192 0.127 
# 6 Ctrl  Neutrophil    253 0.0270


p_props_immune  <-
  dat_props_immune %>%
  ggplot(aes(x = group, y = freq, fill = first_type)) +
  geom_col(width = 0.95) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = all_colors) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(
    x = NULL,
    y = NULL
  )

p_props_immune_flow <- 
  dat_props_immune %>% 
  # mutate(id = 1:n()) %>%
  
  ggplot(aes(x = group, y = freq, stratum = first_type, alluvium = first_type, fill = first_type, label = first_type)) +
  geom_stratum(alpha = 1, color = NA) +
  geom_flow(alpha = 0.5) +
  #geom_text(stat = "stratum") +
  theme_classic2() +
  theme(
    axis.line = element_blank(),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = all_colors) +
  #scale_x_discrete(expand = c(0.01, 0.01)) +
  #scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(
    x = "group",
    y = "cell proportion"
  )

ggsave(file = file.path(result_dir, "props_flow.pdf"), p_props_immune_flow, width = 10,height = 6)
ggsave(file = file.path(result_dir, "props_flow.png"), p_props_immune_flow, width = 10,height = 6)




