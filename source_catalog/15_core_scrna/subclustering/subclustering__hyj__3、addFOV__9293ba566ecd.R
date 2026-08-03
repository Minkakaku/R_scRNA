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

all_data2 <-  JoinLayers(all_data)

###改行名
all_data2@meta.data$cell_id <- gsub("^.+?_(?=[a-z]{8}-\\d+$)", "",
                                 rownames(all_data2@meta.data),
                                 perl = TRUE)
T1<- data_frame(rownames(all_data2@meta.data),all_data2@meta.data[c("orig.ident", "cell_id")])
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
umap <- all_data@reductions$umap@cell.embeddings
rownames(umap) <- rownames(xenium@meta.data)
metadata <- xenium@meta.data
xenium <- Add_meta_umapData_Xenium(object = xenium,umapfile = umap,metedata = metadata)
colnames(xenium@assays$Xenium) <-  colnames(xenium)
row.names(xenium@reductions$umap@cell.embeddings) <- colnames(xenium)
xenium <- RenameCells(xenium,new.names = colnames(xenium))


#所有fov在一个文件夹内-----
#第一张芯片------
has_dt <- requireNamespace("data.table", quietly = TRUE) && requireNamespace("R.utils", quietly = TRUE)
data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/XE2511036121-01/XE2511036121-01/output-XETG00546__0066225__TM1__20260204__084416"
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
samples_TM1 <- rownames(tab)[tab[, "TM1"] > 0]
samples_TMA2 <- rownames(tab)[tab[, "TMA2"] > 0]

print(samples_TM1)
print(samples_TMA2)

#sample_label <- unique(xenium@meta.data$sample)
for (i in samples_TM1) {
  fov <- xenium_fov(object = xenium,label = i,assay = "Xenium")
  xenium@images[[i]] <- fov
}

#第二张芯片--------
has_dt <- requireNamespace("data.table", quietly = TRUE) && requireNamespace("R.utils", quietly = TRUE)
data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/XE2511036121-01/XE2511036121-01/output-XETG00546__0066228__TMA2__20260204__084416"
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
samples_TM1 <- rownames(tab)[tab[, "TM1"] > 0]
samples_TMA2 <- rownames(tab)[tab[, "TMA2"] > 0]

print(samples_TM1)
print(samples_TMA2)

#sample_label <- unique(xenium@meta.data$sample)
for (i in samples_TMA2) {
  fov <- xenium_fov(object = xenium,label = i,assay = "Xenium")
  xenium@images[[i]] <- fov
}


xenium <- NormalizeData(xenium)

saveRDS(xenium,"result/xenium.rds")

xenium <- readRDS("result/xenium.rds")

###测试一个样本----
DefaultBoundary(xenium@images[["P1_POST2"]]) <- "segmentation"
p <- ImageFeaturePlot(
  object = xenium,
  fov = "P1_POST2",
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
ggsave(p,filename = paste0("P1_POST2","_feature_spatial.png"),width = 12 ,height = 12,dpi = 300)


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

#1、明确纤维化时间点：αSMA、是否显著变化，哪个时间点和control比变化明显（结合Masson染色判断）------

xenium$group_raw <- xenium$group

replace_dict <- c(
  "A1" = "H1",
  "A3" = "H12",
  "B1" = "Con",
  "C6" = "D7",
  "D2" = "H3",
  "F1" = "H24",
  "F2" = "H72"
)

xenium$group <- plyr::revalue(xenium$group, replace = replace_dict)
group_order <- c("Con", "H1", "H3", "H12", "H24", "H72", "D7")
xenium$group <- factor(xenium$group, levels = group_order)

table(xenium$group_raw ,xenium$group)

xenium$group %>% table()
xenium$annotation %>% table()

pdf(file.path(result_dir, str_glue("Dotplot group myofib_marker.pdf")), width = 5, height =5)
DotPlot(xenium, group.by = "group", features = c('Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1'))+RotatedAxis()
dev.off()


pdf(file.path(result_dir, str_glue("Vlnplot group myofib_marker111.pdf")), width = 20, height = 10)
VlnPlot(xenium,
        group.by = "group",
        features = c('Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1'),
        pt.size = 0,
        ncol = 4,
)
dev.off()


pdf(file.path(result_dir, str_glue("Vlnplot group myofib_marker pt.pdf")), width = 20, height = 10)
VlnPlot(xenium,group.by = "group",
        features = c('Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1'),
        pt.size = 0.1,
        ncol = 4,
)
dev.off()


##过滤0画小提琴------
Myofibmarker <-c('Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1')

expr_data <- FetchData(
  object = xenium,
  vars = c(Myofibmarker, "sample", "group") 
)

#长格式，适配ggplot2
expr_long <- expr_data %>%
  pivot_longer(
    cols = all_of(Myofibmarker),
    names_to = "gene",         
    values_to = "expression"   
  ) %>%
  filter(expression != 0) %>%
  mutate(
    group = factor(group)
  )

#group_colors <- c("Normal Epi" = "#A23B72", "Tumor Epi" = "#2E86AB")
group_colors <- c("#FF6B6B","#A6CEE3","#00F5A7", "#FDBB69","#feee7d", "#9E8099",
"#00BBF9","#7209B7")


p <- ggplot(expr_long, aes(x = group, y = expression, fill = group)) +
  geom_violin(scale = "width", alpha = 0.8, position = position_dodge(width = 0.8)) +
  #geom_boxplot(width = 0.1, alpha = 1, position = position_dodge(width = 0.8)) +
  facet_wrap(~gene, scales = "free_y", ncol = 4) + 
  scale_fill_manual(values = group_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    panel.grid = element_blank()
  ) +
  labs(
    x = "group",
    y = "Gene Expression",
    fill = "group",
    title = ""
  )


ggsave(
  filename = file.path(result_dir, "Myofibmarker_violin_plot.pdf"),
  plot = p,
  width = 20,  
  height = 10, 
  dpi = 300,   
  device = "pdf"
)

ggsave(
  filename = file.path(result_dir, "Myofibmarker_violin_plot.png"),
  plot = p,
  width = 20,
  height = 10,
  dpi = 300
)



# Fibroblasts --------------------------------------------------------------
all_data$annotation %>% table()
all_data$res_0.8 %>% table()
table(all_data$res_0.8,all_data$annotation)

Fibro_out_dir<-'result/Fibroblast_subtype' 
dir.create(Fibro_out_dir)
Fibroblast<-subset(all_data, subset = annotation %in% c("Fibroblast","Proliferate"))
table(Fibroblast$res_0.8,Fibroblast$annotation)

##修改分组--------
Fibroblast$group_raw <- Fibroblast$group

replace_dict <- c(
  "A1" = "H1",
  "A3" = "H12",
  "B1" = "Con",
  "C6" = "D7",
  "D2" = "H3",
  "F1" = "H24",
  "F2" = "H72"
)

Fibroblast$group <- plyr::revalue(Fibroblast$group, replace = replace_dict)
group_order <- c("Con", "H1", "H3", "H12", "H24", "H72", "D7")
Fibroblast$group <- factor(Fibroblast$group, levels = group_order)

table(Fibroblast$group_raw ,Fibroblast$group)

Fibroblast$group %>% table()
Fibroblast$annotation %>% table()


##修改样本名--------
table(Fibroblast$sample,Fibroblast$group)
Fibroblast$sample %>% table()
Fibroblast$sample_raw <- Fibroblast$sample

replace_rule <- c(
  "A1"  = "H1",
  "A3"  = "H12",
  "B1"  = "Con",
  "C6"  = "D7",
  "D2"  = "H3",
  "F1"  = "H24",
  "F2"  = "H72"
)

Fibroblast$sample <- gsub("A1" , "H1",Fibroblast$sample)  
Fibroblast$sample <- gsub("A3" , "H12",Fibroblast$sample) 
Fibroblast$sample <- gsub("B1" , "Con",Fibroblast$sample) 
Fibroblast$sample <- gsub("C6" , "D7",Fibroblast$sample) 
Fibroblast$sample <- gsub("D2" , "H3",Fibroblast$sample) 
Fibroblast$sample <- gsub("F1" , "H24",Fibroblast$sample) 
Fibroblast$sample <- gsub("F2" , "H72",Fibroblast$sample) 
table(Fibroblast$sample,Fibroblast$group)

sample_order <- Fibroblast[[]] %>%
  arrange(group) %>%
  pull(sample) %>%
  unique()

Fibroblast$sample <- factor(Fibroblast$sample, levels = sample_order)

## PCA --------
DefaultAssay(Fibroblast) <- "RNA"
Fibroblast <- FindVariableFeatures(Fibroblast)
Fibroblast <- ScaleData(Fibroblast, vars.to.regress = 'nCount_RNA')
Fibroblast <- RunPCA(Fibroblast)


### pca inspector ------
pdf(file.path(Fibro_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Fibroblast, ndims = 50)
dev.off()

# pdf(file.path(Fibro_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
# DimHeatmap(Fibroblast, dims = 1:5)
# dev.off()
# pdf(file.path(Fibro_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
# DimHeatmap(Fibroblast, dims = 6:10)
# dev.off()
# pdf(file.path(Fibro_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
# DimHeatmap(Fibroblast, dims = 11:15)
# dev.off()
# pdf(file.path(Fibro_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
# DimHeatmap(Fibroblast, dims = 16:20)
# dev.off()



## dimension reduction --------
dims <- 20

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb

Fibroblast <- RunUMAP(Fibroblast, dims = 1:dims,return.model = TRUE)

pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "sample", label = TRUE) + NoLegend()
dev.off()


pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} group.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "group", label = TRUE)
dev.off()

table(all_data@meta.data$celltype, all_data@meta.data$project)

Fibroblast <- Fibroblast %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

Fibroblast <- RunUMAP(
  Fibroblast, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} group - harmony.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "group", label = TRUE, reduction = "harmonyumap")
dev.off()


reduction.use <- 'harmony' # 'harmony' 'pca'

### cluster -------
Fibroblast  <-
  FindNeighbors(
    Fibroblast ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

Fibroblast <-
  FindClusters(
    Fibroblast,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )


# dimplot
pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} cluster umap.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'umap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

# dimplot
pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} cluster harmonyumap.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

saveRDS(Fibroblast, "result/Fibroblast_subtype/Fibroblast.rds")
#Fibroblast <- readRDS("result/Fibroblast_subtype/Fibroblast.rds")

### dotplot -------
celltype_marker_ls <- c(
  # Epithelial
  'Epithelial' = "Epcam",
  'Epithelial' = "Krt19",
  'Epithelial' = "Cryab",
  
  # POD
  'POD' = "Nphs2",
  'POD' = "Podxl",
  
  # HealthyS1
  'PT' = "Slc34a1",
  'PT' = "Aldob",
  'PT' = "Lrp2",
  
  # LOH
  'LOH' = "Umod",
  'LOH' = "Slc12a1",
  
  # DT
  'DT' = "Slc12a3",
  'DT' = "Calb1",
  
  # PC
  'PC' = "Aqp2",
  'PC' = "Hsd11b2",
  
  # Immune
  'Immune' = "Ptprc",
  'Immune' = "Cd74",
  'Immune' = "Ccl5",
  
  # Fibroblast
  'Fibroblast' = "Dcn",
  'Fibroblast' = "Col1a1",
  
  # SMC
  'SMC' = "Acta2",
  'SMC' = "Tagln",
  
  # Adipocyte
  'Adipocyte' = "Fabp4",
  'Adipocyte' = "Plin4"
)



Fibro_major_marker_ls <- c(
  # S1 fibroblasts
  'S1_fibroblasts' = "Adamdec1", 
  'S1_fibroblasts' = "Cp", 
  'S1_fibroblasts' = "Fabp4",
  
  # S3 fibroblasts
  'S3_fibroblasts' = "Ogn", 
  'S3_fibroblasts' = "Ccdc80", 
  'S3_fibroblasts' = "Grem1",
  
  # Peri-cryptal S2 fibroblasts
  'S2_fibroblasts' = "Sox6",
  'S2_fibroblasts' = "F3",
  
  # S2a fibroblasts
  'S2a_fibroblasts' = "Vstm2",
  
  # S2b fibroblasts
  'S2b_fibroblasts' = "Npy",
  'S2b_fibroblasts' = "Nrg1",
  
  # Ribhi fibroblasts
  'Ribhi_fibroblasts' = "Rpl7a",
  'Ribhi_fibroblasts' = "Rpl28",
  
  # MT fibroblasts
  'MT_fibroblasts' = "Mt-genes",  # 线粒体基因通用标记
  
  # Immediate early-response (IER) fibroblasts
  'IER_fibroblasts' = "Fos", 
  'IER_fibroblasts' = "Fosb", 
  'IER_fibroblasts' = "Irf1",
  
  # Fibroblastic reticular cells (FRCs)
  'FRCs' = "Ccl19", 
  'FRCs' = "Ccl21",
  
  # Smooth muscle cells (SMCs)
  'SMC' = "Tagln",
  'SMC' = "Mylk",
  'SMC' = "Cnn1",
  
  # Fibroblasts (FIB)
  'FIB' = "Lum",
  'FIB' = "Dcn",
  'FIB' = "Col1a1",
  'FIB' = "Col1a2",
  'FIB' = "Fn1",
  
  # Myofibroblasts
  'Myofibroblasts' = "Actg2",
  'Myofibroblasts' = "Sostdc1",
  'Myofibroblasts' = "Myh11",
  
  
  # Inflammatory fibroblasts
  'Inflammatory_fibroblasts' = "Il11",
  'Inflammatory_fibroblasts' = "Cxcl5", 
  'Inflammatory_fibroblasts' = "Fap", 
  'Inflammatory_fibroblasts' = "Inhba",
  'Inflammatory_fibroblasts' = "Il24",
  
  # Endothelial cells
  'Endothelial_cells' = "Vwf", 
  'Endothelial_cells' = "Pecam1", 
  'Endothelial_cells' = "Plvap",
  
  # Activated endothelium
  'Activated_endothelium' = "Ackr1", 
  'Activated_endothelium' = "Sele", 
  'Activated_endothelium' = "Selp",
  
  # Lymphatic endothelium
  'Lymphatic_endothelium' = "Mrc1", 
  'Lymphatic_endothelium' = "Mmrn1", 
  'Lymphatic_endothelium' = "Lyve1",
  
  # Pericytes
  'Pericytes' = "Notch3", 
  'Pericytes' = "Rgs5",
  
  # Glial cells
  'Glial_cells' = "Nrxn1", 
  'Glial_cells' = "S100b", 
  'Glial_cells' = "Cdh19"
)


Fibro_major_marker_ls1 <- c(
  # Smooth muscle cells (SMCs)
  'SMC' = "Tagln",
  'SMC' = "Mylk",
  'SMC' = "Cnn1",#'Acta2'
  
  # Fibroblasts (FIB)
  'FIB' = "Lum",
  'FIB' = "Dcn",
  'FIB' = "Col1a1",
  'FIB' = "Col1a2",
  'FIB' = "Fn1",
  'FIB' = "C7",
  'FIB' = "Ccn2",
  
  # Myofibroblasts
  'Myofibroblasts' = "Actg2",
  'Myofibroblasts' = "Sostdc1",
  'Myofibroblasts' = "Myh11", #'Acta2',
  'Myofibroblasts' = "Notch3"
)


Idents(Fibroblast)<-'RNA_snn_res.0.4'
Fibroblast_markers <-
  FindAllMarkers(
    JoinLayers(Fibroblast),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

Fibroblast_top10_markers <- 
  Fibroblast_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10) %>%
  nest() %>%
  deframe() %>% 
  map(pull, 'gene') %>%
  imap(function(genes, name){
    re <- genes %>%
      map(function(gene){
        # Fibroblastrowser()
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(Fibroblast_top10_markers) <- sub('.*\\.','',names(Fibroblast_top10_markers))

Fibroblast_top10_markers1 <-
  Fibroblast_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10)

# dotplot
pdf(file.path(Fibro_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.4", features = celltype_marker_ls)+RotatedAxis()
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.4", features = Fibro_major_marker_ls)+RotatedAxis()
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.4", features = Fibro_major_marker_ls1)+RotatedAxis()
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.4", features = Fibroblast_top10_markers)+RotatedAxis()
dev.off()

table(Fibroblast$RNA_snn_res.0.4,Fibroblast$group)


pdf(file.path(Fibro_out_dir, str_glue("Dotplot dims{dims} myofib_marker.pdf")), width = 5, height =5)
DotPlot(Fibroblast, 
        group.by = "group", 
        features = c('Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1',"Cnn1","Ccn2"))+RotatedAxis()
dev.off()


pdf(file.path(Fibro_out_dir, str_glue("Vlnplot dims{dims} myofib_marker.pdf")), width = 20, height = 10)
VlnPlot(Fibroblast,group.by = "group",
        feature = c('Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1',"Cnn1","Ccn2"),
        pt.size = 0,
        ncol = 4,
)
dev.off()

pdf(file.path(Fibro_out_dir, str_glue("Vlnplot dims{dims} myofib_marker pt.pdf")), width = 20, height = 10)
VlnPlot(Fibroblast,group.by = "group",
        feature = c('Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1',"Cnn1","Ccn2"),
        pt.size = 0.1,
        ncol = 4,
)
dev.off()

pdf(file.path(Fibro_out_dir, str_glue("Featureplt dims{dims} myofib_marker.pdf")), width = 20, height = 10)
FeaturePlot(Fibroblast,
            #group.by = "group",
            feature = c('Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1',"Cnn1","Ccn2"),
            reduction = 'umap',
            ncol = 4,
)

dev.off()


saveRDS(Fibroblast, file.path(Fibro_out_dir, "Fibroblast.rds"))
#Fibroblast <- readRDS("result/Fibroblast_subtype/Fibroblast.rds")

p_col_proj <-
  Fibroblast[[]] %>%
  group_by(group, RNA_snn_res.0.4) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count),
         RNA_snn_res.0.4 = factor(RNA_snn_res.0.4,levels = c(0:13) )) %>%
  
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.4)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

table(Fibroblast$RNA_snn_res.0.4)

ggsave("result/Fibroblast_subtype/cluster_group_prop.pdf",p_col_proj, height = 6, width = 8)

# dimplot
pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} anno harmonyumap.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "res_0.8", label = TRUE, reduction = 'harmonyumap') 
DimPlot(Fibroblast, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') 
dev.off()

pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} anno umap.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "res_0.8", label = TRUE, reduction = 'umap') 
DimPlot(Fibroblast, group.by = "annotation", label = TRUE, reduction = 'umap') 

dev.off()

# dimplot
pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} cluster _split group harmonyumap.pdf")), height = 6, width = 30)
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap',split.by = "group",pt.size = 1.5) + NoLegend()
#DimPlot(Fibroblast, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()

dev.off()

# dimplot
pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} cluster _split group umap.pdf")), height = 6, width = 30)
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'umap',split.by = "group",pt.size = 1.5) + NoLegend()
#DimPlot(Fibroblast, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()

dev.off()



## annotation ----
cell_cluster <- c(
  "0" = "Myofibroblast",
  "1" = "Mix",
  "2" = "Myofibroblast",
  "3" = "Fibroblast",
  "4" = "Fibroblast",
  "5" = "Pro_Fibroblast",
  "6" = "Myofibroblast",
  "7" = "Mix",
  "8" = "Pro_Fibroblast",
  "9" = "Mix",
  "10" = "Mix",
  "11" = "Myofibroblast"
)

Idents(Fibroblast) <- "RNA_snn_res.0.4"
#Idents(Fibroblast) <- "sub.cluster"

Fibroblast <- RenameIdents(Fibroblast, cell_cluster)
Fibroblast$subtype <- Idents(Fibroblast)

# Fib_marker_ls <- c(
#   "iCAF" = 'CXCL12',
#   "iCAF" = 'CXCL9',
#   "iCAF" = 'CXCL10',
#   "iCAF" = 'CXCL11',
#   
#   "mCAF" = 'POSTN',
#   "mCAF" = 'COL10A1',
#   "mCAF" = 'COL11A1',
#   "mCAF" = 'CTHRC1',
#   "mCAF" = 'MMP11',
#   
#   "Epi_like CAF" = 'EPCAM',
#   "Epi_like CAF" = 'CDH1',
#   "Epi_like CAF" = 'STMN1',
#   
#   "Fib" = 'CCDC80',
#   "Fib" = 'SFRP4',
#   
#   "Mac_like CAF" = 'CD68',
#   
#   "Plasma_like CAF" = 'MZB1',
#   "Plasma_like CAF" = 'XBP1',
#   
#   "meCAF" = "BNIP3" ,
#   "meCAF" ='SLC2A1' ,
#   "meCAF" ='ADM',
#   "meCAF" ='ANGPTL4'
# )
# 
# pdf(file.path(Fibro_out_dir, "dotplot celltype marker.pdf"), width = 20, height = 8)
# DotPlot(Fibroblast, group.by = "subtype", features = Fib_marker_ls, dot.scale = 10) + theme_90angle
# dev.off()

pdf(file.path(Fibro_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
DimPlot(Fibroblast, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 30)
DimPlot(Fibroblast, group.by = "subtype", split.by = 'group',label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
DimPlot(Fibroblast, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

saveRDS(Fibroblast, file.path(Fibro_out_dir,"Fibroblast.rds"))

## cell props barplot  --------
Fibroblast_meta <- Fibroblast[[]]

p_col_proj <-
  Fibroblast_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(Fibro_out_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 16
)


p_col_proj <-
  Fibroblast_meta %>%
  group_by(group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(Fibro_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 6
)


xenium <- readRDS("result/xenium.rds")
#2、炎症发生的时间点：INF、IL等炎症标记在哪个时间点相比control变化明显-----
# inflammatory_genes <- c(
#   # IL 家族
#   "Il1a","Il1b","Il2","Il6","Il8","Il10","Il12a","Il15","Il17a","Il18",
#   # 炎症通路
#   "Tnf","Tnfrsf1a","Nfkb1","Rela","Irf1","Irf7","Ptgs2"
# )

INFLAMMATORY_RESPONSE <- c(
  # IL 白介素家族
  "Il10","Il10ra","Il12b","Il15","Il15ra","Il18","Il18r1","Il18rap","Il1a","Il1b","Il1r1","Il2rb","Il4ra","Il6","Il7r",
  # 趋化因子 Chemokine (Ccl/Cxcl 家族)
  "Ccl17","Ccl20","Ccl22","Ccl24","Ccl5","Ccl7","Ccr7","Ccrl2","Cx3cl1","Cxcl10","Cxcl11","Cxcl15","Cxcl5","Cxcl9","Cxcr6",
  # NF-κB 核心炎症通路分子
  "Nfkb1","Nfkbia","Rela","Irak2","Ripk2","Nod2","Nlrp3",
  # Toll样受体 TLR 家族 & 模式识别受体
  "Tlr1","Tlr2","Tlr3","Olr1","Marco","Adgre1","Csf1","Csf3","Csf3r",
  # 干扰素相关应答基因
  "Ifitm1","Ifnar1","Ifngr2","Irf1","Irf7","Eif2ak2",
  # 细胞黏附分子 ICAM/选择素
  "Icam1","Icam4","Sele","Sell","Itga5","Itgb3","Itgb8","Slamf1",
  # TNF 超家族及受体
  "Tnfaip6","Tnfrsf1b","Tnfrsf9","Tnfsf10","Tnfsf15","Tnfsf9","Lta",
  # 前列腺素、炎症脂质代谢相关
  "Ptger2","Ptger4","Ptgir","Ptafr","Plaur","Pde4b","Sphk1",
  # 缺氧、应激应答
  "Hif1a","Nampt","Inhba","Osm","Osmr","Serpine1","Timp1","Ereg","Ahr",
  # 细胞周期、增殖调控
  "Cdkn1a","Myc","Mxd1","Klf6","Raf1","Rasgrp1","Rnf144b","Tpbg","Met","Axl",
  # 补体、炎症介质受体
  "C3ar1","Calcrl","Bdkrb1","Fpr1","P2rx4","P2rx7","P2ry2","Adora2b","Adm","Aplnr","Hrh1","Tacr1","Tacr3",
  # 转运/离子通道
  "Kcna3","Kcnj2","Kcnmb2","Atp2a2","Atp2b1","Atp2c1","Slc11a2","Slc1a2","Slc28a2","Slc31a1","Slc31a2","Slc4a4","Slc7a1","Slc7a2","Aqp9",
  # 基质金属蛋白酶、组织重塑
  "Mmp14","Pdpn","Stab1","Scarf1","Has2","Sema4d","Vip","Prok2","Psen1","Pvr","Ptpre","Rgs1","Rgs16","Rhog","Ros1","Rtp4","Scn1b","Selenos","Sgms2","Sri","Tapbp","Lamp3","Lck","Lcp2","Ldlr","Lif","Lpar1","Ly6e","Lyn","Mefv","Mep1a","Msr1","Ndp","Nmi","Nmur1","Npffr2","Oprk1","Pcdh7","Pik3r5","Gna15","Gnai3","Gp1ba","Gpc3","Gpr132","Gpr183","Gch1","Gabbr1","Fzd5","Ffar2","F3","Emp3","Edn1","Ebi3","Dcbld2","Cybb","Cd82","Cd70","Cd69","Cd48","Cd40","Cd14","Clec5a","Cmklr1","Chst2","Bt g2","Bst2","Best1","Abi1","Abca1"
)



valid_gene <- INFLAMMATORY_RESPONSE[INFLAMMATORY_RESPONSE %in% rownames(xenium)]
cat("有效基因:", length(valid_gene), "/", length(INFLAMMATORY_RESPONSE))

missing_gene <- setdiff(INFLAMMATORY_RESPONSE, valid_gene)

library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ggsignif)
library(ggpubr)
##整体------
xenium <- AddModuleScore(
  object = xenium,
  features = list(INFLAMMATORY_RESPONSE),
  name = "Inflammatory_Score"
)

xenium$Inflammatory_Score <- xenium$Inflammatory_Score1

xenium$group %>% table()
##小提琴-----
V <- VlnPlot(
  xenium,
  features = "Inflammatory_Score",
  group.by = "group",
  pt.size = 0,
  fill.by = "group"
) +
  geom_boxplot(width=0.1, fill="white") +
  theme_bw() +
  ggtitle("Inflammatory Score Group") +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    ref.group = "Con", # 以Con为对照组，自动和所有其他组比较
    hide.ns = F,
    size = 3,
    vjust = 1.1,
    hjust = 1.1
  )

ggsave("result/Inflammatory_Score/Inflammatory_Score_Violin0701.pdf",V, width=7, height=5)
ggsave("result/Inflammatory_Score/Inflammatory_Score_Violin0701.png",V, width=7, height=5)


group_levels <- levels(xenium$group)
comparisons_list <- combn(group_levels, 2, simplify = FALSE)

##带所有组两两检验的小提琴图------
library(ggplot2)
library(ggpubr)

plot_data <- data.frame(
  Inflammatory_Score = xenium$Inflammatory_Score,
  group = xenium$group
)

plot_data$group %>% levels()

group_levels <- levels(plot_data$group)
comparisons_list <- combn(group_levels, 2, simplify = FALSE)

V <- ggplot(plot_data, aes(x = group, y = Inflammatory_Score, fill = group)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_bw() +
  ggtitle("Inflammatory Score Group") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  stat_compare_means(
    comparisons = comparisons_list,
    method = "wilcox.test",
    label = "p.signif",
    p.adjust.method = "bonferroni",
    hide.ns = FALSE,
    size = 3#,
    # vjust = 0,
    # hjust = 0.5
  )

ggsave("result/Inflammatory_Score/Inflammatory_Score_Violin wilcox 0701.pdf", V, width = 7, height = 10)
ggsave("result/Inflammatory_Score/Inflammatory_Score_Violin wilcox 0701.png", V, width = 7, height = 10)

# df <- data.frame(
#   group = xenium$group,
#   score = xenium$Inflammatory_Score
# )
# 
# B <- ggplot(df, aes(x = group, y = score, fill = group)) +
#   geom_boxplot() +
#   theme_bw() +
#   ggtitle("Inflammatory Score Group") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# ggsave("result/Inflammatory_Score/Inflammatory_Score_box.pdf", B, width=7, height=5)
# ggsave("result/Inflammatory_Score/Inflammatory_Score_box.png",B, width=7, height=5)
# 


##成纤维内------
Fibroblast <- AddModuleScore(
  object = JoinLayers(Fibroblast),
  features = list(INFLAMMATORY_RESPONSE),
  name = "Inflammatory_Score"
)

Fibroblast$Inflammatory_Score <- Fibroblast$Inflammatory_Score1

##小提琴-----
V <- VlnPlot(
  Fibroblast,
  features = "Inflammatory_Score",
  group.by = "group",
  pt.size = 0,
  fill.by = "group"
) +
  geom_boxplot(width=0.1, fill="white") +
  theme_bw() +
  ggtitle("Inflammatory Score Group") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    ref.group = "Con", # 以Con为对照组，自动和所有其他组比较
    hide.ns = F,
    size = 3,
    vjust = 1.1
  )


ggsave("result/Inflammatory_Score/Fibro Inflammatory_Score_Violin0701.pdf",V, width=7, height=5)
ggsave("result/Inflammatory_Score/Fibro Inflammatory_Score_Violin0701.png",V, width=7, height=5)



# df <- data.frame(
#   group = Fibroblast$group,
#   score = Fibroblast$Inflammatory_Score
# )
# 
# B <- ggplot(df, aes(x = group, y = score, fill = group)) +
#   geom_boxplot() +
#   theme_bw() +
#   ggtitle("Inflammatory Score Group") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# ggsave("result/Inflammatory_Score/Fibro Inflammatory_Score_box.pdf", B, width=7, height=5)
# ggsave("result/Inflammatory_Score/Fibro Inflammatory_Score_box.png",B, width=7, height=5)
# 




#3、找到发展为肌成纤维细胞的前体成纤维类型，研究其与其他细胞类型的互作-------
##Fibroblast内拟时序------
library(Seurat)
library(ggplot2)
library(ggrepel)
library(Matrix)
library(dplyr)
library(hdf5r)
library(patchwork)
library(RColorBrewer)
library(monocle)
library(tibble)

Fibroblast <- readRDS(file.path(Fibro_out_dir,'Fibroblast.rds'))
table(Fibroblast@meta.data$subtype)

Fibro_monocle_outdir <- "result/Fibro_monocle_result/"
dir.create(Fibro_monocle_outdir)

unique(Fibroblast@meta.data$subtype)
Idents(Fibroblast)

Fibro_mono <- Fibroblast[,!Fibroblast$subtype %in% c("Mix")]

saveRDS(Fibro_mono,file.path(Fibro_monocle_outdir,"Fibro_mono.rds"))


## 投shell------
#/mnt/NexData_file/home/hyj/shell/FWSP20254307_yurongbing/BMMSC_monocle_result

plan <- c("#ea7070", "#600473",
          "#fdc4b6", 
          "#e59572", 
          "#2694ab", "#96ceb4", 
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")

cds_Fibro <- readRDS("/mnt/NexData_file/home/hyj/shell/XE2511036121-01_mengliqiang/Fibro_monocle_result/cds_disp-1.rds")

pData(cds_Fibro)$RNA_snn_res.0.4 <- factor(pData(cds_Fibro)$RNA_snn_res.0.4, levels = c(0:11)#, levels = c(paste0("Normal_C", 0:6), paste0("Normal_C", 8:11), paste0("Tumor_C", 0:9))
)

pdf(file.path(Fibro_monocle_outdir, "cell trajectory.pdf"), width = 8, height = 6)
plot_cell_trajectory(cds_Fibro,color_by = "Pseudotime", size = 1, show_backbone = TRUE) + theme(legend.position = "right")
plot_cell_trajectory(cds_Fibro,color_by = "subtype", size = 1, show_backbone = TRUE) + 
  scale_color_manual(values = plan) +
  theme(legend.position = "right")
plot_cell_trajectory(cds_Fibro,color_by = "RNA_snn_res.0.4", size = 1, show_backbone = TRUE) + 
  scale_color_manual(values = plan) +
  theme(legend.position = "right")
plot_cell_trajectory(cds_Fibro, color_by = "State", size = 1, show_backbone = TRUE) +
  scale_color_manual(values = plan) +
  theme(legend.position = "right")
plot_cell_trajectory(cds_Fibro, color_by = "group", size = 1, show_backbone = TRUE) +
  scale_color_manual(values = plan) +
  theme(legend.position = "right")

dev.off()

pdf(file.path(Fibro_monocle_outdir, "cell trajectory - - split by subtype.pdf"), width = 10, height = 10)
plot_cell_trajectory(cds_Fibro, color_by = "subtype") + facet_wrap("~subtype", nrow = 2) + theme(legend.position = "right") + scale_color_manual(values = plan) 
dev.off()

pdf(file.path(Fibro_monocle_outdir, "cell trajectory - - split by RNA_snn_res.0.4.pdf"), width = 20, height = 10)
plot_cell_trajectory(cds_Fibro, color_by = "RNA_snn_res.0.4") + facet_wrap("~RNA_snn_res.0.4", nrow = 2) + theme(legend.position = "right") + scale_color_manual(values = plan) 
dev.off()

pdf(file.path(Fibro_monocle_outdir, "cell trajectory - - split by state.pdf"), width = 15, height = 10)
plot_cell_trajectory(cds_Fibro, color_by = "State") + facet_wrap("~State", nrow = 2) + scale_color_manual(values = plan)
dev.off()

pdf(file.path(Fibro_monocle_outdir, "cell trajectory - - split by group.pdf"), width = 20, height = 10)
plot_cell_trajectory(cds_Fibro, color_by = "group") + facet_wrap("~group", nrow = 2) + scale_color_manual(values = plan)
dev.off()



pdf(file.path(Fibro_monocle_outdir, "cell trajectory - subtype - split by group.pdf"), width = 12, height = 10)
local({
  #pData(cds_Fibro)$group <- pData(cds_Fibro)$group %>% factor(levels = c("PA", "PC"))
  plot_cell_trajectory(cds_Fibro, color_by = "subtype", size = 1, show_backbone = TRUE) + 
    scale_color_manual(values = plan) + 
    facet_wrap(~ group) + 
    theme(legend.position = "right")
})
dev.off()

pdf(file.path(Fibro_monocle_outdir, "cell trajectory - cluster - split by group.pdf"), width = 12, height = 10)
local({
  #pData(cds_Fibro)$group <- pData(cds_Fibro)$group %>% factor(levels = c("PA", "PC"))
  plot_cell_trajectory(cds_Fibro, color_by = "RNA_snn_res.0.4", size = 1, show_backbone = TRUE) + 
    scale_color_manual(values = plan) + 
    facet_wrap(~ group) + 
    theme(legend.position = "right")
})
dev.off()




pdf(file.path(Fibro_monocle_outdir, "cell trajectory - tree.pdf"), width = 15, height = 5)
plot_complex_cell_trajectory(cds_Fibro, x = 1, y = 2, color_by = "subtype") + theme(legend.position = "right") + 
  scale_color_manual(values = plan) +
  facet_wrap(~ subtype)
dev.off()

pdf(file.path(Fibro_monocle_outdir, "cell trajectory  group - tree.pdf"), width = 15, height = 15)
plot_complex_cell_trajectory(cds_Fibro, x = 1, y = 2, color_by = "subtype") + theme(legend.position = "right") + 
  scale_color_manual(values = plan) +
  facet_wrap(~ group)
dev.off()



# cell barplot 
p_col_cell_state <- 
  cds_Fibro %>% 
  pData() %>% 
  dplyr::select(RNA_snn_res.0.4, State) %>% 
  group_by(RNA_snn_res.0.4, State) %>% 
  # mutate(RNA_snn_res.0.4 = as.character(RNA_snn_res.0.4)) %>% 
  # mutate(RNA_snn_res.0.4 = factor(RNA_snn_res.0.4#, levels = c(
  #                                 #   paste0("Normal_C", 0:6), paste0("Normal_C", 8:11), paste0("Tumor_C", 0:9)
  #                                 # )
  #)
   #%>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count)) %>% 
  
  ggplot(aes(x = RNA_snn_res.0.4, y = freq)) + 
  geom_col(aes(fill = State)) + 
  scale_fill_manual(values = plan) + 
  scale_y_continuous(labels = scales::percent) + 
  coord_flip() +
  labs(x = "cluster", y = NULL) + 
  ggpubr::theme_pubr() + 
  theme(
    legend.position = "right"
  )
ggsave(p_col_cell_state, filename = file.path(Fibro_monocle_outdir, "cell stat barplot.pdf"), height = 8, width = 8)


# cell barplot 
p_col_cell_state <- 
  cds_Fibro %>% 
  pData() %>% 
  dplyr::select(RNA_snn_res.0.4, group) %>% 
  group_by(RNA_snn_res.0.4, group) %>% 
  # mutate(RNA_snn_res.0.2 = as.character(RNA_snn_res.0.2)) %>% 
  # mutate(RNA_snn_res.0.2 = factor(RNA_snn_res.0.2#, levels = c(
  #                                 #   paste0("Normal_C", 0:6), paste0("Normal_C", 8:11), paste0("Tumor_C", 0:9)
  #                                 # )
  # )
  # ) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count)) %>% 
  
  ggplot(aes(x = RNA_snn_res.0.4, y = freq)) + 
  geom_col(aes(fill = group)) + 
  scale_fill_manual(values = plan) + 
  scale_y_continuous(labels = scales::percent) + 
  coord_flip() +
  labs(x = "cluster", y = "") + 
  ggpubr::theme_pubr() + 
  theme(
    legend.position = "right"
  )
ggsave(p_col_cell_state, filename = file.path(Fibro_monocle_outdir, "cell stat cluster group barplot.pdf"), height = 8, width = 8)


p_col_cell_state <- 
  cds_Fibro %>% 
  pData() %>% 
  dplyr::select(State, group) %>% 
  group_by(State, group) %>% 
  mutate(State = as.character(State)) %>% 
  mutate(State = factor(State#, levels = c(
                        #   paste0("Normal_C", 0:6), paste0("Normal_C", 8:11), paste0("Tumor_C", 0:9)
                        # )
  )
  ) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count)) %>% 
  
  ggplot(aes(x = State, y = freq)) + 
  geom_col(aes(fill = group)) + 
  scale_fill_manual(values = plan) + 
  scale_y_continuous(labels = scales::percent) + 
  coord_flip() +
  labs(x = "State", y = "") + 
  ggpubr::theme_pubr() + 
  theme(
    legend.position = "right"
  )
ggsave(p_col_cell_state, filename = file.path(Fibro_monocle_outdir, "cell stat State group barplot.pdf"), height = 8, width = 8)

p_col_cell_state <- 
  cds_Fibro %>% 
  pData() %>% 
  dplyr::select(State, subtype) %>% 
  group_by(State, subtype) %>% 
  mutate(State = as.character(State)) %>% 
  mutate(State = factor(State#, levels = c(
                        #   paste0("Normal_C", 0:6), paste0("Normal_C", 8:11), paste0("Tumor_C", 0:9)
                        # )
  )
  ) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count)) %>% 
  
  ggplot(aes(x = State, y = freq)) + 
  geom_col(aes(fill = subtype)) + 
  scale_fill_manual(values = plan) + 
  scale_y_continuous(labels = scales::percent) + 
  coord_flip() +
  labs(x = "State", y = "") + 
  ggpubr::theme_pubr() + 
  theme(
    legend.position = "right"
  )
ggsave(p_col_cell_state, filename = file.path(Fibro_monocle_outdir, "cell stat State subtype barplot.pdf"), height = 8, width = 8)



p_col_cell_state <- 
  cds_Fibro %>% 
  pData() %>% 
  dplyr::select(State, group) %>% 
  group_by (State,group) %>% 
  mutate(State = as.character(State)) %>% 
  mutate(State = factor(State#, levels = c(
                        #   paste0("Normal_C", 0:6), paste0("Normal_C", 8:11), paste0("Tumor_C", 0:9)
                        # )
  )
  ) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count)) %>% 
  
  ggplot(aes(x = group, y = freq)) + 
  geom_col(aes(fill = State)) + 
  scale_fill_manual(values = plan) + 
  scale_y_continuous(labels = scales::percent) + 
  coord_flip() +
  labs(x = "group", y = "") + 
  ggpubr::theme_pubr() + 
  theme(
    legend.position = "right"
  )
ggsave(p_col_cell_state, filename = file.path(Fibro_monocle_outdir, "cell stat  group  State 111barplot.pdf"), height = 8, width = 8)



# 沿时间轴的细胞密度图 -----------
colour = c('#C6DADF', '#E26EB8', '#918FBE', '#A2C5A6', '#EED784', '#D4B796', 
           '#F5BDC8', '#F0BDD8', '#77A1EC', '#F8D39C', '#F17591', '#CCC9E6')

df <- pData(cds_Fibro)[, c("Pseudotime", "RNA_snn_res.0.4")]
p_density <- ggplot(df, aes(Pseudotime, colour = RNA_snn_res.0.4, fill = RNA_snn_res.0.4)) + 
  geom_density(bw = 0.5, size = 0.5, alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name = "", values = plan) + 
  scale_color_manual(name = "", values = plan)
ggsave(p_density, filename = file.path(Fibro_monocle_outdir, "cell density cluster plot.pdf"), height = 6, width = 12)


df <- pData(cds_Fibro)[, c("Pseudotime", "group")]
p_density <- ggplot(df, aes(Pseudotime, colour = group, fill = group)) + 
  geom_density(bw = 0.5, size = 0.5, alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name = "", values = plan) + 
  scale_color_manual(name = "", values = plan)
ggsave(p_density, filename = file.path(Fibro_monocle_outdir, "cell density group plot.pdf"), height = 6, width = 12)


df <- pData(cds_Fibro)[, c("Pseudotime", "subtype")]
p_density <- ggplot(df, aes(Pseudotime, colour = subtype, fill = subtype)) + 
  geom_density(bw = 0.5, size = 0.5, alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name = "", values = plan) + 
  scale_color_manual(name = "", values = plan)
ggsave(p_density, filename = file.path(Fibro_monocle_outdir, "cell density subtype plot.pdf"), height = 6, width = 12)


#4、LIF(白血病抑制因子)在空间feature映射，和其他细胞类型共同映射------
Lif
source("ImageFeature-plot.R")

Lif_Imagefeature_dir <- 'result/Lif_Imagefeature' 
dir.create(Lif_Imagefeature_dir)

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


all_colors <- c(
  # 上皮 / 小管上皮 (暖黄 - 橙黄系)
  "Epithelial"    = "#FF9F1C",
  "PT"            = "#FEE440",
  "Podocyte"      = "#DBB466",
  
  # 基质 / 间充质 (绿色系)
  "Fibroblast"    = "#52AF43",
  "SMC"           = "#6F9E4C",
  "Proliferate"   = "#98D277",
  
  # 血管 / 内皮 (蓝青 / 紫色系)
  "Endothelial"   = "#00BBF9",
  "PC"            = "#4361EE",
  "DT"            = "#02c9c9",
  "IC"            = "#2D82AF",
  
  # 免疫 / 炎症 (红色 / 粉色系)
 
  "Macrophage"    = "#F15BB5"
)


####feature--------
for (sample_id in unique(xenium@meta.data$sample)) {

    DefaultBoundary(xenium@images[[sample_id]]) <- "segmentation"
    p <- ImageFeaturePlot(
      object = xenium,
      fov = sample_id,
      nmols = 20000,
      border.color = "black",
      border.size = 0.01,
      features = c("Lif"),
      max.cutoff = "q90",
      axes = FALSE,
      cols = c('#003366', '#FFD700'),
      # cols = c('blue', "red"),
      size = 0.8,
      dark.background = F#,
      #cells = WhichCells(xenium_score, idents = "Macrophage" )
    )

    ggsave(
      filename = file.path(Lif_Imagefeature_dir, paste0(sample_id, "_Lif_feature_spatial.png")),
      plot = p,
      width = 12,
      height = 12,
      dpi = 300
    )
    ggsave(
      filename = file.path(Lif_Imagefeature_dir, paste0(sample_id, "_Lif_feature_spatial.pdf")),
      plot = p,
      width = 12,
      height = 12
    )

    message("成功保存: ", sample_id)
  }


xenium$annotation %>% table()
###feature➕边界--------
for (sample_id in unique(xenium@meta.data$sample)){
    p <- ImageFeaturePlot_wyt(xenium, 
                              features = c("Lif"), 
                              size = 0.5, 
                              cols = c('#003366', '#FFD700'),
                              #cols = c('blue', "red"),
                              fov =sample_id,
                              nmols = 20000,
                              border.color = all_colors, 
                              border.size = 0.2,
                              border.by = "annotation",
                              #cells.use = c("Macrophage"),
                              axes = FALSE,
                              dark.background = F
    )
    
    ggsave(
      filename = file.path(Lif_Imagefeature_dir, paste0(sample_id, "_Lif_spatial.png")),
      plot = p,
      width = 12,
      height = 12,
      dpi = 300
    )
    ggsave(
      filename = file.path(Lif_Imagefeature_dir, paste0(sample_id, "_Lif_spatial.pdf")),
      plot = p,
      width = 12,
      height = 12
    )
    
    message("成功保存: ", sample_id)
}


#5、内皮亚类注释-------
# Endothelial --------------------------------------------------------------
all_data$annotation %>% table()
all_data$res_0.8 %>% table()
table(all_data$res_0.8,all_data$annotation)

Endo_out_dir<-'result/Endothelial_subtype' 
dir.create(Endo_out_dir)
Endothelial<-subset(all_data, subset = annotation %in% c("Endothelial"))
table(Endothelial$res_0.8,Endothelial$annotation)

##修改分组--------
Endothelial$group_raw <- Endothelial$group

replace_dict <- c(
  "A1" = "H1",
  "A3" = "H12",
  "B1" = "Con",
  "C6" = "D7",
  "D2" = "H3",
  "F1" = "H24",
  "F2" = "H72"
)

Endothelial$group <- plyr::revalue(Endothelial$group, replace = replace_dict)
group_order <- c("Con", "H1", "H3", "H12", "H24", "H72", "D7")
Endothelial$group <- factor(Endothelial$group, levels = group_order)

table(Endothelial$group_raw ,Endothelial$group)

Endothelial$group %>% table()
Endothelial$annotation %>% table()


##修改样本名--------
table(Endothelial$sample,Endothelial$group)
Endothelial$sample %>% table()
Endothelial$sample_raw <- Endothelial$sample

Endothelial$sample <- gsub("A1" , "H1",Endothelial$sample)  
Endothelial$sample <- gsub("A3" , "H12",Endothelial$sample) 
Endothelial$sample <- gsub("B1" , "Con",Endothelial$sample) 
Endothelial$sample <- gsub("C6" , "D7",Endothelial$sample) 
Endothelial$sample <- gsub("D2" , "H3",Endothelial$sample) 
Endothelial$sample <- gsub("F1" , "H24",Endothelial$sample) 
Endothelial$sample <- gsub("F2" , "H72",Endothelial$sample) 
table(Endothelial$sample,Endothelial$group)

sample_order <- Endothelial[[]] %>%
  arrange(group) %>%
  pull(sample) %>%
  unique()

Endothelial$sample <- factor(Endothelial$sample, levels = sample_order)

## PCA --------
DefaultAssay(Endothelial) <- "RNA"
Endothelial <- FindVariableFeatures(Endothelial)
Endothelial <- ScaleData(Endothelial, vars.to.regress = 'nCount_RNA')
Endothelial <- RunPCA(Endothelial)


### pca inspector ------
pdf(file.path(Endo_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Endothelial, ndims = 50)
dev.off()

# pdf(file.path(Endo_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
# DimHeatmap(Endothelial, dims = 1:5)
# dev.off()
# pdf(file.path(Endo_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
# DimHeatmap(Endothelial, dims = 6:10)
# dev.off()
# pdf(file.path(Endo_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
# DimHeatmap(Endothelial, dims = 11:15)
# dev.off()
# pdf(file.path(Endo_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
# DimHeatmap(Endothelial, dims = 16:20)
# dev.off()



## dimension reduction --------
dims <- 20

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb

Endothelial <- RunUMAP(Endothelial, dims = 1:dims,return.model = TRUE)

pdf(file.path(Endo_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Endothelial, group.by = "sample", label = TRUE) + NoLegend()
dev.off()


pdf(file.path(Endo_out_dir, str_glue("Dimplot dim{dims} group.pdf")), height = 6, width = 8)
DimPlot(Endothelial, group.by = "group", label = TRUE)
dev.off()

table(all_data@meta.data$celltype, all_data@meta.data$project)

Endothelial <- Endothelial %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

Endothelial <- RunUMAP(
  Endothelial, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Endo_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Endothelial, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



pdf(file.path(Endo_out_dir, str_glue("Dimplot dim{dims} group - harmony.pdf")), height = 6, width = 8)
DimPlot(Endothelial, group.by = "group", label = TRUE, reduction = "harmonyumap")
dev.off()




reduction.use <- 'harmony' # 'harmony' 'pca'

### cluster -------
Endothelial  <-
  FindNeighbors(
    Endothelial ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

Endothelial <-
  FindClusters(
    Endothelial,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )


# dimplot
pdf(file.path(Endo_out_dir, str_glue("Dimplot dim{dims} cluster umap.pdf")), height = 6, width = 8)
DimPlot(Endothelial, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Endothelial, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Endothelial, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Endothelial, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'umap') + NoLegend()
# DimPlot(Endothelial, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Endothelial, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Endothelial, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

# dimplot
pdf(file.path(Endo_out_dir, str_glue("Dimplot dim{dims} cluster harmonyumap.pdf")), height = 6, width = 8)
DimPlot(Endothelial, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Endothelial, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Endothelial, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Endothelial, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Endothelial, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Endothelial, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Endothelial, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

saveRDS(Endothelial, "result/Endothelial_subtype/Endothelial.rds")
#Endothelial <- readRDS("result/Endothelial_subtype/Endothelial.rds")

### dotplot -------
celltype_marker_ls <- c(
  # Epithelial
  'Epithelial' = "Epcam",
  'Epithelial' = "Krt19",
  'Epithelial' = "Cryab",
  
  # POD
  'POD' = "Nphs2",
  'POD' = "Podxl",
  
  # HealthyS1
  'PT' = "Slc34a1",
  'PT' = "Aldob",
  'PT' = "Lrp2",
  
  # LOH
  'LOH' = "Umod",
  'LOH' = "Slc12a1",
  
  # DT
  'DT' = "Slc12a3",
  'DT' = "Calb1",
  
  # PC
  'PC' = "Aqp2",
  'PC' = "Hsd11b2",
  
  # Immune
  'Immune' = "Ptprc",
  'Immune' = "Cd74",
  'Immune' = "Ccl5",
  
  # Fibroblast
  'Fibroblast' = "Dcn",
  'Fibroblast' = "Col1a1",
  
  # SMC
  'SMC' = "Acta2",
  'SMC' = "Tagln",
  
  # Adipocyte
  'Adipocyte' = "Fabp4",
  'Adipocyte' = "Plin4"
)



Fibro_major_marker_ls <- c(
  # S1 fibroblasts
  'S1_fibroblasts' = "Adamdec1", 
  'S1_fibroblasts' = "Cp", 
  'S1_fibroblasts' = "Fabp4",
  
  # S3 fibroblasts
  'S3_fibroblasts' = "Ogn", 
  'S3_fibroblasts' = "Ccdc80", 
  'S3_fibroblasts' = "Grem1",
  
  # Peri-cryptal S2 fibroblasts
  'S2_fibroblasts' = "Sox6",
  'S2_fibroblasts' = "F3",
  
  # S2a fibroblasts
  'S2a_fibroblasts' = "Vstm2",
  
  # S2b fibroblasts
  'S2b_fibroblasts' = "Npy",
  'S2b_fibroblasts' = "Nrg1",
  
  # Ribhi fibroblasts
  'Ribhi_fibroblasts' = "Rpl7a",
  'Ribhi_fibroblasts' = "Rpl28",
  
  # MT fibroblasts
  'MT_fibroblasts' = "Mt-genes",  # 线粒体基因通用标记
  
  # Immediate early-response (IER) fibroblasts
  'IER_fibroblasts' = "Fos", 
  'IER_fibroblasts' = "Fosb", 
  'IER_fibroblasts' = "Irf1",
  
  # Fibroblastic reticular cells (FRCs)
  'FRCs' = "Ccl19", 
  'FRCs' = "Ccl21",
  
  # Smooth muscle cells (SMCs)
  'SMC' = "Tagln",
  'SMC' = "Mylk",
  'SMC' = "Cnn1",
  
  # Fibroblasts (FIB)
  'FIB' = "Lum",
  'FIB' = "Dcn",
  'FIB' = "Col1a1",
  'FIB' = "Col1a2",
  'FIB' = "Fn1",
  
  # Myofibroblasts
  'Myofibroblasts' = "Actg2",
  'Myofibroblasts' = "Sostdc1",
  'Myofibroblasts' = "Myh11",
  
  
  # Inflammatory fibroblasts
  'Inflammatory_fibroblasts' = "Il11",
  'Inflammatory_fibroblasts' = "Cxcl5", 
  'Inflammatory_fibroblasts' = "Fap", 
  'Inflammatory_fibroblasts' = "Inhba",
  'Inflammatory_fibroblasts' = "Il24",
  
  # Endothelial cells
  'Endothelial_cells' = "Vwf", 
  'Endothelial_cells' = "Pecam1", 
  'Endothelial_cells' = "Plvap",
  
  # Activated endothelium
  'Activated_endothelium' = "Ackr1", 
  'Activated_endothelium' = "Sele", 
  'Activated_endothelium' = "Selp",
  
  # Lymphatic endothelium
  'Lymphatic_endothelium' = "Mrc1", 
  'Lymphatic_endothelium' = "Mmrn1", 
  'Lymphatic_endothelium' = "Lyve1",
  
  # Pericytes
  'Pericytes' = "Notch3", 
  'Pericytes' = "Rgs5",
  
  # Glial cells
  'Glial_cells' = "Nrxn1", 
  'Glial_cells' = "S100b", 
  'Glial_cells' = "Cdh19"
)


Fibro_major_marker_ls1 <- c(
  # Smooth muscle cells (SMCs)
  'SMC' = "Tagln",
  'SMC' = "Mylk",
  'SMC' = "Cnn1",#'Acta2'
  
  # Fibroblasts (FIB)
  'FIB' = "Lum",
  'FIB' = "Dcn",
  'FIB' = "Col1a1",
  'FIB' = "Col1a2",
  'FIB' = "Fn1",
  'FIB' = "C7",
  'FIB' = "Ccn2",
  
  # Myofibroblasts
  'Myofibroblasts' = "Actg2",
  'Myofibroblasts' = "Sostdc1",
  'Myofibroblasts' = "Myh11", #'Acta2',
  'Myofibroblasts' = "Notch3"
)


Idents(Endothelial)<-'RNA_snn_res.0.4'
Endothelial_markers <-
  FindAllMarkers(
    JoinLayers(Endothelial),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

Endothelial_top10_markers <- 
  Endothelial_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10) %>%
  nest() %>%
  deframe() %>% 
  map(pull, 'gene') %>%
  imap(function(genes, name){
    re <- genes %>%
      map(function(gene){
        # Endothelialrowser()
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(Endothelial_top10_markers) <- sub('.*\\.','',names(Endothelial_top10_markers))

Endothelial_top10_markers1 <-
  Endothelial_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10)

# dotplot
pdf(file.path(Endo_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Endothelial, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()
DotPlot(Endothelial, group.by = "RNA_snn_res.0.4", features = celltype_marker_ls)+RotatedAxis()
DotPlot(Endothelial, group.by = "RNA_snn_res.0.4", features = Fibro_major_marker_ls)+RotatedAxis()
DotPlot(Endothelial, group.by = "RNA_snn_res.0.4", features = Fibro_major_marker_ls1)+RotatedAxis()
DotPlot(Endothelial, group.by = "RNA_snn_res.0.4", features = Endothelial_top10_markers)+RotatedAxis()
dev.off()

table(Endothelial$RNA_snn_res.0.4,Endothelial$group)


pdf(file.path(Endo_out_dir, str_glue("Dotplot dims{dims} myofib_marker.pdf")), width = 5, height =5)
DotPlot(Endothelial, 
        group.by = "group", 
        features = c('Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1',"Cnn1","Ccn2"))+RotatedAxis()
dev.off()


pdf(file.path(Endo_out_dir, str_glue("Vlnplot dims{dims} myofib_marker.pdf")), width = 20, height = 10)
VlnPlot(Endothelial,group.by = "group",
        feature = c('Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1',"Cnn1","Ccn2"),
        pt.size = 0,
        ncol = 4,
)
dev.off()

pdf(file.path(Endo_out_dir, str_glue("Vlnplot dims{dims} myofib_marker pt.pdf")), width = 20, height = 10)
VlnPlot(Endothelial,group.by = "group",
        feature = c('Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1',"Cnn1","Ccn2"),
        pt.size = 0.1,
        ncol = 4,
)
dev.off()

pdf(file.path(Endo_out_dir, str_glue("Featureplt dims{dims} myofib_marker.pdf")), width = 20, height = 10)
FeaturePlot(Endothelial,
            #group.by = "group",
            feature = c('Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1',"Cnn1","Ccn2"),
            reduction = 'umap',
            ncol = 4,
)

dev.off()


saveRDS(Endothelial, file.path(Endo_out_dir, "Endothelial.rds"))
#Endothelial <- readRDS("result/Endothelial_subtype/Endothelial.rds")

p_col_proj <-
  Endothelial[[]] %>%
  group_by(group, RNA_snn_res.0.4) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count),
         RNA_snn_res.0.4 = factor(RNA_snn_res.0.4,levels = c(0:13) )) %>%
  
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.4)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

table(Endothelial$RNA_snn_res.0.4)

ggsave("result/Endothelial_subtype/cluster_group_prop.pdf",p_col_proj, height = 6, width = 8)

# dimplot
pdf(file.path(Endo_out_dir, str_glue("Dimplot dim{dims} anno harmonyumap.pdf")), height = 6, width = 8)
DimPlot(Endothelial, group.by = "res_0.8", label = TRUE, reduction = 'harmonyumap') 
DimPlot(Endothelial, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') 
dev.off()

pdf(file.path(Endo_out_dir, str_glue("Dimplot dim{dims} anno umap.pdf")), height = 6, width = 8)
DimPlot(Endothelial, group.by = "res_0.8", label = TRUE, reduction = 'umap') 
DimPlot(Endothelial, group.by = "annotation", label = TRUE, reduction = 'umap') 

dev.off()

# dimplot
pdf(file.path(Endo_out_dir, str_glue("Dimplot dim{dims} cluster _split group harmonyumap.pdf")), height = 6, width = 30)
DimPlot(Endothelial, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap',split.by = "group",pt.size = 1.5) + NoLegend()
#DimPlot(Endothelial, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()

dev.off()

# dimplot
pdf(file.path(Endo_out_dir, str_glue("Dimplot dim{dims} cluster _split group umap.pdf")), height = 6, width = 30)
DimPlot(Endothelial, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'umap',split.by = "group",pt.size = 1.5) + NoLegend()
#DimPlot(Endothelial, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()

dev.off()



## annotation ----
cell_cluster <- c(
  "0" = "EndT",
  "1" = "Endothelial",
  "2" = "Endothelial",
  "3" = "Endothelial",
  "4" = "Endothelial",
  "5" = "Endothelial",
  "6" = "Endothelial",
  "7" = "Endothelial",
  "8" = "EndT"
)

Idents(Endothelial) <- "RNA_snn_res.0.4"
#Idents(Endothelial) <- "sub.cluster"

Endothelial <- RenameIdents(Endothelial, cell_cluster)
Endothelial$subtype <- Idents(Endothelial)

# Fib_marker_ls <- c(
#   "iCAF" = 'CXCL12',
#   "iCAF" = 'CXCL9',
#   "iCAF" = 'CXCL10',
#   "iCAF" = 'CXCL11',
#   
#   "mCAF" = 'POSTN',
#   "mCAF" = 'COL10A1',
#   "mCAF" = 'COL11A1',
#   "mCAF" = 'CTHRC1',
#   "mCAF" = 'MMP11',
#   
#   "Epi_like CAF" = 'EPCAM',
#   "Epi_like CAF" = 'CDH1',
#   "Epi_like CAF" = 'STMN1',
#   
#   "Fib" = 'CCDC80',
#   "Fib" = 'SFRP4',
#   
#   "Mac_like CAF" = 'CD68',
#   
#   "Plasma_like CAF" = 'MZB1',
#   "Plasma_like CAF" = 'XBP1',
#   
#   "meCAF" = "BNIP3" ,
#   "meCAF" ='SLC2A1' ,
#   "meCAF" ='ADM',
#   "meCAF" ='ANGPTL4'
# )
# 
# pdf(file.path(Endo_out_dir, "dotplot celltype marker.pdf"), width = 20, height = 8)
# DotPlot(Endothelial, group.by = "subtype", features = Fib_marker_ls, dot.scale = 10) + theme_90angle
# dev.off()

pdf(file.path(Endo_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Endothelial, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
DimPlot(Endothelial, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Endo_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 30)
DimPlot(Endothelial, group.by = "subtype", split.by = 'group',label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
DimPlot(Endothelial, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

saveRDS(Endothelial, file.path(Endo_out_dir,"Endothelial.rds"))

## cell props barplot  --------
Endothelial_meta <- Endothelial[[]]

p_col_proj <-
  Endothelial_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(Endo_out_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 16
)


p_col_proj <-
  Endothelial_meta %>%
  group_by(group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(Endo_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 6
)


#6、肾小管上皮亚类注释-------
# PT--------------------------------------------------------------
all_data$annotation %>% table()
all_data$res_0.8 %>% table()
table(all_data$res_0.8,all_data$annotation)

PT_out_dir<-'result/PT_subtype' 
dir.create(PT_out_dir)
PT <-subset(all_data, subset = annotation %in% c("PT"))
table(PT$res_0.8,PT$annotation)

PT$group %>% table()
##修改分组--------
PT$group_raw <- PT$group

replace_dict <- c(
  "A1" = "H1",
  "A3" = "H12",
  "B1" = "Con",
  "C6" = "D7",
  "D2" = "H3",
  "F1" = "H24",
  "F2" = "H72"
)

PT$group <- plyr::revalue(PT$group, replace = replace_dict)
group_order <- c("Con", "H1", "H3", "H12", "H24", "H72", "D7")
PT$group <- factor(PT$group, levels = group_order)

table(PT$group_raw ,PT$group)

PT$group %>% table()
PT$annotation %>% table()


##修改样本名--------
table(PT$sample,PT$group)
PT$sample %>% table()
PT$sample_raw <- PT$sample

PT$sample <- gsub("A1" , "H1",PT$sample)  
PT$sample <- gsub("A3" , "H12",PT$sample) 
PT$sample <- gsub("B1" , "Con",PT$sample) 
PT$sample <- gsub("C6" , "D7",PT$sample) 
PT$sample <- gsub("D2" , "H3",PT$sample) 
PT$sample <- gsub("F1" , "H24",PT$sample) 
PT$sample <- gsub("F2" , "H72",PT$sample) 
table(PT$sample,PT$group)

sample_order <- PT[[]] %>%
  arrange(group) %>%
  pull(sample) %>%
  unique()

PT$sample <- factor(PT$sample, levels = sample_order)

## PCA --------
DefaultAssay(PT) <- "RNA"
PT <- FindVariableFeatures(PT)

library(future)
# 切换为串行模式，不再拆分导出大对象
plan(sequential)
PT <- ScaleData(PT, vars.to.regress = 'nCount_RNA')
PT <- RunPCA(PT)

### pca inspector ------
pdf(file.path(PT_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(PT, ndims = 50)
dev.off()

# pdf(file.path(PT_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
# DimHeatmap(PT, dims = 1:5)
# dev.off()
# pdf(file.path(PT_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
# DimHeatmap(PT, dims = 6:10)
# dev.off()
# pdf(file.path(PT_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
# DimHeatmap(PT, dims = 11:15)
# dev.off()
# pdf(file.path(PT_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
# DimHeatmap(PT, dims = 16:20)
# dev.off()



## dimension reduction --------
dims <- 20

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb

PT <- RunUMAP(PT, dims = 1:dims,return.model = TRUE)

pdf(file.path(PT_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(PT, group.by = "sample", label = TRUE) + NoLegend()
dev.off()


pdf(file.path(PT_out_dir, str_glue("Dimplot dim{dims} group.pdf")), height = 6, width = 8)
DimPlot(PT, group.by = "group", label = TRUE)
dev.off()

table(all_data@meta.data$celltype, all_data@meta.data$project)

PT <- PT %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

PT <- RunUMAP(
  PT, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(PT_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(PT, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(PT_out_dir, str_glue("Dimplot dim{dims} group - harmony.pdf")), height = 6, width = 8)
DimPlot(PT, group.by = "group", label = TRUE, reduction = "harmonyumap")
dev.off()


reduction.use <- 'harmony' # 'harmony' 'pca'

### cluster -------
PT  <-
  FindNeighbors(
    PT ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

options(future.globals.maxSize = 25 * 1024^3)

PT <-
  FindClusters(
    PT,
    resolution = c(0.2, 0.3, 0.4)
  )


# dimplot
pdf(file.path(PT_out_dir, str_glue("Dimplot dim{dims} cluster umap.pdf")), height = 6, width = 8)
DimPlot(PT, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(PT, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(PT, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'umap') + NoLegend()
# DimPlot(PT, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'umap') + NoLegend()
# DimPlot(PT, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(PT, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(PT, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

# dimplot
pdf(file.path(PT_out_dir, str_glue("Dimplot dim{dims} cluster harmonyumap.pdf")), height = 6, width = 8)
DimPlot(PT, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(PT, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(PT, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(PT, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(PT, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(PT, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(PT, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

saveRDS(PT, "result/PT_subtype/PT.rds")
#PT <- readRDS("result/PT_subtype/PT.rds")

### dotplot -------
celltype_marker_ls <- c(
  # Epithelial
  'Epithelial' = "Epcam",
  'Epithelial' = "Krt19",
  'Epithelial' = "Cryab",
  
  # POD
  'POD' = "Nphs2",
  'POD' = "Podxl",
  
  # HealthyS1
  'PT' = "Slc34a1",
  'PT' = "Aldob",
  'PT' = "Lrp2",
  
  # LOH
  'LOH' = "Umod",
  'LOH' = "Slc12a1",
  
  # DT
  'DT' = "Slc12a3",
  'DT' = "Calb1",
  
  # PC
  'PC' = "Aqp2",
  'PC' = "Hsd11b2",
  
  # Immune
  'Immune' = "Ptprc",
  'Immune' = "Cd74",
  'Immune' = "Ccl5",
  
  # Fibroblast
  'Fibroblast' = "Dcn",
  'Fibroblast' = "Col1a1",
  
  # SMC
  'SMC' = "Acta2",
  'SMC' = "Tagln",
  
  # Adipocyte
  'Adipocyte' = "Fabp4",
  'Adipocyte' = "Plin4"
)



Fibro_major_marker_ls <- c(
  # S1 fibroblasts
  'S1_fibroblasts' = "Adamdec1", 
  'S1_fibroblasts' = "Cp", 
  'S1_fibroblasts' = "Fabp4",
  
  # S3 fibroblasts
  'S3_fibroblasts' = "Ogn", 
  'S3_fibroblasts' = "Ccdc80", 
  'S3_fibroblasts' = "Grem1",
  
  # Peri-cryptal S2 fibroblasts
  'S2_fibroblasts' = "Sox6",
  'S2_fibroblasts' = "F3",
  
  # S2a fibroblasts
  'S2a_fibroblasts' = "Vstm2",
  
  # S2b fibroblasts
  'S2b_fibroblasts' = "Npy",
  'S2b_fibroblasts' = "Nrg1",
  
  # Ribhi fibroblasts
  'Ribhi_fibroblasts' = "Rpl7a",
  'Ribhi_fibroblasts' = "Rpl28",
  
  # MT fibroblasts
  'MT_fibroblasts' = "Mt-genes",  # 线粒体基因通用标记
  
  # Immediate early-response (IER) fibroblasts
  'IER_fibroblasts' = "Fos", 
  'IER_fibroblasts' = "Fosb", 
  'IER_fibroblasts' = "Irf1",
  
  # Fibroblastic reticular cells (FRCs)
  'FRCs' = "Ccl19", 
  'FRCs' = "Ccl21",
  
  # Smooth muscle cells (SMCs)
  'SMC' = "Tagln",
  'SMC' = "Mylk",
  'SMC' = "Cnn1",
  
  # Fibroblasts (FIB)
  'FIB' = "Lum",
  'FIB' = "Dcn",
  'FIB' = "Col1a1",
  'FIB' = "Col1a2",
  'FIB' = "Fn1",
  
  # Myofibroblasts
  'Myofibroblasts' = "Actg2",
  'Myofibroblasts' = "Sostdc1",
  'Myofibroblasts' = "Myh11",
  
  
  # Inflammatory fibroblasts
  'Inflammatory_fibroblasts' = "Il11",
  'Inflammatory_fibroblasts' = "Cxcl5", 
  'Inflammatory_fibroblasts' = "Fap", 
  'Inflammatory_fibroblasts' = "Inhba",
  'Inflammatory_fibroblasts' = "Il24",
  
  # Endothelial cells
  'Endothelial_cells' = "Vwf", 
  'Endothelial_cells' = "Pecam1", 
  'Endothelial_cells' = "Plvap",
  
  # Activated endothelium
  'Activated_endothelium' = "Ackr1", 
  'Activated_endothelium' = "Sele", 
  'Activated_endothelium' = "Selp",
  
  # Lymphatic endothelium
  'Lymphatic_endothelium' = "Mrc1", 
  'Lymphatic_endothelium' = "Mmrn1", 
  'Lymphatic_endothelium' = "Lyve1",
  
  # Pericytes
  'Pericytes' = "Notch3", 
  'Pericytes' = "Rgs5",
  
  # Glial cells
  'Glial_cells' = "Nrxn1", 
  'Glial_cells' = "S100b", 
  'Glial_cells' = "Cdh19"
)


Fibro_major_marker_ls1 <- c(
  # Smooth muscle cells (SMCs)
  'SMC' = "Tagln",
  'SMC' = "Mylk",
  'SMC' = "Cnn1",#'Acta2'
  
  # Fibroblasts (FIB)
  'FIB' = "Lum",
  'FIB' = "Dcn",
  'FIB' = "Col1a1",
  'FIB' = "Col1a2",
  'FIB' = "Fn1",
  'FIB' = "C7",
  'FIB' = "Ccn2",
  
  # Myofibroblasts
  'Myofibroblasts' = "Actg2",
  'Myofibroblasts' = "Sostdc1",
  'Myofibroblasts' = "Myh11", #'Acta2',
  'Myofibroblasts' = "Notch3"
)


Idents(PT)<-'RNA_snn_res.0.3'
PT_markers <-
  FindAllMarkers(
    JoinLayers(PT),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

PT_top10_markers <- 
  PT_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10) %>%
  nest() %>%
  deframe() %>% 
  map(pull, 'gene') %>%
  imap(function(genes, name){
    re <- genes %>%
      map(function(gene){
        # PTrowser()
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(PT_top10_markers) <- sub('.*\\.','',names(PT_top10_markers))

PT_top10_markers1 <-
  PT_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10)

# dotplot
pdf(file.path(PT_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(PT, group.by = "RNA_snn_res.0.3", features = major_marker_ls)+RotatedAxis()
DotPlot(PT, group.by = "RNA_snn_res.0.3", features = celltype_marker_ls)+RotatedAxis()
DotPlot(PT, group.by = "RNA_snn_res.0.3", features = Fibro_major_marker_ls)+RotatedAxis()
DotPlot(PT, group.by = "RNA_snn_res.0.3", features = Fibro_major_marker_ls1)+RotatedAxis()
DotPlot(PT, group.by = "RNA_snn_res.0.3", features = PT_top10_markers)+RotatedAxis()
dev.off()

table(PT$RNA_snn_res.0.3,PT$group)


pdf(file.path(PT_out_dir, str_glue("Dotplot dims{dims} myofib_marker.pdf")), width = 5, height =5)
DotPlot(PT, 
        group.by = "group", 
        features = c('Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1',"Cnn1","Ccn2"))+RotatedAxis()
dev.off()


pdf(file.path(PT_out_dir, str_glue("Vlnplot dims{dims} myofib_marker.pdf")), width = 20, height = 10)
VlnPlot(PT,group.by = "group",
        feature = c('Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1',"Cnn1","Ccn2"),
        pt.size = 0,
        ncol = 4,
)
dev.off()

pdf(file.path(PT_out_dir, str_glue("Vlnplot dims{dims} myofib_marker pt.pdf")), width = 20, height = 10)
VlnPlot(PT,group.by = "group",
        feature = c('Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1',"Cnn1","Ccn2"),
        pt.size = 0.1,
        ncol = 4,
)
dev.off()

pdf(file.path(PT_out_dir, str_glue("Featureplt dims{dims} myofib_marker.pdf")), width = 20, height = 10)
FeaturePlot(PT,
            #group.by = "group",
            feature = c('Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1',"Cnn1","Ccn2"),
            reduction = 'umap',
            ncol = 4,
)

dev.off()


saveRDS(PT, file.path(PT_out_dir, "PT.rds"))
#PT <- readRDS("result/PT_subtype/PT.rds")

p_col_proj <-
  PT[[]] %>%
  group_by(group, RNA_snn_res.0.3) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count),
         RNA_snn_res.0.3 = factor(RNA_snn_res.0.3,levels = c(0:7) )) %>%
  
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.3)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

table(PT$RNA_snn_res.0.3)

ggsave("result/PT_subtype/cluster_group_prop.pdf",p_col_proj, height = 6, width = 8)

# dimplot
pdf(file.path(PT_out_dir, str_glue("Dimplot dim{dims} anno harmonyumap.pdf")), height = 6, width = 8)
DimPlot(PT, group.by = "res_0.8", label = TRUE, reduction = 'harmonyumap') 
DimPlot(PT, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') 
dev.off()

pdf(file.path(PT_out_dir, str_glue("Dimplot dim{dims} anno umap.pdf")), height = 6, width = 8)
DimPlot(PT, group.by = "res_0.8", label = TRUE, reduction = 'umap') 
DimPlot(PT, group.by = "annotation", label = TRUE, reduction = 'umap') 

dev.off()

# dimplot
pdf(file.path(PT_out_dir, str_glue("Dimplot dim{dims} cluster _split group harmonyumap.pdf")), height = 6, width = 30)
DimPlot(PT, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = 'harmonyumap',split.by = "group",pt.size = 1.5) + NoLegend()
#DimPlot(PT, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()

dev.off()

# dimplot
pdf(file.path(PT_out_dir, str_glue("Dimplot dim{dims} cluster _split group umap.pdf")), height = 6, width = 30)
DimPlot(PT, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = 'umap',split.by = "group",pt.size = 1.5) + NoLegend()
#DimPlot(PT, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()

dev.off()

pdf(file.path(PT_out_dir, str_glue("Dimplot dim{dims} cluster _split cluster umap.pdf")), height = 6, width = 30)
DimPlot(PT, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = 'harmonyumap',split.by = "RNA_snn_res.0.3",pt.size = 1.5) + NoLegend()
#DimPlot(PT, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()

dev.off()



## annotation ----
cell_cluster <- c(
  "0" = "PT",
  "1" = "PT",
  "2" = "PT",
  "3" = "tran_PT",
  "4" = "Injured_PT",
  "5" = "tran_PT",
  "6" = "Injured_PT",
  "7" = "De_PT",
  "8" = "De_PT",
  "9" = "MIX",
  "10" = "MIX"
)

Idents(PT) <- "RNA_snn_res.0.3"
#Idents(PT) <- "sub.cluster"

PT <- RenameIdents(PT, cell_cluster)
PT$subtype <- Idents(PT)

# Fib_marker_ls <- c(
#   "iCAF" = 'CXCL12',
#   "iCAF" = 'CXCL9',
#   "iCAF" = 'CXCL10',
#   "iCAF" = 'CXCL11',
#   
#   "mCAF" = 'POSTN',
#   "mCAF" = 'COL10A1',
#   "mCAF" = 'COL11A1',
#   "mCAF" = 'CTHRC1',
#   "mCAF" = 'MMP11',
#   
#   "Epi_like CAF" = 'EPCAM',
#   "Epi_like CAF" = 'CDH1',
#   "Epi_like CAF" = 'STMN1',
#   
#   "Fib" = 'CCDC80',
#   "Fib" = 'SFRP4',
#   
#   "Mac_like CAF" = 'CD68',
#   
#   "Plasma_like CAF" = 'MZB1',
#   "Plasma_like CAF" = 'XBP1',
#   
#   "meCAF" = "BNIP3" ,
#   "meCAF" ='SLC2A1' ,
#   "meCAF" ='ADM',
#   "meCAF" ='ANGPTL4'
# )
# 
# pdf(file.path(PT_out_dir, "dotplot celltype marker.pdf"), width = 20, height = 8)
# DotPlot(PT, group.by = "subtype", features = Fib_marker_ls, dot.scale = 10) + theme_90angle
# dev.off()

pdf(file.path(PT_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(PT, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
DimPlot(PT, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(PT_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 30)
DimPlot(PT, group.by = "subtype", split.by = 'group',label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
DimPlot(PT, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


## cell props barplot  --------
PT_meta <- PT[[]]

p_col_proj <-
  PT_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(PT_out_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 16
)


p_col_proj <-
  PT_meta %>%
  group_by(group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(PT_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 6
)

saveRDS(PT, file.path(PT_out_dir,"PT.rds"))

#7、Lif表达组间趋势，老师补充基因后一起做-------
Experission_dir <- "result/Gene Experission"
dir.create(Experission_dir)

pdf(file.path(Experission_dir, str_glue("Dotplot group myofib_marker.pdf")), width = 5, height =5)
DotPlot(xenium, group.by = "group", features = c('Lif'
  #'Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1'
  ))+RotatedAxis()
dev.off()


pdf(file.path(Experission_dir, str_glue("Vlnplot group myofib_marker111.pdf")), width = 20, height = 10)
VlnPlot(xenium,
        group.by = "group",
        features = c('Lif'
                     #'Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1'
                     ),
        pt.size = 0,
        ncol = 4,
)
dev.off()


pdf(file.path(Experission_dir, str_glue("Vlnplot group myofib_marker pt.pdf")), width = 20, height = 10)
VlnPlot(xenium,group.by = "group",
        features = c('Lif'
                     #'Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1'
                     ),
        pt.size = 0.1,
        ncol = 4,
)
dev.off()


xenium$annotation %>% table()
pdf(file.path(Experission_dir, str_glue("Dotplot annotation marker.pdf")), width = 5, height =5)
DotPlot(xenium, group.by = "annotation", features = c('Lif'
                                                 #'Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1'
                                                     )
        #,split.by = "group"
        )+RotatedAxis()
dev.off()

xenium$annotation %>% table()
pdf(file.path(Experission_dir, str_glue("Dotplot subtype marker.pdf")), width = 5, height =5)
DotPlot(xenium, group.by = "subtype", features = c('Lif'
                                                      #'Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1'
)
#,split.by = "group"
)+RotatedAxis()
dev.off()



#8、组间差异富集分析、损伤到纤维化的关键基因-------
xenium$annotation %>% table()

Fibro_out_dir
Fibroblast <- readRDS(file.path(Fibro_out_dir,'Fibroblast.rds'))
table(Fibroblast@meta.data$subtype)

Endothelial <- readRDS(file.path(Endo_out_dir,'Endothelial.rds'))
table(Endothelial@meta.data$subtype)

PT <- readRDS(file.path(PT_out_dir,'PT.rds'))
table(PT@meta.data$subtype)


##映射—add亚类注释到all_data----
library(tibble)
library(data.table)
subtype <- rbind(Fibroblast[['subtype']] %>% rownames_to_column('cell_id_row'), 
                 Endothelial[['subtype']] %>% rownames_to_column('cell_id_row'),
                 PT[['subtype']]%>% rownames_to_column('cell_id_row'))

xenium@meta.data %>% head()
xenium_meta <- xenium@meta.data

meta.data.t <-
  subtype%>%
  right_join(
    xenium_meta[,c("cell_id_row", "annotation")] %>% rownames_to_column('cells'),
    by = "cell_id_row"
  ) %>% 
  mutate(
    subtype = ifelse((as.character(annotation) %in% c("Fibroblast","Proliferate", "Endothelial", "PT")),
                     as.character(subtype), 
                     as.character(annotation)))%>%
  dplyr::select(cells, subtype) %>%
  column_to_rownames('cells')


xenium<-AddMetaData(xenium,meta.data.t['subtype'])
xenium$annotation %>% table()
xenium$subtype %>% table()

table(xenium$subtype, xenium$annotation, useNA = "ifany")

xenium$subtype <- gsub("MIX","Mix",xenium$subtype)

saveRDS(xenium, 'result/xenium_subtype.rds')

xenium <- readRDS('result/xenium_subtype.rds')
##8.1整体FinderMarker-------
library(tidyverse)
library(Seurat)
library(patchwork)
library(stringr)
library("clusterProfiler")
# library("org.Hs.eg.db")
library("org.Mm.eg.db")
library(RColorBrewer)
library(EnhancedVolcano)
library(dplyr)
library(gtools)

results_FindMarker <- "result/FindMarker/All"
dir.create(results_FindMarker)

# groups <- unique(xenium@meta.data$group)
# all_combns_ls <- gtools::permutations(
#   n = length(groups),
#   r = 2,
#   v = seq_along(groups),
#   repeats.allowed = FALSE
# ) %>% apply(1, c, simplify = FALSE)
# 
# all_combns_ls[[1]]

group_order <- c("Con", "H1", "H3", "H12", "H24", "H72", "D7")
comparisons_list <- combn(group_order, m = 2, simplify = FALSE)

comparisons_list %>% 
  walk(function(grp) {
    print(grp)
    print(grp[1])
    print(grp[2])
    
    control <- grp[1]  
    treat <- grp[2]  
    
    diff_markers <- FindMarkers(
      #JoinLayers(all_data_DG),
      xenium,
      ident.1 = treat, #处理组
      ident.2 = control, #对照组
      test.use = "wilcox",
      logfc.threshold = 0.25,
      min.pct = 0.1,
      group.by = "group"
    )
    if (nrow(diff_markers) == 0) {
      cat("No DEGs found for", celltype, ". Skipping enrichment analysis.\n")
      return(NULL)
    }
    
    diffmarkers_up <- diff_markers %>%
      filter(avg_log2FC > 0) %>%
      arrange(desc(avg_log2FC))
    diffmarkers_down <- diff_markers %>%
      filter(avg_log2FC < 0) %>%
      arrange(desc(avg_log2FC))
    
    diff_outdir <- file.path(results_FindMarker, str_glue("{grp[2]} vs {grp[1]}"))
    dir.create(diff_outdir, recursive = TRUE)
    write.csv(diff_markers, file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]} group diffmarker.csv")))
    write.csv(diffmarkers_up, file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]} group diffmarker up.csv")))
    write.csv(diffmarkers_down, file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]} group diffmarker down.csv")))
    
    p_volcano <- EnhancedVolcano(diff_markers,
                                 lab = rownames(diff_markers),
                                 x = 'avg_log2FC',
                                 y = 'p_val_adj',
                                 title = str_glue('volcano plot'),
                                 pCutoff = 10e-2,
                                 FCcutoff = 0,
                                 pointSize = 3.0,
                                 labSize = 6.0)
    ggsave(p_volcano,
           filename = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]} diffgene valcano.pdf")),
           height = 6,
           width = 8)
  
  #})
    
    # GO up
    enrichGO_up_ls <- enrichGO(
      rownames(diffmarkers_up),
      OrgDb = org.Mm.eg.db,
      ont = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      keyType = "SYMBOL"
    )
    write.csv(enrichGO_up_ls@result, file = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]} Go up enriment table.csv")), row.names = FALSE)
    
    if (nrow(enrichGO_up_ls) > 0) {
      enrichGO_up_dot <- enrichplot::dotplot(
        enrichGO_up_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichGO_up_dot,
        filename = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]} Go up dotplot.pdf")),
        height = 6,
        width = 8
      )
      enrichGO_up_bar <- barplot(
        enrichGO_up_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichGO_up_bar,
        filename = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]} Go up barplot.pdf")),
        height = 6,
        width = 8
      )
    } else {
      cat(paste0(grp[2], " vs ", grp[1], " has No enriched GO up terms found. Skipping plot generation.\n"))
    }
    
    # GO down
    enrichGO_down_ls <- enrichGO(
      rownames(diffmarkers_down),
      OrgDb = org.Mm.eg.db,
      ont = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      keyType = "SYMBOL"
    )
    write.csv(enrichGO_down_ls@result, file = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]} Go down enriment table.csv")), row.names = FALSE)
    
    if (nrow(enrichGO_down_ls) > 0) {
      enrichGO_down_dot <- enrichplot::dotplot(
        enrichGO_down_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichGO_down_dot,
        filename = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]} Go down dotplot.pdf")),
        height = 6,
        width = 8
      )
      enrichGO_down_bar <- barplot(
        enrichGO_down_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichGO_down_bar,
        filename = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]} Go down barplot.pdf")),
        height = 6,
        width = 8
      )
    } else {
      cat(paste0(grp[2], " vs ", grp[1], " has No enriched GO down terms found. Skipping plot generation.\n"))
    }
    
    # KEGG up
    de_up_genes_id_for_enrichmenet <-
      clusterProfiler::bitr(
        rownames(diffmarkers_up),
        OrgDb = "org.Mm.eg.db",
        fromType = "SYMBOL",
        toType = "ENTREZID"
      ) %>%
      .$ENTREZID
    enrichKEGG_up_ls <- enrichKEGG(
      de_up_genes_id_for_enrichmenet,
      organism = "mmu",
      keyType = "ncbi-geneid",
      # keyType = "kegg",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
      # use_internal_data = TRUE
    )
    enrichKEGG_up_ls <- setReadable(enrichKEGG_up_ls, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    write.csv(enrichKEGG_up_ls@result, file = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]} KEGG up enriment table.csv")), row.names = FALSE)
    
    if (nrow(enrichKEGG_up_ls) > 0) {
      enrichKEGG_up_dot <- enrichplot::dotplot(
        enrichKEGG_up_ls,
        showCategory = 10,
        label_format = 50
      )
      enrichKEGG_up_bar <- barplot(
        enrichKEGG_up_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichKEGG_up_dot,
        filename = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]} KEGG up dotplot.pdf")),
        height = 6,
        width = 8
      )
      ggsave(
        enrichKEGG_up_bar,
        filename = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]} KEGG up barplot.pdf")),
        height = 6,
        width = 8
      )
    } else {
      cat(paste0(grp[2], " vs ", grp[1], " has No enriched KEGG up terms found. Skipping plot generation.\n"))
    }
    
    # KEGG down
    de_down_genes_id_for_enrichmenet <-
      clusterProfiler::bitr(
        rownames(diffmarkers_down),
        OrgDb = "org.Mm.eg.db",
        fromType = "SYMBOL",
        toType = "ENTREZID"
      ) %>%
      .$ENTREZID
    enrichKEGG_down_ls <- enrichKEGG(
      de_down_genes_id_for_enrichmenet,
      organism = "mmu",
      keyType = "ncbi-geneid",
      # keyType = "kegg",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
      # use_internal_data = TRUE
    )
    enrichKEGG_down_ls <- setReadable(enrichKEGG_down_ls, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    write.csv(enrichKEGG_down_ls@result, file = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]} KEGG down enriment table.csv")), row.names = FALSE)
    
    if (nrow(enrichKEGG_down_ls) > 0) {
      enrichKEGG_down_dot <- enrichplot::dotplot(
        enrichKEGG_down_ls,
        showCategory = 10,
        label_format = 50
      )
      enrichKEGG_down_bar <- barplot(
        enrichKEGG_down_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichKEGG_down_dot,
        filename = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]} KEGG down dotplot.pdf")),
        height = 6,
        width = 8
      )
      ggsave(
        enrichKEGG_down_bar,
        filename = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]} KEGG down barplot.pdf")),
        height = 6,
        width = 8
      )
    } else {
      cat(paste0(grp[2], " vs ", grp[1], " has No enriched KEGG down terms found. Skipping plot generation.\n"))
    }
  })


##8.2各细胞FinderMarker-----
library(tidyverse)
library(Seurat)
library(patchwork)
library(stringr)
library("clusterProfiler")
# library("org.Hs.eg.db")
library("org.Mm.eg.db")
library(RColorBrewer)
library(EnhancedVolcano)
library(dplyr)
library(gtools)

results_FindMarker <- "result/FindMarker/annotation/"
dir.create(results_FindMarker)

group_order <- c("Con", "H1", "H3", "H12", "H24", "H72", "D7")
comparisons_list <- combn(group_order, m = 2, simplify = FALSE)

#cell <- "Fibroblast"
lapply(unique(xenium$annotation), function(cell){ 
   
   FM_data <- xenium[,xenium$annotation == cell]

   cell_dir <- file.path(results_FindMarker,paste0(cell,"_out"))
   dir.create(cell_dir)
   
   comparisons_list %>% 
     walk(function(grp) {
       print(grp)
       print(grp[1])
       print(grp[2])

       control <- grp[1]  
       treat <- grp[2]  
       
       diff_markers <- FindMarkers(
         #JoinLayers(all_data_DG),
         xenium,
         ident.1 = treat, #处理组
         ident.2 = control, #对照组
         test.use = "wilcox",
         logfc.threshold = 0.25,
         min.pct = 0.1,
         group.by = "group"
       )
       if (nrow(diff_markers) == 0) {
         cat("No DEGs found for", celltype, ". Skipping enrichment analysis.\n")
         return(NULL)
       }
       
       diffmarkers_up <- diff_markers %>%
         filter(avg_log2FC > 0) %>%
         arrange(desc(avg_log2FC))
       diffmarkers_down <- diff_markers %>%
         filter(avg_log2FC < 0) %>%
         arrange(desc(avg_log2FC))
       
       diff_outdir <- file.path(cell_dir, str_glue("{grp[2]} vs {grp[1]}"))
       dir.create(diff_outdir, recursive = TRUE)
       write.csv(diff_markers, file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]} group diffmarker.csv")))
       write.csv(diffmarkers_up, file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]}  group diffmarker up.csv")))
       write.csv(diffmarkers_down, file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]}  group diffmarker down.csv")))
       
       p_volcano <- EnhancedVolcano(diff_markers,
                                    lab = rownames(diff_markers),
                                    x = 'avg_log2FC',
                                    y = 'p_val_adj',
                                    title = str_glue('volcano plot'),
                                    pCutoff = 10e-2,
                                    FCcutoff = 0,
                                    pointSize = 3.0,
                                    labSize = 6.0)
       ggsave(p_volcano,
              filename = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]} diffgene valcano.pdf")),
              height = 6,
              width = 8)
       
#      })
# })
   
       # GO up
       enrichGO_up_ls <- enrichGO(
         rownames(diffmarkers_up),
         OrgDb = org.Mm.eg.db,
         ont = "ALL",
         pAdjustMethod = "BH",
         pvalueCutoff = 0.05,
         qvalueCutoff = 0.05,
         keyType = "SYMBOL"
       )
       write.csv(enrichGO_up_ls@result, file = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]}  Go up enriment table.csv")), row.names = FALSE)
       
       
       if (nrow(enrichGO_up_ls) > 0) {
         enrichGO_up_dot <- enrichplot::dotplot(
           enrichGO_up_ls,
           showCategory = 10,
           label_format = 50
         )
         ggsave(
           enrichGO_up_dot,
           filename = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]}  Go up dotplot.pdf")),
           height = 6,
           width = 8
         )
         enrichGO_up_bar <- barplot(
           enrichGO_up_ls,
           showCategory = 10,
           label_format = 50
         )
         ggsave(
           enrichGO_up_bar,
           filename = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]}  Go up barplot.pdf")),
           height = 6,
           width = 8
         )
       } else {
         cat(paste0(grp[2], " vs ", grp[1], " has No enriched GO up terms found. Skipping plot generation.\n"))
       }
       
       
       # GO down
       enrichGO_down_ls <- enrichGO(
         rownames(diffmarkers_down),
         OrgDb = org.Mm.eg.db,
         ont = "ALL",
         pAdjustMethod = "BH",
         pvalueCutoff = 0.05,
         qvalueCutoff = 0.05,
         keyType = "SYMBOL"
       )
       write.csv(enrichGO_down_ls@result, file = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]}  Go down enriment table.csv")), row.names = FALSE)
       
       
       if (nrow(enrichGO_down_ls) > 0) {
         enrichGO_down_dot <- enrichplot::dotplot(
           enrichGO_down_ls,
           showCategory = 10,
           label_format = 50
         )
         ggsave(
           enrichGO_down_dot,
           filename = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]}  Go down dotplot.pdf")),
           height = 6,
           width = 8
         )
         enrichGO_down_bar <- barplot(
           enrichGO_down_ls,
           showCategory = 10,
           label_format = 50
         )
         ggsave(
           enrichGO_down_bar,
           filename = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]}  Go down barplot.pdf")),
           height = 6,
           width = 8
         )
       } else {
         cat(paste0(grp[2], " vs ", grp[1], " has No enriched GO down terms found. Skipping plot generation.\n"))
       }
       
       
       # KEGG up
       de_up_genes_id_for_enrichmenet <-
         clusterProfiler::bitr(
           rownames(diffmarkers_up),
           OrgDb = "org.Mm.eg.db",
           fromType = "SYMBOL",
           toType = "ENTREZID"
         ) %>%
         .$ENTREZID
       enrichKEGG_up_ls <- enrichKEGG(
         de_up_genes_id_for_enrichmenet,
         organism = "mmu",
         keyType = "ncbi-geneid",
         # keyType = "kegg",
         pvalueCutoff = 0.05,
         pAdjustMethod = "BH",
         minGSSize = 10,
         maxGSSize = 500,
         qvalueCutoff = 0.2,
         use_internal_data = FALSE
         # use_internal_data = TRUE
       )
       enrichKEGG_up_ls <- setReadable(enrichKEGG_up_ls, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
       write.csv(enrichKEGG_up_ls@result, file = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]}  KEGG up enriment table.csv")), row.names = FALSE)
       
       
       if (nrow(enrichKEGG_up_ls) > 0) {
         enrichKEGG_up_dot <- enrichplot::dotplot(
           enrichKEGG_up_ls,
           showCategory = 10,
           label_format = 50
         )
         enrichKEGG_up_bar <- barplot(
           enrichKEGG_up_ls,
           showCategory = 10,
           label_format = 50
         )
         ggsave(
           enrichKEGG_up_dot,
           filename = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]}  KEGG up dotplot.pdf")),
           height = 6,
           width = 8
         )
         ggsave(
           enrichKEGG_up_bar,
           filename = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]}  KEGG up barplot.pdf")),
           height = 6,
           width = 8
         )
       } else {
         cat(paste0(grp[2], " vs ", grp[1], " has No enriched KEGG up terms found. Skipping plot generation.\n"))
       }
       
       
       # KEGG down
       de_down_genes_id_for_enrichmenet <-
         clusterProfiler::bitr(
           rownames(diffmarkers_down),
           OrgDb = "org.Mm.eg.db",
           fromType = "SYMBOL",
           toType = "ENTREZID"
         ) %>%
         .$ENTREZID
       enrichKEGG_down_ls <- enrichKEGG(
         de_down_genes_id_for_enrichmenet,
         organism = "mmu",
         keyType = "ncbi-geneid",
         # keyType = "kegg",
         pvalueCutoff = 0.05,
         pAdjustMethod = "BH",
         minGSSize = 10,
         maxGSSize = 500,
         qvalueCutoff = 0.2,
         use_internal_data = FALSE
         # use_internal_data = TRUE
       )
       enrichKEGG_down_ls <- setReadable(enrichKEGG_down_ls, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
       write.csv(enrichKEGG_down_ls@result, file = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]}  KEGG down enriment table.csv")), row.names = FALSE)
       
       
       if (nrow(enrichKEGG_down_ls) > 0) {
         enrichKEGG_down_dot <- enrichplot::dotplot(
           enrichKEGG_down_ls,
           showCategory = 10,
           label_format = 50
         )
         enrichKEGG_down_bar <- barplot(
           enrichKEGG_down_ls,
           showCategory = 10,
           label_format = 50
         )
         ggsave(
           enrichKEGG_down_dot,
           filename = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]}  KEGG down dotplot.pdf")),
           height = 6,
           width = 8
         )
         ggsave(
           enrichKEGG_down_bar,
           filename = file.path(diff_outdir, str_glue("{grp[2]} vs {grp[1]}  KEGG down barplot.pdf")),
           height = 6,
           width = 8
         )
       } else {
         cat(paste0(grp[2], " vs ", grp[1], " has No enriched KEGG down terms found. Skipping plot generation.\n"))
       }
     })
   
})



##8.3Mfuzz-------

#在徐帅docker里-------
xenium$annotation %>% table()
xenium$subtype %>% table()
#聚类图
library(Mfuzz)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggstatsplot)
library(tidyverse)

Mfuzzdata <- xenium[,xenium$subtype%in% c("Myofibroblast")]

expr_mat <- Mfuzzdata@assays$Xenium@layers$data

rownames(expr_mat) <- rownames(Mfuzzdata)
colnames(expr_mat) <- colnames(Mfuzzdata)

cell_meta <- Mfuzzdata@meta.data[, "group", drop = FALSE]

# 2. 转长数据，按基因+分组汇总均值
expr_df <- as.data.frame(as.matrix(expr_mat)) %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "cell_id",
    values_to = "expr"
  ) %>%
  left_join(rownames_to_column(cell_meta, "cell_id"), by = "cell_id")

expr_df[6000:6010,1:4]

group_mean_df <- expr_df %>%
  group_by(gene, group) %>%
  summarise(expr_mean = mean(expr), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = expr_mean) %>%
  column_to_rownames("gene")

# 1. 数据预处理：去除表达量低或在不同时间点变化太小的基因
# Mfuzz 聚类要求输入数据为 ExpressionSet 类型的对象，因此需先构建此对象。

# 将数据框转换为矩阵
count <- data.matrix(group_mean_df)  # 输入数据转换为矩阵格式，作为构建 ExpressionSet 对象的基础

# 创建 ExpressionSet 对象
eset <- new("ExpressionSet", exprs = count)

eset <- filter.std(eset, min.std = 0) 
eset <- standardise(eset)
