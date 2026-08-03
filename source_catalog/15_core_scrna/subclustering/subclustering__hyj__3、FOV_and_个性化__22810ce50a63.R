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

#all_data2 <- readRDS("R_OUTPUT/all_data.rds") 
all_data <-  JoinLayers(all_data2,layers = 'counts')
colnames(all_data@meta.data)

all_data@meta.data$cell_id %>% head()
rownames(all_data@meta.data) %>% head()

###改行名
all_data@meta.data$cell_id <- gsub("^.+?_(?=[a-z]{8}-\\d+$)", "", 
                                     rownames(all_data@meta.data), 
                                     perl = TRUE)
T1<- data_frame(rownames(all_data@meta.data),all_data@meta.data[c("orig.ident","cell_id")])

all_data@meta.data$cell_id2 <- rownames(all_data@meta.data)

rownames(all_data@meta.data) <- paste0(all_data@meta.data$sample,"_",all_data@meta.data$cell_id)

T1<- data_frame(rownames(all_data@meta.data),all_data@meta.data[c("orig.ident", "cell_id","cell_id2")])

rownames(all_data) %>% head()
colnames(all_data) %>% head()

rownames(all_data@assays$RNA@layers$count) <- rownames(all_data)
colnames(all_data@assays$RNA@layers$count) <- colnames(all_data)
xenium <- CreateSeuratObject(
  all_data@assays$RNA@layers$count,
  assay = "Xenium")

xenium@meta.data <- all_data@meta.data

umap <- all_data@reductions$umap@cell.embeddings
rownames(umap) <- rownames(xenium@meta.data)
metadata <- xenium@meta.data
xenium <- Add_meta_umapData_Xenium(object = xenium,umapfile = umap,metedata = metadata)
colnames(xenium@assays$Xenium) <-  colnames(xenium)
row.names(xenium@reductions$umap@cell.embeddings) <- colnames(xenium)
xenium <- RenameCells(xenium,new.names = colnames(xenium))

#saveRDS(xenium,"xenium_renames.rds")
#xenium_renames <- readRDS("xenium_renames.rds")

has_dt <- requireNamespace("data.table", quietly = TRUE) && requireNamespace("R.utils", quietly = TRUE)
data.dir = "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/XE2510316144-02/XE2510316144-02/output-XETG00305__0086346__TMA__20251226__065149"
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


options(future.globals.maxSize = 20 * 1024 ^ 3)  # 20GB
plan(multisession, workers = 10)

xenium <- NormalizeData(xenium)

saveRDS(xenium,"result/xenium_FOV.rds")

###测试一个样本----
DefaultBoundary(xenium@images[["Shu_yu2"]]) <- "segmentation"
p <- ImageFeaturePlot(
  object = xenium,
  fov = "Shu_yu2",
  nmols = 20000,       
  border.color = "white",  
  border.size = 0.1,    
  features = 'INHBA',
  max.cutoff = "q90",
  axes = FALSE,
  cols = c('#003366', "#E7475E"),
  size = 1.5,
  dark.background = T
)
ggsave(p,filename = paste0("Shu_yu2","_feature_spatial.png"),width = 12 ,height = 12,dpi = 300)


###Imagedimplot----
xenium@meta.data$annotation <- as.factor(xenium@meta.data$annotation)

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
  "#FF6B6B","#A6CEE3",
  "#F15BB5","#00F5A7", "#FDBB69","#feee7d", #"#9E8099",
  "#00BBF9","#7209B7"#,"#9B5DE5"
)

names(all_colors8) <- levels(xenium@meta.data$annotation)

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
  ggsave(p,filename = paste0("result/Imagedimplot/",i,"_anno_spatial.pdf"),width = 12 ,height = 12)
  ggsave(p,filename = paste0("result/Imagedimplot/",i,"_anno_spatial.png"),width = 12 ,height = 12,dpi = 300)
}


#1、炎性非炎性对比------
##柱状图-----
p_col_group <-
  xenium@meta.data %>%
  group_by(group,annotation) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = group)) + 
  facet_grid(~annotation)+
  geom_col() +
  #ggsci::scale_fill_igv() +
  scale_fill_manual(values = c("#4DBBD5", "#E64B35"))+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_group,
  filename = file.path("result/cell_freq_group.pdf"),
  height = 6,
  width  = 12
)
ggsave(
  p_col_group,
  filename = file.path("result/cell_freq_group.png"),
  height = 6,
  width  = 12
)


#2、T、B、Mac、Plasma是否有差异------
## 箱线图------
boxplot_data <-  xenium@meta.data %>%
  group_by(sample, group, annotation) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))  

types<- unique(boxplot_data$annotation)

lapply(types, function(type){
    p <- boxplot_data %>%
      filter(annotation == type) %>% 
      ggplot(aes(x = group, y = freq, fill = group)) +
      geom_boxplot(width = 0.6) +
      #facet_wrap(~annotation, scales = "free_y", ncol = 4) +
      stat_compare_means(
        comparisons = list(c("Infla","Noninfla")),  
        method = "t.test",  
        label = "p.signif",  # p.signif, "p.value"
        p.adjust.method = "none"  
      ) +
      
      theme_bw() +
      labs(
        title = "",
        x = "group",
        y = "freq"
      ) +
      scale_fill_manual(values = c("#4DBBD5", "#E64B35"))  #
    
    print(p)
    
    ggsave(file.path(result_dir,paste0("group_boxplot/",type,"_group_boxplot.pdf")), p, width =5, height = 5)
    ggsave(file.path(result_dir,paste0("group_boxplot/",type,"_group_boxplot.png")), p, width =5, height = 5)

})


#3、T、Mac亚类是否有差异（T暂未注释出来）-------
## 髓系亚类箱线图------
MYELOID_meta <- MYELOID[[]]

boxplot_data <-  MYELOID_meta %>%
  group_by(sample, group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))  

types<- unique(boxplot_data$subtype)

lapply(types, function(type){
  p <- boxplot_data %>%
    filter(subtype == type) %>% 
    ggplot(aes(x = group, y = freq, fill = group)) +
    geom_boxplot(width = 0.6) +
    #facet_wrap(~annotation, scales = "free_y", ncol = 4) +
    stat_compare_means(
      comparisons = list(c("Infla","Noninfla")),  
      method = "t.test",  
      label = "p.signif",  # p.signif, "p.value"
      p.adjust.method = "none"  
    ) +
    
    theme_bw() +
    labs(
      title = "",
      x = "group",
      y = "freq"
    ) +
    scale_fill_manual(values = c("#4DBBD5", "#E64B35"))  #
  
  print(p)
  
  ggsave(file.path(MYELOID_out_dir,paste0("group_boxplot/",type,"_group_boxplot.pdf")), p, width =5, height = 5)
  ggsave(file.path(MYELOID_out_dir,paste0("group_boxplot/",type,"_group_boxplot.png")), p, width =5, height = 5)
  
})


#4、客户genelist表达情况、pathway评分、组间对比------
gene_sets <- c(
  # 补体
  "Complement_producing" = "C1QA",
  "Complement_producing" = "C1QB", 
  "Complement_producing" = "C1QC",
  "Complement_producing" = "C2",
  "Complement_producing" = "C3",
  "Complement_producing" = "C3AR1",
  
  # 损伤
  "DAMP_signaling" = "HMGB1",
  "DAMP_signaling" = "HSPA1A",
  "DAMP_signaling" = "S100A8",
  "DAMP_signaling" = "S100A9", 
  "DAMP_signaling" = "S100A12",
  
  # 模式识别受体
  "PRR_expressing" = "TLR1",
  "PRR_expressing" = "TLR5",
  "PRR_expressing" = "TLR7",
  "PRR_expressing" = "IFIH1",
  "PRR_expressing" = "DHX58",
  
  "IG" = "IGHG1",
  "IG" = "IGHG2",
  "IG" = "IGHG3", 
  "IG" = "IGHG4",
  "IG" = "IGHM",
  "IG" = "IGKC",
  "IG" = "IGLC2",
  "IG" = "IGHA1",
  "IG" = "IGHA2",
  
  "IG_Receptor" = "FCGR1A",
  "IG_Receptor" = "FCGR2A", 
  "IG_Receptor" = "FCGR3A",
  
  "M1_Macrophage" = "IL12",
  "M1_Macrophage" = "IL23",
  "M1_Macrophage" = "TNF",
  "M1_Macrophage" = "IL6",
  "M1_Macrophage" = "CD86",
  "M1_Macrophage" = "MHCII",
  "M1_Macrophage" = "IL1B",
  "M1_Macrophage" = "MARCO",
  "M1_Macrophage" = "NOS2",
  "M1_Macrophage" = "CD64",
  "M1_Macrophage" = "CD80",
  "M1_Macrophage" = "CXCR10",
  "M1_Macrophage" = "CXCL9",
  "M1_Macrophage" = "CXCL10",
  "M1_Macrophage" = "CXCL11",
  "M1_Macrophage" = "IL1A",
  "M1_Macrophage" = "CCL5",
  "M1_Macrophage" = "IRF5",
  "M1_Macrophage" = "IRF1",
  "M1_Macrophage" = "CD40",
  "M1_Macrophage" = "IDO1",
  "M1_Macrophage" = "KYNU",
  "M1_Macrophage" = "CCR7",
  "M1_Macrophage" = "PTPRC",
  "M1_Macrophage" = "CD68",
  "M1_Macrophage" = "CD115",
  "M1_Macrophage" = "HLA-DRA",
  "M1_Macrophage" = "CD205",
  "M1_Macrophage" = "CD14",
  
  "Exhausted" = "PDCD1",
  "Exhausted" = "LAG3",
  "Exhausted" = "HAVCR2",
  "Exhausted" = "CD244",
  "Exhausted" = "CD160",
  "Exhausted" = "IL10R",
  "Exhausted" = "EOMES",
  "Exhausted" = "NR4A2",
  "Exhausted" = "PTGER4",
  "Exhausted" = "TOX",
  "Exhausted" = "TOX2",
  "Exhausted" = "TIGIT",
  "Exhausted" = "CTLA4",
  "Exhausted" = "ENTPD1"
  
)


valid_gene_sets <- gene_sets[gene_sets %in% rownames(xenium)]
cat("有效基因:", length(valid_gene_sets), "/", length(gene_sets))

missing_genes <- setdiff(gene_sets, valid_gene_sets)


###基因集---
modules <- names(valid_gene_sets) |> unique()
modules

colnames(xenium@meta.data)
xenium_score <- AddModuleScore(object = xenium,
                               features = split(valid_gene_sets, modules),
                               name = "Module_",
                               ctrl = 50,             # 对照基因数
                               seed = 2025)           # 确保可重复

# # 批量删除这些列
# module_cols <- c("M1", "M2", "TAM", "Pha", "Met")
# for (col in module_cols) {
#   xenium@meta.data[[col]] <- NULL
# }

# 新旧列名映射关系
new_names <- c("Complement_producing","PRR_expressing","IG_Receptor","M1_Macrophage","Exhausted")
old_names <- paste0("Module_", 1:5)
colnames(xenium_score@meta.data)[match(old_names, colnames(xenium_score@meta.data))] <- new_names
colnames(xenium_score@meta.data)

#测试
Idents(xenium_score) <- "annotation"
ImageDimPlot(xenium_score, fov = "LYS",nmols = 20000,cols =all_colors8,
             cells = WhichCells(xenium, idents = "Macrophage" )
             )

###评分空间featureplt----
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



##评分箱线图（带检验）----
library(Seurat)
library(ggpubr)
library(tidyverse)

ModuleBox_dir <- 'result/ModuleBox' 
dir.create(ModuleBox_dir)

plot_data <- xenium_score@meta.data %>%
  dplyr::select(all_of(new_names), group) %>%
  pivot_longer(cols = all_of(new_names), 
               names_to = "Module", 
               values_to = "Score") %>%
  mutate(Module = factor(Module, levels = new_names))

lapply(new_names, function(module){
  module_data <- plot_data %>%
    filter(Module == module)
  
  y_max <- max(module_data$Score, na.rm = TRUE) * 1.15
  
  p <- ggplot(module_data, aes(x = group, y = Score, fill = group)) +
    geom_boxplot(
      width = 0.5,
      alpha = 0.85,
      outlier.shape = 19,
      outlier.size = 0.8,
      outlier.alpha = 0.6
    ) +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",
      label.x = 1.5,
      size = 5,
      label.y = y_max,
      vjust = -0.5
    ) +
    labs(
      title = module,
      x = "Group",
      y = "Module Score",
      fill = "Group"
    ) +
    theme_bw() +
    scale_fill_manual(values = c("#4DBBD5", "#E64B35"))
  
  
  ggsave(
    filename = file.path(ModuleBox_dir, paste0('group_', module, "_boxplot.png")),
    plot = p,
    width = 6,
    height = 7,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(ModuleBox_dir, paste0('group_', module, "_boxplot.pdf")),
    plot = p,
    width = 6,
    height = 7
  )
  
  message("已保存: ", module, " 箱线图")

  }
)
